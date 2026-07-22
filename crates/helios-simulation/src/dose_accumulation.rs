//! Delivered-dose accumulation: ray-trace every [`DeliveryFrame`] beamlet into a
//! dose [`Volume`].
//!
//! Closes the delivery→dose loop. Each frame carries the machine state (gantry
//! angle, couch position) and the effective per-leaf fluence actually delivered
//! (leakage + tongue-and-groove already applied by the MLC model). This kernel
//! turns that time-ordered fluence into a spatial dose distribution by depositing
//! each leaf's beamlet terma through the attenuation volume, summed over all
//! frames — the input the DVH / gamma analysis consumes.
//!
//! # Beam geometry
//! A helical TomoTherapy fan: at gantry angle `θ` the beam travels along the
//! axial-plane direction `d = (cosθ, sinθ, 0)`; each binary-MLC leaf is a beamlet
//! laterally offset along the in-plane perpendicular `p = (−sinθ, cosθ, 0)` by
//! `(leaf − centre)·leaf_width`, at the couch `z` slice. [`BeamGeometry`] selects
//! whether the beamlets run parallel (small-fan approximation) or diverge from a
//! point source (true fan, with inverse-square fluence falloff).
//!
//! [`accumulate_delivered_dose`] returns the pooled **terma**; a beam-aligned
//! anisotropic collapsed cone is available via
//! [`accumulate_delivered_dose_anisotropic`], which scatters each frame's terma
//! along that frame's own gantry direction (the forward-peaked physics follows
//! the rotating beam). The [`CollapsedCone`] kernel is monoenergetic
//! ([`CollapsedCone::forward_peaked`]) or poly-energetic / beam-hardened
//! ([`CollapsedCone::poly_forward_peaked`]). Per-leaf gaia collimation remains a
//! later increment.

use crate::delivery::DeliveryFrame;
use helios_domain::Volume;
use helios_math::{GeometryScalar, NumericElement, Point3, Ray, Vector3};
use helios_solver::{
    deposit_ray_terma, deposit_ray_terma_diverging, forward_peaked_kernel,
    oriented_forward_scatter, poly_forward_peaked_kernel, symmetric_deposition_kernel,
    SpectralComponent,
};
use hyperion::TransportError;

/// Beam geometry for delivered-dose accumulation — the seam that selects how each
/// MLC leaf's beamlet ray is cast for a gantry angle.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum BeamGeometry<T: GeometryScalar> {
    /// Parallel beamlets (small-fan approximation): every leaf ray runs along the
    /// gantry direction, offset laterally; `standoff_mm` places the origin behind
    /// isocentre. Cheap and exact for a narrow field.
    Parallel {
        /// Distance the beamlet origin stands off behind isocentre (mm).
        standoff_mm: T,
    },
    /// Divergent fan from a point source at `source_axis_mm` from isocentre (SAD):
    /// each beamlet runs from the focal spot through its isocentre-plane offset
    /// point, so beamlets diverge with depth — the true TomoTherapy fan geometry.
    /// Reduces to [`Parallel`](Self::Parallel) as `source_axis_mm → ∞`.
    PointSource {
        /// Source-to-axis distance / SAD (mm).
        source_axis_mm: T,
    },
}

/// Accumulate the delivered dose from a helical-delivery `frames` sequence into a
/// dose [`Volume`] over the same grid as the attenuation volume `mu`.
///
/// `geometry` selects the beam model (parallel vs divergent point-source fan);
/// `leaf_width_mm` is the inter-leaf lateral pitch; `step_mm` is the ray-march
/// sampling step. Dose is linear in the per-leaf fluence, so scaling all fluence
/// scales the dose and independent frames/leaves superpose (the test oracles).
///
/// # Errors
///
/// Returns [`TransportError`] when a sampled attenuation coefficient violates
/// Hyperion's transport contract. Validation is transactional per beamlet.
pub fn accumulate_delivered_dose<T: GeometryScalar>(
    frames: &[DeliveryFrame<T>],
    mu: &Volume<T>,
    geometry: BeamGeometry<T>,
    leaf_width_mm: T,
    step_mm: T,
) -> Result<Volume<T>, TransportError<T>> {
    let grid = *mu.grid();
    let mut dose = Volume::zeros(grid);
    for frame in frames {
        deposit_frame_terma(&mut dose, frame, mu, geometry, leaf_width_mm, step_mm)?;
    }
    Ok(dose)
}

/// Deposit one `frame`'s per-leaf beamlet terma into `dose`, returning the beam's
/// forward unit direction (the gantry central axis) — the SSOT deposition loop
/// shared by the isotropic accumulation and the per-frame anisotropic path, so
/// the beamlet geometry lives in exactly one place.
fn deposit_frame_terma<T: GeometryScalar>(
    dose: &mut Volume<T>,
    frame: &DeliveryFrame<T>,
    mu: &Volume<T>,
    geometry: BeamGeometry<T>,
    leaf_width_mm: T,
    step_mm: T,
) -> Result<Vector3<T>, TransportError<T>> {
    let zero = <T as NumericElement>::ZERO;
    let (centre, dir, perp) = gantry_basis(mu.grid(), frame.gantry_angle_rad);
    for (leaf, &weight) in frame.leaf_fluence.iter().enumerate() {
        if weight <= zero {
            continue; // closed/leak-free leaf deposits nothing.
        }
        let Some(beamlet) = beamlet_ray(centre, dir, perp, frame, leaf, leaf_width_mm, geometry)
        else {
            continue;
        };
        let _deposited = match beamlet.falloff {
            Some((focal, sad)) => {
                deposit_ray_terma_diverging(dose, mu, &beamlet.ray, weight, step_mm, focal, sad)
            }
            None => deposit_ray_terma(dose, mu, &beamlet.ray, weight, step_mm),
        }?;
    }
    Ok(dir)
}

/// Beam-frame collapsed-cone scatter kernel for anisotropic delivered dose: a
/// forward-peaked axial kernel (secondary electrons carried downstream) plus a
/// symmetric lateral kernel across the beam, both sampled at `sample_step_mm`.
///
/// Built once and reused across every frame; the anisotropy is re-oriented to
/// each frame's gantry direction at application time.
#[derive(Debug, Clone, PartialEq)]
pub struct CollapsedCone<T: GeometryScalar> {
    beam_kernel: Vec<T>,
    beam_center: usize,
    lateral: Vec<T>,
    sample_step_mm: T,
}

impl<T: GeometryScalar> CollapsedCone<T> {
    /// Build a forward-peaked cone: exponential deposition with distinct upstream
    /// (`range_up_cm`) and downstream (`range_down_cm`) ranges along the beam and a
    /// symmetric `lateral_range_cm` across it. `voxel_cm` is the sampling pitch in
    /// cm (the trilinear step is `voxel_cm × 10` mm); the radii bound each kernel's
    /// tap support. Equal up/down ranges give an isotropic (direction-independent)
    /// cone.
    #[must_use]
    pub fn forward_peaked(
        range_up_cm: T,
        range_down_cm: T,
        lateral_range_cm: T,
        voxel_cm: T,
        radius_up: usize,
        radius_down: usize,
        lateral_radius: usize,
    ) -> Self {
        let (beam_kernel, beam_center) =
            forward_peaked_kernel(range_up_cm, range_down_cm, voxel_cm, radius_up, radius_down);
        Self::from_beam_kernel(
            beam_kernel,
            beam_center,
            lateral_range_cm,
            voxel_cm,
            lateral_radius,
        )
    }

    /// Build a **poly-energetic** forward-peaked cone: the beam kernel is the
    /// energy-fluence-weighted [`poly_forward_peaked_kernel`] of the spectral
    /// `components` (beam hardening — harder components reach farther downstream),
    /// with the same symmetric lateral spread. A single-component spectrum reduces
    /// to [`forward_peaked`](Self::forward_peaked) exactly.
    #[must_use]
    pub fn poly_forward_peaked(
        components: &[SpectralComponent<T>],
        lateral_range_cm: T,
        voxel_cm: T,
        radius_up: usize,
        radius_down: usize,
        lateral_radius: usize,
    ) -> Self {
        let (beam_kernel, beam_center) =
            poly_forward_peaked_kernel(components, voxel_cm, radius_up, radius_down);
        Self::from_beam_kernel(
            beam_kernel,
            beam_center,
            lateral_range_cm,
            voxel_cm,
            lateral_radius,
        )
    }

    /// Shared assembly: pair a prepared beam kernel with the symmetric lateral
    /// kernel and the trilinear sample step (`voxel_cm × 10` mm). SSOT for the two
    /// public constructors above.
    fn from_beam_kernel(
        beam_kernel: Vec<T>,
        beam_center: usize,
        lateral_range_cm: T,
        voxel_cm: T,
        lateral_radius: usize,
    ) -> Self {
        let lateral = symmetric_deposition_kernel(lateral_range_cm, voxel_cm, lateral_radius);
        let sample_step_mm = voxel_cm * <T as GeometryScalar>::from_f64(10.0);
        Self {
            beam_kernel,
            beam_center,
            lateral,
            sample_step_mm,
        }
    }
}

/// Accumulate delivered dose with a **beam-aligned anisotropic** collapsed cone:
/// each frame's terma is scattered along that frame's own gantry direction before
/// summing, so the forward-peaked physics follows the rotating beam (unlike a
/// single separable scatter applied to the pooled terma, which has no coherent
/// beam axis).
///
/// Identical to [`accumulate_delivered_dose`] in beamlet geometry; it adds the
/// per-frame [`oriented_forward_scatter`] stage. With an isotropic `cone` (equal
/// up/down ranges) and a single frame at gantry angle 0 it reduces to
/// [`scatter_superposition`](helios_solver::scatter_superposition) of the
/// delivered terma (the differential oracle).
///
/// # Errors
///
/// Returns [`TransportError`] under the same attenuation contract as
/// [`accumulate_delivered_dose`].
pub fn accumulate_delivered_dose_anisotropic<T: GeometryScalar>(
    frames: &[DeliveryFrame<T>],
    mu: &Volume<T>,
    geometry: BeamGeometry<T>,
    leaf_width_mm: T,
    step_mm: T,
    cone: &CollapsedCone<T>,
) -> Result<Volume<T>, TransportError<T>> {
    let grid = *mu.grid();
    let mut acc = vec![<T as NumericElement>::ZERO; grid.num_voxels()];
    for frame in frames {
        let mut frame_terma = Volume::zeros(grid);
        let dir = deposit_frame_terma(
            &mut frame_terma,
            frame,
            mu,
            geometry,
            leaf_width_mm,
            step_mm,
        )?;
        let frame_dose = oriented_forward_scatter(
            &frame_terma,
            dir,
            &cone.beam_kernel,
            cone.beam_center,
            &cone.lateral,
            cone.sample_step_mm,
        );
        for (a, d) in acc.iter_mut().zip(frame_dose.as_slice()) {
            *a += *d;
        }
    }
    Ok(Volume::from_shape_vec(grid, acc).expect("output length equals grid voxel count"))
}

/// A single MLC-leaf beamlet: its ray plus, for a divergent fan, the inverse-square
/// falloff `(focal spot, SAD)`.
pub(crate) struct Beamlet<T: GeometryScalar> {
    pub ray: Ray<T>,
    pub falloff: Option<(Point3<T>, T)>,
}

/// Construct the [`Beamlet`] for one MLC `leaf` of a `frame` under the selected
/// [`BeamGeometry`]. Shared by dose accumulation and portal dosimetry so the fan
/// geometry lives in one place.
///
/// `centre` is the grid axial centre; `dir`/`perp` the gantry basis (central-axis
/// and in-plane lateral). Returns `None` if the resulting direction is degenerate.
pub(crate) fn beamlet_ray<T: GeometryScalar>(
    centre: Point3<T>,
    dir: Vector3<T>,
    perp: Vector3<T>,
    frame: &DeliveryFrame<T>,
    leaf: usize,
    leaf_width_mm: T,
    geometry: BeamGeometry<T>,
) -> Option<Beamlet<T>> {
    let zero = <T as NumericElement>::ZERO;
    let leaves = frame.leaf_fluence.len();
    let centre_leaf = <T as GeometryScalar>::from_f64((leaves as f64 - 1.0) * 0.5);
    let offset = (<T as GeometryScalar>::from_f64(leaf as f64) - centre_leaf) * leaf_width_mm;
    // (origin, direction, inverse-square falloff) per the selected geometry. Both
    // branches lie in the couch z-slice (dir.z = perp.z = 0); the ray constructor
    // normalizes the direction vector.
    let (origin, direction, falloff) = match geometry {
        BeamGeometry::Parallel { standoff_mm } => (
            Point3::new(
                centre.x + perp.x * offset - dir.x * standoff_mm,
                centre.y + perp.y * offset - dir.y * standoff_mm,
                frame.couch_mm,
            ),
            dir,
            None,
        ),
        BeamGeometry::PointSource { source_axis_mm } => {
            // Focal spot behind isocentre; ray aims through the leaf's isocentre-
            // plane point `centre + perp·offset`, so direction = (perp·offset +
            // dir·SAD). For offset 0 this is the central axis; off-axis leaves fan
            // out with depth. Fluence falls off inverse-square from the focal spot.
            let focal = Point3::new(
                centre.x - dir.x * source_axis_mm,
                centre.y - dir.y * source_axis_mm,
                frame.couch_mm,
            );
            let aim = Vector3::new(
                perp.x * offset + dir.x * source_axis_mm,
                perp.y * offset + dir.y * source_axis_mm,
                zero,
            );
            (focal, aim, Some((focal, source_axis_mm)))
        }
    };
    Ray::try_new(origin, direction)
        .ok()
        .map(|ray| Beamlet { ray, falloff })
}

/// The gantry basis `(centre, dir, perp)` for a frame's gantry angle over `grid`.
pub(crate) fn gantry_basis<T: GeometryScalar>(
    grid: &helios_domain::VoxelGrid<T>,
    gantry_angle_rad: T,
) -> (Point3<T>, Vector3<T>, Vector3<T>) {
    let zero = <T as NumericElement>::ZERO;
    let [nx, ny, nz] = grid.dims();
    let centre = grid.voxel_center((nx - 1) / 2, (ny - 1) / 2, (nz - 1) / 2);
    let (cos, sin) = (gantry_angle_rad.cos(), gantry_angle_rad.sin());
    (
        centre,
        Vector3::new(cos, sin, zero),
        Vector3::new(-sin, cos, zero),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use eunomia::assert_relative_eq;
    use helios_domain::VoxelGrid;

    // Uniform-μ cube: 9³ voxels, 2 mm spacing → node box [0,16] mm, centre 8 mm,
    // axial chord 16 mm = 1.6 cm.
    fn uniform_cube(mu_val: f64) -> Volume<f64> {
        let grid = VoxelGrid::axis_aligned([9, 9, 9], [2.0, 2.0, 2.0], Point3::new(0.0, 0.0, 0.0))
            .expect("grid");
        Volume::from_shape_fn(grid, move |_| mu_val)
    }

    // A single-leaf frame at gantry angle θ, couch z, fluence w. One leaf → the
    // beamlet is on the central axis (zero lateral offset).
    fn frame(gantry_angle_rad: f64, couch_mm: f64, w: f64) -> DeliveryFrame<f64> {
        DeliveryFrame {
            projection: 0,
            gantry_angle_rad,
            couch_mm,
            leaf_fluence: vec![w],
        }
    }

    // Analytic primary energy removed by one central axial beamlet of weight w.
    fn expected_axial_energy(mu_val: f64, w: f64) -> f64 {
        w * (1.0 - (-mu_val * 1.6_f64).exp())
    }

    #[test]
    fn single_central_beamlet_matches_analytic_energy() {
        // θ=0 → +x beamlet through the centre (couch z = 8 mm). Total delivered
        // dose = w·(1 − e^{−μ·1.6}).
        let mu = uniform_cube(0.05);
        let dose = accumulate_delivered_dose(
            &[frame(0.0, 8.0, 2.0)],
            &mu,
            BeamGeometry::Parallel { standoff_mm: 500.0 },
            2.0,
            0.25,
        )
        .expect("valid attenuation volume");
        assert_relative_eq!(dose.sum(), expected_axial_energy(0.05, 2.0), epsilon = 1e-9);
    }

    #[test]
    fn zero_fluence_delivers_zero_dose() {
        let mu = uniform_cube(0.05);
        let dose = accumulate_delivered_dose(
            &[frame(0.0, 8.0, 0.0)],
            &mu,
            BeamGeometry::Parallel { standoff_mm: 500.0 },
            2.0,
            0.5,
        )
        .expect("valid attenuation volume");
        assert_relative_eq!(dose.sum(), 0.0, epsilon = 1e-15);
    }

    #[test]
    fn dose_is_linear_in_fluence() {
        // Doubling every leaf fluence doubles the dose voxelwise.
        let mu = uniform_cube(0.05);
        let d1 = accumulate_delivered_dose(
            &[frame(0.0, 8.0, 1.0)],
            &mu,
            BeamGeometry::Parallel { standoff_mm: 500.0 },
            2.0,
            0.25,
        )
        .expect("valid attenuation volume");
        let d2 = accumulate_delivered_dose(
            &[frame(0.0, 8.0, 2.0)],
            &mu,
            BeamGeometry::Parallel { standoff_mm: 500.0 },
            2.0,
            0.25,
        )
        .expect("valid attenuation volume");
        assert_relative_eq!(d2.sum(), 2.0 * d1.sum(), epsilon = 1e-12);
        assert_relative_eq!(
            d2.get(0, 4, 4).unwrap(),
            2.0 * d1.get(0, 4, 4).unwrap(),
            epsilon = 1e-12
        );
    }

    #[test]
    fn independent_frames_superpose() {
        // Dose from [A, B] equals dose from [A] plus dose from [B]. Use two
        // different gantry angles so the beams occupy different voxels.
        let mu = uniform_cube(0.05);
        let a = frame(0.0, 8.0, 1.0);
        let b = frame(std::f64::consts::FRAC_PI_2, 8.0, 1.0);
        let together = accumulate_delivered_dose(
            &[a.clone(), b.clone()],
            &mu,
            BeamGeometry::Parallel { standoff_mm: 500.0 },
            2.0,
            0.25,
        )
        .expect("valid attenuation volume");
        let da = accumulate_delivered_dose(
            &[a],
            &mu,
            BeamGeometry::Parallel { standoff_mm: 500.0 },
            2.0,
            0.25,
        )
        .expect("valid attenuation volume");
        let db = accumulate_delivered_dose(
            &[b],
            &mu,
            BeamGeometry::Parallel { standoff_mm: 500.0 },
            2.0,
            0.25,
        )
        .expect("valid attenuation volume");
        assert_relative_eq!(together.sum(), da.sum() + db.sum(), epsilon = 1e-12);
    }

    #[test]
    fn multi_leaf_fan_sums_offset_beamlets() {
        // Three open leaves at θ=0 place parallel +x beamlets at y = 6, 8, 10 mm
        // (perp = +y, 2 mm pitch). Each crosses the full 1.6 cm chord, so the
        // total energy is 3 × the single-beamlet value.
        let mu = uniform_cube(0.05);
        let f = DeliveryFrame {
            projection: 0,
            gantry_angle_rad: 0.0,
            couch_mm: 8.0,
            leaf_fluence: vec![1.0, 1.0, 1.0],
        };
        let dose = accumulate_delivered_dose(
            &[f],
            &mu,
            BeamGeometry::Parallel { standoff_mm: 500.0 },
            2.0,
            0.25,
        )
        .expect("valid attenuation volume");
        assert_relative_eq!(
            dose.sum(),
            3.0 * expected_axial_energy(0.05, 1.0),
            epsilon = 1e-9
        );
        // The three beamlets land in distinct leaf rows (y = 3, 4, 5).
        for j in [3usize, 4, 5] {
            assert!(
                dose.get(0, j, 4).unwrap() > 0.0,
                "row {j} should be irradiated"
            );
        }
        assert_relative_eq!(dose.get(0, 1, 4).unwrap(), 0.0, epsilon = 1e-15);
    }

    #[test]
    fn accumulation_is_generic_over_scalar_f32() {
        let grid =
            VoxelGrid::<f32>::axis_aligned([9, 9, 9], [2.0, 2.0, 2.0], Point3::new(0.0, 0.0, 0.0))
                .unwrap();
        let mu = Volume::from_shape_fn(grid, |_| 0.05_f32);
        let f = DeliveryFrame {
            projection: 0,
            gantry_angle_rad: 0.0_f32,
            couch_mm: 8.0,
            leaf_fluence: vec![2.0_f32],
        };
        let dose = accumulate_delivered_dose(
            &[f],
            &mu,
            BeamGeometry::Parallel { standoff_mm: 500.0 },
            2.0,
            0.25,
        )
        .expect("valid attenuation volume");
        let expected = 2.0_f32 * (1.0 - (-0.05_f32 * 1.6).exp());
        assert_relative_eq!(dose.sum(), expected, epsilon = 1e-5);
    }

    #[test]
    fn point_source_reduces_to_parallel_at_large_sad() {
        // As SAD → ∞ the divergent fan degenerates to parallel: the total dose of
        // a multi-leaf frame matches the parallel geometry.
        let mu = uniform_cube(0.05);
        let f = DeliveryFrame {
            projection: 0,
            gantry_angle_rad: 0.0,
            couch_mm: 8.0,
            leaf_fluence: vec![1.0, 1.0, 1.0],
        };
        let parallel = accumulate_delivered_dose(
            std::slice::from_ref(&f),
            &mu,
            BeamGeometry::Parallel { standoff_mm: 500.0 },
            2.0,
            0.25,
        )
        .expect("valid attenuation volume");
        let far = accumulate_delivered_dose(
            &[f],
            &mu,
            BeamGeometry::PointSource {
                source_axis_mm: 1.0e6,
            },
            2.0,
            0.25,
        )
        .expect("valid attenuation volume");
        assert_relative_eq!(far.sum(), parallel.sum(), max_relative = 1e-4);
    }

    #[test]
    fn point_source_fan_diverges_across_rows() {
        // A far off-axis beamlet stays in a single detector row when parallel, but
        // sweeps several rows under a divergent point-source fan — the defining fan
        // property. Fine 1 mm grid so the divergence resolves past nearest-voxel
        // quantization.
        let grid =
            VoxelGrid::axis_aligned([31, 31, 1], [1.0, 1.0, 1.0], Point3::new(0.0, 0.0, 0.0))
                .unwrap();
        let mu = Volume::from_shape_fn(grid, |_| 0.05);
        // Single lit leaf at +6 mm offset (leaf 2 of 3 at 6 mm pitch).
        let f = DeliveryFrame {
            projection: 0,
            gantry_angle_rad: 0.0,
            couch_mm: 0.0,
            leaf_fluence: vec![0.0, 0.0, 1.0],
        };
        let lit_rows = |dose: &Volume<f64>| -> usize {
            (0..31)
                .filter(|&j| (0..31).any(|i| dose.get(i, j, 0).unwrap() > 0.0))
                .count()
        };
        let par = accumulate_delivered_dose(
            std::slice::from_ref(&f),
            &mu,
            BeamGeometry::Parallel { standoff_mm: 500.0 },
            6.0,
            0.25,
        )
        .expect("valid attenuation volume");
        let pts = accumulate_delivered_dose(
            &[f],
            &mu,
            BeamGeometry::PointSource {
                source_axis_mm: 30.0,
            },
            6.0,
            0.25,
        )
        .expect("valid attenuation volume");
        assert_eq!(lit_rows(&par), 1, "parallel beamlet must stay in one row");
        assert!(
            lit_rows(&pts) >= 3,
            "divergent fan must sweep multiple rows, got {}",
            lit_rows(&pts)
        );
    }

    #[test]
    fn terma_then_scatter_produces_lateral_penumbra() {
        // End-to-end stage 1 → stage 2: a single central +x beamlet deposits terma
        // only on the y = z = centre line; the scatter kernel then spreads dose to
        // laterally-adjacent voxels that received *zero* primary terma (penumbra),
        // while the identity kernel leaves the terma unchanged.
        use helios_solver::{scatter_superposition, symmetric_deposition_kernel};
        let mu = uniform_cube(0.05);
        let terma = accumulate_delivered_dose(
            &[frame(0.0, 8.0, 1.0)],
            &mu,
            BeamGeometry::Parallel { standoff_mm: 500.0 },
            2.0,
            0.25,
        )
        .expect("valid attenuation volume");

        // Off-line voxel (mid-beam x=4, one voxel over in y) gets no primary terma.
        assert_relative_eq!(terma.get(4, 3, 4).unwrap(), 0.0, epsilon = 1e-15);

        // Identity kernel: dose == terma (differential vs the primary reference).
        let identity = scatter_superposition(&terma, &[1.0], &[1.0], &[1.0]);
        assert_relative_eq!(
            identity.get(4, 4, 4).unwrap(),
            terma.get(4, 4, 4).unwrap(),
            epsilon = 1e-15
        );

        // Spread kernel: the off-line voxel now receives scattered dose.
        let k = symmetric_deposition_kernel(0.5_f64, 0.2, 1);
        let dose = scatter_superposition(&terma, &k, &k, &k);
        assert!(
            dose.get(4, 3, 4).unwrap() > 0.0,
            "lateral neighbour must receive scattered penumbra dose"
        );
    }

    // Energy-weighted centroid along the beam axis (x, at θ=0).
    fn beam_axis_centroid_x(dose: &Volume<f64>) -> f64 {
        let [nx, ny, nz] = dose.grid().dims();
        let (mut num, mut den) = (0.0, 0.0);
        for i in 0..nx {
            for j in 0..ny {
                for k in 0..nz {
                    let d = dose.get(i, j, k).unwrap();
                    num += d * dose.grid().voxel_center(i, j, k).x;
                    den += d;
                }
            }
        }
        num / den
    }

    #[test]
    fn anisotropic_isotropic_cone_matches_the_separable_scatter_pipeline() {
        // Single frame at θ=0 (beam = +x), an ISOTROPIC cone (equal up/down
        // ranges) and a 2 mm sample step (= voxel pitch): every sample lands on a
        // node, so the per-frame oriented scatter reduces to scatter_superposition
        // of the delivered terma — the differential oracle tying the new path to
        // the verified isotropic pipeline.
        use helios_solver::scatter_superposition;
        let mu = uniform_cube(0.05);
        let geom = BeamGeometry::Parallel { standoff_mm: 500.0 };
        let frames = [frame(0.0, 8.0, 2.0)];
        let voxel_cm = 0.2;

        let cone = CollapsedCone::forward_peaked(0.4, 0.4, 0.3, voxel_cm, 1, 1, 1);
        let (beam_k, _) = forward_peaked_kernel(0.4, 0.4, voxel_cm, 1, 1);
        let lat = symmetric_deposition_kernel(0.3, voxel_cm, 1);

        let terma = accumulate_delivered_dose(&frames, &mu, geom, 2.0, 0.25)
            .expect("valid attenuation volume");
        let expected = scatter_superposition(&terma, &beam_k, &lat, &lat);
        let got = accumulate_delivered_dose_anisotropic(&frames, &mu, geom, 2.0, 0.25, &cone)
            .expect("valid attenuation volume");

        for i in 0..9 {
            for j in 0..9 {
                for k in 0..9 {
                    assert_relative_eq!(
                        got.get(i, j, k).unwrap(),
                        expected.get(i, j, k).unwrap(),
                        epsilon = 1e-10
                    );
                }
            }
        }
    }

    #[test]
    fn forward_peaked_cone_shifts_delivered_dose_downstream() {
        // Same θ=0 delivery under an isotropic vs a forward-peaked cone: the
        // forward-peaked kernel must move the beam-axis dose centroid downstream
        // (larger x) — proof the anisotropy reaches delivered dose, not just the
        // solver primitive.
        let mu = uniform_cube(0.05);
        let geom = BeamGeometry::Parallel { standoff_mm: 500.0 };
        let frames = [frame(0.0, 8.0, 2.0)];
        let voxel_cm = 0.2;
        let iso = CollapsedCone::forward_peaked(0.4, 0.4, 0.3, voxel_cm, 2, 2, 1);
        let fwd = CollapsedCone::forward_peaked(0.1, 1.0, 0.3, voxel_cm, 2, 2, 1);

        let d_iso = accumulate_delivered_dose_anisotropic(&frames, &mu, geom, 2.0, 0.25, &iso)
            .expect("valid attenuation volume");
        let d_fwd = accumulate_delivered_dose_anisotropic(&frames, &mu, geom, 2.0, 0.25, &fwd)
            .expect("valid attenuation volume");
        let (c_iso, c_fwd) = (beam_axis_centroid_x(&d_iso), beam_axis_centroid_x(&d_fwd));
        assert!(
            c_fwd > c_iso,
            "forward-peaked centroid {c_fwd} must exceed isotropic {c_iso}"
        );
    }

    #[test]
    fn single_component_poly_cone_matches_the_monoenergetic_cone() {
        // A one-component poly-energetic cone delivers exactly the same dose as
        // the monoenergetic cone with the same ranges (the differential oracle
        // tying the poly path to the verified monoenergetic one).
        let mu = uniform_cube(0.05);
        let geom = BeamGeometry::Parallel { standoff_mm: 500.0 };
        let frames = [frame(0.0, 8.0, 2.0)];
        let voxel_cm = 0.2;
        let mono = CollapsedCone::forward_peaked(0.1, 1.0, 0.3, voxel_cm, 2, 3, 1);
        let poly = CollapsedCone::poly_forward_peaked(
            &[SpectralComponent {
                range_up_cm: 0.1,
                range_down_cm: 1.0,
                weight: 5.0,
            }],
            0.3,
            voxel_cm,
            2,
            3,
            1,
        );
        let d_mono = accumulate_delivered_dose_anisotropic(&frames, &mu, geom, 2.0, 0.25, &mono)
            .expect("valid attenuation volume");
        let d_poly = accumulate_delivered_dose_anisotropic(&frames, &mu, geom, 2.0, 0.25, &poly)
            .expect("valid attenuation volume");
        for i in 0..9 {
            for j in 0..9 {
                for k in 0..9 {
                    assert_relative_eq!(
                        d_poly.get(i, j, k).unwrap(),
                        d_mono.get(i, j, k).unwrap(),
                        epsilon = 1e-13
                    );
                }
            }
        }
    }

    #[test]
    fn harder_spectrum_shifts_delivered_dose_further_downstream() {
        // Two-component spectra weighted soft vs hard: the harder-weighted beam
        // pushes the delivered-dose beam-axis centroid further downstream.
        let mu = uniform_cube(0.05);
        let geom = BeamGeometry::Parallel { standoff_mm: 500.0 };
        let frames = [frame(0.0, 8.0, 2.0)];
        let voxel_cm = 0.2;
        let soft = SpectralComponent {
            range_up_cm: 0.2,
            range_down_cm: 0.3,
            weight: 1.0,
        };
        let hard = SpectralComponent {
            range_up_cm: 0.05,
            range_down_cm: 1.5,
            weight: 1.0,
        };
        let mostly_soft = CollapsedCone::poly_forward_peaked(
            &[
                SpectralComponent {
                    weight: 9.0,
                    ..soft
                },
                hard,
            ],
            0.3,
            voxel_cm,
            2,
            2,
            1,
        );
        let mostly_hard = CollapsedCone::poly_forward_peaked(
            &[
                soft,
                SpectralComponent {
                    weight: 9.0,
                    ..hard
                },
            ],
            0.3,
            voxel_cm,
            2,
            2,
            1,
        );
        let d_soft =
            accumulate_delivered_dose_anisotropic(&frames, &mu, geom, 2.0, 0.25, &mostly_soft)
                .expect("valid attenuation volume");
        let d_hard =
            accumulate_delivered_dose_anisotropic(&frames, &mu, geom, 2.0, 0.25, &mostly_hard)
                .expect("valid attenuation volume");
        assert!(
            beam_axis_centroid_x(&d_hard) > beam_axis_centroid_x(&d_soft),
            "harder spectrum centroid {} must exceed softer {}",
            beam_axis_centroid_x(&d_hard),
            beam_axis_centroid_x(&d_soft)
        );
    }

    #[test]
    fn anisotropic_dose_is_linear_in_fluence() {
        let mu = uniform_cube(0.05);
        let geom = BeamGeometry::Parallel { standoff_mm: 500.0 };
        let voxel_cm = 0.2;
        let cone = CollapsedCone::forward_peaked(0.1, 1.0, 0.3, voxel_cm, 2, 2, 1);
        let d1 = accumulate_delivered_dose_anisotropic(
            &[frame(0.0, 8.0, 1.0)],
            &mu,
            geom,
            2.0,
            0.25,
            &cone,
        )
        .expect("valid attenuation volume");
        let d2 = accumulate_delivered_dose_anisotropic(
            &[frame(0.0, 8.0, 2.0)],
            &mu,
            geom,
            2.0,
            0.25,
            &cone,
        )
        .expect("valid attenuation volume");
        assert_relative_eq!(d2.sum(), 2.0 * d1.sum(), epsilon = 1e-12);
    }
}
