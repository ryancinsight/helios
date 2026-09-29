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
//! A helical `TomoTherapy` fan: at gantry angle `θ` the beam travels along the
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
use aequitas::systems::si::{
    quantities::{AbsorbedDose, Angle, Length},
    units::{Millimeter, Radian},
};
use eunomia::UnitScalar;
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
    /// gantry direction, offset laterally; `standoff` places the origin behind
    /// isocentre. Cheap and exact for a narrow field.
    Parallel {
        /// Distance the beamlet origin stands off behind isocentre (mm).
        standoff: Length<T>,
    },
    /// Divergent fan from a point source at `source_axis` from isocentre (SAD):
    /// each beamlet runs from the focal spot through its isocentre-plane offset
    /// point, so beamlets diverge with depth — the true `TomoTherapy` fan geometry.
    /// Reduces to [`Parallel`](Self::Parallel) as `source_axis → ∞`.
    PointSource {
        /// Source-to-axis distance / SAD (mm).
        source_axis: Length<T>,
    },
}

/// Accumulate the delivered dose from a helical-delivery `frames` sequence into a
/// dose [`Volume`] over the same grid as the attenuation volume `mu`.
///
/// `geometry` selects the beam model (parallel vs divergent point-source fan);
/// `leaf_width` is the inter-leaf lateral pitch; `step` is the ray-march
/// sampling step. Dose is linear in the per-leaf fluence, so scaling all fluence
/// scales the dose and independent frames/leaves superpose (the test oracles).
///
/// # Errors
///
/// Returns [`TransportError`] when a sampled attenuation coefficient violates
/// Hyperion's transport contract. Validation is transactional per beamlet.
pub fn accumulate_delivered_dose<T: GeometryScalar + UnitScalar>(
    frames: &[DeliveryFrame<T>],
    mu: &Volume<T>,
    geometry: BeamGeometry<T>,
    leaf_width: Length<T>,
    step: Length<T>,
) -> Result<Volume<T>, TransportError<T>> {
    let grid = *mu.grid();
    let mut dose = Volume::zeros(grid);
    for frame in frames {
        deposit_frame_terma(&mut dose, frame, mu, geometry, leaf_width, step)?;
    }
    Ok(dose)
}

/// Visit each positive-fluence leaf beamlet for one delivery `frame`, sharing the
/// per-frame gantry basis and beamlet geometry construction between dose and
/// portal kernels. The callback receives the leaf index, leaf fluence in base
/// units, and the optional constructed beamlet (degenerate rays map to `None`).
#[inline]
pub(crate) fn for_each_positive_leaf_beamlet<T, E, F>(
    frame: &DeliveryFrame<T>,
    grid: &helios_domain::VoxelGrid<T>,
    geometry: BeamGeometry<T>,
    leaf_width: Length<T>,
    mut visit: F,
) -> Result<(), E>
where
    T: GeometryScalar + UnitScalar,
    F: FnMut(usize, T, Option<Beamlet<T>>) -> Result<(), E>,
{
    let zero = <T as NumericElement>::ZERO;
    let (centre, dir, perp) = gantry_basis(grid, frame.gantry_angle_rad);
    for (leaf, fluence) in frame.leaf_fluence.iter().enumerate() {
        let weight = *fluence.as_base();
        if weight <= zero {
            continue;
        }
        let beamlet = beamlet_ray(centre, dir, perp, frame, leaf, leaf_width, geometry);
        visit(leaf, weight, beamlet)?;
    }
    Ok(())
}

/// Deposit one `frame`'s per-leaf beamlet terma into `dose`, returning the beam's
/// forward unit direction (the gantry central axis) — the SSOT deposition loop
/// shared by the isotropic accumulation and the per-frame anisotropic path, so
/// the beamlet geometry lives in exactly one place.
fn deposit_frame_terma<T: GeometryScalar + UnitScalar>(
    dose: &mut Volume<T>,
    frame: &DeliveryFrame<T>,
    mu: &Volume<T>,
    geometry: BeamGeometry<T>,
    leaf_width: Length<T>,
    step: Length<T>,
) -> Result<Vector3<T>, TransportError<T>> {
    let step_mm = step.in_unit::<Millimeter>();
    let (_, dir, _) = gantry_basis(mu.grid(), frame.gantry_angle_rad);
    for_each_positive_leaf_beamlet(
        frame,
        mu.grid(),
        geometry,
        leaf_width,
        |_, weight, beamlet| {
            let Some(beamlet) = beamlet else {
                return Ok(());
            };
            let _deposited: AbsorbedDose<T> = match beamlet.falloff {
                Some((focal, sad)) => {
                    deposit_ray_terma_diverging(dose, mu, &beamlet.ray, weight, step_mm, focal, sad)
                }
                None => deposit_ray_terma(dose, mu, &beamlet.ray, weight, step_mm),
            }?;
            Ok(())
        },
    )?;
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
    sample_step: Length<T>,
}

impl<T: GeometryScalar> CollapsedCone<T> {
    /// Build a forward-peaked cone: exponential deposition with distinct upstream
    /// (`range_up`) and downstream (`range_down`) ranges along the beam and a
    /// symmetric `lateral_range` across it. `voxel_spacing` is the sampling pitch;
    /// the radii bound each kernel's tap support. Equal up/down ranges give an
    /// isotropic (direction-independent) cone.
    #[must_use]
    pub fn forward_peaked(
        range_up: Length<T>,
        range_down: Length<T>,
        lateral_range: Length<T>,
        voxel_spacing: Length<T>,
        radius_up: usize,
        radius_down: usize,
        lateral_radius: usize,
    ) -> Self {
        let (beam_kernel, beam_center) =
            forward_peaked_kernel(range_up, range_down, voxel_spacing, radius_up, radius_down);
        Self::from_beam_kernel(
            beam_kernel,
            beam_center,
            lateral_range,
            voxel_spacing,
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
        lateral_range: Length<T>,
        voxel_spacing: Length<T>,
        radius_up: usize,
        radius_down: usize,
        lateral_radius: usize,
    ) -> Self {
        let (beam_kernel, beam_center) =
            poly_forward_peaked_kernel(components, voxel_spacing, radius_up, radius_down);
        Self::from_beam_kernel(
            beam_kernel,
            beam_center,
            lateral_range,
            voxel_spacing,
            lateral_radius,
        )
    }

    /// Shared assembly: pair a prepared beam kernel with the symmetric lateral
    /// kernel and the trilinear sample step (`voxel_cm × 10` mm). SSOT for the two
    /// public constructors above.
    fn from_beam_kernel(
        beam_kernel: Vec<T>,
        beam_center: usize,
        lateral_range: Length<T>,
        voxel_spacing: Length<T>,
        lateral_radius: usize,
    ) -> Self {
        let lateral = symmetric_deposition_kernel(lateral_range, voxel_spacing, lateral_radius);
        let sample_step = voxel_spacing;
        Self {
            beam_kernel,
            beam_center,
            lateral,
            sample_step,
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
pub fn accumulate_delivered_dose_anisotropic<T: GeometryScalar + UnitScalar>(
    frames: &[DeliveryFrame<T>],
    mu: &Volume<T>,
    geometry: BeamGeometry<T>,
    leaf_width: Length<T>,
    step: Length<T>,
    cone: &CollapsedCone<T>,
) -> Result<Volume<T>, TransportError<T>> {
    let grid = *mu.grid();
    let mut acc = vec![<T as NumericElement>::ZERO; grid.num_voxels()];
    for frame in frames {
        let mut frame_terma = Volume::zeros(grid);
        let dir = deposit_frame_terma(&mut frame_terma, frame, mu, geometry, leaf_width, step)?;
        let frame_dose = oriented_forward_scatter(
            &frame_terma,
            dir,
            &cone.beam_kernel,
            cone.beam_center,
            &cone.lateral,
            cone.sample_step,
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
pub(crate) fn beamlet_ray<T: GeometryScalar + UnitScalar>(
    centre: Point3<T>,
    dir: Vector3<T>,
    perp: Vector3<T>,
    frame: &DeliveryFrame<T>,
    leaf: usize,
    leaf_width: Length<T>,
    geometry: BeamGeometry<T>,
) -> Option<Beamlet<T>> {
    let zero = <T as NumericElement>::ZERO;
    let leaf_width_mm = leaf_width.in_unit::<Millimeter>();
    let couch_mm = frame.couch.in_unit::<Millimeter>();
    let leaves = frame.leaf_fluence.len();
    let centre_leaf = <T as GeometryScalar>::from_f64((leaves as f64 - 1.0) * 0.5);
    let offset = (<T as GeometryScalar>::from_f64(leaf as f64) - centre_leaf) * leaf_width_mm;
    // (origin, direction, inverse-square falloff) per the selected geometry. Both
    // branches lie in the couch z-slice (dir.z = perp.z = 0); the ray constructor
    // normalizes the direction vector.
    let (origin, direction, falloff) = match geometry {
        BeamGeometry::Parallel { standoff } => {
            let standoff_mm = standoff.in_unit::<Millimeter>();
            (
                Point3::new(
                    centre.x + perp.x * offset - dir.x * standoff_mm,
                    centre.y + perp.y * offset - dir.y * standoff_mm,
                    couch_mm,
                ),
                dir,
                None,
            )
        }
        BeamGeometry::PointSource { source_axis } => {
            let source_axis_mm = source_axis.in_unit::<Millimeter>();
            // Focal spot behind isocentre; ray aims through the leaf's isocentre-
            // plane point `centre + perp·offset`, so direction = (perp·offset +
            // dir·SAD). For offset 0 this is the central axis; off-axis leaves fan
            // out with depth. Fluence falls off inverse-square from the focal spot.
            let focal = Point3::new(
                centre.x - dir.x * source_axis_mm,
                centre.y - dir.y * source_axis_mm,
                couch_mm,
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
pub(crate) fn gantry_basis<T: GeometryScalar + UnitScalar>(
    grid: &helios_domain::VoxelGrid<T>,
    gantry_angle_rad: Angle<T>,
) -> (Point3<T>, Vector3<T>, Vector3<T>) {
    let zero = <T as NumericElement>::ZERO;
    let [nx, ny, nz] = grid.dims();
    let centre = grid.voxel_center((nx - 1) / 2, (ny - 1) / 2, (nz - 1) / 2);
    let angle = gantry_angle_rad.in_unit::<Radian>();
    let (cos, sin) = (angle.cos(), angle.sin());
    (
        centre,
        Vector3::new(cos, sin, zero),
        Vector3::new(-sin, cos, zero),
    )
}

#[cfg(test)]
mod tests;
