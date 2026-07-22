//! Helical MVCT acquisition: rotate the beam per projection and forward-project.

use aequitas::systems::si::quantities::Dimensionless;
use helios_domain::{HelicalDelivery, Volume};
use helios_math::{GeometryScalar, NumericElement, Point3, Ray, Vector3};
use helios_solver::forward_project_ray;
use hyperion::{quantity::OpticalDepth, TransportError};
use moirai_parallel::Adaptive;

/// One projection of a helical acquisition: the delivery state (gantry angle,
/// couch position) and the resulting central-ray measurement.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HelicalProjection<T: GeometryScalar> {
    /// Projection index.
    pub projection: usize,
    /// Gantry angle at this projection (rad).
    pub gantry_angle_rad: T,
    /// Couch position at this projection (mm).
    pub couch_mm: T,
    /// Central-ray optical depth `∫μ dl` (dimensionless).
    pub optical_depth: T,
    /// Central-ray transmitted fraction `exp(−∫μ dl)`.
    pub transmission: T,
}

/// Simulate a helical acquisition of `num_projections` over the attenuation
/// volume `mu`.
///
/// For each projection the beam rotates in the **axial (x–y) plane** by the
/// [`HelicalDelivery`] gantry angle while the couch advances the imaged slice
/// along **z** — together tracing a helix. The central ray, aimed through the
/// grid's axial centre at the projection's couch `z`, starts `source_distance_mm`
/// behind isocentre and is forward-projected with sampling step `step_mm`. Rays
/// that miss the grid (e.g. couch beyond the volume) record zero optical depth
/// (full transmission).
///
/// The independent per-projection forward projections are dispatched through
/// moirai's [`Adaptive`] execution policy (sequential below its threshold, parallel
/// above) — the mandated time-dependent-orchestration seam. The collect is
/// index-ordered, so the result is identical to a sequential run regardless of
/// thread scheduling (each projection is an independent read of `mu`; no reduction).
/// `T: Send + Sync` (satisfied by every real scalar) is required for the dispatch.
///
/// # Errors
///
/// Returns [`TransportError`] if a projected optical depth is negative or
/// non-finite.
pub fn simulate_helical_sinogram<T: GeometryScalar + Send + Sync>(
    delivery: &HelicalDelivery<T>,
    mu: &Volume<T>,
    num_projections: usize,
    source_distance_mm: T,
    step_mm: T,
) -> Result<Vec<HelicalProjection<T>>, TransportError<T>> {
    let zero = <T as NumericElement>::ZERO;
    let grid = *mu.grid();
    let [nx, ny, nz] = grid.dims();
    // Axial centre of the grid (used for the beam's x–y aim point).
    let centre = grid.voxel_center((nx - 1) / 2, (ny - 1) / 2, (nz - 1) / 2);

    let projections =
        moirai_parallel::map_collect_index_with::<Adaptive, _, _>(num_projections, |projection| {
            let gantry_angle_rad = delivery.gantry_angle_rad(projection);
            let couch_mm = delivery.couch_position_mm(projection);

            // Beam direction rotates in the axial plane; z fixed at the couch slice.
            let direction = Vector3::new(gantry_angle_rad.cos(), gantry_angle_rad.sin(), zero);
            // Aim point: axial centre at the couch z; source sits behind isocentre.
            let origin = Point3::new(
                centre.x - direction.x * source_distance_mm,
                centre.y - direction.y * source_distance_mm,
                couch_mm - direction.z * source_distance_mm,
            );

            let optical_depth = Ray::try_new(origin, direction)
                .ok()
                .and_then(|ray| forward_project_ray(mu, &ray, step_mm))
                .unwrap_or(zero);
            let transmission = OpticalDepth::new(Dimensionless::from_base(optical_depth))?
                .transmission()
                .into_quantity()
                .into_base();

            Ok(HelicalProjection {
                projection,
                gantry_angle_rad,
                couch_mm,
                optical_depth,
                transmission,
            })
        });
    projections.into_iter().collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use eunomia::assert_relative_eq;
    use helios_domain::VoxelGrid;
    use helios_math::Point3;

    // Uniform-μ cube: 9³ voxels, 2 mm spacing → node extent 16 mm = 1.6 cm/axis.
    fn uniform_cube(mu_val: f64) -> Volume<f64> {
        let grid = VoxelGrid::axis_aligned([9, 9, 9], [2.0, 2.0, 2.0], Point3::new(0.0, 0.0, 0.0))
            .expect("grid");
        Volume::from_shape_fn(grid, move |_| mu_val)
    }

    // 4 projections/rotation so projection 1 is a clean 90° turn; couch centred.
    fn delivery() -> HelicalDelivery<f64> {
        HelicalDelivery::new(4, 25.0, 0.2, 10.0, 0.0, 8.0).expect("delivery")
    }

    #[test]
    fn sinogram_has_one_entry_per_projection() {
        let sino = simulate_helical_sinogram(&delivery(), &uniform_cube(0.05), 12, 500.0, 0.5)
            .expect("valid attenuation volume");
        assert_eq!(sino.len(), 12);
        assert!(sino.iter().enumerate().all(|(i, p)| p.projection == i));
    }

    #[test]
    fn axial_central_ray_measures_mu_times_chord() {
        // Projection 0: θ=0 → +x ray through the cube centre. Chord = 16 mm =
        // 1.6 cm, μ = 0.05 → τ = 0.08; transmission = exp(-0.08).
        let sino = simulate_helical_sinogram(&delivery(), &uniform_cube(0.05), 4, 500.0, 0.25)
            .expect("valid attenuation volume");
        assert_relative_eq!(sino[0].optical_depth, 0.05 * 1.6, epsilon = 1e-9);
        assert_relative_eq!(
            sino[0].transmission,
            (-0.05 * 1.6_f64).exp(),
            epsilon = 1e-9
        );
    }

    #[test]
    fn rotational_symmetry_of_a_uniform_cube() {
        // For a uniform cube the central-ray line integral is the same at 0° and
        // 90° (equal chords), independent of the couch advance.
        let sino = simulate_helical_sinogram(&delivery(), &uniform_cube(0.05), 4, 500.0, 0.25)
            .expect("valid attenuation volume");
        assert_relative_eq!(
            sino[0].optical_depth,
            sino[1].optical_depth,
            max_relative = 1e-6
        );
    }

    #[test]
    fn couch_advances_monotonically_across_projections() {
        let sino = simulate_helical_sinogram(&delivery(), &uniform_cube(0.05), 20, 500.0, 1.0)
            .expect("valid attenuation volume");
        for pair in sino.windows(2) {
            assert!(pair[1].couch_mm > pair[0].couch_mm, "couch must advance");
        }
    }

    #[test]
    fn empty_region_transmits_fully() {
        // Zero-μ volume → no attenuation → τ=0, transmission=1 everywhere.
        let sino = simulate_helical_sinogram(&delivery(), &uniform_cube(0.0), 6, 500.0, 0.5)
            .expect("valid attenuation volume");
        for p in &sino {
            assert_relative_eq!(p.optical_depth, 0.0, epsilon = 1e-12);
            assert_relative_eq!(p.transmission, 1.0, epsilon = 1e-12);
        }
    }

    #[test]
    fn simulation_is_generic_over_scalar_f32() {
        let grid =
            VoxelGrid::<f32>::axis_aligned([9, 9, 9], [2.0, 2.0, 2.0], Point3::new(0.0, 0.0, 0.0))
                .unwrap();
        let mu = Volume::from_shape_fn(grid, |_| 0.05_f32);
        let del = HelicalDelivery::<f32>::new(4, 25.0, 0.2, 10.0, 0.0, 8.0).unwrap();
        let sino =
            simulate_helical_sinogram(&del, &mu, 4, 500.0, 0.25).expect("valid attenuation volume");
        assert_relative_eq!(sino[0].optical_depth, 0.05_f32 * 1.6, epsilon = 1e-4);
    }

    #[test]
    fn moirai_dispatch_is_deterministic_and_order_preserving() {
        // 256 projections exceed moirai's Adaptive parallel threshold, so this
        // exercises the parallel path. The index-ordered collect makes the result
        // identical run-to-run (no data race; each projection is an independent
        // read) — the differential guarantee vs a sequential run.
        let a = simulate_helical_sinogram(&delivery(), &uniform_cube(0.05), 256, 500.0, 0.5)
            .expect("valid attenuation volume");
        let b = simulate_helical_sinogram(&delivery(), &uniform_cube(0.05), 256, 500.0, 0.5)
            .expect("valid attenuation volume");
        assert_eq!(a, b);
        assert!(a.iter().enumerate().all(|(i, p)| p.projection == i));
    }

    #[test]
    fn negative_projected_optical_depth_is_rejected() {
        let error = simulate_helical_sinogram(&delivery(), &uniform_cube(-0.05), 1, 500.0, 0.25)
            .expect_err("negative optical depth must fail");
        assert!(matches!(
            error,
            TransportError::InvalidValue {
                field: hyperion::ValueKind::OpticalDepth,
                ..
            }
        ));
    }
}
