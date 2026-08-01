//! Helical MVCT acquisition: rotate the beam per projection and forward-project.

use aequitas::systems::si::{
    quantities::{Angle, Dimensionless, Length},
    units::{Millimeter, Radian},
};
use eunomia::UnitScalar;
use helios_domain::{HelicalDelivery, Volume};
use helios_math::{GeometryScalar, NumericElement, Point3, Ray, Vector3};
use helios_solver::forward_project_ray;
use hyperion::{TransportError, quantity::OpticalDepth};
use moirai_parallel::Adaptive;

#[cfg(test)]
use aequitas::systems::si::quantities::Time;

/// One projection of a helical acquisition: the delivery state (gantry angle,
/// couch position) and the resulting central-ray measurement.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HelicalProjection<T: GeometryScalar> {
    /// Projection index.
    pub projection: usize,
    /// Gantry angle at this projection (rad).
    pub gantry_angle_rad: Angle<T>,
    /// Couch position at this projection (mm).
    pub couch_mm: Length<T>,
    /// Central-ray optical depth `∫μ dl`.
    pub optical_depth: Dimensionless<T>,
    /// Central-ray transmitted fraction `exp(−∫μ dl)`.
    pub transmission: Dimensionless<T>,
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
pub fn simulate_helical_sinogram<T: GeometryScalar + UnitScalar + Send + Sync>(
    delivery: &HelicalDelivery<T>,
    mu: &Volume<T>,
    num_projections: usize,
    source_distance_mm: Length<T>,
    step_mm: Length<T>,
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
            let angle = gantry_angle_rad.in_unit::<Radian>();
            let couch = couch_mm.in_unit::<Millimeter>();

            // Beam direction rotates in the axial plane; z fixed at the couch slice.
            let direction = Vector3::new(angle.cos(), angle.sin(), zero);
            // Aim point: axial centre at the couch z; source sits behind isocentre.
            let origin = Point3::new(
                centre.x - direction.x * source_distance_mm.in_unit::<Millimeter>(),
                centre.y - direction.y * source_distance_mm.in_unit::<Millimeter>(),
                couch - direction.z * source_distance_mm.in_unit::<Millimeter>(),
            );

            let optical_depth = Ray::try_new(origin, direction)
                .ok()
                .and_then(|ray| forward_project_ray(mu, &ray, step_mm.in_unit::<Millimeter>()))
                .unwrap_or(zero);
            let optical_depth = Dimensionless::from_base(optical_depth);
            let transmission = OpticalDepth::new(optical_depth)?
                .transmission()
                .into_quantity();

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
    use helios_math::ShippedScalar;

    // Uniform-μ cube: 9³ voxels, 2 mm spacing → node extent 16 mm = 1.6 cm/axis.
    fn uniform_cube(mu_val: f64) -> Volume<f64> {
        let grid = VoxelGrid::axis_aligned([9, 9, 9], [2.0, 2.0, 2.0], Point3::new(0.0, 0.0, 0.0))
            .expect("grid");
        Volume::from_shape_fn(grid, move |_| mu_val)
    }

    // 4 projections/rotation so projection 1 is a clean 90° turn; couch centred.
    fn delivery() -> HelicalDelivery<f64> {
        HelicalDelivery::new(
            4,
            Length::from_unit::<Millimeter>(25.0),
            Dimensionless::from_base(0.2),
            Time::from_unit::<aequitas::systems::si::units::Second>(10.0),
            Angle::from_unit::<Radian>(0.0),
            Length::from_unit::<Millimeter>(8.0),
        )
        .expect("delivery")
    }

    #[test]
    fn sinogram_has_one_entry_per_projection() {
        let sino = simulate_helical_sinogram(
            &delivery(),
            &uniform_cube(0.05),
            12,
            Length::from_unit::<Millimeter>(500.0),
            Length::from_unit::<Millimeter>(0.5),
        )
        .expect("valid attenuation volume");
        assert_eq!(sino.len(), 12);
        assert!(sino.iter().enumerate().all(|(i, p)| p.projection == i));
    }

    #[test]
    fn axial_central_ray_measures_mu_times_chord() {
        // Projection 0: θ=0 → +x ray through the cube centre. Chord = 16 mm =
        // 1.6 cm, μ = 0.05 → τ = 0.08; transmission = exp(-0.08).
        let sino = simulate_helical_sinogram(
            &delivery(),
            &uniform_cube(0.05),
            4,
            Length::from_unit::<Millimeter>(500.0),
            Length::from_unit::<Millimeter>(0.25),
        )
        .expect("valid attenuation volume");
        assert_relative_eq!(
            sino[0].optical_depth.into_base(),
            0.05 * 1.6,
            epsilon = 1e-9
        );
        assert_relative_eq!(
            sino[0].transmission.into_base(),
            (-0.05 * 1.6_f64).exp(),
            epsilon = 1e-9
        );
    }

    #[test]
    fn rotational_symmetry_of_a_uniform_cube() {
        // For a uniform cube the central-ray line integral is the same at 0° and
        // 90° (equal chords), independent of the couch advance.
        let sino = simulate_helical_sinogram(
            &delivery(),
            &uniform_cube(0.05),
            4,
            Length::from_unit::<Millimeter>(500.0),
            Length::from_unit::<Millimeter>(0.25),
        )
        .expect("valid attenuation volume");
        assert_relative_eq!(
            sino[0].optical_depth.into_base(),
            sino[1].optical_depth.into_base(),
            max_relative = 1e-6
        );
    }

    #[test]
    fn couch_advances_monotonically_across_projections() {
        let sino = simulate_helical_sinogram(
            &delivery(),
            &uniform_cube(0.05),
            20,
            Length::from_unit::<Millimeter>(500.0),
            Length::from_unit::<Millimeter>(1.0),
        )
        .expect("valid attenuation volume");
        for pair in sino.windows(2) {
            assert!(
                pair[1].couch_mm.in_unit::<Millimeter>() > pair[0].couch_mm.in_unit::<Millimeter>(),
                "couch must advance"
            );
        }
    }

    #[test]
    fn empty_region_transmits_fully() {
        // Zero-μ volume → no attenuation → τ=0, transmission=1 everywhere.
        let sino = simulate_helical_sinogram(
            &delivery(),
            &uniform_cube(0.0),
            6,
            Length::from_unit::<Millimeter>(500.0),
            Length::from_unit::<Millimeter>(0.5),
        )
        .expect("valid attenuation volume");
        for p in &sino {
            assert_relative_eq!(p.optical_depth.into_base(), 0.0, epsilon = 1e-12);
            assert_relative_eq!(p.transmission.into_base(), 1.0, epsilon = 1e-12);
        }
    }

    /// Sinogram optical depth is a ray-marched accumulation at a 0.25 mm step over
    /// the traversed 1.6 cm, so on the order of 64 additions contribute. Worst-case
    /// growth is `n * T::EPSILON`, and 256 ulps of `T` bounds that. The former
    /// 1e-4 absolute bound was ~1250 ulps of `f32` at this magnitude, so this
    /// tightens the assertion while making its basis explicit.
    const SINOGRAM_ACCUMULATION_ULPS: f64 = 256.0;

    /// Asserts sinogram optical depth through a uniform slab in one scalar width.
    fn helical_sinogram_accumulates_uniform_optical_depth<T>()
    where
        T: ShippedScalar + helios_math::GeometryScalar,
    {
        // GeometryScalar and FloatElement both define `from_f64`; name the
        // one that performs the literal conversion so the call is unambiguous.
        let cast = <T as helios_math::FloatElement>::from_f64;
        let zero = cast(0.0);
        let spacing = cast(2.0);
        let grid =
            VoxelGrid::<T>::axis_aligned([9, 9, 9], [spacing; 3], Point3::new(zero, zero, zero))
                .expect("valid axis-aligned grid");

        let mu = cast(0.05);
        let attenuation = Volume::from_shape_fn(grid, |_| mu);
        let delivery = HelicalDelivery::<T>::new(
            4,
            Length::from_unit::<Millimeter>(cast(25.0)),
            Dimensionless::from_base(cast(0.2)),
            Time::from_unit::<aequitas::systems::si::units::Second>(cast(10.0)),
            Angle::from_unit::<Radian>(zero),
            Length::from_unit::<Millimeter>(cast(8.0)),
        )
        .expect("valid helical geometry");

        let sinogram = simulate_helical_sinogram(
            &delivery,
            &attenuation,
            4,
            Length::from_unit::<Millimeter>(cast(500.0)),
            Length::from_unit::<Millimeter>(cast(0.25)),
        )
        .expect("valid attenuation volume");

        // Uniform mu over the 1.6 cm traversed path.
        assert_relative_eq!(
            sinogram[0].optical_depth.into_base(),
            mu * cast(1.6),
            max_relative = T::EPSILON * cast(SINOGRAM_ACCUMULATION_ULPS)
        );
    }

    #[test]
    fn helical_sinogram_accumulates_uniform_optical_depth_in_single_precision() {
        helical_sinogram_accumulates_uniform_optical_depth::<f32>();
    }

    #[test]
    fn helical_sinogram_accumulates_uniform_optical_depth_in_double_precision() {
        helical_sinogram_accumulates_uniform_optical_depth::<f64>();
    }

    #[test]
    fn moirai_dispatch_is_deterministic_and_order_preserving() {
        // 256 projections exceed moirai's Adaptive parallel threshold, so this
        // exercises the parallel path. The index-ordered collect makes the result
        // identical run-to-run (no data race; each projection is an independent
        // read) — the differential guarantee vs a sequential run.
        let a = simulate_helical_sinogram(
            &delivery(),
            &uniform_cube(0.05),
            256,
            Length::from_unit::<Millimeter>(500.0),
            Length::from_unit::<Millimeter>(0.5),
        )
        .expect("valid attenuation volume");
        let b = simulate_helical_sinogram(
            &delivery(),
            &uniform_cube(0.05),
            256,
            Length::from_unit::<Millimeter>(500.0),
            Length::from_unit::<Millimeter>(0.5),
        )
        .expect("valid attenuation volume");
        assert_eq!(a, b);
        assert!(a.iter().enumerate().all(|(i, p)| p.projection == i));
    }

    #[test]
    fn negative_projected_optical_depth_is_rejected() {
        let error = simulate_helical_sinogram(
            &delivery(),
            &uniform_cube(-0.05),
            1,
            Length::from_unit::<Millimeter>(500.0),
            Length::from_unit::<Millimeter>(0.25),
        )
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
