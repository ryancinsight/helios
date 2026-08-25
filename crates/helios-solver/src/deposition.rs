//! Primary-fluence energy deposition (TERMA) along a beam ray.
//!
//! Companion to the [`forward_project_ray`](crate::forward_project_ray) line
//! integral: where the projector reduces `∫ μ dl` to a single optical depth, this
//! kernel *deposits* the energy the primary beam loses as it attenuates, voxel by
//! voxel, producing the terma (total energy released per unit mass) that a
//! collapsed-cone/convolution dose engine spreads with a scatter kernel.

use crate::projector::ray_grid_interval;
use aequitas::systems::si::{
    quantities::{AbsorbedDose, Length, ReciprocalLength},
    units::{Centimeter, PerCentimeter},
};
use eunomia::UnitScalar;
use helios_core::constants::MM_PER_CM;
use helios_domain::Volume;
use helios_math::{GeometryScalar, NumericElement, Point3, Ray};
use hyperion::{
    coefficient::{InteractionCoefficient, LinearAttenuation},
    quantity::{OpticalDepth, PathLength},
    TransportError,
};

/// Nearest voxel index along one axis for a continuous index `coord`, clamped to
/// `[0, n−1]`. Segment midpoints lie inside the node-centre AABB, so the clamp
/// only guards floating-point boundary rounding.
fn nearest<T: GeometryScalar>(coord: T, n: usize) -> usize {
    let half = <T as GeometryScalar>::from_f64(0.5);
    let r = (coord + half).floor().to_f64();
    if r <= 0.0 {
        0
    } else {
        (r as usize).min(n - 1)
    }
}

/// Deposit primary-beam energy along `ray` into `dose`, returning the total
/// energy removed from the primary beam.
///
/// # Model
/// The primary energy fluence attenuates as `Ψ(s) = weight · e^{−τ(s)}`. The
/// energy lost in a path segment `[s_i, s_{i+1}]` is
/// `weight · (e^{−τ_i} − e^{−τ_{i+1}})`; it is scattered into the voxel nearest
/// the segment midpoint. Because the per-segment losses telescope, the returned
/// total is **exactly** `weight · (1 − e^{−τ_total})` — independent of `step_mm`
/// — and equals the sum of the deposited voxel values (energy conservation).
/// This is the terma along the ray; lateral scatter is a later increment.
///
/// # Units
/// `dose` and `mu` must share the same grid. `mu` is in cm⁻¹ and the grid / `ray`
/// in mm, so segment lengths are converted mm→cm (matching the projector). A ray
/// that misses the grid deposits nothing and returns zero. The returned total is
/// an Aequitas [`AbsorbedDose`] quantity; the voxel field remains the established
/// scalar `Volume` storage boundary.
///
/// # Errors
///
/// Returns [`TransportError`] when a sampled attenuation coefficient is negative
/// or non-finite, or when an optical-depth product or partial sum is non-finite.
/// Validation completes before `dose` is mutated, so an error leaves the output
/// unchanged.
pub fn deposit_ray_terma<T: GeometryScalar + UnitScalar>(
    dose: &mut Volume<T>,
    mu: &Volume<T>,
    ray: &Ray<T>,
    weight: T,
    step_mm: T,
) -> Result<AbsorbedDose<T>, TransportError<T>> {
    deposit_terma_impl(dose, mu, ray, weight, step_mm, None)
}

/// Divergent-fan variant of [`deposit_ray_terma`]: the per-segment terma is
/// additionally scaled by the inverse-square fluence falloff `(sad_mm / r)²` from
/// the point source at `focal`, with `r` the focal-to-segment distance.
///
/// The factor is 1 at isocentre (`r = sad_mm`), `> 1` nearer the source, and `< 1`
/// beyond — the geometric divergence of a real fan beam. It reduces to
/// [`deposit_ray_terma`] as `sad_mm → ∞`. The returned total is no longer the
/// closed-form `weight·(1 − e^{−τ})` (the falloff breaks the telescoping) but still
/// equals the summed deposited voxel dose.
///
/// # Errors
///
/// Returns the same typed transport failures as [`deposit_ray_terma`]. The
/// returned total is an Aequitas [`AbsorbedDose`] quantity; voxel storage keeps
/// the established scalar `Volume` boundary.
pub fn deposit_ray_terma_diverging<T: GeometryScalar + UnitScalar>(
    dose: &mut Volume<T>,
    mu: &Volume<T>,
    ray: &Ray<T>,
    weight: T,
    step_mm: T,
    focal: Point3<T>,
    sad_mm: T,
) -> Result<AbsorbedDose<T>, TransportError<T>> {
    deposit_terma_impl(dose, mu, ray, weight, step_mm, Some((focal, sad_mm)))
}

/// Shared ray-march for [`deposit_ray_terma`] and [`deposit_ray_terma_diverging`];
/// `falloff = Some((focal, sad))` applies the inverse-square divergence factor.
fn deposit_terma_impl<T: GeometryScalar + UnitScalar>(
    dose: &mut Volume<T>,
    mu: &Volume<T>,
    ray: &Ray<T>,
    weight: T,
    step_mm: T,
    falloff: Option<(Point3<T>, T)>,
) -> Result<AbsorbedDose<T>, TransportError<T>> {
    let grid = *mu.grid();
    debug_assert_eq!(
        grid.dims(),
        dose.grid().dims(),
        "dose and mu must share the same grid"
    );
    let Some((t_enter, t_exit)) = ray_grid_interval(&grid, ray) else {
        return Ok(AbsorbedDose::from_base(T::ZERO));
    };
    let length = t_exit - t_enter;
    if length <= T::ZERO {
        return Ok(AbsorbedDose::from_base(T::ZERO));
    }

    // Substeps so the step divides the traversed length exactly (>= 1).
    let n = ((length * step_mm.recip()).ceil().to_f64() as usize).max(1);
    let actual_step = length * <T as GeometryScalar>::from_f64(n as f64).recip();
    let step_cm = actual_step * <T as GeometryScalar>::from_f64(MM_PER_CM).recip();
    let path = PathLength::new(Length::from_unit::<Centimeter>(step_cm))?;
    let half = <T as GeometryScalar>::from_f64(0.5);
    let [nx, ny, nz] = grid.dims();
    let sample = |i: usize| {
        let t_mid = t_enter + (<T as GeometryScalar>::from_f64(i as f64) + half) * actual_step;
        let world_pt: Point3<T> = ray.point_at(t_mid);
        let index = grid.world_to_index(world_pt);
        let mu_sample = mu.sample_trilinear(index).unwrap_or(T::ZERO);
        (world_pt, index, mu_sample)
    };

    // Validate the complete ray before mutating the output. Non-negative
    // segment depths make every partial sum bounded by this checked total.
    let _validated_total = (0..n).try_fold(OpticalDepth::zero(), |total, i| {
        let (_, _, mu_sample) = sample(i);
        let coefficient =
            InteractionCoefficient::<T, LinearAttenuation>::new(ReciprocalLength::from_unit::<
                PerCentimeter,
            >(mu_sample))?;
        total.checked_add(coefficient.optical_depth(path)?)
    })?;

    let mut optical_depth = OpticalDepth::zero();
    let mut trans_before = <T as NumericElement>::ONE; // e^{−τ} at τ = 0.
    let mut total = T::ZERO;
    for i in 0..n {
        let (world_pt, index, mu_sample) = sample(i);
        let coefficient =
            InteractionCoefficient::<T, LinearAttenuation>::new(ReciprocalLength::from_unit::<
                PerCentimeter,
            >(mu_sample))?;
        optical_depth = optical_depth.checked_add(coefficient.optical_depth(path)?)?;
        let trans_after = optical_depth.transmission().into_quantity().into_base();
        let mut absorbed = weight * (trans_before - trans_after);
        if let Some((focal, sad)) = falloff {
            let dx = world_pt.x - focal.x;
            let dy = world_pt.y - focal.y;
            let dz = world_pt.z - focal.z;
            let r2 = dx * dx + dy * dy + dz * dz;
            if r2 > T::ZERO {
                absorbed *= sad * sad * r2.recip();
            }
        }
        dose.add_at(
            nearest(index.x, nx),
            nearest(index.y, ny),
            nearest(index.z, nz),
            absorbed,
        );
        total += absorbed;
        trans_before = trans_after;
    }
    Ok(AbsorbedDose::from_base(total))
}

#[cfg(test)]
mod tests {
    #![expect(
        clippy::unwrap_used,
        reason = "ratchet HELIOS-UNWRAP-1: pre-existing debt"
    )]
    use super::*;
    use eunomia::assert_relative_eq;
    use helios_domain::VoxelGrid;
    use helios_math::ShippedScalar;
    use helios_math::{Point3, Vector3};

    // Uniform-μ cube: 9³ voxels, 2 mm spacing → node box [0,16] mm = 1.6 cm/axis.
    fn uniform_cube(mu_val: f64) -> Volume<f64> {
        let grid = VoxelGrid::axis_aligned([9, 9, 9], [2.0, 2.0, 2.0], Point3::new(0.0, 0.0, 0.0))
            .expect("grid");
        Volume::from_shape_fn(grid, move |_| mu_val)
    }

    // +x ray through the cube centre (y = z = 8 mm), starting outside the box.
    fn central_x_ray() -> Ray<f64> {
        Ray::try_new(Point3::new(-50.0, 8.0, 8.0), Vector3::new(1.0, 0.0, 0.0))
            .expect("unit +x ray")
    }

    fn oriented_cube(mu_val: f64) -> Volume<f64> {
        let rotation = helios_math::UnitQuaternion::try_from_rotation_columns(
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(-1.0, 0.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
            1.0e-12,
        )
        .expect("right-handed quarter-turn basis");
        let grid = VoxelGrid::oriented(
            [9, 9, 9],
            [2.0, 2.0, 2.0],
            Point3::new(10.0, 20.0, 30.0),
            rotation,
        )
        .expect("grid");
        Volume::from_shape_fn(grid, move |_| mu_val)
    }

    #[test]
    fn total_deposited_equals_primary_energy_lost() {
        // Uniform μ = 0.05 cm⁻¹, chord 1.6 cm → τ = 0.08. Energy removed from a
        // unit-weight beam = 1 − e^{−0.08}. Exact (telescoping), any step.
        let mu = uniform_cube(0.05);
        let mut dose = Volume::zeros(*mu.grid());
        let total = deposit_ray_terma(&mut dose, &mu, &central_x_ray(), 1.0, 0.5)
            .expect("valid attenuation volume")
            .into_base();
        let expected = 1.0 - (-0.05 * 1.6_f64).exp();
        assert_relative_eq!(total, expected, epsilon = 1e-12);
    }

    #[test]
    fn oriented_grid_deposition_preserves_primary_energy_loss() {
        // The local x-span is 16 mm; rotation maps it to world +y. The same
        // Beer–Lambert oracle must hold after index-space clipping.
        let mu = oriented_cube(0.05);
        let mut dose = Volume::zeros(*mu.grid());
        let ray = Ray::try_new(Point3::new(2.0, -20.0, 38.0), Vector3::new(0.0, 1.0, 0.0))
            .expect("unit +y ray");
        let total = deposit_ray_terma(&mut dose, &mu, &ray, 1.0, 0.5)
            .expect("valid attenuation volume")
            .into_base();
        let expected = 1.0 - (-0.05 * 1.6_f64).exp();
        assert_relative_eq!(total, expected, epsilon = 1e-12);
        assert_relative_eq!(dose.sum(), expected, epsilon = 1e-12);
    }

    #[test]
    fn deposited_voxels_conserve_the_returned_total() {
        // Sum of scattered voxel dose must equal the returned total exactly.
        let mu = uniform_cube(0.08);
        let mut dose = Volume::zeros(*mu.grid());
        let total = deposit_ray_terma(&mut dose, &mu, &central_x_ray(), 3.0, 0.3)
            .expect("valid attenuation volume")
            .into_base();
        assert_relative_eq!(dose.sum(), total, epsilon = 1e-12);
    }

    #[test]
    fn step_size_does_not_change_the_total() {
        // Telescoping makes the total independent of the sampling step.
        let mu = uniform_cube(0.05);
        let (mut d_coarse, mut d_fine) = (Volume::zeros(*mu.grid()), Volume::zeros(*mu.grid()));
        let coarse = deposit_ray_terma(&mut d_coarse, &mu, &central_x_ray(), 1.0, 4.0)
            .expect("valid attenuation volume")
            .into_base();
        let fine = deposit_ray_terma(&mut d_fine, &mu, &central_x_ray(), 1.0, 0.05)
            .expect("valid attenuation volume")
            .into_base();
        assert_relative_eq!(coarse, fine, epsilon = 1e-12);
    }

    #[test]
    fn energy_is_front_loaded_by_attenuation() {
        // Primary fluence decays with depth, so the entry voxel absorbs more than
        // the exit voxel along the beam.
        let mu = uniform_cube(0.3); // strong attenuation to make the gradient clear
        let mut dose = Volume::zeros(*mu.grid());
        let _ = deposit_ray_terma(&mut dose, &mu, &central_x_ray(), 1.0, 0.1)
            .expect("valid attenuation volume");
        let entry = dose.get(0, 4, 4).unwrap();
        let exit = dose.get(8, 4, 4).unwrap();
        assert!(entry > exit, "entry {entry} should exceed exit {exit}");
    }

    #[test]
    fn zero_attenuation_and_zero_weight_deposit_nothing() {
        let empty = uniform_cube(0.0);
        let mut d0 = Volume::zeros(*empty.grid());
        assert_relative_eq!(
            deposit_ray_terma(&mut d0, &empty, &central_x_ray(), 1.0, 0.5)
                .expect("valid attenuation volume")
                .into_base(),
            0.0,
            epsilon = 1e-15
        );
        assert_relative_eq!(d0.sum(), 0.0, epsilon = 1e-15);

        let mu = uniform_cube(0.05);
        let mut dw = Volume::zeros(*mu.grid());
        assert_relative_eq!(
            deposit_ray_terma(&mut dw, &mu, &central_x_ray(), 0.0, 0.5)
                .expect("valid attenuation volume")
                .into_base(),
            0.0,
            epsilon = 1e-15
        );
    }

    #[test]
    fn missing_ray_deposits_nothing() {
        let mu = uniform_cube(0.05);
        let mut dose = Volume::zeros(*mu.grid());
        let miss =
            Ray::try_new(Point3::new(-50.0, 500.0, 8.0), Vector3::new(1.0, 0.0, 0.0)).unwrap();
        assert_relative_eq!(
            deposit_ray_terma(&mut dose, &mu, &miss, 1.0, 0.5)
                .expect("valid attenuation volume")
                .into_base(),
            0.0,
            epsilon = 1e-15
        );
        assert_relative_eq!(dose.sum(), 0.0, epsilon = 1e-15);
    }

    #[test]
    fn invalid_attenuation_leaves_dose_unchanged() {
        let mu = uniform_cube(-0.05);
        let mut dose = Volume::from_shape_fn(*mu.grid(), |_| 2.0);
        let before = dose.as_slice().to_vec();
        let error = deposit_ray_terma(&mut dose, &mu, &central_x_ray(), 1.0, 0.5)
            .expect_err("negative attenuation must fail");
        assert!(matches!(
            error,
            TransportError::InvalidValue {
                field: hyperion::ValueKind::LinearAttenuation,
                ..
            }
        ));
        assert_eq!(dose.as_slice(), before);
    }

    /// Absorbed fraction `1 - exp(-mu*L)` is a subtractive cancellation: at
    /// `mu*L = 0.08` the result is ~0.0769, so the relative error of `exp` is
    /// amplified by `1/0.0769`, roughly 13x. Budgeting `exp` at a few ulps
    /// (IEEE 754 does not require it correctly rounded) and applying that
    /// amplification gives a few tens of ulps; 64 bounds it.
    const ABSORBED_FRACTION_ULPS: f64 = 64.0;

    /// The energy-conservation check compares the returned total against a sum
    /// over the 9x9x9 dose grid, so it carries both the cancellation above and
    /// the accumulation over the deposited voxels. Twice the single-value budget
    /// covers the independent second computation.
    const ENERGY_BALANCE_ULPS: f64 = 128.0;

    /// Asserts ray TERMA deposition and its energy balance in one scalar width.
    fn ray_deposition_conserves_absorbed_energy<T: ShippedScalar + helios_math::GeometryScalar>() {
        // Both FloatElement and GeometryScalar define `from_f64`; name the one
        // that performs the literal conversion so the call is unambiguous.
        let cast = <T as helios_math::FloatElement>::from_f64;
        let zero = cast(0.0);
        let spacing = cast(2.0);
        let grid =
            VoxelGrid::<T>::axis_aligned([9, 9, 9], [spacing; 3], Point3::new(zero, zero, zero))
                .expect("valid axis-aligned grid");

        let mu = cast(0.05);
        let attenuation = Volume::from_shape_fn(grid, |_| mu);
        let mut dose = Volume::zeros(*attenuation.grid());
        let ray = Ray::try_new(
            Point3::new(cast(-50.0), cast(8.0), cast(8.0)),
            Vector3::new(cast(1.0), zero, zero),
        )
        .expect("non-degenerate ray direction");

        let total = deposit_ray_terma(&mut dose, &attenuation, &ray, cast(1.0), cast(0.25))
            .expect("valid attenuation volume")
            .into_base();

        // Absorbed fraction over the 1.6 cm traversed path.
        let path_cm = cast(1.6);
        let expected = cast(1.0) - (-(mu * path_cm)).exp();
        assert_relative_eq!(
            total,
            expected,
            max_relative = T::EPSILON * cast(ABSORBED_FRACTION_ULPS)
        );

        // Energy conservation: what the kernel reports equals what it deposited.
        assert_relative_eq!(
            dose.sum(),
            total,
            max_relative = T::EPSILON * cast(ENERGY_BALANCE_ULPS)
        );
    }

    #[test]
    fn ray_deposition_conserves_absorbed_energy_in_single_precision() {
        ray_deposition_conserves_absorbed_energy::<f32>();
    }

    #[test]
    fn ray_deposition_conserves_absorbed_energy_in_double_precision() {
        ray_deposition_conserves_absorbed_energy::<f64>();
    }

    #[test]
    fn diverging_reduces_to_no_falloff_at_large_sad() {
        // As SAD → ∞ the inverse-square factor → 1 everywhere → the diverging
        // deposition matches the plain energy-conserving one.
        let mu = uniform_cube(0.05);
        let (mut plain_d, mut div_d) = (Volume::zeros(*mu.grid()), Volume::zeros(*mu.grid()));
        let plain = deposit_ray_terma(&mut plain_d, &mu, &central_x_ray(), 1.0, 0.25)
            .expect("valid attenuation volume")
            .into_base();
        let focal = Point3::new(-1.0e9, 8.0, 8.0);
        let div =
            deposit_ray_terma_diverging(&mut div_d, &mu, &central_x_ray(), 1.0, 0.25, focal, 1.0e9)
                .expect("valid attenuation volume")
                .into_base();
        assert_relative_eq!(div, plain, max_relative = 1e-6);
    }

    /// Relative tolerance, in ulps of `T`, for the inverse-square ratio.
    ///
    /// `deposit_ray_terma` and `deposit_ray_terma_diverging` share one march, so
    /// with identical inputs the pre-falloff per-segment loss is bit-identical
    /// and cancels in the per-voxel ratio: what is compared is
    /// `fl(loss · isf) / loss`, two roundings. The kernel's
    /// `isf = fl(sad²) · fl(1/fl(r²))` adds two more (here `sad = 28` and
    /// `r ∈ {21,…,35}` are integers, so the squares are exact), and the test's
    /// own `sad²/r²` one. Six roundings; eight bounds them — 9.5e-7 for `f32`,
    /// 1.8e-15 for `f64`.
    const INVERSE_SQUARE_ULPS: f64 = 8.0;

    /// Asserts the point-source inverse-square falloff in one scalar width.
    ///
    /// The oracle is geometric and independent of the transport: the step is set
    /// to the voxel pitch so each voxel receives exactly one segment, at a known
    /// midpoint, and the divergent deposition must differ from the parallel one
    /// by exactly `(SAD/r)²` there. The attenuation physics cancels in the
    /// ratio, so only the divergence law is under test.
    fn point_source_falloff_matches_inverse_square<
        T: ShippedScalar + helios_math::GeometryScalar,
    >() {
        let cast = <T as helios_math::FloatElement>::from_f64;
        let zero = cast(0.0);
        let spacing = cast(2.0);
        let grid =
            VoxelGrid::<T>::axis_aligned([9, 9, 9], [spacing; 3], Point3::new(zero, zero, zero))
                .expect("valid axis-aligned grid");
        let mu = Volume::from_shape_fn(grid, |_| cast(0.05));
        let ray = Ray::try_new(
            Point3::new(cast(-50.0), cast(8.0), cast(8.0)),
            Vector3::new(cast(1.0), zero, zero),
        )
        .expect("non-degenerate ray direction");

        // Step == pitch over the 16 mm node box: 8 segments with midpoints at
        // x = 1, 3, …, 15 mm, one per voxel i = 1..=8 (voxel 0 receives none).
        let (weight, step) = (cast(1.0), spacing);
        // Source 20 mm before the entry face; SAD 28 mm reaches the grid centre.
        let focal = Point3::new(cast(-20.0), cast(8.0), cast(8.0));
        let sad = cast(28.0);

        let mut parallel = Volume::zeros(grid);
        deposit_ray_terma(&mut parallel, &mu, &ray, weight, step)
            .expect("valid attenuation volume");
        let mut divergent = Volume::zeros(grid);
        deposit_ray_terma_diverging(&mut divergent, &mu, &ray, weight, step, focal, sad)
            .expect("valid attenuation volume");

        for i in 1..=8usize {
            let midpoint_mm = cast(2.0 * i as f64 - 1.0);
            let r = midpoint_mm - focal.x;
            let expected = sad * sad * (r * r).recip();
            let ratio = divergent.get(i, 4, 4).expect("in-grid voxel")
                * parallel.get(i, 4, 4).expect("in-grid voxel").recip();
            assert_relative_eq!(
                ratio,
                expected,
                max_relative = T::EPSILON * cast(INVERSE_SQUARE_ULPS)
            );
        }

        // The factor spans 784/441 = 1.78 at the entry voxel down to
        // 784/1225 = 0.64 at the exit voxel, so no constant factor — in
        // particular none at all — satisfies the loop above.
        let voxel = |v: &Volume<T>, i: usize| v.get(i, 4, 4).expect("in-grid voxel");
        assert!(
            voxel(&divergent, 1) > voxel(&parallel, 1),
            "near-source voxel must gain fluence"
        );
        assert!(
            voxel(&divergent, 8) < voxel(&parallel, 8),
            "far voxel must lose fluence"
        );
    }

    #[test]
    fn point_source_falloff_matches_inverse_square_in_single_precision() {
        point_source_falloff_matches_inverse_square::<f32>();
    }

    #[test]
    fn point_source_falloff_matches_inverse_square_in_double_precision() {
        point_source_falloff_matches_inverse_square::<f64>();
    }

    /// Water linear attenuation at the `TomoTherapy` 6 MV working point, `cm⁻¹`:
    /// the `(μ/ρ) = 0.06 cm²/g` used across the Helios fixtures at unit density.
    const WATER_MU_PER_CM: f64 = 0.06;

    /// Source-to-surface distance for the depth-dose phantom, mm.
    const PDD_SSD_MM: f64 = 800.0;

    /// Relative tolerance, in ulps of `T`, for the primary depth-dose ratio.
    ///
    /// The per-voxel deposit is a subtractive cancellation:
    /// `e^{−τ} − e^{−(τ+μΔ)} = e^{−τ}(1 − e^{−0.12}) = 0.1131·e^{−τ}`, so the
    /// relative error of the two transmissions is amplified by `1/0.1131 ≈ 8.8`.
    /// Each transmission carries the accumulated optical depth (up to seven
    /// `checked_add`s summing to 0.84, ≲4.2 ulps absolute) plus `exp`, which
    /// IEEE 754 does not require correctly rounded (≈2 ulps): ≈6.2 ulps, doubled
    /// for the pair and amplified gives ≈110 ulps at the deepest voxel and
    /// ≈66 at the reference (whose leading transmission is exactly 1). With the
    /// inverse-square factor, the ratio, and the oracle's own `exp` and squaring,
    /// ≈191 ulps; 256 bounds it — 3.0e-5 for `f32`, 5.7e-14 for `f64`.
    const DEPTH_DOSE_ULPS: f64 = 256.0;

    /// Asserts percentage depth dose in water in one scalar width.
    ///
    /// In a uniform medium with one ray segment per voxel, the primary terma at
    /// depth `d` relative to a reference depth `d₀` obeys the textbook law
    ///
    /// ```text
    /// PDD(d) = exp(−μ·(d − d₀)) · ((SSD + d₀)/(SSD + d))²
    /// ```
    ///
    /// — exponential attenuation of the primary fluence times inverse-square
    /// divergence. The `weight·(1 − e^{−μΔ})` factor common to every segment
    /// cancels in the ratio, so the oracle is that analytical law rather than a
    /// restatement of the kernel's recursion. Buildup and phantom scatter are
    /// absent by construction: this kernel transports the primary beam only, so
    /// the curve is monotonically decreasing from the surface.
    fn depth_dose_in_water_matches_the_primary_law<
        T: ShippedScalar + helios_math::GeometryScalar,
    >() {
        let cast = <T as helios_math::FloatElement>::from_f64;
        let zero = cast(0.0);
        // 16 cm of water along x, 2 cm voxels; the beam runs down the y = z = 20 mm axis.
        let pitch_mm = 20.0;
        let spacing = cast(pitch_mm);
        let grid =
            VoxelGrid::<T>::axis_aligned([9, 3, 3], [spacing; 3], Point3::new(zero, zero, zero))
                .expect("valid axis-aligned grid");
        let mu = Volume::from_shape_fn(grid, |_| cast(WATER_MU_PER_CM));

        let focal = Point3::new(cast(-PDD_SSD_MM), cast(20.0), cast(20.0));
        let ray = Ray::try_new(focal, Vector3::new(cast(1.0), zero, zero))
            .expect("non-degenerate ray direction");
        let mut dose = Volume::zeros(grid);
        deposit_ray_terma_diverging(
            &mut dose,
            &mu,
            &ray,
            cast(1.0),
            spacing,
            focal,
            cast(PDD_SSD_MM + 100.0),
        )
        .expect("valid attenuation volume");

        // Step == pitch: segment midpoints at depth 10, 30, …, 150 mm land one
        // per voxel i = 1..=8. Voxel 1 (d₀ = 10 mm) is the reference depth.
        let depth_mm = |i: usize| pitch_mm * i as f64 - pitch_mm / 2.0;
        let reference_depth = depth_mm(1);
        let reference = dose.get(1, 1, 1).expect("in-grid voxel");

        for i in 2..=8usize {
            let depth = depth_mm(i);
            // exp(−μ·Δdepth) with μ in cm⁻¹ and the depth difference in cm.
            let attenuation = cast(-WATER_MU_PER_CM * (depth - reference_depth) / MM_PER_CM).exp();
            let divergence = cast(((PDD_SSD_MM + reference_depth) / (PDD_SSD_MM + depth)).powi(2));
            let expected = attenuation * divergence;
            let pdd = dose.get(i, 1, 1).expect("in-grid voxel") * reference.recip();
            assert_relative_eq!(
                pdd,
                expected,
                max_relative = T::EPSILON * cast(DEPTH_DOSE_ULPS)
            );
            // Primary-only transport: no buildup, so the curve only descends.
            assert!(
                pdd < <T as helios_math::NumericElement>::ONE,
                "primary depth dose must fall below the reference depth"
            );
        }

        // The divergence term is not decorative: at 15 cm the full law sits
        // strictly below pure attenuation, by the (810/950)² = 0.727 factor.
        let deepest = dose.get(8, 1, 1).expect("in-grid voxel") * reference.recip();
        let attenuation_only =
            cast(-WATER_MU_PER_CM * (depth_mm(8) - reference_depth) / MM_PER_CM).exp();
        assert!(
            deepest < attenuation_only,
            "inverse-square divergence must steepen the depth-dose curve: \
             {deepest:?} !< {attenuation_only:?}"
        );
    }

    #[test]
    fn depth_dose_in_water_matches_the_primary_law_in_single_precision() {
        depth_dose_in_water_matches_the_primary_law::<f32>();
    }

    #[test]
    fn depth_dose_in_water_matches_the_primary_law_in_double_precision() {
        depth_dose_in_water_matches_the_primary_law::<f64>();
    }

    #[test]
    fn inverse_square_steepens_the_entry_to_exit_ratio() {
        // Point source 20 mm before the entry face (SAD = 28 mm to the centre):
        // the entry voxel (near source, isf > 1) gains dose relative to the exit
        // voxel (far, isf < 1) beyond the pure-attenuation ratio.
        let mu = uniform_cube(0.1);
        let mut plain_d = Volume::zeros(*mu.grid());
        let _ = deposit_ray_terma(&mut plain_d, &mu, &central_x_ray(), 1.0, 0.1)
            .expect("valid attenuation volume");
        let mut div_d = Volume::zeros(*mu.grid());
        let focal = Point3::new(-20.0, 8.0, 8.0);
        let _ =
            deposit_ray_terma_diverging(&mut div_d, &mu, &central_x_ray(), 1.0, 0.1, focal, 28.0)
                .expect("valid attenuation volume");

        let ratio_plain = plain_d.get(0, 4, 4).unwrap() / plain_d.get(8, 4, 4).unwrap();
        let ratio_div = div_d.get(0, 4, 4).unwrap() / div_d.get(8, 4, 4).unwrap();
        assert!(
            ratio_div > ratio_plain,
            "inverse-square should steepen entry/exit ratio: {ratio_div} !> {ratio_plain}"
        );
    }
}
