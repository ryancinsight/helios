//! Primary-fluence energy deposition (TERMA) along a beam ray.
//!
//! Companion to the [`forward_project_ray`](crate::forward_project_ray) line
//! integral: where the projector reduces `∫ μ dl` to a single optical depth, this
//! kernel *deposits* the energy the primary beam loses as it attenuates, voxel by
//! voxel, producing the terma (total energy released per unit mass) that a
//! collapsed-cone/convolution dose engine spreads with a scatter kernel.

use crate::projector::{ray_grid_interval, sample_volume_along_ray, RayMarchPlan};
use aequitas::systems::si::{
    quantities::{AbsorbedDose, Length, ReciprocalLength},
    units::{Centimeter, PerCentimeter},
};
use eunomia::UnitScalar;
use helios_domain::Volume;
use helios_math::{GeometryScalar, NumericElement, Point3, Ray};
use hyperion::{
    coefficient::{InteractionCoefficient, LinearAttenuation},
    quantity::{OpticalDepth, PathLength},
    TransportError,
};

#[inline]
fn segment_optical_depth<T: GeometryScalar + UnitScalar>(
    mu_sample: T,
    path: PathLength<T>,
) -> Result<OpticalDepth<T>, TransportError<T>> {
    let coefficient =
        InteractionCoefficient::<T, LinearAttenuation>::new(ReciprocalLength::from_unit::<
            PerCentimeter,
        >(mu_sample))?;
    coefficient.optical_depth(path)
}

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
/// that misses the grid deposits nothing and returns zero. Deposition covers the
/// half-line `t ≥ 0`: an origin inside the grid deposits only ahead of itself,
/// and the voxels behind the origin are untouched. The returned total is
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
    let Some(plan) = RayMarchPlan::from_interval(t_enter, t_exit, step_mm) else {
        return Ok(AbsorbedDose::from_base(T::ZERO));
    };
    let path = PathLength::new(Length::from_unit::<Centimeter>(plan.step_cm()))?;
    let [nx, ny, nz] = grid.dims();

    // Validate the complete ray before mutating the output. Non-negative
    // segment depths make every partial sum bounded by this checked total.
    let _validated_total = (0..plan.steps()).try_fold(OpticalDepth::zero(), |total, i| {
        let (_, _, mu_sample) = sample_volume_along_ray(mu, ray, plan, i);
        total.checked_add(segment_optical_depth(mu_sample, path)?)
    })?;

    let mut optical_depth = OpticalDepth::zero();
    let mut trans_before = <T as NumericElement>::ONE; // e^{−τ} at τ = 0.
    let mut total = T::ZERO;
    for i in 0..plan.steps() {
        let (world_pt, index, mu_sample) = sample_volume_along_ray(mu, ray, plan, i);
        optical_depth = optical_depth.checked_add(segment_optical_depth(mu_sample, path)?)?;
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
    use crate::test_support::{
        quarter_turn_rotation, ORIENTED_FRAME_ROUNDINGS, RAY_DIRECTION_ROUNDINGS,
        SLAB_ENDPOINT_ROUNDINGS,
    };
    use aequitas::systems::si::units::Millimeter;
    use eunomia::assert_relative_eq;
    use helios_domain::VoxelGrid;
    use helios_math::ShippedScalar;
    use helios_math::{Point3, Vector3};

    /// Node-box extent of the 9-voxel, 2 mm cube along each axis.
    const CUBE_EXTENT_MM: f64 = 16.0;

    /// Distance from the central beam's origin to the cube's entry face.
    const BEAM_STANDOFF_MM: f64 = 50.0;

    // Uniform-μ cube: 9³ voxels, 2 mm spacing → node box [0,16] mm = 1.6 cm/axis.
    fn uniform_cube(mu_val: f64) -> Volume<f64> {
        let grid = VoxelGrid::axis_aligned(
            [9, 9, 9],
            [Length::from_unit::<Millimeter>(2.0); 3],
            Point3::new(0.0, 0.0, 0.0),
        )
        .expect("grid");
        Volume::from_shape_fn(grid, move |_| mu_val)
    }

    // +x ray through the cube centre (y = z = 8 mm), starting outside the box.
    fn central_x_ray() -> Ray<f64> {
        Ray::try_new(
            Point3::new(-BEAM_STANDOFF_MM, 8.0, 8.0),
            Vector3::new(1.0, 0.0, 0.0),
        )
        .expect("unit +x ray")
    }

    fn oriented_cube(mu_val: f64) -> Volume<f64> {
        let grid = VoxelGrid::oriented(
            [9, 9, 9],
            [Length::from_unit::<Millimeter>(2.0); 3],
            Point3::new(10.0, 20.0, 30.0),
            quarter_turn_rotation(),
        )
        .expect("grid");
        Volume::from_shape_fn(grid, move |_| mu_val)
    }

    /// Unit roundoff `u = EPSILON / 2` of `T`.
    fn unit_roundoff<T: ShippedScalar>() -> T {
        T::EPSILON / T::from_f64(2.0)
    }

    /// Assumed relative error of `exp`, in `u`: "a few ulps", since IEEE 754
    /// does not require it correctly rounded.
    const EXP_ROUNDINGS: f64 = 4.0;

    /// Roundings of one segment's optical depth that do not scale with the
    /// segment count, 9 in all:
    ///
    /// - `mu` to base units (1): `from_unit::<PerCentimeter>` multiplies by
    ///   the exactly representable `1e2`, so only the product rounds.
    /// - The step to base units (2): `Length::from_unit::<Centimeter>` is
    ///   `value * from_f64(1e-2)` (aequitas `quantity/construction.rs`
    ///   `from_unit`; eunomia `impls/unit.rs` `scale_by_factor`;
    ///   `Centimeter::SCALE = 1.0e-2` in aequitas
    ///   `systems/si/units/scaled.rs`, at the locked revisions). `0.01` is
    ///   not representable, so the constant rounds (relative error up to `u`)
    ///   and then the product does.
    /// - The product of the two base-unit values (1).
    /// - `step_cm = step_mm * recip(10)` (2): the reciprocal and the product.
    /// - `step_mm = length * recip(n)` (2): the reciprocal and the product.
    /// - The subtraction forming `length` (1).
    const SEGMENT_DEPTH_ROUNDINGS: f64 = 9.0;

    /// Conditioning `(|t_enter| + |t_exit|) / length` of the chord length
    /// against errors in the two slab parameters.
    fn slab_conditioning<T: ShippedScalar>(t_enter: T, t_exit: T) -> T {
        (t_enter.abs() + t_exit.abs()) / (t_exit - t_enter)
    }

    /// Conditioning of the central beam: entry at the standoff, exit one cube
    /// extent beyond it.
    fn central_beam_conditioning<T: ShippedScalar>() -> T {
        slab_conditioning(
            T::from_f64(BEAM_STANDOFF_MM),
            T::from_f64(BEAM_STANDOFF_MM + CUBE_EXTENT_MM),
        )
    }

    /// Cancellation factor `T / (1 - T) = 1 / expm1(tau)` that carries a
    /// relative error of `T = exp(-tau)` into the absorbed fraction `1 - T`.
    fn cancellation<T: ShippedScalar>(tau: T) -> T {
        tau.exp_m1().recip()
    }

    /// First-order worst-case relative error of a kernel total over `segments`
    /// uniform segments of total optical depth `tau`, with every rounding at its
    /// maximum and all signs aligned. `pose_roundings` is
    /// `ORIENTED_FRAME_ROUNDINGS` for an oriented grid and zero otherwise.
    ///
    /// The terms, in `u`: the relative error of `tau` is `segments - 1`
    /// (accumulation), `SEGMENT_DEPTH_ROUNDINGS`, the two slab endpoints
    /// `SLAB_ENDPOINT_ROUNDINGS * conditioning`, and the direction terms
    /// `RAY_DIRECTION_ROUNDINGS` and `pose_roundings`, which scale both
    /// endpoints alike and are not amplified (see their docs). It reaches the
    /// absorbed fraction `1 - T` through `tau * T / (1 - T)`, as does the
    /// `exp` error of `T` through `T / (1 - T)`. Each difference of
    /// neighbouring transmissions is exact (Sterbenz: the ratio lies in
    /// `[1/2, 2]`), so the total telescopes to `weight * (1 - T)` up to its
    /// `segments`-term positive sum and the per-term weight product:
    /// `segments`.
    ///
    /// Valid for `tau > 0` with `tau / segments <= ln 2` (the Sterbenz
    /// condition) and `k * u << 1` for the `k` roundings counted.
    ///
    /// The model is first order and carries no margin for the neglected
    /// terms: `k` roundings of at most `u` compound to at most
    /// `(1 + u)^k - 1 <= k u / (1 - k u)`, so the bound falls short of that
    /// worst case by a relative `k u`, second order in `u`. Here `k` is the
    /// coefficient of `u` in this function's result. At the central beam of
    /// `ray_deposition_conserves_absorbed_energy` (64 segments, `tau = 0.08`,
    /// conditioning 7.25, no pose) `k` is 217.4, so the shortfall is about
    /// `1.30e-5` of the bound in `f32` and `2.4e-14` in `f64`
    /// (`kernel_error_note_figures_reproduce` checks them). The bound is not
    /// a rigorous upper bound.
    fn kernel_total_error<T: ShippedScalar>(
        segments: u32,
        tau: T,
        conditioning: T,
        pose_roundings: T,
    ) -> T {
        let n = T::from_f64(f64::from(segments));
        let depth = n - T::from_f64(1.0)
            + T::from_f64(SEGMENT_DEPTH_ROUNDINGS)
            + T::from_f64(SLAB_ENDPOINT_ROUNDINGS) * conditioning
            + T::from_f64(RAY_DIRECTION_ROUNDINGS)
            + pose_roundings;
        unit_roundoff::<T>() * ((T::from_f64(EXP_ROUNDINGS) + depth * tau) * cancellation(tau) + n)
    }

    /// First-order worst-case relative error of the closed-form reference
    /// `1 - exp(-mu * path)`: the rounded `path` literal and the product (2
    /// in `tau`), the `exp` error, and the final subtraction. Valid for
    /// `tau > 0` and `k * u << 1`, as [`kernel_total_error`].
    fn closed_form_error<T: ShippedScalar>(tau: T) -> T {
        unit_roundoff::<T>()
            * ((T::from_f64(EXP_ROUNDINGS) + T::from_f64(2.0) * tau) * cancellation(tau)
                + T::from_f64(1.0))
    }

    /// Worst-case relative disagreement of two summations of the same `terms`
    /// positive values in different orders: each is within `(terms - 1) u`.
    /// Valid for `terms >= 1` and `terms * u << 1`.
    fn reordered_sum_error<T: ShippedScalar>(terms: u32) -> T {
        T::from_f64(2.0 * f64::from(terms - 1)) * unit_roundoff::<T>()
    }

    /// Worst-case relative disagreement of a kernel total with the closed-form
    /// absorbed fraction `1 - exp(-tau)`: the kernel error and the reference
    /// error add, since each carries its own roundings. The bound for the
    /// uniform-field fixtures of a grid-aligned ray, which has no pose.
    fn absorbed_fraction_bound<T: ShippedScalar>(segments: u32, tau: T, conditioning: T) -> T {
        kernel_total_error(segments, tau, conditioning, <T as NumericElement>::ZERO)
            + closed_form_error(tau)
    }

    /// Half unit of the last printed digit of `k = 217.4` in the
    /// [`kernel_total_error`] note, relative to `k`.
    const K_NOTE_HALF_UNIT: f64 = 0.05 / 217.4;

    /// Asserts the note's figures on `kernel_total_error` for the central beam
    /// of `ray_deposition_conserves_absorbed_energy`: `k = 217.4` and the
    /// bound `k u` equal to `figure`, each within the half unit of its last
    /// printed digit (`figure_half_unit`, relative). The bound is evaluated
    /// in `T` with about a dozen roundings, at most `12 u` (`7e-7` in `f32`),
    /// below the margins of the printed digits (`2.5e-5` for `k`, `8e-4` for
    /// the `f32` figure).
    fn kernel_error_note_figures_reproduce<T: ShippedScalar>(figure: f64, figure_half_unit: f64) {
        let zero = <T as NumericElement>::ZERO;
        let bound = kernel_total_error(64, T::from_f64(0.08), central_beam_conditioning(), zero);
        assert_relative_eq!(
            bound / unit_roundoff::<T>(),
            T::from_f64(217.4),
            epsilon = zero,
            max_relative = T::from_f64(K_NOTE_HALF_UNIT)
        );
        assert_relative_eq!(
            bound,
            T::from_f64(figure),
            epsilon = zero,
            max_relative = T::from_f64(figure_half_unit)
        );
    }

    #[test]
    fn kernel_error_note_figures_reproduce_in_single_precision() {
        kernel_error_note_figures_reproduce::<f32>(1.30e-5, 0.005 / 1.30);
    }

    #[test]
    fn kernel_error_note_figures_reproduce_in_double_precision() {
        kernel_error_note_figures_reproduce::<f64>(2.4e-14, 0.05 / 2.4);
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
        // 32 segments at the 0.5 mm step over the 16 mm chord.
        let bound = absorbed_fraction_bound(32, 0.08, central_beam_conditioning());
        assert_relative_eq!(total, expected, max_relative = bound);
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
        // 32 segments; the ray enters at 40 mm and exits at 56 mm along +y.
        let total_bound = kernel_total_error(
            32,
            0.08,
            slab_conditioning(40.0, 56.0),
            ORIENTED_FRAME_ROUNDINGS,
        ) + closed_form_error(0.08);
        assert_relative_eq!(total, expected, max_relative = total_bound);
        // The voxel sum adds one reordering of the 32 positive segment terms.
        assert_relative_eq!(
            dose.sum(),
            expected,
            max_relative = total_bound + reordered_sum_error::<f64>(32)
        );
    }

    #[test]
    fn deposited_voxels_conserve_the_returned_total() {
        // The voxel sum and the returned total add the same 54 positive terms
        // in different orders, so they agree to the reordering error.
        let mu = uniform_cube(0.08);
        let mut dose = Volume::zeros(*mu.grid());
        let total = deposit_ray_terma(&mut dose, &mu, &central_x_ray(), 3.0, 0.3)
            .expect("valid attenuation volume")
            .into_base();
        // 16 mm / 0.3 mm rounds up to 54 segments.
        assert_relative_eq!(dose.sum(), total, max_relative = reordered_sum_error(54));
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
        // 16 mm / 4 mm = 4 segments and 16 mm / 0.05 mm = 320 segments; two
        // independent totals differ by at most the sum of their errors.
        let conditioning = central_beam_conditioning();
        let bound = kernel_total_error(4, 0.08, conditioning, 0.0)
            + kernel_total_error(320, 0.08, conditioning, 0.0);
        assert_relative_eq!(coarse, fine, max_relative = bound);
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
        // tau = 0 exactly, so T = exp(0) = 1 and every difference is zero.
        assert_eq!(
            deposit_ray_terma(&mut d0, &empty, &central_x_ray(), 1.0, 0.5)
                .expect("valid attenuation volume")
                .into_base(),
            0.0
        );
        assert_eq!(d0.sum(), 0.0);

        let mu = uniform_cube(0.05);
        let mut dw = Volume::zeros(*mu.grid());
        // Zero weight scales every finite loss to exactly zero.
        assert_eq!(
            deposit_ray_terma(&mut dw, &mu, &central_x_ray(), 0.0, 0.5)
                .expect("valid attenuation volume")
                .into_base(),
            0.0
        );
    }

    #[test]
    fn missing_ray_deposits_nothing() {
        let mu = uniform_cube(0.05);
        let mut dose = Volume::zeros(*mu.grid());
        let miss =
            Ray::try_new(Point3::new(-50.0, 500.0, 8.0), Vector3::new(1.0, 0.0, 0.0)).unwrap();
        // A miss returns before any arithmetic and leaves the zeroed dose as is.
        assert_eq!(
            deposit_ray_terma(&mut dose, &mu, &miss, 1.0, 0.5)
                .expect("valid attenuation volume")
                .into_base(),
            0.0
        );
        assert_eq!(dose.sum(), 0.0);
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

    /// Asserts ray TERMA deposition and its energy balance in one scalar width.
    ///
    /// The central beam crosses the 16 mm chord in 64 segments of 0.25 mm at
    /// `tau = 0.08`, so `1 - exp(-tau)` is about 0.077 and the cancellation
    /// factor `1 / expm1(tau)` amplifies each relative error by about 12. The
    /// total agrees with the closed form within `absorbed_fraction_bound` of
    /// those 64 segments at the beam's conditioning. The voxel sum and the
    /// total add the same 64 positive terms in different orders, so they agree
    /// within `reordered_sum_error(64)`.
    fn ray_deposition_conserves_absorbed_energy<T: ShippedScalar + helios_math::GeometryScalar>() {
        // Both FloatElement and GeometryScalar define `from_f64`; name the one
        // that performs the literal conversion so the call is unambiguous.
        let cast = <T as helios_math::FloatElement>::from_f64;
        let zero = cast(0.0);
        let spacing = cast(2.0);
        let grid = VoxelGrid::<T>::axis_aligned(
            [9, 9, 9],
            [Length::from_unit::<Millimeter>(spacing); 3],
            Point3::new(zero, zero, zero),
        )
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
        let tau = mu * path_cm;
        let expected = cast(1.0) - (-tau).exp();
        // `epsilon = zero` removes the macro's default absolute floor
        // (`EPSILON`, about 26u relative for values near 0.077). It does not
        // bind at the derived bound, which is larger, but it would accept any
        // bound or regression that shrinks below it, so the derived bound is
        // made the only tolerance.
        assert_relative_eq!(
            total,
            expected,
            epsilon = zero,
            max_relative = absorbed_fraction_bound(64, tau, central_beam_conditioning())
        );

        // Energy conservation: what the kernel reports equals what it deposited.
        assert_relative_eq!(
            dose.sum(),
            total,
            epsilon = zero,
            max_relative = reordered_sum_error(64)
        );
    }

    /// Asserts TERMA deposition from an origin inside the grid in one scalar
    /// width: only the forward half-line deposits, the total follows
    /// Beer-Lambert over the forward chord, and the voxels behind the origin
    /// stay exactly zero.
    ///
    /// `mu = 0.5 cm^-1` over the 0.8 cm forward chord gives `tau = 0.4` in 32
    /// segments of 0.25 mm. The entry is clamped to the origin, an exact
    /// `t = 0`, so only the exit endpoint at `t = 8 mm` carries slab roundings
    /// and the conditioning is `slab_conditioning(0, 8) = 1`.
    ///
    /// A clamped entry puts the origin error `delta o` in `t_exit` alone,
    /// adding `|delta o| / (bound - o)` (see `SLAB_ENDPOINT_ROUNDINGS`). Here
    /// `delta o = 0` in both widths: the pose is the identity, whose Rodrigues
    /// form returns a point unchanged, the grid origin is zero, and the 2 mm
    /// spacing is a power of two, so `world_to_index`, which the kernel
    /// applies to the ray origin, maps `(8, 8, 8)` mm to `(4, 4, 4)` exactly.
    /// The fixture asserts that mapping, so the term is absent from the bound
    /// by construction and not by omission.
    ///
    /// The total agrees with the closed form within
    /// `absorbed_fraction_bound`, and the voxel sum within
    /// `reordered_sum_error(32)` of the total.
    fn inside_origin_deposits_only_the_forward_half_line<
        T: ShippedScalar + helios_math::GeometryScalar,
    >() {
        let cast = <T as helios_math::FloatElement>::from_f64;
        let zero = cast(0.0);
        let grid = VoxelGrid::<T>::axis_aligned(
            [9, 9, 9],
            [Length::from_unit::<Millimeter>(cast(2.0)); 3],
            Point3::new(zero, zero, zero),
        )
        .expect("valid axis-aligned grid");
        let centre = cast(8.0);
        let origin_index = grid.world_to_index(Point3::new(centre, centre, centre));
        let index_centre = cast(4.0);
        assert_eq!(
            (origin_index.x, origin_index.y, origin_index.z),
            (index_centre, index_centre, index_centre),
            "the ray origin maps exactly, so the origin error is zero"
        );
        let mu = cast(0.5);
        let attenuation = Volume::from_shape_fn(grid, |_| mu);
        let mut dose = Volume::zeros(*attenuation.grid());
        let ray = Ray::try_new(
            Point3::new(centre, centre, centre),
            Vector3::new(cast(1.0), zero, zero),
        )
        .expect("non-degenerate ray direction");

        let total = deposit_ray_terma(&mut dose, &attenuation, &ray, cast(1.0), cast(0.25))
            .expect("valid attenuation volume")
            .into_base();

        let tau = mu * cast(0.8);
        let expected = cast(1.0) - (-tau).exp();
        // The ray runs along +x with a unit direction, so the exit parameter
        // is the distance from the origin to the node box's far face.
        let exit_mm = cast(CUBE_EXTENT_MM) - centre;
        assert_relative_eq!(
            total,
            expected,
            epsilon = zero,
            max_relative = absorbed_fraction_bound(32, tau, slab_conditioning(zero, exit_mm))
        );
        assert_relative_eq!(
            dose.sum(),
            total,
            epsilon = zero,
            max_relative = reordered_sum_error(32)
        );
        for ix in 0..4 {
            for iy in 0..9 {
                for iz in 0..9 {
                    assert_eq!(
                        dose.get(ix, iy, iz).expect("voxel in range"),
                        zero,
                        "voxel ({ix}, {iy}, {iz}) lies behind the ray origin"
                    );
                }
            }
        }
    }

    #[test]
    fn inside_origin_deposits_only_the_forward_half_line_in_single_precision() {
        inside_origin_deposits_only_the_forward_half_line::<f32>();
    }

    #[test]
    fn inside_origin_deposits_only_the_forward_half_line_in_double_precision() {
        inside_origin_deposits_only_the_forward_half_line::<f64>();
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
        let focal = Point3::new(-LARGE_SAD_MM, 8.0, 8.0);
        let div = deposit_ray_terma_diverging(
            &mut div_d,
            &mu,
            &central_x_ray(),
            1.0,
            0.25,
            focal,
            LARGE_SAD_MM,
        )
        .expect("valid attenuation volume")
        .into_base();
        assert_relative_eq!(div, plain, max_relative = large_sad_bound());
    }

    /// Source-to-axis distance, and focal-point offset, that stand in for
    /// "infinite" in the large-SAD limit.
    const LARGE_SAD_MM: f64 = 1.0e9;

    /// Roundings of the falloff factor `fl(sad^2) * fl(1 / fl(dx^2))` applied to
    /// a segment loss, in `u`: the focal offset `dx` (1, doubled to 2 in `dx^2`),
    /// the square (1), the reciprocal (1), the product by `sad^2` (1) and the
    /// product into the loss (1). `dy` and `dz` are exactly zero on the axis, so
    /// the sum of squares adds none, and `sad^2 = 1e18` is exact in `f64`.
    const FALLOFF_ROUNDINGS: f64 = 6.0;

    /// Bound on the relative difference of the diverging and plain totals at
    /// `LARGE_SAD_MM`.
    ///
    /// Segment midpoints lie at `x` in `[0, CUBE_EXTENT_MM]` ahead of a focus at
    /// `-R`, so the exact falloff is `R^2 / (R + x)^2` in `[1 - 2x / R, 1]`.
    /// Both totals sum the same 64 positive segment losses (16 mm at the 0.25 mm
    /// step) and the falloff weights them by factors inside that interval, so
    /// the totals differ by at most `2 * CUBE_EXTENT_MM / R`: a truncation of
    /// the limit, not a rounding. Rounding adds the falloff roundings above and
    /// the two 64-term sums, `2 * 63` in `u`.
    fn large_sad_bound() -> f64 {
        2.0 * CUBE_EXTENT_MM / LARGE_SAD_MM
            + unit_roundoff::<f64>() * (FALLOFF_ROUNDINGS + 2.0 * 63.0)
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
        let grid = VoxelGrid::<T>::axis_aligned(
            [9, 9, 9],
            [Length::from_unit::<Millimeter>(spacing); 3],
            Point3::new(zero, zero, zero),
        )
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
        let grid = VoxelGrid::<T>::axis_aligned(
            [9, 3, 3],
            [Length::from_unit::<Millimeter>(spacing); 3],
            Point3::new(zero, zero, zero),
        )
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
            let attenuation = cast(
                -WATER_MU_PER_CM * (depth - reference_depth)
                    / helios_quantities::constants::MM_PER_CM,
            )
            .exp();
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
        let attenuation_only = cast(
            -WATER_MU_PER_CM * (depth_mm(8) - reference_depth)
                / helios_quantities::constants::MM_PER_CM,
        )
        .exp();
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
