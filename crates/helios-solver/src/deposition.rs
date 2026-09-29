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
mod tests;
