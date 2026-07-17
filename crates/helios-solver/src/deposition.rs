//! Primary-fluence energy deposition (TERMA) along a beam ray.
//!
//! Companion to the [`forward_project_ray`](crate::forward_project_ray) line
//! integral: where the projector reduces `∫ μ dl` to a single optical depth, this
//! kernel *deposits* the energy the primary beam loses as it attenuates, voxel by
//! voxel, producing the terma (total energy released per unit mass) that a
//! collapsed-cone/convolution dose engine spreads with a scatter kernel.

use crate::projector::ray_grid_interval;
use helios_core::constants::MM_PER_CM;
use helios_domain::Volume;
use helios_math::{GeometryScalar, NumericElement, Point3, Ray};

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
/// that misses the grid deposits nothing and returns zero.
#[must_use]
pub fn deposit_ray_terma<T: GeometryScalar>(
    dose: &mut Volume<T>,
    mu: &Volume<T>,
    ray: &Ray<T>,
    weight: T,
    step_mm: T,
) -> T {
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
#[must_use]
pub fn deposit_ray_terma_diverging<T: GeometryScalar>(
    dose: &mut Volume<T>,
    mu: &Volume<T>,
    ray: &Ray<T>,
    weight: T,
    step_mm: T,
    focal: Point3<T>,
    sad_mm: T,
) -> T {
    deposit_terma_impl(dose, mu, ray, weight, step_mm, Some((focal, sad_mm)))
}

/// Shared ray-march for [`deposit_ray_terma`] and [`deposit_ray_terma_diverging`];
/// `falloff = Some((focal, sad))` applies the inverse-square divergence factor.
fn deposit_terma_impl<T: GeometryScalar>(
    dose: &mut Volume<T>,
    mu: &Volume<T>,
    ray: &Ray<T>,
    weight: T,
    step_mm: T,
    falloff: Option<(Point3<T>, T)>,
) -> T {
    let grid = *mu.grid();
    debug_assert_eq!(
        grid.dims(),
        dose.grid().dims(),
        "dose and mu must share the same grid"
    );
    let (t_enter, t_exit) = match ray_grid_interval(&grid, ray) {
        Some(v) => v,
        None => return T::ZERO,
    };
    let length = t_exit - t_enter;
    if length <= T::ZERO {
        return T::ZERO;
    }

    // Substeps so the step divides the traversed length exactly (>= 1).
    let n = ((length * step_mm.recip()).ceil().to_f64() as usize).max(1);
    let actual_step = length * <T as GeometryScalar>::from_f64(n as f64).recip();
    let step_cm = actual_step * <T as GeometryScalar>::from_f64(MM_PER_CM).recip();
    let half = <T as GeometryScalar>::from_f64(0.5);
    let [nx, ny, nz] = grid.dims();

    let mut tau = T::ZERO;
    let mut trans_before = <T as NumericElement>::ONE; // e^{−τ} at τ = 0.
    let mut total = T::ZERO;
    for i in 0..n {
        let t_mid = t_enter + (<T as GeometryScalar>::from_f64(i as f64) + half) * actual_step;
        let world_pt: Point3<T> = ray.point_at(t_mid);
        let index = grid.world_to_index(world_pt);
        let mu_sample = mu.sample_trilinear(index).unwrap_or(T::ZERO);
        tau += mu_sample * step_cm;
        let trans_after = (-tau).exp();
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
    total
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use helios_domain::VoxelGrid;
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
        let total = deposit_ray_terma(&mut dose, &mu, &central_x_ray(), 1.0, 0.5);
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
        let total = deposit_ray_terma(&mut dose, &mu, &ray, 1.0, 0.5);
        let expected = 1.0 - (-0.05 * 1.6_f64).exp();
        assert_relative_eq!(total, expected, epsilon = 1e-12);
        assert_relative_eq!(dose.sum(), expected, epsilon = 1e-12);
    }

    #[test]
    fn deposited_voxels_conserve_the_returned_total() {
        // Sum of scattered voxel dose must equal the returned total exactly.
        let mu = uniform_cube(0.08);
        let mut dose = Volume::zeros(*mu.grid());
        let total = deposit_ray_terma(&mut dose, &mu, &central_x_ray(), 3.0, 0.3);
        assert_relative_eq!(dose.sum(), total, epsilon = 1e-12);
    }

    #[test]
    fn step_size_does_not_change_the_total() {
        // Telescoping makes the total independent of the sampling step.
        let mu = uniform_cube(0.05);
        let (mut d_coarse, mut d_fine) = (Volume::zeros(*mu.grid()), Volume::zeros(*mu.grid()));
        let coarse = deposit_ray_terma(&mut d_coarse, &mu, &central_x_ray(), 1.0, 4.0);
        let fine = deposit_ray_terma(&mut d_fine, &mu, &central_x_ray(), 1.0, 0.05);
        assert_relative_eq!(coarse, fine, epsilon = 1e-12);
    }

    #[test]
    fn energy_is_front_loaded_by_attenuation() {
        // Primary fluence decays with depth, so the entry voxel absorbs more than
        // the exit voxel along the beam.
        let mu = uniform_cube(0.3); // strong attenuation to make the gradient clear
        let mut dose = Volume::zeros(*mu.grid());
        let _ = deposit_ray_terma(&mut dose, &mu, &central_x_ray(), 1.0, 0.1);
        let entry = dose.get(0, 4, 4).unwrap();
        let exit = dose.get(8, 4, 4).unwrap();
        assert!(entry > exit, "entry {entry} should exceed exit {exit}");
    }

    #[test]
    fn zero_attenuation_and_zero_weight_deposit_nothing() {
        let empty = uniform_cube(0.0);
        let mut d0 = Volume::zeros(*empty.grid());
        assert_relative_eq!(
            deposit_ray_terma(&mut d0, &empty, &central_x_ray(), 1.0, 0.5),
            0.0,
            epsilon = 1e-15
        );
        assert_relative_eq!(d0.sum(), 0.0, epsilon = 1e-15);

        let mu = uniform_cube(0.05);
        let mut dw = Volume::zeros(*mu.grid());
        assert_relative_eq!(
            deposit_ray_terma(&mut dw, &mu, &central_x_ray(), 0.0, 0.5),
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
            deposit_ray_terma(&mut dose, &mu, &miss, 1.0, 0.5),
            0.0,
            epsilon = 1e-15
        );
        assert_relative_eq!(dose.sum(), 0.0, epsilon = 1e-15);
    }

    #[test]
    fn deposition_is_generic_over_scalar_f32() {
        let grid =
            VoxelGrid::<f32>::axis_aligned([9, 9, 9], [2.0, 2.0, 2.0], Point3::new(0.0, 0.0, 0.0))
                .unwrap();
        let mu = Volume::from_shape_fn(grid, |_| 0.05_f32);
        let mut dose = Volume::zeros(*mu.grid());
        let ray = Ray::try_new(
            Point3::new(-50.0_f32, 8.0, 8.0),
            Vector3::new(1.0_f32, 0.0, 0.0),
        )
        .unwrap();
        let total = deposit_ray_terma(&mut dose, &mu, &ray, 1.0_f32, 0.25);
        let expected = 1.0_f32 - (-0.05_f32 * 1.6).exp();
        assert_relative_eq!(total, expected, epsilon = 1e-6);
        assert_relative_eq!(dose.sum(), total, epsilon = 1e-5);
    }

    #[test]
    fn diverging_reduces_to_no_falloff_at_large_sad() {
        // As SAD → ∞ the inverse-square factor → 1 everywhere → the diverging
        // deposition matches the plain energy-conserving one.
        let mu = uniform_cube(0.05);
        let (mut plain_d, mut div_d) = (Volume::zeros(*mu.grid()), Volume::zeros(*mu.grid()));
        let plain = deposit_ray_terma(&mut plain_d, &mu, &central_x_ray(), 1.0, 0.25);
        let focal = Point3::new(-1.0e9, 8.0, 8.0);
        let div =
            deposit_ray_terma_diverging(&mut div_d, &mu, &central_x_ray(), 1.0, 0.25, focal, 1.0e9);
        assert_relative_eq!(div, plain, max_relative = 1e-6);
    }

    #[test]
    fn inverse_square_steepens_the_entry_to_exit_ratio() {
        // Point source 20 mm before the entry face (SAD = 28 mm to the centre):
        // the entry voxel (near source, isf > 1) gains dose relative to the exit
        // voxel (far, isf < 1) beyond the pure-attenuation ratio.
        let mu = uniform_cube(0.1);
        let mut plain_d = Volume::zeros(*mu.grid());
        let _ = deposit_ray_terma(&mut plain_d, &mu, &central_x_ray(), 1.0, 0.1);
        let mut div_d = Volume::zeros(*mu.grid());
        let focal = Point3::new(-20.0, 8.0, 8.0);
        let _ =
            deposit_ray_terma_diverging(&mut div_d, &mu, &central_x_ray(), 1.0, 0.1, focal, 28.0);

        let ratio_plain = plain_d.get(0, 4, 4).unwrap() / plain_d.get(8, 4, 4).unwrap();
        let ratio_div = div_d.get(0, 4, 4).unwrap() / div_d.get(8, 4, 4).unwrap();
        assert!(
            ratio_div > ratio_plain,
            "inverse-square should steepen entry/exit ratio: {ratio_div} !> {ratio_plain}"
        );
    }
}
