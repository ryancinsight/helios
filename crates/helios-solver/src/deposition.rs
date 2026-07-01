//! Primary-fluence energy deposition (TERMA) along a beam ray.
//!
//! Companion to the [`forward_project_ray`](crate::forward_project_ray) line
//! integral: where the projector reduces `∫ μ dl` to a single optical depth, this
//! kernel *deposits* the energy the primary beam loses as it attenuates, voxel by
//! voxel, producing the terma (total energy released per unit mass) that a
//! collapsed-cone/convolution dose engine spreads with a scatter kernel.

use crate::projector::world_aabb;
use helios_domain::Volume;
use helios_math::{GeometryScalar, NumericElement, Point3, Ray};

/// Millimetres per centimetre — the world grid is in mm, `μ` in cm⁻¹.
const MM_PER_CM: f64 = 10.0;

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
    let grid = *mu.grid();
    debug_assert_eq!(
        grid.dims(),
        dose.grid().dims(),
        "dose and mu must share the same grid"
    );
    let aabb = world_aabb(&grid);
    let (t_enter, t_exit) = match ray.intersect_aabb(&aabb) {
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
        let absorbed = weight * (trans_before - trans_after);
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
        Ray::try_from_direction(Point3::new(-50.0, 8.0, 8.0), Vector3::new(1.0, 0.0, 0.0))
            .expect("unit +x ray")
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
            Ray::try_from_direction(Point3::new(-50.0, 500.0, 8.0), Vector3::new(1.0, 0.0, 0.0))
                .unwrap();
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
        let ray = Ray::try_from_direction(
            Point3::new(-50.0_f32, 8.0, 8.0),
            Vector3::new(1.0_f32, 0.0, 0.0),
        )
        .unwrap();
        let total = deposit_ray_terma(&mut dose, &mu, &ray, 1.0_f32, 0.25);
        let expected = 1.0_f32 - (-0.05_f32 * 1.6).exp();
        assert_relative_eq!(total, expected, epsilon = 1e-6);
        assert_relative_eq!(dose.sum(), total, epsilon = 1e-5);
    }
}
