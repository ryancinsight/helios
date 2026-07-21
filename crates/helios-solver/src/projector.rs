//! MVCT forward projector: ray-marched line integral of attenuation.
//!
//! Computes the optical depth `τ = ∫ μ dl` of a [`Ray`] through a linear-
//! attenuation [`Volume`] — the core of MVCT forward projection and dose
//! ray-tracing. Geometry comes from the gaia kernel ([`Ray`]/[`Aabb`], consumed
//! via `helios-math`): the ray is clipped to the grid's world bounding box, then
//! marched with midpoint sampling of the trilinearly-interpolated `μ` field.
//!
//! Midpoint integration is exact for `μ` fields affine along the ray (the
//! analytical oracle in the tests); for a homogeneous slab it reproduces
//! `τ = μ·L` exactly.
//!
//! The ray is transformed into the grid's scaled-index frame for clipping, then
//! sampled along its original world-space parameter. Consequently the same
//! integration contract holds for every rigid [`VoxelGrid`] pose, including
//! DICOM orientation-cosine grids.

use helios_core::constants::MM_PER_CM;
use helios_domain::{Volume, VoxelGrid};
use helios_math::{Aabb, GeometryScalar, Point3, Ray, Vector3};

/// Intersect a world-space unit ray with the grid's node-centre box.
///
/// Clipping occurs in continuous index coordinates, where the node-centre box
/// is axis-aligned. `local_speed` converts the normalized local-ray parameter
/// back to the original world-space millimetre parameter, so midpoint samples
/// preserve the physical path-length convention of [`forward_project_ray`].
pub(crate) fn ray_grid_interval<T: GeometryScalar>(
    grid: &VoxelGrid<T>,
    ray: &Ray<T>,
) -> Option<(T, T)> {
    let [nx, ny, nz] = grid.dims();
    let local_origin = grid.world_to_index(ray.origin());
    let local_direction_mm = grid.pose().inverse().transform_vector(ray.direction());
    let spacing = grid.spacing();
    let local_direction = Vector3::new(
        local_direction_mm.x * spacing[0].recip(),
        local_direction_mm.y * spacing[1].recip(),
        local_direction_mm.z * spacing[2].recip(),
    );
    let local_speed = local_direction.norm();
    if !local_speed.is_finite() || local_speed <= T::ZERO {
        return None;
    }
    let local_ray = Ray::try_new(local_origin, local_direction).ok()?;
    let local_aabb = Aabb::new(
        Point3::new(T::ZERO, T::ZERO, T::ZERO),
        Point3::new(
            <T as GeometryScalar>::from_f64((nx - 1) as f64),
            <T as GeometryScalar>::from_f64((ny - 1) as f64),
            <T as GeometryScalar>::from_f64((nz - 1) as f64),
        ),
    );
    let (enter, exit) = local_ray.intersect_aabb(&local_aabb)?;
    Some((enter * local_speed.recip(), exit * local_speed.recip()))
}

/// Ray-march the optical depth `τ = ∫ μ dl` of `ray` through the `mu` volume.
///
/// The `mu` volume holds the linear attenuation coefficient in **cm⁻¹** (physics
/// convention, matching `helios_physics::LinearAttenuation`), while the grid /
/// `ray` are in **mm** (DICOM convention). The path length is converted mm→cm so
/// `τ` is a true dimensionless optical depth.
///
/// `step_mm` is the nominal sampling step (the actual step is `L / ceil(L/step)`
/// so it divides the traversed length exactly). Returns `None` if the ray misses
/// the grid.
#[must_use]
pub fn forward_project_ray<T: GeometryScalar>(
    mu: &Volume<T>,
    ray: &Ray<T>,
    step_mm: T,
) -> Option<T> {
    let grid = *mu.grid();
    let (t_enter, t_exit) = ray_grid_interval(&grid, ray)?;

    let length = t_exit - t_enter;
    if length <= T::ZERO {
        return Some(T::ZERO);
    }
    // Number of substeps so the step divides the length exactly (>= 1).
    let n_f = (length * step_mm.recip()).ceil();
    let n = (n_f.to_f64() as usize).max(1);
    let actual_step = length * <T as GeometryScalar>::from_f64(n as f64).recip();
    // Segment length in cm so cm⁻¹ · cm is dimensionless.
    let step_cm = actual_step * <T as GeometryScalar>::from_f64(MM_PER_CM).recip();
    let half = <T as GeometryScalar>::from_f64(0.5);

    let mut tau = T::ZERO;
    for i in 0..n {
        let t_mid = t_enter + (<T as GeometryScalar>::from_f64(i as f64) + half) * actual_step;
        let world_pt: Point3<T> = ray.point_at(t_mid);
        let index = grid.world_to_index(world_pt);
        let mu_sample = mu.sample_trilinear(index).unwrap_or(T::ZERO);
        tau += mu_sample * step_cm;
    }
    Some(tau)
}

#[cfg(test)]
mod tests {
    use super::*;
    use eunomia::assert_relative_eq;
    use helios_math::{Ray, Vector3};

    fn axis_grid() -> VoxelGrid<f64> {
        // Node box: x∈[0,20], y∈[0,4], z∈[0,4]; 2 mm spacing.
        VoxelGrid::axis_aligned([11, 3, 3], [2.0, 2.0, 2.0], Point3::new(0.0, 0.0, 0.0))
            .expect("grid")
    }

    fn ray_along_x(origin_x: f64) -> Ray<f64> {
        Ray::try_new(Point3::new(origin_x, 2.0, 2.0), Vector3::new(1.0, 0.0, 0.0))
            .expect("unit +x ray")
    }

    fn oriented_grid() -> VoxelGrid<f64> {
        let rotation = helios_math::UnitQuaternion::try_from_rotation_columns(
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(-1.0, 0.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
            1.0e-12,
        )
        .expect("right-handed quarter-turn basis");
        VoxelGrid::oriented(
            [11, 3, 3],
            [2.0, 2.0, 2.0],
            Point3::new(10.0, 20.0, 30.0),
            rotation,
        )
        .expect("grid")
    }

    #[test]
    fn homogeneous_slab_gives_mu_times_length() {
        // Uniform μ=0.06 cm⁻¹ over a 20 mm = 2 cm crossing → τ = 0.06·2 = 0.12.
        let mu = Volume::from_shape_fn(axis_grid(), |_| 0.06);
        let tau = forward_project_ray(&mu, &ray_along_x(-5.0), 0.5).expect("hit");
        assert_relative_eq!(tau, 0.06 * 2.0, epsilon = 1e-12);
    }

    #[test]
    fn linear_field_is_integrated_exactly_by_midpoint() {
        // μ(index i) = 0.01·i + 0.02 cm⁻¹ → along x, index i = x_mm/2, so
        // μ(x) = 0.005·x_mm + 0.02. τ = ∫₀²⁰ μ dx_cm = 0.1·∫₀²⁰(0.005x+0.02)dx_mm
        //      = 0.1·(0.005·200 + 0.02·20) = 0.1·1.4 = 0.14. Midpoint is exact for
        // an affine field.
        let mu = Volume::from_shape_fn(axis_grid(), |idx| 0.01 * idx[0] as f64 + 0.02);
        let tau = forward_project_ray(&mu, &ray_along_x(-5.0), 1.0).expect("hit");
        assert_relative_eq!(tau, 0.14, epsilon = 1e-10);
    }

    #[test]
    fn ray_missing_the_grid_returns_none() {
        let mu = Volume::from_shape_fn(axis_grid(), |_| 0.06);
        // Offset in y beyond [0,4].
        let miss =
            Ray::try_new(Point3::new(-5.0, 100.0, 2.0), Vector3::new(1.0, 0.0, 0.0)).unwrap();
        assert_eq!(forward_project_ray(&mu, &miss, 0.5), None);
    }

    #[test]
    fn step_size_does_not_change_homogeneous_result() {
        let mu = Volume::from_shape_fn(axis_grid(), |_| 0.06);
        let coarse = forward_project_ray(&mu, &ray_along_x(-5.0), 5.0).expect("hit");
        let fine = forward_project_ray(&mu, &ray_along_x(-5.0), 0.1).expect("hit");
        assert_relative_eq!(coarse, fine, epsilon = 1e-12);
    }

    #[test]
    fn oriented_homogeneous_slab_gives_mu_times_length() {
        // The local +x node span is 20 mm. The grid rotates it onto world +y,
        // so this world-space ray must still integrate τ = 0.06 cm⁻¹ · 2 cm.
        let mu = Volume::from_shape_fn(oriented_grid(), |_| 0.06);
        let ray = Ray::try_new(Point3::new(10.0, 15.0, 30.0), Vector3::new(0.0, 1.0, 0.0))
            .expect("unit +y ray");
        let tau = forward_project_ray(&mu, &ray, 0.5).expect("hit");
        assert_relative_eq!(tau, 0.12, epsilon = 1e-12);
    }

    #[test]
    fn projector_is_generic_over_scalar_f32() {
        let g =
            VoxelGrid::<f32>::axis_aligned([11, 3, 3], [2.0, 2.0, 2.0], Point3::new(0.0, 0.0, 0.0))
                .unwrap();
        let mu = Volume::from_shape_fn(g, |_| 0.06_f32);
        let ray = Ray::try_new(
            Point3::new(-5.0_f32, 2.0, 2.0),
            Vector3::new(1.0_f32, 0.0, 0.0),
        )
        .unwrap();
        let tau = forward_project_ray(&mu, &ray, 0.5_f32).expect("hit");
        assert_relative_eq!(tau, 0.06_f32 * 2.0, epsilon = 1e-5);
    }
}
