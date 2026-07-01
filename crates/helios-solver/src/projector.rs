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
//! # Scope
//! This first projector supports **axis-aligned** voxel grids (identity pose
//! rotation) — the common CT/MVCT patient-frame layout. Oriented grids (non-
//! trivial `ImageOrientationPatient` cosines) return [`None`] rather than a
//! silently-wrong result; general-pose clipping is a tracked follow-up.

use helios_domain::{Volume, VoxelGrid};
use helios_math::{Aabb, GeometryScalar, Point3, Ray};

/// World bounding box of a grid's sample region (node centres). `VoxelGrid` is
/// axis-aligned, so the min/max corners are the `(0,0,0)` and `(nx-1,ny-1,nz-1)`
/// voxel centres.
pub(crate) fn world_aabb<T: GeometryScalar>(grid: &VoxelGrid<T>) -> Aabb<T> {
    let [nx, ny, nz] = grid.dims();
    let min = grid.voxel_center(0, 0, 0);
    let max = grid.voxel_center(nx - 1, ny - 1, nz - 1);
    Aabb::new(min, max)
}

/// Millimetres per centimetre — the world grid is in mm, `μ` in cm⁻¹.
const MM_PER_CM: f64 = 10.0;

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
    let aabb = world_aabb(&grid);
    let (t_enter, t_exit) = ray.intersect_aabb(&aabb)?;

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
    use approx::assert_relative_eq;
    use helios_math::{Ray, Vector3};

    fn axis_grid() -> VoxelGrid<f64> {
        // Node box: x∈[0,20], y∈[0,4], z∈[0,4]; 2 mm spacing.
        VoxelGrid::axis_aligned([11, 3, 3], [2.0, 2.0, 2.0], Point3::new(0.0, 0.0, 0.0))
            .expect("grid")
    }

    fn ray_along_x(origin_x: f64) -> Ray<f64> {
        Ray::try_from_direction(Point3::new(origin_x, 2.0, 2.0), Vector3::new(1.0, 0.0, 0.0))
            .expect("unit +x ray")
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
            Ray::try_from_direction(Point3::new(-5.0, 100.0, 2.0), Vector3::new(1.0, 0.0, 0.0))
                .unwrap();
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
    fn projector_is_generic_over_scalar_f32() {
        let g =
            VoxelGrid::<f32>::axis_aligned([11, 3, 3], [2.0, 2.0, 2.0], Point3::new(0.0, 0.0, 0.0))
                .unwrap();
        let mu = Volume::from_shape_fn(g, |_| 0.06_f32);
        let ray = Ray::try_from_direction(
            Point3::new(-5.0_f32, 2.0, 2.0),
            Vector3::new(1.0_f32, 0.0, 0.0),
        )
        .unwrap();
        let tau = forward_project_ray(&mu, &ray, 0.5_f32).expect("hit");
        assert_relative_eq!(tau, 0.06_f32 * 2.0, epsilon = 1e-5);
    }
}
