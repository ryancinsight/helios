//! Ray and axis-aligned bounding box with slab intersection.
//!
//! These are the Helios-owned geometry primitives layered on the leto
//! substrate: [`leto::Vector3`] supplies the vector algebra, and this module adds
//! the ray/AABB traversal used by imaging projectors (Siddon/Joseph forward
//! projection) and voxel-grid dose transport.

use eunomia::RealField;
use leto::Vector3;

/// A parametric ray `p(t) = origin + t · direction`, `t ≥ 0`.
///
/// `direction` need not be unit length; intersection parameters are expressed in
/// units of `direction`, so a unit direction yields parameters in world-space
/// distance.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Ray<T: RealField> {
    /// Ray origin.
    pub origin: Vector3<T>,
    /// Ray direction (not required to be normalized).
    pub direction: Vector3<T>,
}

impl<T: RealField> Ray<T> {
    /// Construct a ray from an origin and a direction.
    #[must_use]
    pub const fn new(origin: Vector3<T>, direction: Vector3<T>) -> Self {
        Self { origin, direction }
    }

    /// The point `origin + t · direction`.
    #[must_use]
    pub fn at(&self, t: T) -> Vector3<T> {
        Vector3::new(
            self.origin.data[0] + t * self.direction.data[0],
            self.origin.data[1] + t * self.direction.data[1],
            self.origin.data[2] + t * self.direction.data[2],
        )
    }
}

/// An axis-aligned bounding box defined by its minimum and maximum corners.
///
/// The invariant `min[i] ≤ max[i]` for each axis is the caller's responsibility;
/// [`Aabb::new`] does not reorder corners.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Aabb<T: RealField> {
    /// Minimum corner (smallest coordinate on each axis).
    pub min: Vector3<T>,
    /// Maximum corner (largest coordinate on each axis).
    pub max: Vector3<T>,
}

/// The entry/exit parameters of a ray–AABB intersection.
///
/// `t_near` is the parameter at which the ray enters the box (clamped to `0` when
/// the origin is inside), `t_far` the parameter at which it exits. Both are in
/// units of the ray's `direction`, and `t_near ≤ t_far`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RayHit<T: RealField> {
    /// Entry parameter (≥ 0).
    pub t_near: T,
    /// Exit parameter (≥ `t_near`).
    pub t_far: T,
}

impl<T: RealField> Aabb<T> {
    /// Construct an AABB from its minimum and maximum corners.
    #[must_use]
    pub const fn new(min: Vector3<T>, max: Vector3<T>) -> Self {
        Self { min, max }
    }

    /// Intersect a forward ray (`t ≥ 0`) with the box using the slab method.
    ///
    /// Returns the entry/exit parameters when the ray enters the box at some
    /// `t ≥ 0`, or `None` when it misses or the box lies entirely behind the
    /// origin. When the origin is inside the box, `t_near` is `0`.
    ///
    /// # Boundary contract
    /// A ray exactly parallel to a face *and* coincident with it (a
    /// measure-zero grazing case that can produce `0 · ∞`) has unspecified
    /// membership and may report either outcome; projector callers perturb such
    /// rays off-axis, so this does not affect dose/projection accuracy.
    #[must_use]
    pub fn intersect_ray(&self, ray: &Ray<T>) -> Option<RayHit<T>> {
        let mut t_near = T::ZERO;
        let mut t_far = T::infinity();

        for axis in 0..3 {
            let inv = ray.direction.data[axis].recip();
            let mut t0 = (self.min.data[axis] - ray.origin.data[axis]) * inv;
            let mut t1 = (self.max.data[axis] - ray.origin.data[axis]) * inv;
            if t1 < t0 {
                core::mem::swap(&mut t0, &mut t1);
            }
            t_near = t_near.max_scalar(t0);
            t_far = t_far.min_scalar(t1);
            if t_far < t_near {
                return None;
            }
        }

        Some(RayHit { t_near, t_far })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    fn unit_offset_box() -> Aabb<f64> {
        Aabb::new(Vector3::new(2.0, -1.0, -1.0), Vector3::new(4.0, 1.0, 1.0))
    }

    #[test]
    fn axis_aligned_ray_enters_and_exits_at_slab_faces() {
        let ray = Ray::new(Vector3::new(0.0, 0.0, 0.0), Vector3::new(1.0, 0.0, 0.0));
        let hit = unit_offset_box()
            .intersect_ray(&ray)
            .expect("ray along +x through the box hits");
        // Direction is unit, so parameters equal world-space distances: enters
        // the near x-face at x=2, exits the far x-face at x=4.
        assert_eq!(hit.t_near, 2.0);
        assert_eq!(hit.t_far, 4.0);
    }

    #[test]
    fn ray_pointing_away_misses() {
        let ray = Ray::new(Vector3::new(0.0, 0.0, 0.0), Vector3::new(-1.0, 0.0, 0.0));
        assert_eq!(unit_offset_box().intersect_ray(&ray), None);
    }

    #[test]
    fn parallel_ray_outside_slab_misses() {
        // Travels along +x but offset in y beyond the box's y-extent.
        let ray = Ray::new(Vector3::new(0.0, 5.0, 0.0), Vector3::new(1.0, 0.0, 0.0));
        assert_eq!(unit_offset_box().intersect_ray(&ray), None);
    }

    #[test]
    fn origin_inside_box_clamps_t_near_to_zero() {
        let ray = Ray::new(Vector3::new(3.0, 0.0, 0.0), Vector3::new(1.0, 0.0, 0.0));
        let hit = unit_offset_box()
            .intersect_ray(&ray)
            .expect("origin inside the box always hits");
        assert_eq!(hit.t_near, 0.0);
        assert_eq!(hit.t_far, 1.0); // exits at x=4 from x=3.
    }

    #[test]
    fn diagonal_ray_enters_at_near_corner_exits_at_far_corner() {
        // Unnormalized (1,1,1) direction: t is in units of that direction, so a
        // box [1,2]³ is entered at t=1 (corner (1,1,1)) and exited at t=2.
        let ray = Ray::new(Vector3::new(0.0, 0.0, 0.0), Vector3::new(1.0, 1.0, 1.0));
        let bbox = Aabb::new(Vector3::new(1.0, 1.0, 1.0), Vector3::new(2.0, 2.0, 2.0));
        let hit = bbox.intersect_ray(&ray).expect("diagonal ray hits");
        assert_eq!(hit.t_near, 1.0);
        assert_eq!(hit.t_far, 2.0);
        // `at(t_near)` must land on the near corner.
        let entry = ray.at(hit.t_near);
        assert_relative_eq!(entry.data[0], 1.0, epsilon = 1e-12);
        assert_relative_eq!(entry.data[1], 1.0, epsilon = 1e-12);
        assert_relative_eq!(entry.data[2], 1.0, epsilon = 1e-12);
    }

    #[test]
    fn intersection_is_generic_over_scalar_f32() {
        // Exercising the same code path at f32 confirms the kernel is genuinely
        // generic over the Scalar seam (native precision, no widen-narrow).
        let ray = Ray::new(
            Vector3::new(0.0_f32, 0.0, 0.0),
            Vector3::new(1.0_f32, 0.0, 0.0),
        );
        let bbox = Aabb::new(
            Vector3::new(2.0_f32, -1.0, -1.0),
            Vector3::new(4.0_f32, 1.0, 1.0),
        );
        let hit = bbox.intersect_ray(&ray).expect("f32 ray hits");
        assert_eq!(hit.t_near, 2.0_f32);
        assert_eq!(hit.t_far, 4.0_f32);
    }
}
