//! Geometry predicates evaluated over the Atlas geometry primitives.
//!
//! One shared definition per predicate, so callers do not each hand-roll the
//! same field (and drift). The `Aabb`/`Point3`/`Ray` types themselves are owned
//! by gaia and leto and only re-exported by this crate.

use crate::{Aabb, GeometryScalar, NumericElement, Point3};

/// Signed distance from `p` to the surface of the axis-aligned box `aabb` (the
/// standard AABB signed-distance field): negative inside, exactly `0` on the
/// boundary, positive outside, in the box's own length units.
///
/// This is the single shared evaluation of the AABB SDF used by collimator
/// penumbra and any other geometry predicate that needs "how far inside/outside"
/// rather than the boolean [`Aabb::contains_point`].
#[must_use]
pub fn aabb_signed_distance<T: GeometryScalar>(aabb: &Aabb<T>, p: &Point3<T>) -> T {
    let zero = <T as NumericElement>::ZERO;
    let centre = aabb.center();
    let half = <T as GeometryScalar>::from_f64(0.5);
    // q_i = |p_i − c_i| − h_i : >0 outside that axis's slab, <0 inside.
    let q = [
        (p.x - centre.x).abs() - (aabb.max.x - aabb.min.x) * half,
        (p.y - centre.y).abs() - (aabb.max.y - aabb.min.y) * half,
        (p.z - centre.z).abs() - (aabb.max.z - aabb.min.z) * half,
    ];
    let outside_sq = q[0].max_scalar(zero) * q[0].max_scalar(zero)
        + q[1].max_scalar(zero) * q[1].max_scalar(zero)
        + q[2].max_scalar(zero) * q[2].max_scalar(zero);
    // Inside distance: the least-negative q (nearest face), clamped ≤ 0.
    let inside = q[0].max_scalar(q[1]).max_scalar(q[2]).min_scalar(zero);
    outside_sq.sqrt() + inside
}
