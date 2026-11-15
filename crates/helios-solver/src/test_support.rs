//! Shared fixtures for the solver kernel unit tests.
//!
//! The dose and scatter kernels operate on the same axis-aligned grid geometry,
//! so the test grid and the centimetre-length constructor live here once rather
//! than being copied (and drifting) in each `#[cfg(test)]` module.

use aequitas::systems::si::{quantities::Length, units::Millimeter};
use helios_domain::VoxelGrid;
use helios_math::{Point3, ShippedScalar, UnitQuaternion, Vector3};

/// A `dims`-shaped axis-aligned `f64` test grid: 2 mm spacing, origin at the
/// world origin.
pub(crate) fn grid(dims: [usize; 3]) -> VoxelGrid<f64> {
    VoxelGrid::axis_aligned(
        dims,
        [Length::from_unit::<Millimeter>(2.0); 3],
        Point3::new(0.0, 0.0, 0.0),
    )
    .expect("grid")
}

/// A `value`-centimetre `Length`, spelled in base SI metres (`0.01` m per cm).
pub(crate) fn length_cm<T: ShippedScalar>(value: T) -> Length<T> {
    Length::from_base(value * T::from_f64(0.01))
}

/// Acceptance tolerance for [`quarter_turn_rotation`]'s basis columns.
///
/// The columns `(0, 1, 0)`, `(-1, 0, 0)` and `(0, 0, 1)` have components in
/// `{-1, 0, 1}`, so every norm, dot product and triple product the
/// constructor evaluates is exact and every residual is exactly zero. The
/// constructor accepts any tolerance in `(0, 1)`; `EPSILON` is the unit
/// roundoff scale of the unit-length columns, the smallest tolerance that
/// still admits a basis carrying rounding in its last place.
pub(crate) const QUARTER_TURN_BASIS_TOLERANCE: f64 = f64::EPSILON;

/// Rotation taking local `+x` onto world `+y` (a quarter turn about `z`).
pub(crate) fn quarter_turn_rotation() -> UnitQuaternion<f64> {
    UnitQuaternion::try_from_rotation_columns(
        Vector3::new(0.0, 1.0, 0.0),
        Vector3::new(-1.0, 0.0, 0.0),
        Vector3::new(0.0, 0.0, 1.0),
        QUARTER_TURN_BASIS_TOLERANCE,
    )
    .expect("right-handed quarter-turn basis")
}

/// Roundings per slab endpoint of a ray-grid chord that are independent of
/// the other endpoint, in units of the unit roundoff `u = EPSILON / 2`.
///
/// `ray_grid_interval` forms `length = (t_exit - t_enter) * speed_recip` with
/// `t = (bound - o) * recip(d)`. The first-order worst case (every rounding at
/// its maximum, all signs aligned) holds when the entry and the exit lie on
/// one slab axis, where the origin error cancels in their difference. A
/// clamped entry is exact `t = 0` on no slab, so the origin error enters
/// `t_exit` alone and adds `|delta o| / (bound - o)`. Each endpoint rounds
/// `bound - o`, the product with `recip(d)` and the product with
/// `speed_recip`: 3. Each is a relative error of one `t`, so it reaches the
/// length amplified by the conditioning `c = (|t_enter| + |t_exit|) / length`.
pub(crate) const SLAB_ENDPOINT_ROUNDINGS: f64 = 3.0;

/// Roundings of the ray-grid chord length that scale both endpoints by one
/// factor, so conditioning does not amplify them, in units of `u`.
///
/// - `d * recip(spacing)` for the on-axis component: 2.
/// - The speed `|d|` and its reciprocal: the norm divides by the largest
///   component, accumulates three fused multiply-adds, takes the root and
///   rescales, `(2 + 3) / 2 + 1 + 1 = 4.5`, plus 1.
/// - `Ray::try_new` rescales by the largest component (1 while the chord
///   axis carries it, as in every fixture; up to 2 otherwise) and
///   `new_normalize` takes the norm (4.5), its reciprocal (1) and one
///   product (1).
/// - `recip(d)` of the unit component in `intersect_aabb`: 1.
///
/// The total is 16.
pub(crate) const RAY_DIRECTION_ROUNDINGS: f64 = 16.0;

/// Relative error the quarter-turn pose adds to the rotated direction, in
/// units of `u`: one `transform_vector` of at most 10 roundings per
/// component, bounded by inspection of its Rodrigues form, not measured.
///
/// On one slab axis the rotated origin enters both slab distances through the
/// same value `o`, so its error shifts `t_enter` and `t_exit` equally and
/// cancels in their difference; only the roundings of `bound - o` remain,
/// counted in [`SLAB_ENDPOINT_ROUNDINGS`]. The rotated direction scales both
/// distances by one factor, so its error is added to the length unamplified.
pub(crate) const ORIENTED_FRAME_ROUNDINGS: f64 = 10.0;
