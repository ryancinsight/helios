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
