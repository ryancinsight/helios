//! Shared fixtures for the solver kernel unit tests.
//!
//! The dose and scatter kernels operate on the same axis-aligned grid geometry,
//! so the test grid and the centimetre-length constructor live here once rather
//! than being copied (and drifting) in each `#[cfg(test)]` module.

use aequitas::systems::si::quantities::Length;
use helios_domain::VoxelGrid;
use helios_math::{Point3, ShippedScalar};

/// A `dims`-shaped axis-aligned `f64` test grid: 2 mm spacing, origin at the
/// world origin.
pub(crate) fn grid(dims: [usize; 3]) -> VoxelGrid<f64> {
    VoxelGrid::axis_aligned(dims, [2.0, 2.0, 2.0], Point3::new(0.0, 0.0, 0.0)).expect("grid")
}

/// A `value`-centimetre `Length`, spelled in base SI metres (`0.01` m per cm).
pub(crate) fn length_cm<T: ShippedScalar>(value: T) -> Length<T> {
    Length::from_base(value * T::from_f64(0.01))
}
