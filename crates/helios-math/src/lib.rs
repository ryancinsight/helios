//! Helios numeric and geometry layer.
//!
//! This crate is a thin, additive layer over the Atlas numeric SSOT
//! ([`eunomia`]) and array/geometry substrate ([`leto`]). It does not reinvent
//! scalars or vectors; it re-exports them as the Helios vocabulary and adds only
//! the primitives Helios needs that do not exist upstream.
//!
//! ## The `Scalar` seam
//!
//! [`Scalar`] is exactly [`eunomia::RealField`] — the single real-field trait
//! that all Helios compute is generic over (`f32`, `f64`, and the reduced-
//! precision `eunomia` types implement it natively). It is re-exported under the
//! canonical Helios name so higher layers write `T: Scalar`; there is no parallel
//! trait, only one documented name for the upstream SSOT.
//!
//! ## Geometry
//!
//! Vectors, points, and rigid transforms come from [`leto`] ([`Vector3`],
//! [`Point3`], [`Isometry3`], …). Helios adds [`Ray`] and [`Aabb`] with a robust
//! slab intersection ([`Aabb::intersect_ray`]) — the ray/voxel-grid traversal
//! primitive used by the imaging projectors and dose engines, which leto (an
//! array library) does not provide.
#![forbid(unsafe_code)]
#![deny(missing_docs)]

mod geometry;

pub use geometry::{Aabb, Ray, RayHit};

/// The Helios scalar seam: the real-field trait every compute kernel is generic
/// over. Identical to [`eunomia::RealField`]; re-exported here as the canonical
/// Helios name.
pub use eunomia::RealField as Scalar;

// Re-export the rest of the numeric SSOT surface so downstream crates depend on
// one vocabulary source.
pub use eunomia::{CastFrom, CastTo, FloatElement, NumericElement};

/// Geometry primitives from the leto substrate, re-exported as the Helios
/// geometry vocabulary.
pub use leto::{Isometry3, Point3, Quaternion, Translation3, UnitQuaternion, Vector2, Vector3};
