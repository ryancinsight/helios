//! Helios numeric and geometry-vocabulary layer.
//!
//! This crate is a thin, additive layer over the Atlas numeric SSOT
//! ([`eunomia`]) and the array/linear-algebra substrate ([`leto`]). It does not
//! reinvent scalars, vectors, or geometry primitives; it re-exports them as the
//! Helios vocabulary so higher layers depend on one source.
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
//! The linear-algebra substrate — vectors, points, and rigid transforms — comes
//! from [`leto`] ([`Vector3`], [`Point3`], [`Isometry3`], …) and is re-exported
//! here. Higher-level geometry *primitives* (axis-aligned boxes, rays and their
//! intersection, meshes, CSG solids) are owned by the **gaia** geometry kernel,
//! not by Helios: `gaia::Aabb`/`gaia::Ray` are the canonical types, consumed once
//! gaia's leto-based geometry lands on its default branch (see `gap_audit.md`
//! G-11 / `backlog.md` H-003b). Helios never re-implements them.
#![forbid(unsafe_code)]
#![deny(missing_docs)]

/// The Helios scalar seam: the real-field trait every compute kernel is generic
/// over. Identical to [`eunomia::RealField`]; re-exported here as the canonical
/// Helios name.
pub use eunomia::RealField as Scalar;

// Re-export the rest of the numeric SSOT surface so downstream crates depend on
// one vocabulary source.
pub use eunomia::{CastFrom, CastTo, FloatElement, NumericElement};

/// Linear-algebra substrate from leto, re-exported as the Helios vocabulary.
/// (Geometry *primitives* — `Aabb`, `Ray`, meshes — are owned by gaia.)
pub use leto::geometry::{UnitVector2, UnitVector3};
pub use leto::{
    Isometry3, Point3, Quaternion, Translation3, Unit, UnitQuaternion, Vector2, Vector3,
};
