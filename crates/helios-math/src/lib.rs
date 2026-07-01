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
pub use leto::geometry::{UnitVector2, UnitVector3};
pub use leto::{
    Isometry3, Point3, Quaternion, Translation3, Unit, UnitQuaternion, Vector2, Vector3,
};

/// Geometry *primitives* from the gaia geometry kernel, re-exported as the Helios
/// geometry vocabulary (upstream ownership — Helios never re-implements these).
/// `gaia::Ray` carries a validated unit direction; `Ray::intersect_aabb` provides
/// the ray/voxel traversal the imaging projectors and dose ray-trace use.
pub use gaia::{Aabb, Ray};

#[cfg(test)]
mod gaia_geometry_bridge_tests {
    //! Proves Helios consumes gaia's migrated (leto/eunomia) geometry through the
    //! synchronized-checkout wiring (H-050): a gaia `Ray` intersects a gaia `Aabb`
    //! built from leto `Point3`/`Vector3`.
    use super::{Aabb, Point3, Ray, Vector3};

    #[test]
    fn gaia_ray_intersects_gaia_aabb_through_helios() {
        let ray =
            Ray::try_from_direction(Point3::new(-2.0_f64, 0.5, 0.5), Vector3::new(1.0, 0.0, 0.0))
                .expect("non-zero direction");
        let aabb = Aabb::new(Point3::new(0.0_f64, 0.0, 0.0), Point3::new(1.0, 1.0, 1.0));
        // Unit +x ray from x=-2 enters the unit box at t=2, exits at t=3.
        let (t_enter, t_exit) = ray.intersect_aabb(&aabb).expect("ray hits the box");
        assert!((t_enter - 2.0).abs() < 1e-12, "t_enter = {t_enter}");
        assert!((t_exit - 3.0).abs() < 1e-12, "t_exit = {t_exit}");
    }
}
