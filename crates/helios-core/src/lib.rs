//! Helios foundation layer.
//!
//! `helios-core` is the innermost crate of the Helios workspace: it depends on
//! nothing project-local and everything else depends inward on it. It owns the
//! cross-cutting vocabulary shared by every higher layer — the typed error
//! surface, physical constants for radiation transport and dosimetry, the
//! validating domain newtypes that make invalid states unrepresentable at the
//! parse/deserialize boundary, and the memory seam ([`memory`]) that names the
//! workspace allocator once for everyone else.
//!
//! Domain compute is deliberately absent here: the generic numeric seam
//! (`Scalar`, backed by `hermes`/`leto`) lands in `helios-math`, and physics,
//! geometry, GPU, and orchestration layers build outward from that. Constants in
//! this crate are `f64` literals at their definition boundary; higher layers
//! convert into `T: Scalar` at call sites.
#![forbid(unsafe_code)]
#![deny(missing_docs)]

pub mod constants;
mod error;
/// The workspace memory seam: the one place Helios names its allocator.
#[cfg(feature = "mnemosyne-memory")]
pub mod memory;
mod units;

pub use error::{HeliosError, Result};
pub use units::{EnergyMeV, HounsfieldUnit, VoxelSpacingMm};
