//! Helios computational engines.
//!
//! Deterministic kernels over voxel volumes: material-property maps (this
//! module's [`attenuation_map`]), and — as they land — dose engines and imaging
//! projectors. Each engine is authored as a CPU reference first; the GPU-
//! accelerated path (`helios-gpu`, on the `hephaestus` device seam) is
//! differentially validated against that reference.
//!
//! All kernels are generic over the [`Scalar`](helios_math::Scalar) seam.
#![forbid(unsafe_code)]
#![deny(missing_docs)]

mod attenuation_map;

pub use attenuation_map::attenuation_map;
