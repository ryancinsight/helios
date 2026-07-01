//! Helios radiation-imaging layer.
//!
//! MVCT acquisition modeling and reconstruction. This module provides the
//! parallel-beam **Radon forward transform** ([`parallel_beam_radon`]) — the
//! sinogram of line integrals `p(θ, s) = ∫ μ dl` over projection angles `θ` and
//! signed detector offsets `s` — built on the `helios-solver` ray-march projector.
//! Reconstruction (filtered back-projection) consumes the same [`Sinogram`].
//!
//! Operates on an axis-aligned [`Volume`](helios_domain::Volume) slice; all types
//! are generic over the geometry scalar seam.
#![forbid(unsafe_code)]
#![deny(missing_docs)]

mod fbp;
mod noise;
mod radon;

pub use fbp::filtered_back_projection;
pub use noise::add_quantum_noise;
pub use radon::{parallel_beam_radon, Sinogram};
