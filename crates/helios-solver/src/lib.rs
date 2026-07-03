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
mod deposition;
mod dose;
mod oriented_scatter;
mod projector;
mod scatter;

pub use attenuation_map::attenuation_map;
pub use deposition::{deposit_ray_terma, deposit_ray_terma_diverging};
pub use dose::{dose_convolution_x, exponential_deposition_kernel, primary_fluence_parallel_x};
pub use oriented_scatter::{directional_convolve, oriented_forward_scatter};
pub use projector::forward_project_ray;
pub use scatter::{
    anisotropic_scatter_superposition, forward_peaked_kernel, poly_forward_peaked_kernel,
    scatter_superposition, symmetric_deposition_kernel, SpectralComponent,
};
