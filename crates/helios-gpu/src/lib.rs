//! Helios GPU compute layer.
//!
//! Dispatches Helios kernels onto the Atlas accelerator substrate — the
//! [`hephaestus_core::ComputeDevice`] seam and its `hephaestus-wgpu` backend.
//! Helios programs against that seam directly (it does not reinvent a device
//! abstraction); this crate adds the Helios-specific GPU operations and a device
//! accessor, keeping GPU dependencies out of the pure domain/physics layers.
//!
//! GPU buffers are `f32` (the wgpu compute precision); callers stage `f32` data
//! at this boundary. Every GPU kernel here has a CPU reference it is
//! differentially validated against (see the tests), per the differential-
//! verification discipline.
#![forbid(unsafe_code)]
#![deny(missing_docs)]

mod projection;
mod transmission;

pub use hephaestus_core::{HephaestusError, Result};
pub use hephaestus_wgpu::WgpuDevice;

pub use projection::GpuProjector;
pub use transmission::beam_transmission_into;

/// Acquire the default wgpu compute device (highest-power adapter available).
///
/// # Errors
/// Returns [`HephaestusError`] if no compatible GPU adapter/device is available.
pub fn default_device() -> Result<WgpuDevice> {
    WgpuDevice::try_default("helios-gpu")
}
