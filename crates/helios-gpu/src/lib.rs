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
//!
//! # Running the tests
//!
//! Every test that acquires a device is `#[ignore]`d: a host with no wgpu
//! adapter reports them as *skipped* instead of passing a test that dispatched
//! nothing. On a machine with an adapter, opt in with
//! `cargo nextest run -p helios-gpu --run-ignored all`; a missing adapter is
//! then a hard failure.
#![forbid(unsafe_code)]
#![deny(missing_docs)]

mod attenuation;
mod projection;
mod transmission;

pub use hephaestus_core::{HephaestusError, Result};
pub use hephaestus_wgpu::WgpuDevice;

pub use attenuation::GpuAttenuationMapper;
pub use projection::GpuProjector;
pub use transmission::beam_transmission_into;

/// Acquire the default wgpu compute device (highest-power adapter available).
///
/// # Errors
/// Returns [`HephaestusError`] if no compatible GPU adapter/device is available.
pub fn default_device() -> Result<WgpuDevice> {
    WgpuDevice::try_default("helios-gpu")
}

/// Device accessor for the adapter-requiring tests.
///
/// Those tests are `#[ignore]`d (see the module tests), so reaching this helper
/// means the run explicitly opted in with `--run-ignored all`. A missing adapter
/// is therefore a hard failure: the alternative — returning early — reports a
/// test that executed no kernel as *passed*, which is a false green rather than
/// an honest skip.
#[cfg(test)]
fn required_device() -> WgpuDevice {
    default_device().expect(
        "GPU tests are opt-in (`cargo nextest run -p helios-gpu --run-ignored all`) and require a \
         wgpu adapter; none was found",
    )
}
