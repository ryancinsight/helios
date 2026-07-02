//! GPU Beer–Lambert transmission: `exp(−τ)` over a projection buffer.
//!
//! Converts a buffer of optical depths `τ = ∫μ dl` (from the forward projector)
//! into MVCT detector transmission `exp(−τ)` on the GPU with hephaestus-wgpu's
//! fused [`ExpNegOp`] kernel — one dispatch, no intermediate device buffer
//! (previously a `NegOp` → `ExpOp` chain; the fused op was upstreamed for this
//! path, H-043b). Differentially validated against the CPU `f32::exp` reference.

use hephaestus_core::{BlockWidth, ComputeDevice, Result};
use hephaestus_wgpu::{unary_elementwise_strided, ExpNegOp, StridedOperand, WgpuDevice};
use leto::Layout;

/// Compute `out[i] = exp(-optical_depth[i])` on the GPU.
///
/// `out` must have the same length as `optical_depth`.
///
/// # Errors
/// Returns [`HephaestusError`](crate::HephaestusError) on device transfer/dispatch
/// failure or if `out.len()` does not match `optical_depth.len()`.
pub fn beam_transmission_into(
    device: &WgpuDevice,
    optical_depth: &[f32],
    out: &mut [f32],
) -> Result<()> {
    if optical_depth.is_empty() && out.is_empty() {
        return Ok(());
    }
    let input = device.upload(optical_depth)?;
    // A rank-1 contiguous layout over `n` elements is always valid.
    let layout = Layout::c_contiguous([optical_depth.len()])
        .expect("invariant: rank-1 contiguous layout of a slice length is valid");

    let transmitted = unary_elementwise_strided::<ExpNegOp, f32, 1>(
        device,
        StridedOperand {
            buffer: &input,
            layout: &layout,
        },
        [optical_depth.len()],
        BlockWidth::DEFAULT,
    )?;
    device.download(&transmitted, out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gpu_transmission_matches_cpu_exp() {
        // Requires a GPU adapter; skip cleanly if none is present.
        let Ok(device) = crate::default_device() else {
            eprintln!("no GPU adapter — skipping GPU transmission test");
            return;
        };
        let tau = [0.0_f32, 0.1, 0.5, 1.0, 2.0, 3.5, 7.0];
        let mut got = [0.0_f32; 7];
        beam_transmission_into(&device, &tau, &mut got).expect("gpu transmission");
        for (&t, &g) in tau.iter().zip(got.iter()) {
            let cpu = (-t).exp();
            // f32 GPU exp vs CPU exp: agree to a few ULP → relative bound 1e-5.
            assert!(
                (g - cpu).abs() <= 1e-5 * (1.0 + cpu.abs()),
                "τ={t}: gpu={g} vs cpu={cpu}"
            );
        }
    }

    #[test]
    fn empty_input_is_ok() {
        let Ok(device) = crate::default_device() else {
            return;
        };
        let mut out: [f32; 0] = [];
        beam_transmission_into(&device, &[], &mut out).expect("empty is a no-op");
    }
}
