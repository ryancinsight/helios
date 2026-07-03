//! GPU HU→μ attenuation map: fused affine-clamp over a CT HU buffer (H-010b).
//!
//! Implements `μ = max(scale·HU + offset, 0)` in one dispatch, where
//! `scale = (μ/ρ)·ρ_water/1000` and `offset = (μ/ρ)·ρ_water` — the GPU form of
//! [`helios_solver::attenuation_map`]'s per-voxel math
//! (`μ = (μ/ρ)·ρ_water·max(1 + HU/1000, 0)`, the Compton-dominated MV
//! approximation). Authored as a consumer-side kernel over the hephaestus
//! ADR-0004 seam ([`KernelInterface`] + [`KernelSource<Wgsl>`]) — no substrate
//! changes were needed to add it, which is the seam's acceptance criterion.
//!
//! Differentially validated against the CPU closed form and against
//! `helios_solver::attenuation_map` (see tests). Inputs are CT Hounsfield
//! units, finite by construction upstream; non-finite HU is outside the
//! contract (the clamp uses WGSL `max`, whose NaN ordering is
//! implementation-defined).

use bytemuck::{Pod, Zeroable};
use hephaestus_core::{
    Binding, BindingDecl, DispatchGrid, HephaestusError, KernelDevice, KernelInterface,
    KernelSource, Result, Wgsl,
};
use std::borrow::Cow;

/// Uniform parameter block of the fused affine-clamp kernel.
///
/// `_pad` keeps the block a 16-byte multiple (WGSL uniform layout).
#[repr(C)]
#[derive(Clone, Copy, Debug, Pod, Zeroable)]
struct AffineClampParams {
    scale: f32,
    offset: f32,
    len: u32,
    _pad: u32,
}

/// The fused `max(fma(scale, hu, offset), 0)` elementwise kernel.
#[derive(Clone, Copy, Debug, Default)]
struct HuToMuKernel;

const WORKGROUP_WIDTH: u32 = 256;

impl KernelInterface for HuToMuKernel {
    type Params = AffineClampParams;
    const LABEL: &'static str = "helios-hu-to-mu";
    const BINDINGS: &'static [BindingDecl] = &[
        BindingDecl::read_only::<f32>(),
        BindingDecl::read_write::<f32>(),
    ];
    const WORKGROUP: [u32; 3] = [WORKGROUP_WIDTH, 1, 1];
}

impl KernelSource<Wgsl> for HuToMuKernel {
    const ENTRY: &'static str = "main";

    fn source(&self) -> Cow<'static, str> {
        // Bindings follow the seam ABI: storage bindings 0..N-1 in
        // KernelInterface::BINDINGS order, uniform params at binding N.
        Cow::Borrowed(
            r#"struct Params {
    scale: f32,
    offset: f32,
    len: u32,
    _pad: u32,
}

@group(0) @binding(0) var<storage, read> hu: array<f32>;
@group(0) @binding(1) var<storage, read_write> mu: array<f32>;
@group(0) @binding(2) var<uniform> params: Params;

@compute @workgroup_size(256)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    let i = gid.x;
    if (i >= params.len) {
        return;
    }
    mu[i] = max(fma(params.scale, hu[i], params.offset), 0.0);
}
"#,
        )
    }
}

/// Resident GPU HU→μ mapper: the pipeline is prepared once at construction
/// and reused across dispatches (mirrors [`GpuProjector`](crate::GpuProjector)'s
/// residency pattern).
///
/// Generic over any WGSL-dialect [`KernelDevice`] so the same mapper runs on
/// the native wgpu backend or a Metal-pinned device without change.
pub struct GpuAttenuationMapper<D: KernelDevice<Dialect = Wgsl>> {
    device: D,
    prepared: D::Prepared<HuToMuKernel>,
    scale: f32,
    offset: f32,
}

impl<D: KernelDevice<Dialect = Wgsl>> GpuAttenuationMapper<D> {
    /// Prepare the mapper for a beam energy's water mass-attenuation
    /// coefficient `mu_over_rho_cm2_g` (cm²/g) and reference water density
    /// `water_density_g_cm3` (g/cm³).
    ///
    /// # Errors
    /// Returns [`HephaestusError::DispatchFailed`] when a coefficient is
    /// non-finite or negative, or the backend's typed failure when pipeline
    /// preparation fails.
    pub fn new(device: D, mu_over_rho_cm2_g: f32, water_density_g_cm3: f32) -> Result<Self> {
        if !mu_over_rho_cm2_g.is_finite() || mu_over_rho_cm2_g < 0.0 {
            return Err(HephaestusError::DispatchFailed {
                message: format!(
                    "mass attenuation μ/ρ must be finite and ≥ 0, got {mu_over_rho_cm2_g}"
                ),
            });
        }
        if !water_density_g_cm3.is_finite() || water_density_g_cm3 <= 0.0 {
            return Err(HephaestusError::DispatchFailed {
                message: format!("water density must be finite and > 0, got {water_density_g_cm3}"),
            });
        }
        let offset = mu_over_rho_cm2_g * water_density_g_cm3;
        let scale = offset / 1000.0;
        let prepared = device.prepare(&HuToMuKernel)?;
        Ok(Self {
            device,
            prepared,
            scale,
            offset,
        })
    }

    /// Compute `μ[i] = max(scale·HU[i] + offset, 0)` on the GPU.
    ///
    /// # Errors
    /// Returns [`HephaestusError::LengthMismatch`] when `out.len()` differs
    /// from `ct_hu.len()`, [`HephaestusError::DispatchFailed`] when the
    /// element count exceeds `u32` range, or the backend's typed
    /// transfer/dispatch failure.
    pub fn map_into(&self, ct_hu: &[f32], out: &mut [f32]) -> Result<()> {
        if ct_hu.len() != out.len() {
            return Err(HephaestusError::LengthMismatch {
                host_len: ct_hu.len(),
                device_len: out.len(),
            });
        }
        if ct_hu.is_empty() {
            return Ok(());
        }
        let len = u32::try_from(ct_hu.len()).map_err(|_| HephaestusError::DispatchFailed {
            message: format!("HU buffer length {} exceeds u32 range", ct_hu.len()),
        })?;

        let input = self.device.upload(ct_hu)?;
        let output = self.device.alloc_zeroed::<f32>(ct_hu.len())?;
        let params = AffineClampParams {
            scale: self.scale,
            offset: self.offset,
            len,
            _pad: 0,
        };
        let grid =
            DispatchGrid::covering_domain([ct_hu.len(), 1, 1], [WORKGROUP_WIDTH as usize, 1, 1])?;
        self.device.dispatch(
            &self.prepared,
            &[Binding::read(&input), Binding::read_write(&output)],
            &params,
            grid,
        )?;
        self.device.download(&output, out)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::default_device;

    // NIST-representative water μ/ρ at ~1 MeV and unit water density; tests
    // verify the defining relation, not these specific magnitudes.
    const MU_OVER_RHO: f32 = 0.0707;
    const WATER_DENSITY: f32 = 1.0;

    fn cpu_reference(hu: f32) -> f32 {
        // μ = (μ/ρ)·ρ_water·max(1 + HU/1000, 0), computed in f32 like the GPU.
        MU_OVER_RHO * WATER_DENSITY * (1.0_f32 + hu / 1000.0).max(0.0)
    }

    #[test]
    fn gpu_matches_cpu_closed_form_including_clamp() {
        let Ok(device) = default_device() else {
            eprintln!("no GPU adapter — skipping HU→μ closed-form test");
            return;
        };
        let mapper = GpuAttenuationMapper::new(device, MU_OVER_RHO, WATER_DENSITY).expect("mapper");
        // Air (−1000, clamps to exactly 0), below-air (clamped), lung, water,
        // soft tissue, trabecular and cortical bone, metal.
        let hu = [
            -2000.0_f32,
            -1000.0,
            -700.0,
            -100.0,
            0.0,
            40.0,
            300.0,
            1500.0,
            3000.0,
        ];
        let mut got = [0.0_f32; 9];
        mapper.map_into(&hu, &mut got).expect("gpu HU→μ");
        for (&h, &g) in hu.iter().zip(got.iter()) {
            let cpu = cpu_reference(h);
            // fma vs mul-add differ by ≤1 ulp; bound relative 1e-6 + tiny abs.
            assert!(
                (g - cpu).abs() <= 1e-6 * (1.0 + cpu.abs()),
                "HU={h}: gpu={g} vs cpu={cpu}"
            );
            assert!(
                g >= 0.0,
                "attenuation must be non-negative, got {g} at HU={h}"
            );
        }
        // The clamp region must be exactly zero, not merely close.
        assert_eq!(got[0], 0.0, "below-air HU must clamp to exactly 0");
    }

    #[test]
    fn gpu_matches_helios_solver_attenuation_map() {
        use helios_math::Point3;

        let Ok(device) = default_device() else {
            eprintln!("no GPU adapter — skipping HU→μ solver-differential test");
            return;
        };
        let grid = helios_domain::VoxelGrid::axis_aligned(
            [4, 3, 2],
            [1.0_f32, 1.0, 1.0],
            Point3::new(0.0_f32, 0.0, 0.0),
        )
        .expect("valid grid");
        // Deterministic HU pattern spanning clamp and bone regions.
        let ct = helios_domain::Volume::from_shape_fn(grid, |idx| {
            -1500.0_f32 + 260.0 * (idx[0] + 4 * idx[1] + 12 * idx[2]) as f32
        });

        let mass_atten = helios_physics::MassAttenuation::new(MU_OVER_RHO).expect("μ/ρ ≥ 0");
        let reference = helios_solver::attenuation_map(&ct, mass_atten, WATER_DENSITY);

        let mut hu_flat = Vec::with_capacity(4 * 3 * 2);
        for k in 0..2 {
            for j in 0..3 {
                for i in 0..4 {
                    hu_flat.push(ct.get(i, j, k).expect("in-grid index"));
                }
            }
        }
        let mapper = GpuAttenuationMapper::new(device, MU_OVER_RHO, WATER_DENSITY).expect("mapper");
        let mut got = vec![0.0_f32; hu_flat.len()];
        mapper.map_into(&hu_flat, &mut got).expect("gpu HU→μ");

        let mut flat_idx = 0;
        for k in 0..2 {
            for j in 0..3 {
                for i in 0..4 {
                    let want = reference.get(i, j, k).expect("in-grid index");
                    let g = got[flat_idx];
                    assert!(
                        (g - want).abs() <= 1e-6 * (1.0 + want.abs()),
                        "voxel ({i},{j},{k}): gpu={g} vs solver={want}"
                    );
                    flat_idx += 1;
                }
            }
        }
    }

    #[test]
    fn length_mismatch_is_a_typed_error() {
        let Ok(device) = default_device() else {
            return;
        };
        let mapper = GpuAttenuationMapper::new(device, MU_OVER_RHO, WATER_DENSITY).expect("mapper");
        let mut out = [0.0_f32; 3];
        let err = mapper.map_into(&[0.0; 4], &mut out).unwrap_err();
        assert!(matches!(err, HephaestusError::LengthMismatch { .. }));
    }

    #[test]
    fn invalid_coefficients_are_rejected() {
        let Ok(device) = default_device() else {
            return;
        };
        assert!(GpuAttenuationMapper::new(device.clone(), f32::NAN, 1.0).is_err());
        assert!(GpuAttenuationMapper::new(device.clone(), -0.1, 1.0).is_err());
        assert!(GpuAttenuationMapper::new(device, 0.0707, 0.0).is_err());
    }

    #[test]
    fn empty_input_is_ok() {
        let Ok(device) = default_device() else {
            return;
        };
        let mapper = GpuAttenuationMapper::new(device, MU_OVER_RHO, WATER_DENSITY).expect("mapper");
        let mut out: [f32; 0] = [];
        mapper.map_into(&[], &mut out).expect("empty is a no-op");
    }
}
