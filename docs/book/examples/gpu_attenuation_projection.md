# Example: GPU Attenuation Map and Forward Projection

**Crate**: `helios-gpu`  
**Run**: `cargo run -p helios-gpu --example gpu_attenuation_projection`  
**Source**: [`crates/helios-gpu/examples/gpu_attenuation_projection.rs`](../../../crates/helios-gpu/examples/gpu_attenuation_projection.rs)

## What This Example Demonstrates

Two-stage GPU pipeline over the Atlas `hephaestus` compute substrate:

| Stage | API | GPU Kernel |
|---|---|---|
| HU → μ | `GpuAttenuationMapper::new(device, mu_over_rho, water_rho)` | Fused affine-clamp WGSL: `max(scale·HU + offset, 0)` |
| GPU dispatch | `mapper.map_into(&hu_f32, &mut mu_gpu)` | One wgpu compute pass |
| CPU reference | `helios_solver::attenuation_map(&ct_hu, mu_over_rho, water_rho)` | Differential validation |

## Key Code Snippet

```rust
use helios_gpu::{default_device, GpuAttenuationMapper};

let device = default_device()?;

// Prepare the kernel once (pipeline compiled and cached on-device)
let mapper = GpuAttenuationMapper::new(device, 0.0636_f32, 1.0_f32)?;

// Dispatch: one GPU compute pass for all voxels
mapper.map_into(&hu_f32, &mut mu_out)?;
// mu_out[i] = max(0.0636/1000 * hu[i] + 0.0636, 0)
```

## Physics Background

The GPU kernel implements the Compton-dominated 6 MV attenuation approximation:

```
μ = (μ/ρ) · ρ_water · max(1 + HU/1000, 0)
  = scale · HU + offset   where  scale = (μ/ρ)·ρ/1000, offset = (μ/ρ)·ρ
```

The `max(…, 0)` clamp encodes `ρ ≥ 0` (air is non-negative attenuation).
The fused multiply-add (`fma`) in WGSL eliminates a floating-point rounding
step that would separate the multiply and add.

## GPU Availability

The example detects the GPU at runtime via `default_device()`. If no
compatible adapter is found (CI, embedded systems) it falls back gracefully
and the CPU physics assertions still pass.

## Differential Validation

GPU and CPU must agree to `|μ_gpu − μ_cpu| < 1e-4 cm⁻¹` (f32 vs f64
precision). This matches the `helios-gpu` differential-verification contract
(see `crates/helios-gpu/src/attenuation.rs` tests).

## Book Chapter

[← GPU-Accelerated Dose Kernels](../gpu_dose.md)
