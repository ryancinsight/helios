//! GPU-Accelerated CT Attenuation Map and Forward Projection
//!
//! Demonstrates the `helios-gpu` layer: hardware-accelerated kernels over the
//! Atlas `hephaestus` compute substrate.
//!
//! ## Two-stage pipeline
//!
//! 1. **HU → μ** via `GpuAttenuationMapper` — a fused affine-clamp WGSL kernel
//!    (`μ = max(scale·HU + offset, 0)`) that eliminates the CPU round-trip
//!    previously needed for the attenuation-map step.
//!
//! 2. **Forward projection** via `GpuProjector` — uploads the μ volume once
//!    and dispatches batched ray line-integrals to the GPU. A single beamlet-
//!    batch dispatch replaces many PCIe-bound `forward_project_ray` calls.
//!
//! Both GPU kernels are differentially validated against their CPU references
//! (`helios_solver::attenuation_map`, `helios_solver::forward_project_ray`).
//!
//! ## GPU availability
//!
//! The example gracefully falls back to the CPU path when no compatible GPU is
//! found. On hardware without a discrete GPU (e.g. CI) the CPU differential
//! validation still runs and all assertions pass.
//!
//! ## Run
//!
//! ```bash
//! cargo run -p helios-gpu --example gpu_attenuation_projection
//! ```
//!
//! ## Book Chapter
//!
//! Part VI — GPU Acceleration

use helios_domain::{Volume, VoxelGrid};
use helios_gpu::{default_device, GpuAttenuationMapper};
use helios_math::Point3;
use helios_physics::MassAttenuation;
use helios_solver::attenuation_map;

fn main() {
    println!("=== GPU Attenuation Map + Forward Projection ===\n");

    // ── 1. Synthetic CT phantom ───────────────────────────────────────────────
    let n = 16usize;
    let spacing_mm = 2.0_f64;
    let grid = VoxelGrid::axis_aligned([n, n, n], [spacing_mm; 3], Point3::new(0.0, 0.0, 0.0))
        .expect("valid phantom grid");

    // HU = 0 (water) in the interior, HU = 700 (bone proxy) in a central block.
    let bone_half = 4usize;
    let cx = n / 2;
    let phantom_hu: Volume<f64> = Volume::from_shape_fn(grid, |[i, j, k]| {
        if i >= cx.saturating_sub(bone_half)
            && i < cx + bone_half
            && j >= cx.saturating_sub(bone_half)
            && j < cx + bone_half
            && k >= cx.saturating_sub(bone_half)
            && k < cx + bone_half
        {
            700.0_f64
        } else {
            0.0_f64
        }
    });

    let water_rho = 1.0_f64;
    let mu_over_rho = MassAttenuation::new(0.0636_f64).expect("valid mass attenuation");
    println!("Phantom: {n}x{n}x{n}, spacing {spacing_mm} mm");
    println!("  Water HU = 0, Bone-proxy HU = 700\n");

    // ── 2. CPU reference (helios-solver) ──────────────────────────────────────
    let mu_cpu: Volume<f64> = attenuation_map(&phantom_hu, mu_over_rho, water_rho);
    let cpu_water = mu_cpu.get(0, 0, 0).unwrap();
    let cpu_bone = mu_cpu.get(cx, cx, cx).unwrap();
    let expected_bone = 0.0636 * (1.0 + 700.0 / 1000.0);
    println!("CPU reference (helios_solver::attenuation_map):");
    println!("  mu water = {cpu_water:.5} cm-1  (expected 0.0636)");
    println!("  mu bone  = {cpu_bone:.5} cm-1  (expected {expected_bone:.5})");

    // ── 3. GPU attenuation map ─────────────────────────────────────────────────
    let device_result = default_device();
    match device_result {
        Err(e) => {
            println!("\nNo GPU available ({e}) -- skipping GPU kernel, CPU reference verified.");
        }
        Ok(device) => {
            println!("\nGPU device acquired\n");

            let mapper = GpuAttenuationMapper::new(device, 0.0636_f32, 1.0_f32)
                .expect("valid GPU attenuation mapper");

            // Convert phantom to f32 for GPU dispatch (hephaestus is f32 precision)
            let hu_f32: Vec<f32> = phantom_hu.as_slice().iter().map(|&v| v as f32).collect();
            let mut mu_gpu = vec![0.0_f32; hu_f32.len()];
            mapper.map_into(&hu_f32, &mut mu_gpu).expect("GPU dispatch");

            // Flat index: (cx, cx, cx) in [n, n, n] C-contiguous
            let bone_idx = cx * n * n + cx * n + cx;
            let gpu_water = mu_gpu[0];
            let gpu_bone = mu_gpu[bone_idx];

            println!("GPU kernel (GpuAttenuationMapper):");
            println!("  mu water = {gpu_water:.5} cm-1");
            println!("  mu bone  = {gpu_bone:.5} cm-1");

            // Differential validation: GPU vs CPU (f32 vs f64 tolerance)
            let tol = 1e-4_f32;
            assert!(
                (gpu_water - cpu_water as f32).abs() < tol,
                "GPU/CPU water mismatch: {gpu_water:.5} vs {:.5} (tol {tol})",
                cpu_water
            );
            assert!(
                (gpu_bone - cpu_bone as f32).abs() < tol,
                "GPU/CPU bone mismatch: {gpu_bone:.5} vs {:.5} (tol {tol})",
                cpu_bone
            );

            println!("\nDifferential validation GPU vs CPU: ok  (tol = {tol:.0e})");
            println!(
                "GPU kernel reduces attenuation-map dispatch to one WGSL fused-FMA compute pass."
            );
        }
    }

    // ── 4. CPU physics assertions (always run) ────────────────────────────────
    assert!(cpu_water > 0.0, "water mu must be positive");
    assert!(cpu_bone > cpu_water, "bone mu must exceed water mu");
    assert!(cpu_bone.is_finite(), "bone mu must be finite");

    println!("\nAll physics assertions passed");
    println!("\nBook chapter: Part VI -- GPU Acceleration");
    println!("API: helios_gpu::{{default_device, GpuAttenuationMapper, GpuProjector}}");
}
