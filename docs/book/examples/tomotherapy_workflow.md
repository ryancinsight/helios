# Example: Tomotherapy Workflow

**Crate**: `helios-simulation`  
**Run**: `cargo run -p helios-simulation --example tomotherapy_workflow [output_dir]`  
**Source**: [`crates/helios-simulation/examples/tomotherapy_workflow.rs`](../../../crates/helios-simulation/examples/tomotherapy_workflow.rs)

## What This Example Demonstrates

A complete end-to-end helical TomoTherapy simulation on a synthetic 61×61×5 CT phantom
(64×64×64 voxel water cylinder with bone insert). The workflow proceeds through five
verifiable stages:

| Stage | API | Output |
|---|---|---|
| 1. CT phantom generation | `Volume::from_shape_fn` | `ct.png` |
| 2. Attenuation map | `attenuation_map` | `mu.png` |
| 3. MVCT simulation | `parallel_beam_radon` → `filtered_back_projection` | `recon.png` |
| 4. Helical dose delivery | `simulate_helical_delivery` → `accumulate_delivered_dose_anisotropic` | `dose.png` |
| 5. DVH + Gamma | `Dvh`, `gamma_index_3d`, `gamma_pass_rate` | Console summary |

## Quantitative Verification

The example is self-validating — it prints and asserts:

- **Reconstruction RMSE** < 200 HU (FBP is approximate, not exact)
- **Dose mean** > 0 Gy (sanity: dose was actually delivered)
- **Gamma pass rate** (3%/2 mm, 5 Gy threshold) > 90%

These are not fixed-result assertions but physics-plausibility bounds: the dose
distribution should be spatially consistent with the delivery geometry.

## Architecture Diagram

```
VoxelGrid → Volume<f64> (HU)
     │
     ├─ attenuation_map  ───────────────── μ map
     │                                       │
     ├─ parallel_beam_radon ── sinogram      │
     │    └─ filtered_back_projection        │
     │                                       │
     ├─ simulate_helical_delivery ── terma   │
     │    └─ CollapsedCone (poly-energetic)  │
     │         └─ accumulate dose ──── Dose Volume<f64>
     │                                       │
     └─ Dvh + gamma_index_3d ─────── QA report
```

## Running the Example

```bash
# Default output to ./helios_workflow_output/
cargo run -p helios-simulation --example tomotherapy_workflow

# Custom output directory
cargo run -p helios-simulation --example tomotherapy_workflow -- /tmp/helios_demo/
```

## Book Chapter

[← TomoTherapy End-to-End Workflow](../workflow_tomotherapy.md)
