# Example: SIRT Iterative CT Reconstruction

**Crate**: `helios-imaging`  
**Run**: `cargo run -p helios-imaging --example sirt_reconstruction`  
**Source**: [`crates/helios-imaging/examples/sirt_reconstruction.rs`](../../../crates/helios-imaging/examples/sirt_reconstruction.rs)

## What This Example Demonstrates

Algebraic iterative reconstruction of a 64×64 water/bone phantom from a 32-angle
sinogram, comparing SIRT convergence against the FBP analytic baseline.

| Stage | API | Physics |
|---|---|---|
| Phantom | `Volume::from_shape_fn` | Water μ = 0.02 cm⁻¹, Bone μ = 0.04 cm⁻¹ |
| Sinogram | `parallel_beam_radon(&phantom, &angles, &offsets, sad_mm, step_mm)` | 32-angle parallel beam |
| FBP (baseline) | `filtered_back_projection(&sinogram, &grid, sad_mm, step_mm)` | Ramp-filtered back-projection |
| SIRT | `sirt_reconstruction(&sinogram, &grid, sad_mm, step_mm, iters, λ)` | 10 iterations, λ = 1.0 |

## Key Code Snippet

```rust
use helios_imaging::{parallel_beam_radon, filtered_back_projection, sirt_reconstruction};

let sinogram = parallel_beam_radon(&phantom, &angles, &offsets, 500.0, 1.0);

// Analytic (one-pass) baseline
let fbp = filtered_back_projection(&sinogram, &grid, 500.0, 1.0);

// Iterative — robust to undersampling
let sirt = sirt_reconstruction(&sinogram, &grid, 500.0, 1.0, 10, 1.0);
```

## Physics Background

SIRT solves the linear system `Ax = b` (sinogram = projection matrix × attenuation)
iteratively via normalized back-projection updates:

```
x ← max(0, x + λ · C⁻¹ ⊙ Aᵀ( R⁻¹ ⊙ (b − Ax) ))
```

where `R⁻¹` normalizes by ray chord length and `C⁻¹` by voxel hit weight.
The non-negativity projection encodes `μ ≥ 0`. FBP is O(N² log N) but
creates streak artefacts under sparse sampling; SIRT is O(N² × iterations)
and converges to the least-squares solution.

## Book Chapter

[← MVCT and Correction Workflows](../imaging_mvct.md)
