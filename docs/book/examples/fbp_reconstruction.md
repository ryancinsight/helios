# Example: FBP Reconstruction

**Crate**: `helios-imaging`  
**Run**: `cargo run -p helios-imaging --example fbp_reconstruction`  
**Source**: [`crates/helios-imaging/examples/fbp_reconstruction.rs`](../../crates/helios-imaging/examples/fbp_reconstruction.rs)

## What This Example Demonstrates

Full round-trip: disk phantom → parallel-beam Radon sinogram → Ram-Lak filtered
back-projection → reconstructed μ volume. Validates that the centre voxel
recovers the true attenuation `μ₀` to within 15% FBP discretization tolerance.

| Step | API |
|---|---|
| Forward project | `parallel_beam_radon(phantom, angles, offsets, src_dist, step)` |
| Reconstruct | `filtered_back_projection(sinogram, recon_grid)` |
| Sample result | `recon.get(i, j, k)` |

## Key Code Snippet

```rust
use helios_imaging::{fbp::filtered_back_projection, radon::parallel_beam_radon};

// 180 angles, 181 detector offsets
let sino = parallel_beam_radon(&phantom, &angles, &offsets, 400.0, 0.25);

// 41×41×1 reconstruction at 2 mm
let recon_grid = VoxelGrid::axis_aligned([41, 41, 1], [2.0; 3], origin)?;
let recon = filtered_back_projection(&sino, &recon_grid);

// Centre pixel should recover μ₀ = 0.04 cm⁻¹
let mu_centre = recon.get(20, 20, 0)?;
```

## Physics Background

Filtered back-projection inverts the Radon transform:

1. **Ramp filtering** — each projection `p(θ, s)` is convolved with the
   Ram-Lak kernel `h[n] = 1/(4Δs²)` at `n=0`, `−1/(π²n²Δs²)` for odd `n`.
2. **Back-projection** — the filtered projections are smeared back across the
   reconstruction grid, weighted by the angular step `Δθ`.

FBP is the analytical inverse for ideal continuous data; the 15% tolerance
accounts for voxelization, limited angular sampling, and ray-march discretization.

## Book Chapter

[← Filtered Back Projection](../imaging_fbp.md)
