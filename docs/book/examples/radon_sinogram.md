# Example: Radon Sinogram

**Crate**: `helios-imaging`  
**Run**: `cargo run -p helios-imaging --example radon_sinogram`  
**Source**: [`crates/helios-imaging/examples/radon_sinogram.rs`](../../crates/helios-imaging/examples/radon_sinogram.rs)

## What This Example Demonstrates

Validates the parallel-beam Radon transform by projecting a uniform-density
disk phantom and comparing the sinogram against the analytical chord-length
integral.

| Step | API |
|---|---|
| Build disk phantom | `Volume::from_shape_fn(grid, ...)` |
| Project at each angle | `parallel_beam_radon(mu, angles, offsets, src_dist, step)` |
| Inspect sinogram | `sinogram.get(angle_idx, offset_idx)` |

## Key Code Snippet

```rust
use helios_imaging::radon::parallel_beam_radon;
use helios_domain::{Volume, VoxelGrid};
use helios_math::Point3;

let grid = VoxelGrid::axis_aligned([161, 161, 1], [0.5; 3], Point3::new(0., 0., 0.))?;
let phantom = Volume::from_shape_fn(grid, |[i, j, _]| /* disk μ */ 0.2_f64);

let sino = parallel_beam_radon(&phantom, &angles, &offsets, 200.0, 0.5);
// sino.get(angle_idx, offset_idx) → line integral (mm·cm⁻¹)
```

## Physics Background

The **Radon transform** maps a 2-D attenuation distribution `μ(x, y)` to its
line integrals `p(θ, s)`:

```
p(θ, s) = ∫ μ(s·cosθ − t·sinθ,  s·sinθ + t·cosθ) dt
```

For a uniform disk of radius `R` and attenuation `μ₀`, the analytical chord at
offset `s` is `2√(R² − s²)`, so `p(θ, s) = μ₀ · 2√(R² − s²) / 10` (mm → cm
conversion). This provides a closed-form oracle for validating the numerical
ray-march implementation.

## Book Chapter

[← Parallel-Beam Radon Transform](../imaging_radon.md)
