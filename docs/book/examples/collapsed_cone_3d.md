# Example: Collapsed-Cone 3-D Dose Engine

**Crate**: `helios-solver`  
**Run**: `cargo run -p helios-solver --example collapsed_cone_3d`  
**Source**: [`crates/helios-solver/examples/collapsed_cone_3d.rs`](../../crates/helios-solver/examples/collapsed_cone_3d.rs)

## What This Example Demonstrates

Two-stage deterministic photon dose calculation on a 10×10×10 synthetic water phantom:

| Stage | API | Physics |
|---|---|---|
| 1. CT → μ map | `attenuation_map(&ct_hu, mass_attenuation, water_rho)` | Compton-dominated MV μ = (μ/ρ)·ρ |
| 2. Primary transport | `primary_fluence_parallel_x(&mu, Ψ₀)` | Beer–Lambert attenuation along +x |
| 3. TERMA | `Ψ(i,j,k) · μ(i,j,k) · Δx_cm` | Energy released per voxel |
| 4. 1-D convolution | `dose_convolution_x(&terma, &kernel)` | Depth build-up along beam axis |
| 5. 3-D scatter | `scatter_superposition(&terma, kx, ky, kz)` | Lateral penumbra via separable kernels |

## Key Code Snippet

```rust
use helios_solver::{
    attenuation_map, dose_convolution_x, exponential_deposition_kernel,
    primary_fluence_parallel_x, scatter_superposition, symmetric_deposition_kernel,
};

// Stage 1: CT → μ
let mu = attenuation_map(&ct_hu, MassAttenuation::new(0.0636)?, 1.0);

// Stage 2: primary fluence (Beer–Lambert)
let primary = primary_fluence_parallel_x(&mu, 1.0);

// Stage 4: 1-D dose (collapsed cone along beam +x)
let kernel = exponential_deposition_kernel(0.5_f64, voxel_cm, 8);
let dose_1d = dose_convolution_x(&terma, &kernel);

// Stage 5: 3-D scatter superposition
let kx = symmetric_deposition_kernel(0.5, voxel_cm, 3);
let ky = symmetric_deposition_kernel(0.3, voxel_cm, 2);
let kz = symmetric_deposition_kernel(0.3, voxel_cm, 2);
let dose_3d = scatter_superposition(&terma, &kx, &ky, &kz);
```

## Physics Background

The **collapsed-cone convolution** dose engine is a deterministic (non-Monte Carlo)
method for photon dose calculation:

1. **Beer–Lambert**: `Ψ(x) = Ψ₀ · exp(−∫₀ˣ μ dl)` attenuates the primary beam.
2. **TERMA**: the energy deposited per unit mass at each voxel is `Ψ · μ · Δx`.
3. **Kernel superposition**: TERMA is spread to dose by convolution with a
   dose-deposition kernel derived from measured depth-dose data or Monte Carlo
   point-kernel calculations.

The separable 3-D scatter kernel approximates isotropic scatter with independent
per-axis Laplacian (exponential) shapes — fast and accurate for homogeneous tissue
at 6 MV where Compton scatter is dominant.

**Boundary truncation** (~22 % energy loss on a 30 mm phantom) is expected physics:
the lateral scatter tails of a 10-voxel phantom reach outside the grid. A clinical
calculation on a full CT volume has negligible boundary effects.

## Book Chapter

[← Collapsed-Cone Convolution](../dose_collapsed_cone.md)
