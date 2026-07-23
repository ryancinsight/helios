# Example: Regression and Analytical Validation

> **Source:** `crates/helios-analysis/examples/validation_regression.rs`
>
> **Run:** `cargo run -p helios-analysis --example validation_regression`

## Overview

This example verifies three fundamental mathematical properties of `helios-analysis`
against exact oracles, establishing the quantitative accuracy floor for the entire
dose-engine and imaging stack.

| Case | What is checked | Expected result |
|------|-----------------|-----------------|
| 1 — Gamma self-consistency | γ of any volume against itself | γ = 0 everywhere; 100 % pass rate |
| 2 — Radon / FBP round-trip | RMSE of water-cylinder FBP reconstruction | < 1.5 × 10⁻² cm⁻¹ (coarse 32×32 grid) |
| 3 — DVH monotonicity | Cumulative DVH of a ramp dose | Strictly non-increasing D(v) |

## Case 1 — Gamma Self-Consistency

```rust
let gamma = gamma_index_3d(&dose, &dose, 0.03, 2.0, 1.0, 6.0)?;
let pass_rate = gamma_pass_rate(&gamma, &dose, 0.0);
```

`gamma_index_3d` takes `(reference, evaluated, dose_diff_criterion, dta_mm,
normalization_dose, search_radius_mm)`. Comparing a distribution against itself
yields γ = 0 at every voxel → 100 % pass rate.

## Case 2 — Radon / FBP Round-Trip

```rust
let sinogram = parallel_beam_radon(&phantom, &angles, &offsets, source_mm, step_mm);
let recon    = filtered_back_projection(&sinogram, &recon_grid);
let rmse     = volume_rmse(&recon, &phantom)?;
```

A uniform water cylinder (μ = 0.0636 cm⁻¹, r = 25 mm) is projected through 180
angles and reconstructed with the Ram-Lak FBP kernel. On a coarse 32 × 32 / 2 mm
grid the RMSE is ≈ 7 × 10⁻³ cm⁻¹ (≈ 11 % of μ_water), well within the 1.5 × 10⁻²
tolerance.  Higher-resolution grids approach clinical MVCT accuracy (≤ 1 %).

## Case 3 — DVH Monotonicity

```rust
let dvh = Dvh::from_volume(&ramp);
let d   = dvh.dose_at_volume_fraction(v).into_base();
```

The cumulative DVH of any physical dose must be non-increasing: every higher volume
fraction receives *at most* the dose of the smaller fraction.  A linearly-ramping
field provides a known monotone ground truth.

## Key APIs

| Function / type | Crate |
|-----------------|-------|
| `gamma_index_3d`, `gamma_pass_rate` | `helios-analysis` |
| `volume_rmse`, `volume_relative_l2_error` | `helios-analysis` |
| `Dvh::from_volume`, `dose_at_volume_fraction` | `helios-analysis` |
| `parallel_beam_radon`, `filtered_back_projection` | `helios-imaging` |
| `Volume`, `VoxelGrid` | `helios-domain` |
