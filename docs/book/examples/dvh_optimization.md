# Example: DVH-Constrained Beam-Weight Optimization

**Crate**: `helios-planning`  
**Run**: `cargo run -p helios-planning --example dvh_optimization`  
**Source**: [`crates/helios-planning/examples/dvh_optimization.rs`](../../../crates/helios-planning/examples/dvh_optimization.rs)

## What This Example Demonstrates

Inverse treatment planning on a synthetic 3-field head-and-neck geometry:
4 PTV voxels and 2 OAR voxels with 3 beamlets.

| Stage | API | Purpose |
|---|---|---|
| Build dose matrix | `DoseInfluence::from_rows(voxels, beamlets, data)` | Linear dose model `dose = A·x` |
| Optimize weights | `optimize_beam_weights(&A, &prescription, iters, step)` | Projected gradient descent `x ← max(0, x − step·Aᵀ(Ax−d))` |
| Evaluate objective | `objective_value(&A, &x, &prescription)` | `½‖Ax − d‖²` |
| DVH metrics | Manual D95, D_mean, D_max | Coverage and OAR constraints |

## Key Code Snippet

```rust
use helios_planning::{objective_value, optimize_beam_weights, DoseInfluence};

// 6 voxels (4 PTV + 2 OAR) × 3 beamlets
let influence = DoseInfluence::from_rows(6, 3, matrix)?;
let prescription = vec![2.0, 2.0, 2.0, 2.0, 0.0, 0.0]; // Gy

let weights = optimize_beam_weights(&influence, &prescription, 2000, 0.1);
// → Beam 1: 1.32, Beam 2: 0.87, Beam 3: 1.44

let final_obj = objective_value(&influence, &weights, &prescription);
```

## Physics Background

The quadratic objective `½‖Ax − d‖²` is convex and differentiable. The gradient is
`Aᵀ(Ax − d)`. Projected gradient descent enforces the non-negativity constraint
`x ≥ 0` (beamlet weights cannot be negative) by clamping after each step.

The **DVH** (Dose-Volume Histogram) D95 metric is the dose received by ≥ 95 % of the
target volume — the standard ICRU clinical coverage criterion.

## Book Chapter

[← Dose-Volume Histograms](../planning_dvh.md)
