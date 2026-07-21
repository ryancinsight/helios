# Example: LINAC Dose Accumulation

> **Source:** `crates/helios-simulation/examples/linac_dose_accumulation.rs`
>
> **Run:** `cargo run -p helios-simulation --example linac_dose_accumulation`

## Overview

Demonstrates a 4-field box LINAC step-and-shoot delivery on a uniform water phantom using `helios_simulation::accumulate_delivered_dose`. Constructs four `DeliveryFrame`s (0°/90°/180°/270°), applies parallel-beam ray-tracing through the attenuation map, and verifies the resulting dose DVH.

| Case | Verified |
|------|----------|
| Non-zero total dose | `D_max > 0` |
| DVH monotonicity | `D(v)` non-increasing |
| Physical uniformity | `D_min / D_max ≥ 0` |

## Key APIs

```rust
use helios_simulation::{accumulate_delivered_dose, BeamGeometry, DeliveryFrame};

let frames = vec![
    DeliveryFrame { projection: 0, gantry_angle_rad: 0.0, couch_mm: 0.0,
                    leaf_fluence: vec![1.0; 16] },
    // ... 90°, 180°, 270°
];

let dose = accumulate_delivered_dose(
    &frames, &mu_volume,
    BeamGeometry::Parallel { standoff_mm: 500.0 },
    2.0,  // leaf_width_mm
    0.5,  // step_mm
);
```

## Delivery Model

Each `DeliveryFrame` carries:
- `gantry_angle_rad` — the beam direction in the axial plane
- `leaf_fluence` — per-leaf monitor-unit weight (16 leaves × 2 mm = 32 mm field)

`BeamGeometry::Parallel` uses a small-fan approximation: all beamlets run parallel along the gantry direction, offset laterally by `(leaf - centre) × leaf_width_mm`. Switch to `BeamGeometry::PointSource { source_axis_mm: 850.0 }` for a true divergent SAD geometry.

## Book Chapter

[← LINAC-Based Step-and-Shoot Delivery](../workflow_linac.md)
