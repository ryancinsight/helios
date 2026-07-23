# Example: Adaptive Radiotherapy Workflow

> **Source:** `crates/helios-simulation/examples/adaptive_rt_workflow.rs`
>
> **Run:** `cargo run -p helios-simulation --example adaptive_rt_workflow`

## Overview

Demonstrates the daily IGRT–ART decision loop on a synthetic 32 × 32 phantom:
plan → daily MVCT setup correction → dose recomputation → adaptive decision gate.

| Phase | Action | Key API |
|-------|--------|---------|
| 1 — Planning | Build CT phantom, compute 4-field box dose | `attenuation_map`, `accumulate_delivered_dose` |
| 2 — Daily MVCT | Simulate 2-voxel setup error; register | `register_translation` |
| 3 — Corrected delivery | Apply couch shift; recompute dose | `accumulate_delivered_dose` |
| 4 — Adaptive gate | Gamma + DVH comparison; proceed/replan | `gamma_index_3d`, `Dvh` |

## Expected Output

```
Phase 2: Daily MVCT (simulated setup error)
  Applied shift:    [2, 1]
  Recovered shift:  [2, 1]
  ✓ Registration exact
Phase 4: Adaptive decision gate
  Gamma pass-rate (3%/2mm):    100.0%
  Decision: PROCEED — couch-corrected delivery within clinical tolerance
```

## Clinical Acceptance Criteria

| Metric | Threshold | API |
|--------|-----------|-----|
| 3 %/2 mm gamma pass-rate | ≥ 95 % | `gamma_index_3d`, `gamma_pass_rate` |
| Mean dose deviation | < 5 % | `Dvh::mean().into_base()` at the scalar gamma boundary |

When either criterion fails the workflow returns `REPLAN`, triggering an online
or offline replanning cycle.

## Key APIs

```rust
use helios_imaging::register_translation;
use helios_simulation::{accumulate_delivered_dose, BeamGeometry, DeliveryFrame};
use aequitas::systems::si::{quantities::{AbsorbedDose, Length}, units::Millimeter};
use helios_analysis::{gamma_index_3d, gamma_pass_rate, roi_statistics, Dvh};

// Register daily to planning
let shift = register_translation(&plan_ct, &daily_ct, [4, 4, 0]);

// Recompute dose on corrected anatomy
let corrected = shift_phantom(&daily_ct, -shift[0], -shift[1]);
let corrected_dose = accumulate_delivered_dose(&frames, &mu, geometry, lw, step);

// Adaptive gate
let gamma = gamma_index_3d(
    &plan_dose,
    &corrected_dose,
    0.03,
    Length::from_unit::<Millimeter>(2.0),
    d_max,
    Length::from_unit::<Millimeter>(6.0),
)?;
let pass = gamma_pass_rate(&gamma, &plan_dose, AbsorbedDose::from_base(0.0));
```

## Book Chapter

[← Adaptive Radiotherapy with MVCT](../workflow_adaptive.md)
