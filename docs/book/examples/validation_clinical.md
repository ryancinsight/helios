# Example: Clinical Protocol Validation

**Crate**: `helios-analysis`  
**Run**: `cargo run --example validation_clinical -p helios-analysis`  
**Source**: [`crates/helios-analysis/examples/validation_clinical.rs`](../../../crates/helios-analysis/examples/validation_clinical.rs)

## Overview

Clinical plan validation for a head-and-neck treatment scenario with three structures:

- **PTV** (planning target volume) — ~60 Gy prescription
- **Parotid L** (organ at risk) — lower dose, sparing objective
- **Spinal cord** (serial OAR, strict dose limit)

For each structure the example computes:

1. **DVH coverage**: D₉₅, D_mean, homogeneity index (ICRU-83)
2. **Biological outcome**: gEUD, TCP (logistic), NTCP (Lyman-Kutcher-Burman)
3. **Dose ROI analysis**: Aequitas dose-valued ROI statistics and
   dimensionless contrast-to-noise ratio between PTV and parotid

## Atlas Integration

Uses `helios-analysis` directly: `Dvh`, `gamma_index_3d`, `gamma_pass_rate`, `dose_roi_statistics`, `contrast_to_noise_ratio`, `michelson_contrast`.

## Part Reference

Part VII — Clinical Protocol Compliance


| Case | Protocol reference | Metric |
|------|--------------------|--------|
| 1 — TRS-398 reference point | IAEA TRS-398 §7 | Absorbed dose in water at 10 cm, 6 MV: ± 2 % |
| 2 — TG-119 C-shape coverage | AAPM TG-119 | PTV D95 ≥ 95 % prescribed; OAR D5 ≤ 50 % |
| 3 — TomoTherapy self-gamma | 3 %/2 mm global | Pass rate ≥ 95 % |

## Planned APIs

```rust
use helios_analysis::{Dvh, gamma_index_3d, gamma_pass_rate};
use helios_solver::dose_engine;
```

See [Regression and Analytical Validation](validation_regression.md) for the
foundational metrics this chapter builds on.
