# Example: Clinical Protocol Validation

> **Status:** Planned — implementation forthcoming.
>
> **Source:** `crates/helios-analysis/examples/validation_clinical.rs` *(not yet created)*

## Overview

This chapter will verify helios dose engines against standard clinical acceptance
protocols:

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
