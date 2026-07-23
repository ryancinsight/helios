# Example: Gamma Index Comparison

**Crate**: `helios-analysis`  
**Run**: `cargo run -p helios-analysis --example gamma_index`  
**Source**: [`crates/helios-analysis/examples/gamma_index.rs`](../../../crates/helios-analysis/examples/gamma_index.rs)

## What This Example Demonstrates

Computes the 3%/2 mm global-normalization gamma index — the standard clinical
pass/fail criterion for comparing a calculated or delivered dose against a
reference. Tests three scenarios:

| Scenario | Shift | Expected Pass Rate |
|---|---|---|
| Identical plan | 0 mm | 100% (γ ≈ 0) |
| Small shift | 1.5 mm | High (within DTA) |
| Large shift | 4.0 mm | Lower (exceeds DTA) |

## Key Code Snippet

```rust
use aequitas::systems::si::{
    quantities::{AbsorbedDose, Length},
    units::{Gray, Millimeter},
};
use helios_analysis::{gamma_index_3d, gamma_pass_rate};

// 3%/2 mm global normalization
let gamma = gamma_index_3d(
    &reference, &evaluated,
    0.03,   // dose-difference criterion (3%)
    Length::from_unit::<Millimeter>(2.0),
    AbsorbedDose::from_unit::<Gray>(60.0),
    Length::from_unit::<Millimeter>(6.0),
)?;

let pass_rate = gamma_pass_rate(
    &gamma,
    &reference,
    AbsorbedDose::from_unit::<Gray>(6.0),
);
// pass_rate >= 0.95 → clinical acceptance (95% pass rate)
```

## Clinical Context

Low's gamma index (Med. Phys. 25, 1998):

```
γ(r) = min_e √( |Δx|²/Δd²  +  ΔD²/ΔD_crit² )
```

- `Δd`: distance-to-agreement criterion (typically 2 or 3 mm)
- `ΔD_crit`: dose-difference criterion (typically 3% of global max)
- `γ ≤ 1` at a point → that point **passes**
- Clinical acceptance: **≥ 95%** of reference points above the 10% dose
  threshold must pass the 3%/3 mm or 3%/2 mm criterion

## Book Chapter

[← Gamma Index and Plan Verification](../planning_gamma.md)
