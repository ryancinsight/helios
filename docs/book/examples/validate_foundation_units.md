# Example: Validating Foundation Units

**Crate**: `helios-core`  
**Run**: `cargo run -p helios-core --example validate_foundation_units`  
**Source**: [`crates/helios-core/examples/validate_foundation_units.rs`](../../../crates/helios-core/examples/validate_foundation_units.rs)

## What This Example Demonstrates

This example exercises the three typestate boundary-guards in `helios-core`:

| Type | Accepts | Rejects |
|---|---|---|
| `EnergyMeV` | `(0, ∞)` finite | zero, negative, `NaN`, `±Inf` |
| `HounsfieldUnit` | `[-1024, 31743]` | out-of-range, non-finite |
| `VoxelSpacingMm` | `(0, ∞)` finite | zero, negative, non-finite |

It confirms that:

1. **Happy-path construction** round-trips `get()` with identity.
2. **`Display` carries the unit suffix** — `6.0 MeV` → `"6 MeV"`.
3. **Every invalid value** (non-finite, out-of-range, zero) produces
   `HeliosError::InvalidDomainValue`.
4. **Boundary values are inclusive** — `HounsfieldUnit::MIN` and `MAX` succeed.

## Key Code Snippet

```rust
use helios_core::{EnergyMeV, HeliosError, HounsfieldUnit, VoxelSpacingMm};

// Valid construction
let energy = EnergyMeV::try_from(6.0).expect("6 MeV is a clinically valid beam energy");
assert_eq!(energy.to_string(), "6 MeV");

// Rejection contract
assert!(matches!(
    EnergyMeV::try_from(0.0),
    Err(HeliosError::InvalidDomainValue { .. })
));
```

## Book Chapter

[← Physics Domain Types and Safety Boundaries](../foundations.md)
