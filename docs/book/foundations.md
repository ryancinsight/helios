# Chapter 1 — Physics Domain Types and Safety Boundaries

Helios enforces physical validity at the type level through **validating newtypes** in
`helios-core`.  Each domain quantity is a distinct type that cannot be constructed
from an arbitrary float — it must pass a domain-specific validation contract.

## The Three Foundation Slots

| Type | Unit | Valid range | Rejects |
|---|---|---|---|
| `EnergyMeV` | MeV | `(0, ∞)` finite | zero, negative, `NaN`, `±Inf` |
| `HounsfieldUnit` | HU | `[-1024, 31743]` | out-of-range, non-finite |
| `VoxelSpacingMm` | mm | `(0, ∞)` finite | zero, negative, non-finite |

All three are created via `TryFrom<f64>`, returning `HeliosError::InvalidDomainValue`
on failure.

```rust
use helios_core::{EnergyMeV, HounsfieldUnit, VoxelSpacingMm};

let energy  = EnergyMeV::try_from(6.0).unwrap();   // 6 MV photon beam
let water   = HounsfieldUnit::try_from(0.0).unwrap(); // water = 0 HU
let spacing = VoxelSpacingMm::try_from(1.0).unwrap(); // 1 mm isotropic voxels

assert_eq!(energy.get(), 6.0);
assert_eq!(energy.to_string(), "6 MeV");
```

## Design Rationale

The typestate pattern ensures that:

1. **No domain violation at construction** — out-of-range values fail early with
   a structured `HeliosError`, not a panic deep inside a solver.
2. **Units are self-documenting** — `fn dose(energy: EnergyMeV, hu: HounsfieldUnit)`
   is unambiguous; raw floats are not accepted.
3. **Zero-cost** — each type wraps a single `f64`; `get()` is a transparent accessor.

## Further Reading

- [`helios-core` source](../../crates/helios-core/src/)
- [Example: Validating Foundation Units](examples/validate_foundation_units.md)
