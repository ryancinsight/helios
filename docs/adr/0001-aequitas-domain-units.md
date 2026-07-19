# ADR 0001: Aequitas-backed clinical quantities

- Status: Accepted
- Date: 2026-07-19
- Class: [arch]

## Context

`EnergyMeV` and `VoxelSpacingMm` validate Helios-specific clinical invariants
but also duplicate storage and conversion rules for physical energy and length.
Hounsfield units are different: HU is a calibrated imaging scale with a
Helios-owned representable range, not an SI linear unit.

## Decision

- Keep the public Helios validating newtypes as the domain boundary.
- Store `EnergyMeV` as Aequitas `Energy<f64>` and `VoxelSpacingMm` as Aequitas
  `Length<f64>`.
- Validate raw values before constructing the Aequitas quantity.
- Preserve `HounsfieldUnit` as a Helios-owned transparent `f64` newtype.
- Pin size and alignment equivalence with compile-time assertions.

## Alternatives rejected

- Replace the Helios types with Aequitas aliases: rejected because dimensional
  correctness does not encode strict positivity or the clinical error surface.
- Keep raw `f64` storage: rejected because it duplicates physical-unit
  ownership and cannot participate in Aequitas dimensional arithmetic.
- Model HU as an Aequitas unit: rejected because HU conversion depends on a
  calibration reference and is not a context-free linear SI conversion.

## Consequences

The public validation contract is unchanged, but `EnergyMeV::get` and
`VoxelSpacingMm::get` are no longer `const` because Aequitas resolves unit
conversion through trait methods. Physical conversion law has one provider.
Both domain wrappers remain one-word transparent values with no allocation or
runtime dimension metadata; conversion arithmetic occurs only when values
cross the public unit boundary.
