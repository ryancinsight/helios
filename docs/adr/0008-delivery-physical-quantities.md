# ADR 0008: Typed delivery physical quantities

- Status: accepted
- Date: 2026-07-23
- Class: [arch] [major]

## Context

The helical delivery boundary represented leaf fluence and delivery geometry as
unvalidated scalars. Portal fluence was typed only inside the transmission
calculation and then converted back, while dose accumulation accepted
millimetre-suffixed scalar arguments. The names documented units but did not
make unit mismatches impossible for callers.

## Decision

- Represent per-leaf and total delivered fluence as Aequitas
  `EnergyPerArea<T>`.
- Represent couch position, leaf width, ray sampling step, parallel standoff,
  and point-source axis as Aequitas `Length<T>`.
- Keep gantry angles dimensionless and keep dense `Volume<T>` dose and
  attenuation storage scalar for the existing field-kernel contract.
- Convert typed distances to millimetres once at the ray/voxel numerical
  boundary. `VoxelGrid` coordinates and Hyperion ray deposition remain in their
  established millimetre representation.
- Retain typed fluence through delivery, collimation, portal transmission, and
  the public dose-deposition boundary; convert only where the existing scalar
  provider kernel requires it.

All in-tree callers, examples, tests, and book snippets migrate in the same
change. No scalar overloads, suffix-renamed compatibility fields, or forwarding
wrappers remain.

## Alternatives rejected

- Keep `_mm` scalar parameters: rejected because unit names do not enforce a
  common length unit or distinguish length from another scalar.
- Type only portal internals: rejected because the public delivery and portal
  boundaries would still erase fluence semantics between operations.
- Make dense voxel storage a quantity-valued array: rejected because it would
  change the established field-kernel layout without a field-descriptor
  contract; the dimensional boundary is preserved at the API instead.

## Consequences

This is a pre-1.0 public breaking change. Callers construct physical delivery
criteria with Aequitas quantities, and the numerical kernel receives an
explicitly documented millimetre conversion. The delivery, portal, and dose
APIs now prevent accidental unit mixing while retaining the existing numerical
algorithms and scalar field layout.
