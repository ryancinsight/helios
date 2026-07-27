# ADR 0012: Type Helical Delivery Kinematic Metrics

- Status: Accepted
- Date: 2026-07-24
- Scope: `helios-domain::HelicalDelivery` and `helios-simulation`

## Context

The prior delivery quantity migration typed fluence and distance geometry but
left gantry angle, couch position, rotation time, and couch velocity as `T`
values with unit-bearing names. This left the public helical synchronization
contract open to dimension and unit confusion. Aequitas now provides a distinct
dimensionless `Angle` semantic with a canonical `Radian` unit.

## Decision

`HelicalDelivery` stores and returns Aequitas `Angle`, `Length`, `Time`,
`Velocity`, and `Dimensionless` values. `HelicalProjection` and `DeliveryFrame`
carry typed angle and couch state. The dose, portal, and acquisition kernels
convert typed values to the existing scalar millimetre/radian coordinate seam
once, immediately before trigonometry or voxel traversal. The constructor and
all in-tree callers migrate together; no raw compatibility overload remains.

## Alternatives rejected

- Keep angle as `Dimensionless`: rejected because machine rotation and plane
  wave angles are semantically distinct from scores and fractions.
- Type only the simulation output: rejected because `HelicalDelivery` would
  remain an untyped producer and reintroduce raw values at every caller.
- Store typed values in parallel with raw fields: rejected because duplicated
  representations create synchronization and ownership drift.

## Verification

The domain tests preserve full-rotation angle, wrapped-angle, couch-travel,
projection/time equivalence, monotonic couch advance, invalid-input, and f32
oracles. Simulation tests preserve Beer–Lambert central-ray values, delivery
fluence, dose/portal geometry, generic f32 paths, examples, and the end-to-end
workflow. Cargo verification is blocked before source compilation by the peer
Coeus checkout.
