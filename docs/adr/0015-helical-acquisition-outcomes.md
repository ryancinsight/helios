# ADR 0015: Type helical acquisition outcomes

- Status: Accepted
- Date: 2026-07-31
- Scope: `helios-simulation::HelicalProjection`

## Context

The public helical-acquisition result carried optical depth and transmitted
fraction as raw `T` values even though both outcomes are dimensionless and the
same simulation already uses Aequitas for geometry and Hyperion for the
Beer–Lambert law. This left the result contract weaker than its inputs and
allowed dimensionless values to be confused with arbitrary scalar outputs.

## Decision

Return Aequitas `Dimensionless<T>` for `HelicalProjection::optical_depth` and
`HelicalProjection::transmission`. Construct optical depth at the ray-marching
scalar boundary, pass it through Hyperion `OpticalDepth`, and retain the typed
transmission result through the public projection record. Dense attenuation
fields and ray-kernel coordinates remain scalar representations.

The contract is real-valued under Eunomia `RealField`. No imaginary-unit
quantity applies to the real Beer–Lambert acquisition law; Eunomia complex
values remain available to numerical domains that genuinely require them.

## Alternatives rejected

- Retain raw `T` result fields: rejected because the public dimensionless
  contract remains implicit.
- Add parallel typed accessors: rejected because two public representations
  duplicate the result contract and preserve the weaker path.
- Introduce a complex-valued optical-depth unit: rejected because optical depth
  and transmission are real-valued in this acquisition model.

## Consequences

This is a pre-1.0 breaking change for callers that compare or serialize the
projection outcomes. Callers extract `into_base()` only at scalar assertions,
serialization, or numerical boundaries.

## Verification

`cargo check -p helios-simulation --tests --offline` passes. Nextest run
`50ccc92b-b632-4e66-8800-1f2fee63fe77` passes 42/42; warning-denied Clippy,
doctests, Rustdoc, and targeted Rustfmt pass. The repository `gap_audit.md`
contains peer-owned dirty changes and is not modified in this increment; the
H-096 backlog entry records that reconciliation requirement.
