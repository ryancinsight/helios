# ADR 0017: Type inverse-planning dose objectives

- Status: Accepted
- Date: 2026-08-05
- Scope: `helios-planning` autodiff and `helios-analysis::Dvh` gEUD contracts

## Context

The inverse-planning autodiff API accepted clinical dose floors, ceilings, and
gEUD reference values as `f64`. The shared `helios-analysis::Dvh` gEUD, TCP,
and NTCP entry points also accepted the volume-effect parameter `a` as an
untyped scalar. The implementations computed real dose values, but the public
contracts allowed a caller to pass an unrelated scalar metric.

## Decision

Carry DVH penalty bands and gEUD references as Aequitas
`AbsorbedDose<f64>` values. Carry the gEUD volume-effect parameter as
`Dimensionless<f64>`. Extract the real base values once at the Coeus or
Asclepius formula boundary, where the numerical representation is scalar
storage. Keep beamlet weights, penalty coefficients, response slopes, and dense
influence data scalar because their units are defined by the optimization or
response model rather than a fixed SI quantity.

The contract remains real-valued under Eunomia. No phasor crosses this planning
boundary, so no imaginary dose unit or complex physical wrapper is introduced.

## Alternatives rejected

- Keep raw dose scalars: rejected because dose semantics remain unenforced.
- Add parallel typed and raw constructors: rejected because it preserves two
  contracts and requires compatibility plumbing.
- Assign an imaginary or complex dose unit: rejected because the objective is a
  real dose-response law.

## Consequences

This is a pre-1.0 breaking change. All in-tree callers construct typed dose and
dimensionless values directly. The numerical autodiff tensor and Asclepius
response constructors remain unchanged and own the scalar extraction
boundaries.

## Verification

The focused planning and analysis package checks with all features pass against
the live Eunomia 0.8/Aequitas graph. Value-semantic tests cover the independent
gEUD oracle, finite-difference gradient, zero-violation band, DVH optimizer,
and the shared DVH/clinical/end-to-end callers. The repository's locked hosted
matrix remains the final source-identity gate; the local Atlas overlay cannot
execute the standalone lock because its path patches intentionally replace the
committed Git sources.
