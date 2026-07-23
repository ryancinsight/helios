# ADR 0007: Typed gamma physical criteria

- Status: accepted
- Date: 2026-07-23
- Class: [arch] [major]

## Context

`helios-analysis::gamma_index_3d` and its local-normalization variant accepted
distance-to-agreement, normalization dose, low-dose cutoff, and search radius
as raw scalar values. The gamma result is dimensionless, but those criteria are
physical quantities and the scalar boundary allowed incompatible units to be
combined. `gamma_pass_rate` likewise accepted a raw scalar dose threshold.

## Decision

- Represent distance-to-agreement and search radius as Aequitas `Length<T>`.
- Represent global normalization dose, local low-dose cutoff, and pass-rate
  dose thresholds as Aequitas `AbsorbedDose<T>`.
- Keep the fractional dose-difference criterion, gamma field, and pass rate
  scalar because they are dimensionless.
- Convert typed criteria to base scalars once at the numerical-kernel boundary;
  retain the existing Low gamma equations, search neighborhood, validation, and
  scalar `Volume<T>` result storage.

All in-tree callers, examples, tests, and book snippets migrate in the same
change. No scalar overloads or forwarding wrappers remain.

## Alternatives rejected

- Keep suffix-named scalar parameters such as `dta_mm`: rejected because the
  unit remains a naming convention rather than a type invariant.
- Return a typed gamma field: rejected because gamma is a dimensionless metric
  and dense scalar volume storage is the existing analysis boundary.
- Introduce a second typed API beside the scalar API: rejected because it keeps
  two contracts and preserves the obsolete boundary.

## Consequences

This is a pre-1.0 public breaking change. Callers must construct physical
criteria with Aequitas quantities, while the numerical result and dimensionless
acceptance statistics retain their existing scalar contracts. The public
analysis boundary now prevents distance/dose unit confusion at compile time.
