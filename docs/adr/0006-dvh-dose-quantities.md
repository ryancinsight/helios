# ADR 0006: Typed dose-valued DVH metrics

- Status: accepted
- Date: 2026-07-23
- Class: [arch] [major]

## Context

`helios-analysis::Dvh` already stores its sampled dose as Aequitas
`AbsorbedDose<T>`, but its public extrema, mean, quantile, and radiobiology
threshold APIs converted those values back to raw `T`. That boundary allowed a
dose to be compared with an unrelated scalar and forced every caller to
reconstruct the physical meaning of a result.

## Decision

- `Dvh::min`, `max`, `mean`, and `dose_at_volume_fraction` return
  `AbsorbedDose<T>`.
- `Dvh::volume_fraction_at_dose` accepts `AbsorbedDose<T>`.
- `Dvh::generalized_eud` returns `AbsorbedDose<T>`.
- `Dvh::tcp_logistic` and `ntcp_lkb` accept typed TCD50/TD50 doses.
- Dimensionless volume fractions, homogeneity indices, probabilities, and
  scalar gamma-engine boundaries remain raw values because they are not dose
  quantities.
- Voxel storage remains the existing scalar `Volume<T>` boundary; conversion
  to `AbsorbedDose<T>` occurs once while constructing the borrowed DVH sample.

All in-tree callers migrate in the same change. Scalar conversion is explicit
with `into_base()` or a unit conversion such as `in_unit::<Gray>()`; no
forwarding overloads or compatibility wrappers remain.

## Alternatives rejected

- Keep scalar accessors and rely on naming/documentation: rejected because the
  type system would still permit dimensionally invalid dose operations.
- Add parallel typed methods beside the scalar API: rejected because it keeps
  two public contracts and preserves the obsolete boundary.
- Replace the scalar `Volume<T>` dose field: rejected because dense field
  storage is a domain boundary used by transport and imaging kernels; the DVH
  owns the typed analysis view.

## Consequences

This is a pre-1.0 public breaking change. Clinical dose metrics now compose
directly with Aequitas and Asclepius response quantities, while fractions,
indices, and probabilities retain their scalar contracts. Examples, tests, and
book snippets show explicit unit conversion at presentation boundaries.
