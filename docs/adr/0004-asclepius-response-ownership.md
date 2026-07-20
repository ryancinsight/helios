# ADR 0004: Asclepius owns biological-response laws

- Status: accepted
- Date: 2026-07-20
- Class: architectural, public breaking

## Context

Helios analysis contained local generalized equivalent uniform dose, logistic
tumour-control probability, and Lyman complication probability equations.
Helios planning independently rebuilt generalized equivalent uniform dose from
Coeus primitives. Asclepius now owns these biological-response laws and its
Coeus tape construction over Aequitas physical quantities.

Keeping either Helios implementation would create two authorities for the same
equations, validation domains, numerical stabilization, and proofs.

## Decision

`helios-analysis::Dvh` retains structure sampling and histogram behavior. Its
sorted response sample is stored once as `AbsorbedDose<T>` and borrowed directly
by Asclepius. Response methods validate their parameters through Asclepius and
return its typed failures.

`helios-planning` retains dose-influence application, penalty selection, and
objective construction. It passes the Coeus dose variable to
`asclepius-coeus` for differentiable gEUD construction.

Proteus, Asclepius, and Helios resolve the same merged Aequitas revision. The
dependency direction is:

```text
aequitas -> asclepius -> helios-analysis
                   \-> asclepius-coeus -> helios-planning
aequitas -> proteus -> helios-physics
```

No forwarding free functions or compatibility module remain.

## Migration

- Replace `helios_analysis::generalized_eud`, `tcp_logistic`, and `ntcp_lkb`
  calls with the corresponding Asclepius response model when no DVH exists.
- For structure analysis, call the natural `Dvh` receiver method and propagate
  or handle its `ResponseError` result.
- Treat `Dvh::dose_sample` as `&[AbsorbedDose<T>]`; use `as_base()` only at an
  external scalar boundary.
- Existing `EudPenalty` callers keep their Helios API. Enabling `autodiff` now
  selects the Asclepius Coeus adapter internally.

## Rejected alternative

Wrapping the Helios equations around Asclepius would preserve obsolete names
and duplicate the validation boundary. Converting `&[T]` to a temporary
`Vec<AbsorbedDose<T>>` on each response evaluation would add allocation and
copying to a repeated analysis path.

## Consequences

The three `Dvh` response methods become fallible. The removed public free
functions are a pre-1.0 breaking change. Dose histogram scalar accessors retain
their existing signatures, while the response sample accessor exposes the
dimensionally typed borrowed slice.

Asclepius theorems and proofs are the durable law specification. Helios tests
provide consumer evidence: direct-law equality over the identical borrowed
sample, typed failure propagation, analytical gEUD value, finite-difference
gradient agreement, and end-to-end PTV/OAR outcome ordering.
