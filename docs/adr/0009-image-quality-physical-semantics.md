# ADR 0009: Partition image-quality physical semantics

- Status: Accepted
- Date: 2026-07-23
- Class: [arch] [minor]

## Context

`helios-analysis::image_quality` was authored for MVCT intensity volumes, but
the clinical validation example also used its raw ROI statistics for dose
volumes. The dense `Volume<T>` representation must remain scalar for imaging
and solver kernels, so the physical semantic has to be expressed at the
analysis boundary rather than by changing storage.

## Decision

Keep `roi_statistics` and `volume_rmse` as raw scalar APIs for MVCT intensity
metrics. Add `dose_roi_statistics` and `dose_volume_rmse`, returning
`AbsorbedDose<T>` for dose-valued mean, standard deviation, and RMSE results.
Both paths share the same private value kernels. Michelson contrast and
contrast-to-noise remain scalar because their results are dimensionless; dose
callers convert typed dose values to the scalar boundary only for those
dimensionless equations.

The clinical validation example uses the dose-semantic ROI API and displays
values through the Gray unit. No compatibility wrappers or duplicated numeric
kernels are retained.

## Rejected alternative

Changing `Volume<T>` to store `AbsorbedDose<T>` would add unit metadata to every
dense voxel and couple imaging, solver, and storage code to one field semantic.
Keeping raw ROI results for dose callers would leave the same physical-unit
ambiguity that this audit identified.

## Verification

- Shared ROI statistics and RMSE value kernels preserve existing MVCT results.
- Dose ROI mean/std and dose RMSE tests assert `Gray` values.
- f32 instantiation covers the typed dose RMSE path.
- The clinical example and book API index use the dose-semantic surface.
