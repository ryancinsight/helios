# ADR 0013: Typed GPU attenuation inputs

- Status: accepted
- Date: 2026-07-27

## Context

The CPU CT-number attenuation path already accepted typed mass attenuation and
water-density quantities, but `helios-gpu::GpuAttenuationMapper::new` accepted
the same physical values as raw `f32` parameters named in cm-based units. This
allowed a coefficient or density from another unit system to cross the public
GPU boundary without dimensional checking.

## Decision

Accept Aequitas `AreaPerMass<f32>` and `MassDensity<f32>` in
`GpuAttenuationMapper::new`. Convert those quantities to
`SquareCentimeterPerGram` and `GramPerCubicCentimeter` once at the explicit
formula/GPU staging boundary because the fused kernel's scale and offset are
defined in cm-based units. Keep the HU buffer and attenuation output as dense
`f32` numerical storage, and do not retain scalar compatibility overloads.

## Consequences

This is a pre-1.0 public breaking change. Unit mismatches are rejected at the
constructor boundary while the GPU ABI and closed-form attenuation law remain
unchanged.

## Verification

The `helios-gpu` package check with tests/examples passes. Nextest passes 10/10,
including CPU closed-form and solver differential tests. Warning-denied
Clippy, doctests, Rustdoc, rustfmt, and diff checks pass. Existing shared-graph
unused-patch warnings are outside this change.
