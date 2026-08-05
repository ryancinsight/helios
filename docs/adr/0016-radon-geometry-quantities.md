# ADR 0016: Type Radon imaging geometry quantities

- Status: Accepted
- Date: 2026-08-02
- Scope: `helios-imaging::Sinogram`, Radon, FBP, and SIRT geometry boundaries

## Context

The public parallel-beam imaging API represented projection angles, detector
offsets, source distance, and ray-march step as the scalar type `T`, although
the contract documented radians and millimetres. This allowed angle, length,
and reconstruction values to be interchanged at call sites and left the
imaging boundary weaker than the typed delivery and attenuation boundaries.

## Decision

Carry projection angles as Aequitas `Angle<T>` and detector offsets, source
distance, and ray-march step as `Length<T>` through `Sinogram`, Radon, FBP, SIRT,
noise, examples, and workflow validation. Convert to radians only for scalar
trigonometry and to millimetres or centimetres only at the existing mesh and
filter formula boundaries. The dense voxel grid and line-integral readings
remain scalar storage because their representations are kernel data, not
public physical metadata.

The imaging contract is real-valued under Eunomia's `UnitScalar` seam. It has
no phasor or Fourier-valued physical field, so no imaginary unit or complex
length extension is introduced. Eunomia complex values remain available for
numerical domains whose public law genuinely carries a complex phasor under one
physical unit.

## Alternatives rejected

- Retain raw `T` geometry: rejected because documented units would remain
  unenforced at the public boundary.
- Add typed accessors beside raw constructors: rejected because it preserves
  the weaker representation and creates two geometry contracts.
- Add a complex or imaginary length unit: rejected because Radon geometry and
  the reconstruction laws are real-valued.

## Consequences

This is a pre-1.0 breaking change. All in-tree callers construct Aequitas
quantities directly. The `eunomia::UnitScalar` bound is explicit on kernels
that convert quantities to formula units; no compatibility adapter or scalar
constructor remains.

## Verification

The standalone locked imaging check, warning-denied Clippy, doctests, Rustdoc,
touched-file Rustfmt, and analysis validation example check pass. Nextest run
`26905c71-03f1-447a-be77-df0c84278c3c` passes 33/33 imaging tests. The
downstream simulation gate remains blocked by provider-owned Mnemosyne/Moirai
diagnostics recorded in `gap_audit.md`; that blocker does not alter this
contract's value-semantic evidence.
