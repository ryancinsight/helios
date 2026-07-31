# ADR 0014: Type scatter-kernel physical inputs

- Status: Accepted
- Date: 2026-07-31
- Scope: `helios-solver` deposition/scatter kernels and
  `helios-simulation::CollapsedCone`

## Context

The dose-deposition kernels already depended on Aequitas elsewhere in Helios,
but their public scatter ranges, voxel sampling pitches, and spectral weights
were raw `T` values. Parameter names such as `range_cm` and `voxel_cm` carried
unit meaning only in prose, so a caller could pass a value expressed in a
different unit system without a type error. The oriented scatter seam repeated
the same problem for its physical sampling step.

## Decision

Accept Aequitas `Length<T>` for all public deposition ranges and sampling
pitches, and `Dimensionless<T>` for poly-energetic relative spectral weights.
Rename the public fields and parameters to semantic names without unit suffixes.
Convert lengths to the solver's centimetre or millimetre scalar coordinate only
at the exponential formula or voxel-mesh boundary. Store the collapsed-cone
sampling step as `Length<T>` until it reaches oriented voxel sampling.

Keep kernel taps, dense TERMA/dose volumes, indices, and trigonometric/vector
coordinates as solver-native scalar representations. These are representation
or formula boundaries, not missing physical metrics. The kernels use Eunomia's
real scalar traits; no imaginary or complex-valued unit is appropriate for this
real transport law.

## Alternatives rejected

- Retain unit-bearing raw scalar parameters: rejected because names do not
  enforce dimensional correctness.
- Add parallel raw and typed overloads: rejected because a compatibility path
  would preserve the unsafe public contract and duplicate the API.
- Store both typed and scalar lengths: rejected because duplicated
  representations can diverge; scalar extraction remains local to the formula
  or mesh boundary.

## Consequences

This is a pre-1.0 public breaking change. Existing callers construct Aequitas
length and dimensionless quantities, while the numerical kernels retain their
native precision and existing normalized tap/value semantics.

## Verification

`cargo check -p helios-solver -p helios-simulation --all-targets --offline`
passes after migrating package tests, examples, and the scatter benchmark.
Focused Nextest, Clippy, doctest, Rustdoc, and full locked-workspace gates are
still required for closure; the shared dirty Cargo manifest/lockfile and peer
workspace changes are external to this slice.
