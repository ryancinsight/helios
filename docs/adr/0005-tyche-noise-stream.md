# ADR 0005: Tyche-owned quantum-noise stream

- Status: Accepted
- Date: 2026-07-21
- Class: [arch] [major]

## Context

`helios-imaging` delegates standard-normal sampling to Tyche, but the initial
consumer used Tyche's retired implicit counter algorithm and a local path patch
masked its declared Git revision. A persisted seed therefore did not identify
the complete replay schedule.

## Decision

- Resolve Tyche's merged default-branch revision directly, without a local
  patch override.
- Select `SplitMix64` in `StandardNormal<f64, SplitMix64>` so the counter
  algorithm is part of the monomorphized type.
- Treat Tyche's stream version, algorithm, seed, row-major sample index, and
  draw index as the complete noise-field replay identity.
- Keep photon statistics and the generic sinogram mapping in Helios; Tyche owns
  only counter scheduling and standard-normal transformation.

## Alternatives rejected

- Preserve the retired values with a Helios compatibility sampler: rejected
  because it would restore a second PRNG and Box-Muller owner.
- Retain the local patch: rejected because it makes `Cargo.toml` and
  `Cargo.lock` non-authoritative outside the Atlas checkout.
- Add dynamic algorithm dispatch: rejected because this path has one required
  schedule and static selection has no runtime storage or vtable cost.

## Consequences

Seeded quantum-noise readings intentionally change once. The exact consumer
vector pins the new schedule, while the analytical variance, attenuation
monotonicity, high-flux limit, and deterministic-replay tests verify physical
behavior independently of that vector. The sampling call remains allocation
free and counter-addressed; the existing sinogram result allocation is the
owned output required by the public API.
