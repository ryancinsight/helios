# Chapter 37 — Migration Validation: TG-119 and Atlas Parity

Helios's Atlas migration is governed by a **TG-119/I IAEA-style parity
test harness** that runs both the legacy and Atlas ports of each solver
on the same input and demands agreement within a documented tolerance.
Until a solver passes parity, it does not graduate from the legacy crate
to the Atlas-only path.

## Validation Reference Set

| Phantom / Scenario | Tolerance | Atlas path |
|---|---|---|
| IAEA 2-D benchmark | 1.5 Gy dose ±2% | `helios-solver::iaea_2d` |
| TG-119 C-shape | DVH-of-target ±5% | `helios-analysis::gamma_index` |
| AAPM TG-244 heterogeneity | 3%/3 mm γ pass-rate 95% | `helios-analysis::gamma` |
| Chest-wall Monte Carlo | cross-section within 5% | `helios-physics::compton.rs` |
| TomoTherapy helical standard | gamma pass-rate ≥ 95% | `helios-simulation` |
| MVCT registration | mutual-info ≥ 95% (modal) | `helios-imaging` |

## The Parity Test Anatomy

A TG-119 parity test has three components:

1. **Reference input**. A canonical phantom (TG-119 C-shape, IAEA 2-D,
   AAPM heterogeneity).  Inputs are version-pinned in
   [`helios-validation`] data fixtures.
2. **Legacy baseline**. The pre-migration sweep, run via the legacy
   `ndarray`/`nalgebra`-based crate.  Outputs are checked into
   `validation_reports/` with cryptographic hashes.
3. **Atlas candidate**. The Atlas port running the **exact same** input
   through the new Atlas stack.  Outputs are byte-compared against the
   baseline hash under a tolerance budget.

## Gamma-Index Pass Rate

The gamma index is the canonical helios lossy-comparison metric:

  γ(r_e, r_r) = sqrt((Δr/DTA)^2 + (ΔD/DD)^2)

A pass rate above 95% for `(3 mm, 3%)` is the gating metric for the
Plan-Verification workflow.

## Per-Crate Migration Reports

Each migration produces a `validation_reports/migration_<crate>.md`
document recording:

- The parity-test scenarios used.
- The legacy vs. Atlas timings.
- The gamma pass-rate and observed error.
- The Atlas crates involved.
- The legacy crates that can now be removed.

These reports are filed in
[`validation_reports`](../../validation_reports/) and are referenced from
the Atlas Dependency Map.

## Cleanup Gate (Helios)

Markers that the migration is complete and legacy crates can be removed:

- [ ] `cargo tree | grep -E '(ndarray|nalgebra|rayon|tokio|burn|ch)'` is empty.
- [ ] TG-119 C-shape passes gamma ≥ 95%.
- [ ] IAEA 2-D passes 2% dose tolerance.
- [ ] AAPM TG-244 passes 95%.
- [ ] TomoTherapy helical standard passes gamma ≥ 95%.
- [ ] No `unsafe` SIMD intrinsics remain outside `hermes-simd`.
- [ ] No `unsafe { Vec::from_raw_parts ... }` remains outside
      `mnemosyne`.

## Validation Examples

- [`validation_regression`](examples/validation_regression.md) —
  analytical dose-volume regression tests.
- [`validation_clinical`](examples/validation_clinical.md) — TG-119
  clinical protocol validation.
- [`gamma_index`](examples/gamma_index.md) — gamma-index dose comparison.
- [`tomotherapy_workflow`](examples/tomotherapy_workflow.md) — helical
  end-to-end.

## Further Reading

- [`helios-validation`](../../crates/) — parity harness source.
- [Atlas Dependency Map](appendix_dependencies.md)
- [Migration Overview](migration_overview.md)
