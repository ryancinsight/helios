# Helios Checklist (tactical)

**Sprint target version:** `0.0.1` (Foundation)
**Current phase:** Phase 1 → Phase 2 boundary (Foundation established; entering
Execution). Sprint 1 goal: workspace skeleton + `helios-core` + ritk/gaia domain
integration.

## Owner: claude-helios

### In-flight item: H-003 `helios-math` Scalar seam + geometry primitives — `todo`

Decomposed plan (each step has an observable completion condition):

1. [ ] Add `crates/helios-math` to workspace `members`; manifest depends on
   `hermes-simd`, `leto`, `num-traits`. — *builds empty lib.*
2. [ ] Define sealed `Scalar` trait (assoc `Accumulator`, native ops, `From`
   bridges for `f32`/`f64`) re-exporting/adapting hermes' scalar surface. — *impls
   for `f32`,`f64` compile; value-semantic ops test green.*
3. [ ] `Vec3<T: Scalar>`, `AffineTransform<T>` (patient/beam frames), ray/AABB
   intersection (Smits' slab method). — *analytical intersection tests
   (known ray/box) green, epsilon derived.*
4. [ ] clippy `-D warnings`, fmt, nextest, doctests green; sync artifacts. — *gate
   clean; CHANGELOG + this file updated.*

### Completed this sprint

- [x] **H-001** Workspace skeleton (Cargo.toml edition 2021/resolver 2,
  rust-toolchain, `.config/nextest.toml` 30s/60s budget, `.gitignore`) + Foundation
  artifacts (README, ARCHITECTURE with Atlas dependency map, backlog, gap_audit,
  CHANGELOG, SPRINT_1).
- [x] **H-002** `helios-core`: `HeliosError` (thiserror, `#[non_exhaustive]`),
  CODATA/ICRU physical constants with derivation tests, validating newtypes
  (`EnergyMeV`, `HounsfieldUnit`, `VoxelSpacingMm`). 13 tests pass; build + clippy
  `-D warnings` + fmt + nextest green.

## Gate status (last run, H-001/H-002)

| Gate | Result |
|------|--------|
| `cargo build` | pass (18.1 s cold) |
| `cargo clippy --all-targets --all-features -D warnings` | pass, 0 warnings |
| `cargo fmt --check` | pass |
| `cargo nextest run` | 13 passed / 0 failed (0.19 s) |

## Decision log (this sprint)

- **Edition 2021 / resolver 2** chosen over the edition-2024 default heuristic:
  explicit user directive in the goal + "exact kwavers/cfdrs process" (kwavers uses
  resolver 2). Recorded override of the standards default.
- **`helios-core` constants are `f64`** at their definition boundary (not generic
  over `Scalar`): the generic numeric seam lives in `helios-math` (H-003); constants
  are literals converted by callers. Avoids a premature `Scalar` dependency in the
  foundation crate.
- **No speculative empty crates:** only `helios-core` is a workspace member; the
  remaining 10 crates are added when their layer is built (architecture_scoping
  growth triggers). `workspace.dependencies` declares the full Atlas set now as the
  integration SSOT.
