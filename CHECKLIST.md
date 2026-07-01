# Helios Checklist (tactical)

**Sprint target version:** `0.0.1` (Foundation)
**Current phase:** Phase 1 → Phase 2 boundary (Foundation established; entering
Execution). Sprint 1 goal: workspace skeleton + `helios-core` + ritk/gaia domain
integration.

## Owner: claude-helios

### In-flight item: H-004 `helios-domain` CT/MVCT volume + voxel-grid geometry — `todo`

Decomposed plan (each step has an observable completion condition):

1. [ ] Add `crates/helios-domain` to workspace `members`; depends on `helios-core`,
   `helios-math`, and (later) `ritk-io` for DICOM. — *builds empty lib.*
2. [ ] `VoxelGrid` geometry: dimensions (const/runtime), `VoxelSpacingMm` per axis,
   origin (`Point3`), index↔world affine mapping via `helios-math` `Isometry3`.
   — *round-trip index→world→index tests exact; out-of-bounds rejected.*
3. [ ] `Volume<T: Scalar>`: dense scalar field over a `VoxelGrid` (leto `Array3`),
   trilinear sample; `CtVolume` = `Volume` of `HounsfieldUnit`-validated data.
   — *trilinear sample matches analytical linear field; boundary clamped.*
4. [ ] Verify `ritk-io` DICOM surface (anti-hallucination) and load path in a
   feature-gated module; CPU-only test with a synthetic in-memory volume.
   — *symbols confirmed via source/doc before use.*
5. [ ] clippy `-D warnings`, fmt, nextest, doctests green; sync artifacts.

### Completed this sprint

- [x] **H-003** `helios-math`: `Scalar = eunomia::RealField` seam re-export; leto
  geometry (`Vector3`/`Point3`/`Isometry3`) re-export; Helios-owned `Ray`/`Aabb` +
  slab `intersect_ray` (voxel-traversal primitive absent upstream). 6 analytical
  tests (axis-aligned, miss-behind, parallel-miss, origin-inside, diagonal, f32
  generic). Worked around leto→mnemosyne→themis skew (G-10) via
  `default-features=false`.

- [x] **H-001** Workspace skeleton (Cargo.toml edition 2021/resolver 2,
  rust-toolchain, `.config/nextest.toml` 30s/60s budget, `.gitignore`) + Foundation
  artifacts (README, ARCHITECTURE with Atlas dependency map, backlog, gap_audit,
  CHANGELOG, SPRINT_1).
- [x] **H-002** `helios-core`: `HeliosError` (thiserror, `#[non_exhaustive]`),
  CODATA/ICRU physical constants with derivation tests, validating newtypes
  (`EnergyMeV`, `HounsfieldUnit`, `VoxelSpacingMm`). 13 tests pass; build + clippy
  `-D warnings` + fmt + nextest green.

## Gate status (last run, H-003)

| Gate | Result |
|------|--------|
| `cargo build` | pass (leto/eunomia git deps compiled) |
| `cargo clippy --all-targets --all-features -D warnings` | pass, 0 warnings |
| `cargo fmt --check` | pass |
| `cargo nextest run` | 19 passed / 0 failed (0.19 s) |
| `cargo test --doc` | pass (0 doctests) |

## Decision log (this sprint)

- **Scalar seam = `eunomia::RealField`, geometry from `leto`** (H-003): eunomia is
  the Atlas datatype SSOT (`RealField`/`FloatElement`/`NumericElement`) and leto
  owns `Vector3`/`Point3`/`Isometry3`. `helios-math` re-exports them rather than
  reinventing (consolidation/subtractive bias). Only `Ray`/`Aabb` + slab
  intersection are Helios-owned — absent upstream (leto is an array library).
- **leto `default-features = false`** (G-10): leto's default `mnemosyne-memory`
  pulls an mnemosyne rev bound to `themis ^0.8`, conflicting with themis HEAD 0.9.x.
  Consuming leto with only `std` sidesteps the skew; mnemosyne placement is opted
  into at the layer that needs it. Upstream fix filed as G-10.

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
