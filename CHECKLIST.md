# Helios Checklist (tactical)

**Sprint target version:** `0.0.1` (Foundation → Sprint 2 physics/GPU)
**Current phase:** Phase 2 (Execution). Sprint 1 domain core complete
(`helios-core`/`math`/`domain`); Sprint 2 opened with `helios-physics`.

## Owner: claude-helios

### In-flight item: H-011b `helios-physics` NIST μ/ρ data tables — `todo`

Decomposed plan (each step has an observable completion condition):

1. [ ] Define an energy-indexed `MassAttenuation` table type + interpolation
   (log-log in energy, the NIST-standard scheme). — *interpolation reproduces
   tabulated node values exactly; monotone between nodes.*
2. [ ] Embed a small citable NIST XCOM μ/ρ dataset for water (and air) over the
   MV/kV energy range, sourced from the NIST database (verified, not memorized).
   — *values traceable to the cited table.*
3. [ ] `Material` → `MassAttenuation(E)` lookup; mixture rule (mass-weighted
   Σ wᵢ (μ/ρ)ᵢ). — *water vs mixture value-semantic checks.*
4. [ ] clippy `-D warnings`, fmt, nextest, doctests green; sync artifacts.

*Parallel-unblock note:* H-011c (ray-marched ∫μ dl) and H-004b (ritk DICOM) remain
sequenced — H-011c behind gaia geometry (G-11), H-004b behind a heavy ritk build.

### Completed

- [x] **H-011** `helios-physics`: `LinearAttenuation`/`MassAttenuation` newtypes,
  Beer–Lambert `transmission`, `half_value_layer`, `μ=(μ/ρ)·ρ`, HU→density
  calibration. 9 analytical tests. Closes G-1 (attenuation relations).
- [x] **H-004** `helios-domain`: `VoxelGrid` + `Volume` trilinear (see SPRINT_1).

### Completed (Sprint 1)

- [x] **H-004** `helios-domain`: `VoxelGrid<T>` (dims, per-axis spacing, leto
  `Isometry3` pose; `index_to_world`/`world_to_index`/`voxel_center`) + `Volume<T>`
  backed by leto `Array3` with `sample_trilinear`/`sample_world`. 11 tests: affine-
  field exact-reproduction oracle, C-contiguous layout lock, identity + 90°-rotated
  pose round-trips, out-of-bounds/NaN → None, f32 genericity.
- [x] **H-003** `helios-math`: `Scalar = eunomia::RealField` seam + leto substrate
  re-export (geometry primitives corrected to gaia ownership; local `Ray`/`Aabb`
  removed — see decision log). Worked around leto→mnemosyne→themis skew (G-10) via
  `default-features=false`.

- [x] **H-001** Workspace skeleton (Cargo.toml edition 2021/resolver 2,
  rust-toolchain, `.config/nextest.toml` 30s/60s budget, `.gitignore`) + Foundation
  artifacts (README, ARCHITECTURE with Atlas dependency map, backlog, gap_audit,
  CHANGELOG, SPRINT_1).
- [x] **H-002** `helios-core`: `HeliosError` (thiserror, `#[non_exhaustive]`),
  CODATA/ICRU physical constants with derivation tests, validating newtypes
  (`EnergyMeV`, `HounsfieldUnit`, `VoxelSpacingMm`). 13 tests pass; build + clippy
  `-D warnings` + fmt + nextest green.

## Gate status (last run, H-011)

| Gate | Result |
|------|--------|
| `cargo build` | pass |
| `cargo clippy --all-targets --all-features -D warnings` | pass, 0 warnings |
| `cargo fmt --check` | pass |
| `cargo nextest run` | 34 passed / 0 failed (0.56 s) |
| `cargo test --doc` | pass |

## Decision log (this sprint)

- **Scalar seam = `eunomia::RealField`; substrate from `leto`** (H-003): eunomia is
  the Atlas datatype SSOT (`RealField`/`FloatElement`/`NumericElement`) and leto
  owns `Vector3`/`Point3`/`Isometry3`. `helios-math` re-exports them rather than
  reinventing (consolidation/subtractive bias).
- **Geometry primitives belong to gaia, not Helios** (correction, user directive):
  the initial `helios-math` `Ray`/`Aabb` were a downstream duplication and were
  **removed**. gaia already owns `Aabb` (default branch) and a validated `Ray` +
  `intersect_aabb` (leto-migration branch). Helios will re-export gaia's types once
  that migration lands on gaia's default branch (H-003b, blocked; G-11). Do not
  re-implement geometry in Helios.
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
