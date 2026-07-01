# Helios Checklist (tactical)

**Sprint target version:** `0.0.1` (Foundation → Sprint 2 physics/GPU)
**Current phase:** Phase 2 (Execution). Sprint 1 domain core complete
(`helios-core`/`math`/`domain`); Sprint 2 opened with `helios-physics`.

## Owner: claude-helios

### In-flight item: H-010 `helios-gpu` — hephaestus `ComputeDevice` foundation — `todo (DoR)`

Decomposed plan (each step has an observable completion condition):

1. [ ] Create `crates/helios-gpu`; depend on `hephaestus-core` (+ `hephaestus-wgpu`
   behind a `wgpu` feature). Program against `hephaestus_core::ComputeDevice` — do
   **not** reinvent a device trait. — *builds; feature matrix compiles.*
2. [ ] Runtime backend selection (wgpu if an adapter is present, else CPU
   reference), surfacing the choice + reason (capability detection, not a silent
   stub). — *selection logged; no `todo!()` behind a feature.*
3. [ ] A first real kernel (e.g. per-voxel `μ` map: HU→density→μ) with a CPU
   reference; differential test CPU vs GPU (epsilon by reduction order) when an
   adapter exists, else CPU-reference value-semantic test. — *bitwise/epsilon match.*
4. [ ] clippy `-D warnings`, fmt, nextest, doctests green; sync artifacts.

*Risk:* wgpu is a heavy first compile and GPU-adapter availability in the test
env is unverified — the CPU-reference path keeps the increment green regardless.

*Sequenced/blocked:* H-011b (NIST μ/ρ data), H-011c segment-generation (gaia G-11),
H-004b (ritk DICOM, heavy build).

### Completed

- [x] **H-011c (reduction)** `helios-physics::projection`: `optical_depth`
  (τ=Σμᵢ·Lᵢ) + `beam_transmission` (exp(−τ)) over `(μ,length)` segments. 5 tests
  (homogeneous=μ·L oracle, additivity, multiplicative composition, empty, f32).
- [x] **H-011** `helios-physics`: attenuation relations + HU→density (9 tests).
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

## Gate status (last run, H-011c reduction)

| Gate | Result |
|------|--------|
| `cargo build` | pass |
| `cargo clippy --all-targets --all-features -D warnings` | pass, 0 warnings |
| `cargo fmt --check` | pass |
| `cargo nextest run` | 39 passed / 0 failed (0.34 s) |
| `cargo test --doc` | pass |

## Decision log (Sprint 2)

- **Program against `hephaestus_core::ComputeDevice`, don't reinvent a Helios
  `ComputeBackend`.** hephaestus-core is the Atlas GPU-agnostic device seam that
  apollo/coeus already target; Helios adds GPU *utilities* (buffer helpers, kernel
  cache, backend selection) on top, per DIP + upstream ownership. Avoids a parallel
  device abstraction (anti-duplication).
- **Projector split (physics vs geometry).** The line-integral *reduction*
  (`optical_depth`/`beam_transmission`) is geometry-free and implemented/tested now;
  *segment generation* (voxel DDA/Siddon) needs gaia `Ray`/`Aabb` and is sequenced
  behind G-11. Same reduction will run over CPU- or GPU-generated segments unchanged.

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
