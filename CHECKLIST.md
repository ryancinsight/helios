# Helios Checklist (tactical)

**Sprint target version:** `0.0.1` (Foundation → Sprint 2 physics/GPU)
**Current phase:** Phase 2 (Execution). Sprint 1 domain core complete
(`helios-core`/`math`/`domain`); Sprint 2 opened with `helios-physics`.

## Owner: claude-helios

### In-flight item: H-011b `helios-physics` NIST μ/ρ data tables — `todo`

Next unblocked increment (GPU H-010 blocked on stack convergence, G-12).

1. [ ] Energy-indexed `MassAttenuation` table + log-log interpolation. —
   *interpolation reproduces node values exactly; monotone between nodes.*
2. [ ] Embed a small citable NIST XCOM μ/ρ dataset (water, air) over the MV/kV
   range, sourced from the NIST database (verified, not memorized). — *values
   traceable to the cited table.*
3. [ ] `Material`→`MassAttenuation(E)` lookup + mass-weighted mixture rule. —
   *water vs mixture value-semantic checks.*
4. [ ] clippy `-D warnings`, fmt, nextest, doctests green; sync artifacts.

*Blocked:* H-010 GPU kernel (G-12: hephaestus/leto stack convergence + wgpu build),
H-011c segment-generation (gaia G-11), H-004b ritk DICOM (heavy build).

### Completed

- [x] **H-032** `helios-analysis`: cumulative `Dvh` (Dx/Vx/mean) + `gamma_index_3d`
  (Low, global norm) + `gamma_pass_rate`. 8 tests (identical→γ=0, γ scales with
  dose ratio, 2×criterion→fail, uniform-DVH step, ramp quantiles, f32). Builds the
  3%/2 mm + DVH quality-gate machinery (G-3 partial). *(Sprint-4 crate pulled
  forward — it is unblocked and gate-relevant, unlike the GPU/geometry work.)*
- [x] **H-012b** `helios-solver::attenuation_map`: deterministic per-voxel HU→μ
  engine (CT `Volume` → μ `Volume`, Compton-MV approximation). CPU reference / GPU
  differential oracle. 5 tests (uniform water, air/bone, closed-form differential,
  grid preservation, f32).
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

## Gate status (last run, H-032)

| Gate | Result |
|------|--------|
| `cargo build` | pass |
| `cargo clippy --all-targets --all-features -D warnings` | pass, 0 warnings |
| `cargo fmt --check` | pass |
| `cargo nextest run` | 52 passed / 0 failed (1.4 s) |
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
