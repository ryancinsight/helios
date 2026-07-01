# Helios Checklist (tactical)

**Sprint target version:** `0.0.1` (Foundation → Sprint 2 physics/GPU)
**Current phase:** Phase 2 (Execution). Sprint 1 domain core complete
(`helios-core`/`math`/`domain`); Sprint 2 opened with `helios-physics`.

## Owner: claude-helios

### In-flight item: H-004b `helios-domain` ritk DICOM load path — `todo`

**H-021 done:** `helios-simulation::simulate_helical_sinogram` integrates
`HelicalDelivery` + forward projector into a helical MVCT acquisition (8/11 crates,
103 tests). Next: real inputs (ritk DICOM, H-004b) and reconstruction
(`helios-imaging`, H-042) to close the imaging/therapy validation loop; moirai
parallel projection dispatch (H-021b); planning (`helios-planning`, coeus).

### (done) H-021 helical delivery simulation

**G-14 RESOLVED (H-003c):** the concurrent leto geometry rewrite settled; leto+gaia
build against the new `leto::geometry` API. Helios adapted — `helios-math` re-exports
the new leto types + gaia `Aabb`/`Ray`; `VoxelGrid` simplified to axis-aligned
(origin+spacing); projector pose-check removed. **Full workspace builds; 97 tests
pass**; dose kernel superposition (H-013b) verified. The H-055 geometry-feature split
remains (physics still builds standalone).

Next: H-021 (moirai-orchestrated helical delivery simulation combining
`HelicalDelivery` kinematics + per-projection forward projection / dose over time),
then end-to-end dose→gamma/DVH validation.

*Also queued:* H-020b (binary-MLC sinogram), H-003d (oriented grids when leto
`Isometry3` gains transforms), H-011d (exact Siddon), H-010b (GPU HU→μ + throughput),
H-004b (ritk DICOM), H-011b (NIST μ/ρ tables).

## Gate status (last run, H-021 — helical simulation)

| Gate | Result |
|------|--------|
| `cargo build` (whole workspace) | pass (all 8 crates) |
| `cargo nextest run` | 103 passed / 0 failed (incl. live GPU) |
| `cargo clippy --all-targets --all-features -D warnings` | 0 code warnings |
| `cargo test --doc` | pass |
| `cargo fmt --check` | pass |

### Completed

- [x] **H-013a** `helios-solver::primary_fluence_parallel_x`: dose primary-transport
  stage (Ψ=Ψ₀·exp(−∫μ dl), +x parallel beam, O(N) cumulative). 4 oracles
  (homogeneous exponential, unattenuated entry, heterogeneous accumulation, f32).
  Also fixed projector optical-depth units (mm→cm; G-13).
- [x] **H-010** `helios-gpu`: real GPU kernel — `beam_transmission_into` computes
  `exp(-τ)` on the GPU (hephaestus-wgpu `NegOp`+`ExpOp`), differentially validated
  vs CPU `f32::exp` on a live adapter; `default_device`. Replicated hephaestus's
  mnemosyne/moirai/hermes `[patch]` set (resolves the G-12 cluster skew). Closes G-12.
- [x] **H-011c** `helios-solver::forward_project_ray`: MVCT forward projector —
  clip gaia `Ray` to grid `Aabb`, midpoint ray-march trilinear μ `Volume` → ∫μ dl.
  5 oracles (uniform slab τ=μ·L, affine midpoint-exact, step-invariance, miss, f32).
- [x] **H-050 / H-003b** Wired Helios to the local synchronized Atlas checkout
  (`[patch]` leto/eunomia/gaia → local paths); `helios-math` re-exports
  `gaia::{Aabb, Ray}`; bridge test green. Consumes gaia's migrated geometry;
  unblocks the projector.

### Deferred (still blocked / sequenced)

- **H-020b** binary-MLC leaf-open-time sinogram (unblocked; queued after projector).
- **H-010** GPU HU→μ kernel (add hephaestus-wgpu patch; adapter verified available).
- **H-004b** ritk DICOM (heavy build).

### Superseded in-flight plan: H-020b binary-MLC sinogram — `todo`

Unblocked (timing/modulation model, not spatial MLC geometry which needs gaia).

1. [ ] `LeafOpenTimeSinogram`: per-projection × per-leaf open-time fractions in
   `[0,1]`; validated bounds; indexed by (projection, leaf). — *round-trip
   set/get; out-of-range rejected.*
2. [ ] Effective per-leaf fluence weight = open_fraction·(1−leakage) + leakage
   (binary-MLC transmission/leakage model). — *closed vs open vs partial value
   checks; leakage floor honored.*
3. [ ] Tongue-and-groove inter-leaf underdose factor between adjacent open leaves.
   — *analytical factor vs hand-computed.*
4. [ ] clippy `-D warnings`, fmt, nextest, doctests green; sync artifacts.

*Blocked:* H-010 GPU kernel (G-12), H-011c segment-generation + spatial MLC
geometry (gaia G-11), H-004b ritk DICOM (heavy build). *Unblocked queue:* H-011b
NIST μ/ρ tables, H-021 delivery simulation stepping.

### Completed

- [x] **H-020** `helios-domain::HelicalDelivery`: helical delivery kinematics
  (gantry rotation + couch translation + pitch/time synchronization). 7 tests
  (one-rotation angle 2π + couch travel, projection↔time agreement, half-rotation π,
  monotonic couch, f32). Clinical-realism "helical synchronization" capability.
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

## Gate status (last run, H-011c)

| Gate | Result |
|------|--------|
| `cargo build` | pass (local gaia/leto/eunomia via `[patch]`) |
| `cargo clippy --all-targets --all-features -D warnings` | pass, 0 code warnings |
| `cargo fmt --check` | pass |
| `cargo nextest run` | 71 passed / 0 failed (incl. live GPU test) |
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
