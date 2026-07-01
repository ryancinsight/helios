# Helios Backlog (strategic board)

Single source of cross-session strategy. Each item carries a stable ID, a
change-class tag, a status, an owner, and a claimed scope. Triage order: correctness
gaps → architecture drift → missing tests → docs → PM cleanup.

Status: `todo` · `in-progress` · `review` · `done`

## Sprint 1 — Foundation

| ID | Item | Class | Status | Owner | Scope |
|----|------|-------|--------|-------|-------|
| H-001 | Workspace skeleton + Foundation artifacts (README, ARCHITECTURE, PM files) | [arch] | done | claude-helios | `Cargo.toml`, root `*.md`, `.config/` |
| H-002 | `helios-core`: typed errors, physical constants, validating newtypes | [minor] | done | claude-helios | `crates/helios-core/**` |
| H-003 | `helios-math`: `Scalar` seam (= `eunomia::RealField`) + leto linear-algebra substrate re-export | [minor] | done | claude-helios | `crates/helios-math/**` |
| H-003b | Consume gaia `Aabb`/`Ray` in `helios-math` (re-export as Helios geometry). Blocked on gaia's leto-geometry migration landing on its default branch (G-11). | [minor] | blocked | — | `crates/helios-math/**` |
| H-004 | `helios-domain`: `VoxelGrid` (index↔world via leto `Isometry3`) + `Volume<T>` trilinear sampling over leto `Array3` | [minor] | done | claude-helios | `crates/helios-domain/**` |
| H-004b | `helios-domain`: `ritk-io` DICOM load path (CT/MVCT → `Volume`); `CtVolume`/`MvctVolume` HU-semantic newtypes; DICOM `ImageOrientationPatient` → grid pose | [minor] | todo | — | `crates/helios-domain/**` |
| H-005 | `helios-domain`: gaia-backed binary-MLC + collimator/jaw geometry model | [minor] | todo | — | `crates/helios-domain/**` |
| H-006 | ~~Shared `CARGO_TARGET_DIR`~~ — resolved: inherited from `repos/.cargo/config.toml` (shared `D:/atlas/target`) | [patch] | done | claude-helios | — |

## Integration unblock (gaia/hephaestus now green)

| ID | Item | Class | Status | Owner | Scope |
|----|------|-------|--------|-------|-------|
| H-050 | Wire Helios to the synchronized local Atlas checkout: `[patch]` redirecting `leto`/`eunomia`/`gaia` git sources to local paths (one consistent source). **Done** for geometry; hephaestus-wgpu patch added when GPU kernel lands (H-010). | [arch] | done | claude-helios | `Cargo.toml` |
| H-003b | `helios-math` re-exports `gaia::{Aabb, Ray}` (consumed the migrated geometry); bridge test green. | [minor] | done | claude-helios | `crates/helios-math/**` |

Context: as of this session gaia's leto migration is finalized + green (927 tests)
and hephaestus builds with wgpu GPU tests passing (130 tests, adapter available).
The remaining blocker is purely dependency *wiring* (git-dep version skew across the
leto/mnemosyne/themis cluster → use local path/patch). Merging gaia's `refactor!`
migration to its default branch + pushing is a separate co-evolution step (breaking
for kwavers) requiring consumer coordination.

## Sprint 2 — GPU foundation

| ID | Item | Class | Status | Owner | Scope |
|----|------|-------|--------|-------|-------|
| H-010 | `helios-gpu`: program against `hephaestus_core::ComputeDevice`; runtime backend selection (wgpu/cuda); GPU HU→μ kernel differentially validated vs `helios-solver::attenuation_map`. **Blocked (G-12)** on Atlas leto/hephaestus stack convergence + wgpu build | [minor] | blocked | — | `crates/helios-gpu/**` |
| H-011 | `helios-physics`: photon attenuation relations — `LinearAttenuation`/`MassAttenuation`, Beer–Lambert, HVL, HU→density calibration | [minor] | done | claude-helios | `crates/helios-physics/**` |
| H-011b | `helios-physics`: NIST XCOM μ/ρ data tables (energy-indexed, per material) loaded into `MassAttenuation` | [minor] | todo | — | `crates/helios-physics/**` |
| H-011c | `helios-solver::forward_project_ray`: ray-march optical depth ∫μ dl of a gaia `Ray` through a μ `Volume` (clip to grid `Aabb`, midpoint trilinear sampling). MVCT forward-projection / dose ray-trace core. | [minor] | done | claude-helios | `crates/helios-solver/**` |
| H-011d | `helios-solver`: exact Siddon voxel-DDA + oriented-grid (rotated pose) projection; full parallel/fan sinogram over a detector | [minor] | todo | — | `crates/helios-solver/**` |
| H-012 | `helios-solver`: GPU MVCT forward projector (Siddon/Joseph); CPU reference | [minor] | todo | — | `crates/helios-solver/**` |
| H-012b | `helios-solver`: HU→μ attenuation-map engine (CPU reference; differential oracle for the GPU kernel) | [minor] | done | claude-helios | `crates/helios-solver/**` |
| H-013 | `helios-solver`: collapsed-cone / convolution-superposition dose engine (CPU ref first) | [major] | todo | — | `crates/helios-solver/**` |

## Sprint 3 — Delivery

| ID | Item | Class | Status | Owner | Scope |
|----|------|-------|--------|-------|-------|
| H-020 | `helios-domain`: helical delivery kinematics (gantry rotation + couch velocity + pitch/time synchronization) — `HelicalDelivery<T>` | [minor] | done | claude-helios | `crates/helios-domain/**` |
| H-020b | `helios-domain`: binary-MLC leaf-open-time sinogram model (per-projection leaf pattern) + leakage/transmission/tongue-and-groove factors | [minor] | todo | — | `crates/helios-domain/**` |
| H-021 | `helios-simulation`: moirai-orchestrated time-dependent helical delivery + synchronized MVCT | [major] | todo | — | `crates/helios-simulation/**` |
| H-022 | Binary-MLC leakage/transmission/tongue-and-groove model | [minor] | todo | — | `crates/helios-domain/**` |

## Sprint 4 — Planning & imaging

| ID | Item | Class | Status | Owner | Scope |
|----|------|-------|--------|-------|-------|
| H-030 | `helios-imaging`: MVCT reconstruction (FBP + iterative) | [major] | todo | — | `crates/helios-imaging/**` |
| H-031 | `helios-planning`: coeus-autodiff inverse planning (gradient-based) | [major] | todo | — | `crates/helios-planning/**` |
| H-032 | `helios-analysis`: DVH (cumulative, Dx/Vx/mean) + 3D gamma index (Low, global norm) + pass rate | [minor] | done | claude-helios | `crates/helios-analysis/**` |
| H-032b | `helios-analysis`: structure-masked DVH (RT-struct ROIs via ritk) + local-normalization gamma + low-dose threshold cutoff | [minor] | todo | — | `crates/helios-analysis/**` |

## Sprint 5 — End-to-end

| ID | Item | Class | Status | Owner | Scope |
|----|------|-------|--------|-------|-------|
| H-040 | `helios-python`: PyO3 high-level API (maturin, pytest equivalence) | [minor] | todo | — | `crates/helios-python/**` |
| H-041 | End-to-end helical TomoTherapy workflow example (Rust + Python) | [minor] | todo | — | `examples/**` |
| H-042 | Validation report: gamma/DVH vs reference; MVCT image metrics | [minor] | todo | — | `validation_reports/**` |
| H-043 | Performance: GPU scaling study + criterion baselines | [minor] | todo | — | `benches/**` |
