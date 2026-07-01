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

## Sprint 2 — GPU foundation

| ID | Item | Class | Status | Owner | Scope |
|----|------|-------|--------|-------|-------|
| H-010 | `helios-gpu`: program against `hephaestus_core::ComputeDevice`; runtime backend selection (wgpu/cuda); GPU HU→μ kernel differentially validated vs `helios-solver::attenuation_map`. **Blocked (G-12)** on Atlas leto/hephaestus stack convergence + wgpu build | [minor] | blocked | — | `crates/helios-gpu/**` |
| H-011 | `helios-physics`: photon attenuation relations — `LinearAttenuation`/`MassAttenuation`, Beer–Lambert, HVL, HU→density calibration | [minor] | done | claude-helios | `crates/helios-physics/**` |
| H-011b | `helios-physics`: NIST XCOM μ/ρ data tables (energy-indexed, per material) loaded into `MassAttenuation` | [minor] | todo | — | `crates/helios-physics/**` |
| H-011c | `helios-physics`: **reduction done** — `optical_depth`/`beam_transmission` over `(μ, length)` segments. Remaining: voxel-DDA/Siddon *segment generation* through a `Volume` needs gaia `Ray`/`Aabb` (blocked on G-11) | [minor] | partial | claude-helios | `crates/helios-physics/**` |
| H-012 | `helios-solver`: GPU MVCT forward projector (Siddon/Joseph); CPU reference | [minor] | todo | — | `crates/helios-solver/**` |
| H-012b | `helios-solver`: HU→μ attenuation-map engine (CPU reference; differential oracle for the GPU kernel) | [minor] | done | claude-helios | `crates/helios-solver/**` |
| H-013 | `helios-solver`: collapsed-cone / convolution-superposition dose engine (CPU ref first) | [major] | todo | — | `crates/helios-solver/**` |

## Sprint 3 — Delivery

| ID | Item | Class | Status | Owner | Scope |
|----|------|-------|--------|-------|-------|
| H-020 | `helios-domain`: helical delivery kinematics (gantry rotation + couch velocity + leaf sequencing) | [minor] | todo | — | `crates/helios-domain/**` |
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
