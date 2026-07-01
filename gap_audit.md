# Helios Gap Audit

Physics, numerics, accuracy, architecture, and integration gaps. Closed by
evidence, not silence. Each gap: ID, description, class, current evidence tier,
target closure.

## Open gaps

### Physics / numerics

- **G-1 (physics):** *Partially closed (H-011).* Photon attenuation **relations**
  implemented and analytically verified in `helios-physics`: Beer–Lambert
  transmission, half-value layer, `μ = (μ/ρ)·ρ`, and first-order HU→density CT
  calibration (property/value-semantic tests: `T(HVL)=½`, `T(0)=1`, water/air/bone
  calibration points, f32 genericity). **Still open:** concrete NIST XCOM μ/ρ data
  tables (H-011b) and an electron-transport model; these are data/algorithm gaps,
  not framework gaps. *Evidence tier: analytical (relations) — reference cross
  sections pending.*
- **G-2 (numerics):** ~~No `Scalar` seam.~~ **CLOSED (H-003).** `helios-math`
  establishes `Scalar = eunomia::RealField` (the Atlas numeric SSOT) as the Helios
  compute seam and re-exports the leto linear-algebra substrate. `helios-core`
  constants remain `f64` literals by design and are converted by callers. The seam
  is exercised natively (`f32`/`f64`) by the first compute kernels as they land.
- **G-3 (accuracy):** *Partially closed (H-032).* The **validation machinery** now
  exists: `helios-analysis` implements the cumulative DVH (Dx/Vx/mean) and the 3D
  gamma index (Low, global normalization) + pass rate, with analytical oracles
  (identical→γ=0, criterion-scaled γ, uniform-DVH step, ramp quantiles). **Still
  open:** the dose-engine/projector *reference solutions* to validate (need
  H-013/H-011c) and clinical comparison vs VoLO/TOPAS/GATE/EGSnrc (H-042). The
  gates are implemented; the distributions to feed them are the remaining work.
  *Evidence tier: analytical (metrics verified) — clinical inputs pending.*
- **G-4 (numerics):** Reduction-order sensitivity for future GPU vs CPU differential
  tests not yet characterized; epsilon bounds must be derived per reduction depth
  when the projector/dose kernels land. → H-012.

### Architecture / integration

- **G-5 (integration):** Atlas crate *APIs* partially exercised. **eunomia**
  (`RealField`/`FloatElement`/`NumericElement`) and **leto** (`Vector3`, `Point3`,
  `Isometry3`, `Translation3`, `UnitQuaternion`/`Unit`, `Array3` C-contiguous +
  `as_slice`) verified against real usage and building in-tree (H-003, H-004).
  `ritk-io` (DICOM/MVCT), `gaia` (MLC geometry), moirai, coeus, consus surfaces
  remain unverified; **hephaestus** `ComputeDevice` seam (GAT `Buffer<T: Pod>`,
  `alloc_zeroed`/`upload`/`download`/`write_buffer` with `themis::PlacementHint`)
  read and scoped for H-010 (not yet built — heavy wgpu compile + GPU-device
  availability are the gating risks). Symbol existence must be confirmed via
  `cargo doc`/source before each first use (anti-hallucination). → H-004b, H-005,
  H-010+.
- **G-10 (integration, upstream co-evolution):** leto's **default** features pull
  `mnemosyne` at a rev pinned to `themis ^0.8`, which conflicts with themis HEAD
  `0.9.17` — a version skew in the Atlas stack's transitive git graph. *Workaround
  applied:* Helios consumes leto with `default-features = false, features=["std"]`,
  deferring mnemosyne placement to the layer that needs it (themis/mnemosyne
  integration, later sprint). *Upstream item:* leto's pinned mnemosyne rev (or
  mnemosyne's themis bound) should be advanced to themis 0.9.x so the default
  feature set resolves. File against the leto/mnemosyne repos when that layer is
  built. *Evidence tier: reproduced (cargo resolution error), worked around.*
- **G-6 (build hygiene):** ~~Helios target-dir sharing.~~ **CLOSED.** Helios
  automatically routes its build through the shared `D:/atlas/target` via the
  inherited `repos/.cargo/config.toml` (`[build] target-dir`); Cargo discovers it by
  walking up from the package dir. Verified: `cargo doc` emitted to
  `D:/atlas/target/doc` and no per-`helios` `target/` exists. No action needed;
  backlog H-006 closed.
- **G-7 (toolchain):** `rust-toolchain.toml` pins `channel = "stable"` (currently
  1.95) but does not pin an exact version; MSRV floor declared as 1.85 in
  `Cargo.toml` (`rust-version`) but not yet CI-verified. → revisit at first CI.

- **G-11 (integration, geometry ownership):** Geometry primitives (`Aabb`, `Ray`,
  intersection, meshes, CSG) are owned by **gaia**, not Helios. gaia already has
  `Aabb<T: Scalar>` (over `leto::Point3`) and a validated-`UnitVector3` `Ray` with
  `intersect_aabb`. *Update (this session):* gaia's leto/eunomia migration is now
  **finalized and green** — gaia builds across all targets, **927 tests pass**,
  doctests pass, fmt clean; `Ray`/`Aabb` are committed and re-exported from `gaia`'s
  crate root (commits `b058eb0`, `ecd4060`). The source blocker is **resolved**.
  *Action taken earlier:* removed the duplicate `Ray`/`Aabb` from `helios-math`
  (upstream ownership). **Remaining (consumption wiring, H-003b):** the migration
  lives on gaia's `refactor/migrate-to-leto-geometry` branch, not yet merged to
  gaia's default branch — merging is a `refactor!` breaking change that also affects
  kwavers (co-evolution). **Update: consumption wired (H-050).** Helios now `[patch]`-
  redirects `leto`/`eunomia`/`gaia` git sources to the local synchronized checkout
  (one consistent source) and `helios-math` re-exports `gaia::{Aabb, Ray}`; a bridge
  test (gaia `Ray` ∩ gaia `Aabb` through Helios) passes. **G-11 is effectively
  closed** for local development — the projector (H-011c) is unblocked. Remaining
  release step: merge gaia migration to its default branch + update kwavers, then
  drop the patch. *Evidence tier: verified — Helios builds + 60 tests with local
  gaia geometry.*

- **G-12 (integration, GPU backend blocked):** `helios-gpu` on `hephaestus-wgpu` is
  blocked on the Atlas stack's leto/hephaestus dependency convergence — the same
  migration the goal flags ("gaia will need to move to leto/hephaestus"). Evidence:
  hephaestus's workspace consumes `leto`/`mnemosyne`/`themis` via **local path deps**
  with the `mnemosyne-memory` feature and a pinned `themis` rev, i.e. the same
  leto→mnemosyne→themis cluster that failed resolution in G-10, now compounded by a
  heavy `wgpu` build. Consuming `hephaestus-wgpu` as a git dep would not resolve
  cleanly against the current stack. *Decision:* do not force the GPU backend now;
  author every engine as a CPU reference first (`helios-solver`) so the GPU path
  (H-010) is a differential drop-in once the stack stabilizes. The
  `hephaestus_core::ComputeDevice` seam and `hephaestus-wgpu` op surface
  (`WgpuDevice::try_default`, `unary/scalar_elementwise_strided`, `reduction`) are
  already scoped for that increment.
  *Update (this session):* hephaestus is **verified green locally** — the workspace
  builds, `hephaestus-core` (21 tests) and `hephaestus-wgpu` (109 tests) pass, fmt
  clean, 0 code clippy warnings. Crucially the **wgpu GPU contract tests pass, so a
  usable GPU adapter exists in this environment** (upload/download round-trips,
  strided-elementwise-vs-CPU, sparse spmv/spmm all green). The source repo is not
  broken. **Remaining:** the git-dep *version-alignment* skew (hephaestus uses local
  path deps to the leto/mnemosyne/themis cluster) means Helios must consume it via a
  local `[patch]`/path (synchronized checkout), same wiring as G-11. *Evidence tier:
  verified — hephaestus builds + 130 tests + GPU adapter working locally.*

### Testing / tooling

- **G-8 (coverage):** No `cargo-llvm-cov` run yet; >80% core-logic coverage target
  unmeasured. Applies from first `[minor]`. → measure after H-003.
- **G-9 (CI):** No CI pipeline wired (fmt/clippy/nextest/doc/audit/deny). Gates are
  run locally only. → file when the workspace has ≥2 crates.

## Closed gaps

- **G-6 (build hygiene):** Helios inherits the shared `D:/atlas/target` build dir
  via `repos/.cargo/config.toml`; no per-repo target. Verified this session.

## Residual risk register

- Atlas upstream APIs may drift (multi-repo co-evolution); Helios must pin commits
  in `Cargo.lock` and add cross-repo contract tests as it consumes each crate
  (G-5). Currently no lockfile committed for git deps because none are used yet.
- Physical constants (G-2) are CODATA-2018/ICRU-90 values verified by inter-constant
  derivation tests, not by an external authoritative fetch this session; values are
  standard and cross-checked, but a future audit should confirm against the live
  NIST database.
