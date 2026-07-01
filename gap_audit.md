# Helios Gap Audit

Physics, numerics, accuracy, architecture, and integration gaps. Closed by
evidence, not silence. Each gap: ID, description, class, current evidence tier,
target closure.

## Open gaps

### Recently closed

- **G-14 — RESOLVED (H-003c).** The concurrent leto geometry rewrite settled: leto
  and gaia now build against the new `leto::geometry` API (Vector3/Point3 with
  `.x/.y/.z` fields; `Isometry3` reduced to `{rotation, translation}`). Helios was
  adapted: `helios-math` re-exports the new leto types (`Point2/Point3/Vector3/
  UnitVector3`) + gaia `Aabb/Ray`; `VoxelGrid` simplified to **axis-aligned**
  (origin + spacing), dropping the now-incomplete `Isometry3` pose (oriented grids
  are a follow-up pending a rigid-transform primitive with `transform`/`inverse`);
  the projector's pose-rotation check was removed. **Full workspace builds; 97 tests
  pass** (all crates incl. live GPU), clippy `-D warnings` clean, fmt clean. The
  previously-blocked dose kernel-superposition engine (H-013b) is now built and
  verified. *Evidence tier: verified — whole-workspace build + 97 tests green.*

### (historical) BLOCKER — concurrent foundation refactor

- **G-14 (integration, BLOCKED — concurrent leto geometry relocation):** Mid-session
  the shared **leto** submodule advanced (peer/concurrent work) and its `geometry`
  module (`Vector3`/`Point3`/`Isometry3`/`UnitQuaternion`/…) is no longer present at
  leto's current HEAD (`git ls-files crates/leto/src/geometry` is empty; earlier this
  session gaia built 927 tests against `leto::geometry`). The types are not yet
  relocated to a discoverable home, so **gaia fails to compile** (86 errors,
  `unresolved import leto::geometry`), and every Helios crate that depends on
  `helios-math` (which re-exports gaia geometry) transitively fails to build —
  including at the last green commit `2ce36787` (the foundation shifted under it).
  *Interpretation:* the Atlas stack is mid-migration to **gaia-native geometry**
  (geometry moving out of leto into gaia), the end-state the earlier feedback
  intended. *Action (discipline):* do NOT fix leto/gaia's in-flight relocation
  (peer's active, cross-stack work; unknown target), do NOT revert the shared
  submodule, do NOT commit broken/unverified Helios code. **Deferred item H-013b
  (dose kernel superposition, `dose_convolution_x` + `exponential_deposition_kernel`)
  is written in `crates/helios-solver/src/dose.rs` with exact analytical oracles
  (delta-kernel identity, normalized-kernel interior conservation, physical build-up)
  but is UNVERIFIED — it cannot build until the geometry foundation settles.**
  *Unblock:* when gaia's native geometry lands, update `helios-math` to re-export all
  geometry from gaia (H-003c), then verify + commit H-013b. *Evidence tier:
  reproduced (leto HEAD has no geometry; gaia 86-error build failure).*
  *Mitigation (H-055):* `helios-math`'s geometry vocabulary is now behind a default
  `geometry` feature and `helios-physics` consumes it with `default-features=false`,
  so `helios-core`, `helios-math` (scalar seam) and `helios-physics` **build/test
  independently** of the churning geometry stack (`cargo nextest run -p helios-core
  -p helios-physics`). Only geometry-dependent crates (`helios-domain`/`-solver`,
  whole-workspace `cargo build`) remain blocked until the foundation settles.

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
  H-013) and clinical comparison vs VoLO/TOPAS/GATE/EGSnrc (H-042). *Update:* the
  **MVCT forward projector** (`helios-solver::forward_project_ray`, H-011c) now
  produces line-integral projections (∫μ dl), analytically verified (uniform slab
  τ=μ·L, affine-field midpoint-exact). The dose engine's **primary-transport stage**
  (H-013a, `primary_fluence_parallel_x`) now produces the analytical exponential
  depth curve `Ψ₀·exp(−μx)`; the remaining stage is kernel superposition → dose
  (H-013b). Clinical comparison vs VoLO/TOPAS/GATE/EGSnrc (H-042) still pending real
  inputs.
  *Evidence tier: analytical (metrics + projector + primary transport verified) —
  full dose engine + clinical inputs pending.*
- **G-13 (numerics, units):** ~~Projector optical depth was 10× too large.~~
  **CLOSED (H-013a):** `forward_project_ray` mixed `μ` [cm⁻¹] with mm path length;
  now converts mm→cm for a dimensionless `τ`. Units contract documented (μ volumes
  in cm⁻¹, grid in mm). *Evidence tier: analytical (τ = μ·L_cm verified).*
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
  local `[patch]`/path (synchronized checkout), same wiring as G-11.
  *Update (H-010): CLOSED.* `helios-gpu` dispatches a real GPU kernel —
  `beam_transmission_into` computes `exp(-τ)` on the GPU (hephaestus-wgpu
  `NegOp`+`ExpOp`); a differential test vs CPU `f32::exp` passes on the live adapter.
  Wiring: replicated hephaestus's mnemosyne/moirai/hermes `[patch]` set so the
  leto→mnemosyne(git 1e014d25)→themis ^0.8 skew resolves to the local consistent
  cluster; hephaestus-wgpu consumed with default features (its `linalg` uses
  `leto-ops` ungated). *Evidence tier: verified — Helios GPU kernel runs + matches
  CPU (67 tests).* Remaining: fused HU→μ GPU kernel needs a custom affine-clamp
  `UnaryWgslOp` (H-010b); throughput benchmark vs VoLO pending.

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
