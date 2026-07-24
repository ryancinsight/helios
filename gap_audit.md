# Helios Gap Audit

Physics, numerics, accuracy, architecture, and integration gaps. Closed by
evidence, not silence. Each gap: ID, description, class, current evidence tier,
target closure.

## Aequitas metric gap audit (2026-07-23)

The dose field itself remains `helios_domain::Volume<T>` storage. This audit
targets values that cross a public analysis, delivery, or geometry API. It does
not count dense field storage, gamma values, volume fractions, homogeneity
indices, or response probabilities as dimensional gaps.

### Existing coverage

Helios already uses Aequitas `EnergyPerArea` for portal fluence, `AbsorbedDose`
for deposition totals and DVH samples, `Length`/`ReciprocalLength` for geometry
and attenuation, and `AreaPerMass` for mass attenuation. `EnergyMeV` and
`VoxelSpacingMm` are validated consumer newtypes backed by Aequitas quantities.

### Open implementation ledger

| ID | Evidence | Remaining implementation | Owner | Status / acceptance oracle |
|---|---|---|---|---|
| `HELIOS-AEQ-MET-01` | `helios-analysis/src/dvh.rs` stores `Vec<AbsorbedDose<T>>`, but `min`, `max`, `mean`, `dose_at_volume_fraction`, and `generalized_eud` return raw `T`; dose criteria and TD50/TCD50 parameters also enter as `T`. | Return `AbsorbedDose<T>` for dose-valued results and parameters; keep Vx, HI, TCP, and NTCP dimensionless or probability-typed. | Helios | **RESOLVED.** `Dvh` extrema, mean, Dx, gEUD, and TCP/NTCP dose parameters now use `AbsorbedDose<T>`; nearest-rank, masked, NaN, Asclepius-law, and end-to-end PTV/OAR value semantics remain covered. |
| `HELIOS-AEQ-MET-02` | `helios-analysis/src/gamma.rs` accepted `dta_mm`, normalization dose, low-dose cutoff, and search radius as raw `T`; gamma volume and pass rate are dimensionless. | Type distances as `Length`, dose thresholds as `AbsorbedDose`, and keep the result storage scalar/dimensionless. | Helios | **RESOLVED.** `gamma_index_3d`, `gamma_index_3d_local`, and `gamma_pass_rate` now type physical criteria with Aequitas while retaining the Low gamma kernel, local/global normalization, grid checks, scalar gamma field, and scalar pass rate. Focused value-semantic gamma tests and all in-tree callers migrate; ADR 0007 records the breaking boundary. |
| `HELIOS-AEQ-MET-03` | `helios-simulation/src/delivery.rs` stored leaf fluence as `T`; `total_delivered_fluence` returned `T`. `portal.rs` constructed `EnergyPerArea` internally and converted it back. `dose_accumulation.rs` accepted `*_mm` geometry and sampling values as `T`. | Carry fluence as `EnergyPerArea` and geometry distances as `Length` through delivery, portal, and dose accumulation. | Helios | **RESOLVED.** `DeliveryFrame`, collimation, portal transmission, total fluence, and dose geometry now use Aequitas quantities. Typed values convert once to the existing millimetre ray/voxel kernel; closed-leaf zero, Beer–Lambert darkening, fluence linearity, geometry-limit, f32, example, and end-to-end regressions pass. ADR 0008 records the breaking boundary. |
| `HELIOS-AEQ-MET-04` | `helios-analysis/src/image_quality.rs` returns raw intensity/RMSE values, while the same analysis can operate on dose volumes. | Decide the semantic input at the analysis boundary: retain raw image intensity for MVCT, but return `AbsorbedDose` statistics when the API contract is dose-specific. | Helios | **RESOLVED.** Shared ROI/RMSE value kernels now back raw MVCT `roi_statistics`/`volume_rmse` and typed-dose `dose_roi_statistics`/`dose_volume_rmse`; the clinical validation example uses typed dose means/stddev and converts only for dimensionless contrast/CNR. Value tests cover f64/f32 and Gray outputs; ADR 0009 records the partition. |
| `HELIOS-AEQ-MET-05` | `helios-physics` Compton APIs accepted raw scalars documented as MeV for photon energy, so the public boundary carried no energy dimension. | Accept Aequitas `Energy<T>` for Klein–Nishina, energy-transfer, and Compton mass-coefficient APIs; convert to MeV only at the dimensionless kernel boundary. | Helios | **RESOLVED.** Rust APIs, examples, tests, and the Python conversion boundary now use typed Aequitas energy; the 1 MeV/1,000,000 eV equivalence test preserves the analytical result. ADR 0010 records the breaking boundary. |

### Explicit non-gaps and constraints

- TERMA/dose arrays remain scalar storage until a field-descriptor contract can
  carry dimensions without changing dense-kernel representation.
- Beam angles, gamma values, fractions, homogeneity indices, and TCP/NTCP are
  dimensionless; they must not be wrapped as length or dose merely because they
  are reported beside physical quantities.
- `HELIOS-AEQ-MET-04` is closed. Raw MVCT intensity remains scalar, while
  dose-specific ROI and RMSE results carry `AbsorbedDose`; no Aequitas provider
  extension was required. Future public signature changes must update their
  examples, Python surface, and focused tests in the same change.

## H-003d oriented-grid provider convergence (closed)

- Leto 0.38 now owns checked conversion from world-space rotation columns to a
  `UnitQuaternion`; Helios uses that provider contract to restore the
  local-index-to-world `Isometry3` grid pose. CPU ray clipping and terma
  deposition run in scaled-index space and retain the world-space millimetre
  parameter; the HDF5 boundary records three validated rotation columns. The
  grid core does not duplicate DICOM tags or matrix-to-quaternion logic.
- The current Hephaestus `FieldGeometry` has no rigid pose. `GpuProjector`
  therefore returns a typed dispatch error for a non-identity rotation before
  data upload (correctness evidence: type/contract plus a value-semantic test),
  rather than silently omitting orientation. A pose-bearing GPU field geometry
  remains an upstream Hephaestus capability gap.
- Evidence tier: type-level rigid pose and checked Leto basis construction;
  analytical/differential validation through the 104/104 focused nextest run
  (including oriented Beer–Lambert and HDF5 pose round trips, live GPU checks);
  warning-denied Clippy, doctest/rustdoc, workspace example build, workspace
  format check, and four 196/196-package SemVer checks are clean.
- H-004d remains externally sequenced: RITK's public DICOM tags currently omit
  `ImageOrientationPatient`, and both permitted RITK worktree lanes carry
  active peer migrations. Helios will consume the named provider tag once that
  owner lane is available.

## Open gaps

### H-088 — deterministic book-figure SSOT gate (implemented, PR #32 merged)

- `xtask` now owns `FIGURE_SPECS`, deterministic SHA-256 manifest generation,
  and the `check-figures` command. The command validates the seven committed
  SVGs, scans `docs/book/SUMMARY.md` and `docs/book/README.md`, and fails on
  either an unlisted asset or a docs/spec mismatch.
- The workflow runs the gate after Rustdoc. The Python and benchmark lanes
  allow Cargo to refresh the lock after Atlas path-dependency materialization;
  the Rust workspace remains locked. The previous CI failure was the absent
  subcommand, not a figure mismatch.
- Local evidence: `cargo check -p xtask --offline`, formatter check, `mdbook
  build docs/book`, and `cargo run -p xtask --offline --locked -- check-figures`
  pass with `SSOT_IN_SYNC` and 7/7 references. Hosted PR #32 run
  `30070400660` passes build, Rust workspace, Python bindings, and the
  replicated benchmark gate; the PR merges as `02d7a775`.

### H-087 — portal fluence quantity boundary (implemented, PR #32 merged)

- `helios-simulation::frame_portal_fluence` now carries the transmitted portal
  fluence as Aequitas `EnergyPerArea<T>` through Hyperion's dimensionless
  transmission product before converting at the established scalar frame API.
  The direct Aequitas pin is `3ae0b6b`; implementation commit `b2a9ebe`.
- This closes the remaining quantity-conversion seam in the portal workflow;
  the dense fluence frame remains representation storage, not a second metric
  owner. The direct provider pin is now the merged Aequitas revision
  `e0fc5f3`. Existing full-transmission, Beer–Lambert, closed-leaf, f32, and
  invalid-optical-depth regressions remain the behavioral oracle.
- PR #32 merges as `02d7a775` from implementation head `31147f0` and PM
  follow-up `5832ffa`. The hosted build, Rust workspace, Python, and
  replicated benchmark checks pass. A local focused `helios-simulation`
  Nextest attempt could not start because active peer CFDrs/Leto builds held
  the shared Atlas lock; no local package-gate result is claimed.

### G-29 — DICOM charset dependency (externally blocked)

- `dicom-encoding` 0.10.0 declares `encoding` 0.2.33 unconditionally and uses
  it for the DICOM Specific Character Set codecs. The current release exposes
  no feature that can remove the dependency. A 2026-07-20 registry and resolved
  manifest audit confirms 0.10.0 remains the latest release and still carries
  the unconditional edge.
- RUSTSEC-2021-0153 reports maintenance status, not a known vulnerability.
  CI quarantines only that advisory ID while continuing to deny every other
  warning and vulnerability. Reopen H-073 when `dicom-rs` publishes a release
  backed by a maintained charset implementation; a consumer fork or reduced
  character-set implementation would duplicate provider ownership.

### Recently closed

- **G-30 — RESOLVED (H-083).** A stale unclaimed mdBook expansion described
  numerous types and methods absent from the current Helios source, while four
  older appendix/numeric pages contained control characters created by escaped
  Markdown delimiters. The expansion was not published. The recovered book now
  names APIs verified against crate re-exports, delegates release history and
  architecture to their root SSOT documents, contains no control characters,
  resolves every relative Markdown link, and builds without warnings under
  mdBook 0.5.4. The legacy scanner independently proves that its emptied
  migration allowlist matches zero current legacy dependency or source surfaces.

- **G-27 — RESOLVED (H-071).** The copied same-run classifier and bare native
  test invocation are deleted. Implementation head `44fb2768d` uses Atlas gate
  `9bfb722`, holds the candidate Criterion harness constant, measures
  phase-reversed ABBA and BAAB replications, runs native tests through the
  committed Nextest budget, and passes Rust, Python, and benchmark jobs in
  hosted run `29784712768`.

- **G-28 — RESOLVED (H-072).** The isolated Python binding crate previously
  resolved PyO3 0.23.5, which is affected by RUSTSEC-2025-0020 and
  RUSTSEC-2026-0177. PyO3 0.29.0 closes both vulnerabilities; the thin boundary
  uses `Python::detach` around the existing Rust planning call and retains no
  domain logic. The built extension is covered by the value-semantic Python
  contract suite. G-29 records the sole exact unmaintained dependency
  quarantine.

- **G-26 — RESOLVED (H-068).** `EnergyMeV` and `VoxelSpacingMm` previously
  stored dimensionless `f64` values despite representing physical quantities.
  Their validated newtypes now store Aequitas `Energy<f64>` and `Length<f64>`,
  preserve MeV and millimetre at the public boundary, and retain their
  zero-overhead scalar layout through compile-time size/alignment assertions.
  Round-trip properties use a bound derived from four machine-epsilon
  roundings. `HounsfieldUnit` remains Helios-owned because it is a calibrated
  non-SI scale. Evidence: warning-denied all-target Clippy, 17/17 configured
  Nextest tests, doctests, and warning-clean rustdoc for `helios-core`.

- **G-25 — RESOLVED (H-067).** The stale local lock changed only the
  `apollo-fft` version field, which was not a complete Cargo resolution and
  failed the warning-denied locked gate. Regenerating the Apollo package
  closure selects `apollo-fft` 0.25.0, Eunomia 0.4.0, Leto 0.38.2, and
  Hephaestus 0.17.0, removes Eunomia's obsolete `num-traits` edge, and removes
  Hephaestus WGPU's `num-complex` edge plus the package itself. The root
  manifest already follows the Apollo default branch; no source, manifest,
  compatibility wrapper, or fallback change is required. Locked metadata and
  format pass; warning-denied all-target/all-feature workspace Clippy passes;
  configured Nextest is 272/272; all ten Rust library doctest targets pass
  with zero examples; workspace rustdoc is warning-clean. Evidence tier:
  compiler-checked dependency resolution, warning-denied diagnostics, and
  value-semantic workspace regression execution.

- **G-24 — RESOLVED (H-066).** The workspace declared `num-traits` directly
  even though no Helios manifest or source consumed it. The direct declaration
  is removed. `cargo check --workspace --locked` passes, and
  `cargo tree -i num-traits --locked --edges normal` shows only transitive
  provider paths through Eunomia, Gaia, Half, WGPU, and their dependencies.
  Evidence tier: compiler-checked dependency resolution plus an inverse
  dependency-tree ownership audit.

- **G-23 — RESOLVED (H-005 reconciliation).** The foundation roadmap still
  listed a binary-MLC plus collimator/jaw model as todo after H-020b delivered
  `LeafOpenTimeSinogram`/`MlcModel` and H-020k delivered gaia-backed
  `FieldAperture` plus `collimate_frames`. The board marks duplicate H-022
  done and the README states the delivered ownership without creating a
  duplicate implementation track.
  The same reconciliation check found and restored the exact workspace formatter
  output in `helios-math::lib` and the solver deposition/projector tests.
  Evidence tier: source/API and board reconciliation plus formatter verification.

- **G-19 — RESOLVED (H-062).** Repeated `Dvh::volume_fraction_at_dose` queries
  previously scanned the complete sorted sample for every threshold, which made
  a plan with `q` DVH queries O(q·n). The query now uses a zero-allocation
  `partition_point` lower bound over the existing sorted slice, reducing the
  finite/infinite path to O(q·log n). A `contains_nan` marker preserves the
  previous `>=` filter semantics for unordered samples, and a NaN threshold
  returns zero as before. The fixed 64³/1,024-query Criterion workload measured
  30.090 ms [29.717, 30.472] for the scan reference and 29.229 μs [28.426,
  30.023] for production after the change; the paired median ratio is 1,029×.
  Value-semantic focused nextest passes 34/34. Evidence tier: empirical
  Criterion comparison plus value-semantic/differential boundary tests; see
  `validation_reports/2026-07-15-dvh-query-optimization.md`.

- **G-20 — RESOLVED (H-063).** Helios's direct DICOM dependency declared 0.8
  while the local `ritk-dicom` provider supplied `dicom-object` 0.10. The
  resulting duplicate `dicom_core::Tag` types caused four E0308 errors in
  `helios-domain/src/dicom.rs` during the workspace example check. Helios now
  declares DICOM 0.10, and the lockfile resolves one DICOM 0.10.0 graph across
  Helios and `ritk-dicom`. Locked workspace examples, all-target all-feature
  Clippy, 261/261 workspace nextest tests, doctests, and rustdoc pass. Evidence
  tier: compile-time dependency/type verification plus value-semantic tests.

- **G-21 — RESOLVED (H-064).** Helios previously used `dicom::core` and
  `dicom::object` directly for typed attributes and synthetic test fixtures,
  despite `ritk-dicom` owning the DICOM boundary. Helios now consumes only the
  `ritk-dicom` public tags, attribute-read trait, parser, transfer-syntax, and
  decoder contracts. Production and test scans contain no direct `dicom::`
  imports; the focused provider-backed domain suite remains 41/41.

- **G-22 — RESOLVED (H-065).** Helios's lockfile still selected Moirai 0.2.0
  after the upstream 0.3.0 release retired its unowned NUMA iterator and
  benchmark. The regenerated lockfile selects 0.3.0 for every Moirai package;
  `cargo check --workspace --examples --all-features` compiles the complete
  example graph. Evidence tier: compiler-checked dependency resolution and
  example compilation.
  Evidence tier: dependency/identifier scan plus value-semantic nextest.

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

- **G-15 (imaging accuracy):** *Partially addressed (H-030, H-033).* MVCT
  reconstruction (parallel-beam FBP) validated by a forward→reconstruct round-trip on a
  disk phantom, now *quantified* with `helios-analysis::image_quality` metrics
  (interior-ROI accuracy within 15% of μ₀, background suppression, disk/air contrast
  >0.85, CNR >1), and *quantum noise* (H-033b: `Var(τ')≈e^{τ}/N₀` validated vs
  analytical photon statistics; end-to-end noisy-recon noise scales with flux).
  **Remaining:** statistical reconstruction (OS-SEM/MLEM) and validation vs *published
  TomoTherapy MVCT data*; SIRT iterative reconstruction landed (H-030c, converges to its
  forward model, robust to noise/sparse-angle). The DICOM real-input path now ingests
  both a single slice
  (H-004b) and a full multi-slice **series** → 3-D HU `Volume` (H-004c:
  `load_ct_series` via `ritk-dicom`), so a real CT/MVCT study can drive the pipeline —
  clinical *dataset* validation still needs a licensed reference dataset. *Evidence tier: analytical/round-trip + synthetic-phantom
  metrics + real DICOM parse (synthetic round-trip through the ritk-dicom provider) — published-data
  comparison pending.*

### Physics / numerics

- **G-1 (physics):** *Closed through H-011 and H-011b.* Photon attenuation **relations**
  implemented and analytically verified in `helios-physics`: Beer–Lambert
  transmission, half-value layer, `μ = (μ/ρ)·ρ`, and first-order HU→density CT
  calibration (property/value-semantic tests: `T(HVL)=½`, `T(0)=1`,
  water/air/bone calibration points, f32 genericity). H-011b adds the selected
  NIST dry-air, liquid-water, and cortical-bone mass-attenuation knots over
  10 keV–20 MeV with an explicitly bounded log-linear interpolation contract.
  The five table value, boundary, interpolation, and invalid-domain tests pass
  on the final head. An electron-transport model remains a separate algorithm
  item, not a mass-attenuation table residual. *Evidence tier: analytical
  relations plus primary-source table-value tests.*
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
- **G-7 (toolchain):** Helios 0.1.0 declares Rust 1.95, matching the merged
  provider graph. `rust-toolchain.toml` remains `stable` rather than an exact
  channel pin; the configured Rust 1.95 package gates are the current evidence.

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
  CPU (67 tests).* Remaining: throughput benchmark vs VoLO pending. (H-010b fused HU→μ
  kernel delivered 2026-07-02 via the hephaestus ADR-0004 authored-kernel
  seam — consumer-side `GpuAttenuationMapper`, differential vs solver 9/9.)

### Testing / tooling

- **G-8 (coverage):** No `cargo-llvm-cov` run yet; >80% core-logic coverage target
  unmeasured. Applies from first `[minor]`. → measure after H-003.

## Closed gaps

- **G-9 (CI):** CI now runs format, warning-denied Clippy, configured Nextest,
  doctests, warning-clean rustdoc, RustSec audit, and cargo-deny
  license/source policy. Benchmark classification is separately owned by the
  exact Atlas gate recorded under G-27.
- **G-6 (build hygiene):** Helios inherits the shared `D:/atlas/target` build dir
  via `repos/.cargo/config.toml`; no per-repo target. Verified this session.

## Integrated-workflow status

- **Integrated imaging-delivery workflow (clinical-realism gate): demonstrated end-to-end
  on synthetic/self-consistent data (H-041).** `helios-simulation/tests/end_to_end.rs`
  runs a shared μ through both branches (Radon→FBP→registration; helical MLC delivery→
  divergent-fan dose→scatter→DVH/gamma) with self-consistency oracles. What remains for
  the *clinical* therapy gate is orthogonal and environment-blocked: a licensed CT/plan
  dataset and an external MC reference engine (G-16), plus the anisotropic CC kernel
  (H-020g). The workflow *plumbing* across all layers is verified.

## Concurrent-agent status

- **RESOLVED (next cycle).** The peer `mnemosyne-arena` breakage that blocked `helios-gpu`
  and the moirai consumption last cycle has been reconciled by the peer (new mnemosyne
  commits landed; the crate compiles). The full workspace builds green again, and the
  moirai consumption (H-021b) was re-applied and verified — `simulate_helical_sinogram`
  dispatches per-projection work through moirai's `Adaptive` policy. Kept per the
  concurrent-agent discipline: the peer crate was never touched; the change was designed,
  reverted to stay green while blocked, then re-landed once upstream compiled.

## Residual risk register

- **G-20 (H-011b table interpolation).** The embedded dry-air, liquid-water,
  and cortical-bone `μ/ρ` values are transcribed from the NIST X-Ray Mass
  Attenuation Coefficients tables at their common 10 keV–20 MeV knots. The
  range excludes the selected tables' absorption-edge rows, and Helios defines
  between-knot behavior as native-precision log-linear interpolation. This is
  deliberately not represented as XCOM output: NIST documents log-log cubic
  spline fitting and explicit edge handling for XCOM. Exact table knots,
  boundaries, and the interpolation identity are value-semantic tests; an
  independent clinical-spectrum or Monte-Carlo validation remains outside this
  data-loading slice. *Evidence tier: source-value and analytical-contract
  tests.*

- Atlas upstream APIs may drift (multi-repo co-evolution); Helios pins the local
  synchronized checkout via `[patch]` and commits `Cargo.lock`. `ritk-dicom` is now
  consumed (H-004b) and is **skew-free** (no leto/mnemosyne/themis/eunomia cluster —
  only anyhow/arrayvec plus the RITK DICOM provider and `ritk-codecs`), so it needed no
  patch-cluster work. Helios has no direct dicom-rs dependency; the provider's
  parser implementation remains upstream-owned.
  Remaining ritk surfaces (`ritk-registration`) pull the burn stack and are heavier
  (G-5); add cross-repo contract tests as each is consumed.
- **G-18 — RESOLVED (H-043b).** The residency step landed: hephaestus gained a volume
  ray-integral kernel (`ray_line_integrals`) and `helios_gpu::GpuProjector` keeps μ
  on-device, projecting whole sinograms per dispatch. Measured 171×/371× vs the
  single-thread CPU projector (report `2026-07-02-gpu-projection-throughput.md`);
  differential per-ray agreement within a derived 1e-3 f32 bound. The *elementwise*
  `exp(−τ)` path remains transfer-bound by physics (documented; use the resident
  pipeline instead). "VoLO-competitive" is still unclaimable (no reference engine).
- **(historical) G-18 (performance, GPU transfer-bound).** The GPU-vs-CPU study (H-043,
  `validation_reports/2026-07-01-gpu-transmission-throughput.md`) shows the isolated
  `beam_transmission_into` kernel is memory-/transfer-bound: even on an RTX 5080 it
  reaches only ~0.5–0.72× a single-threaded CPU loop, because every call round-trips the
  buffer over PCIe for ~1 flop/element. This is a correct roofline result, not a defect.
  GPU throughput requires an **on-device fused pipeline** (H-043b) that keeps τ resident
  across HU→μ / projection / transmission so one CT upload + one sinogram download
  amortize many kernels. Until then the GPU path is a differentially-correct reference,
  not a speedup. The "competitive with VoLO-class throughput" gate additionally needs an
  external VoLO reference not available here. *Evidence tier: empirical (criterion, this
  machine).*
- **G-17 (tooling, coverage gate — link unblocked, attribution still empty).** Refined:
  the original blocker (the mingw `ld` bfd linker failing on `__llvm_profile_runtime` /
  profiler builtins) **is fixable** — `RUSTFLAGS="-Clink-arg=-fuse-ld=lld"` (LLVM `lld`,
  present in the MSYS2 ucrt64 toolchain) links the instrumented binaries, and the full
  suite runs under instrumentation (183 tests pass, 356 `.profraw` generated;
  `LLVM_COV`/`LLVM_PROFDATA` point at the MSYS2 llvm-cov 22.1.4 ≈ rustc-LLVM 22.1.3).
  A *distinct* secondary issue remains and is now **conclusively diagnosed**: source/
  region *attribution* is broken on this GNU target — `cargo llvm-cov report` gives
  **0 regions** and `grcov 0.10.5` (which uses its own profraw parser) gives an empty
  file table / **NaN%**, from the *same* 145 profraw. Two independent tools failing
  identically confirms the coverage-map is not read from the mingw (`x86_64-pc-windows-
  gnu`) instrumented binaries — a toolchain-level limitation, not a tool bug. Coverage %
  is therefore **not obtainable on this host**; it requires `x86_64-pc-windows-msvc` or a
  Linux CI container (H-060 re-scoped to CI). Test breadth is high (189 value-semantic
  tests across the CPU crates) but the coverage number is unquantified — not fabricated.
- Physical constants (G-2) are CODATA-2018/ICRU-90 values verified by inter-constant
  derivation tests, not by an external authoritative fetch this session; values are
  standard and cross-checked, but a future audit should confirm against the live
  NIST database.
- **G-16 (dose model fidelity, H-020d/H-020e).** *Partially addressed.* Stage 1
  (`deposit_ray_terma`/`accumulate_delivered_dose`) deposits primary terma along
  **parallel** beamlets; stage 2 (`scatter_superposition`, H-020e) now spreads it with
  a **separable-isotropic** deposition kernel, so lateral penumbra and depth build-up
  are present and energy-conserving (verified). The beam geometry now supports a
  divergent point-source fan (H-020f, `BeamGeometry::PointSource`; verified parallel
  limit + multi-row divergence). Still approximate vs a clinical collapsed-cone dose:
  inverse-square fluence falloff along the divergent fan is now modelled (H-020g,
  `deposit_ray_terma_diverging`; verified SAD→∞ limit + entry/exit steepening). The
  remaining approximation is the scatter kernel: separable-isotropic, not the
  anisotropic forward-peaked beam-aligned CC kernel, tracked as H-020h. Sufficient to exercise DVH/gamma
  on self-consistent phantoms; the therapy gamma/DVH clinical-agreement gate still
  needs the H-020g kernel upgrade AND a licensed real CT dataset AND an external
  Monte-Carlo/reference dose engine (VoLO/TOPAS/GATE/EGSnrc) — the last of which is
  **not runnable in this environment**, so that specific gate cannot be closed here and
  will not be fabricated. Evidence tier: analytical oracles (conservation, identity
  differential, symmetry); NOT validated against a reference dose engine.
