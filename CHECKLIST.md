# Helios Checklist (tactical)

**Sprint target version:** `0.1.0` (oriented-grid provider convergence)
**Current phase:** Phase 2 (Execution). Sprint 1 domain core complete
(`helios-core`/`math`/`domain`); Sprint 2 opened with `helios-physics`.

## Owner: claude-helios

## Owner: Codex

## Codex — H-072 Python binding security [patch] — review

- [x] Upgrade the isolated Python binding boundary to PyO3 0.29.0, the first
      release closing both applicable RustSec advisories.
- [x] Replace the removed GIL-release API with `Python::detach` while
      preserving the existing Rust-core call and value-semantic Python API.
- [x] Add a hosted job that builds the abi3 extension and runs the Python
      contract suite against the installed wheel.
- [x] Pass exact-head local binding verification: warning-denied all-feature
      Clippy, 273/273 Nextest cases, doctests, rustdoc, and 13/13 Pytest cases.
- [ ] Pass exact-head hosted binding verification.

## Codex — H-071 benchmark gate repair [arch] [patch] — review

- [x] Replace the copied same-run Python comparison with the exact Atlas
      `9bfb722` phase-replicated Criterion gate.
- [x] Hold the candidate benchmark harness constant across baseline and
      candidate revisions, run ABBA followed by BAAB on one runner, and fail
      closed across all four report universes.
- [x] Invoke each declared Criterion binary directly so measurement arguments
      cannot be consumed by workspace library test harnesses.
- [x] Recover the linked candidate's replicated scatter regression in the
      production loop: const-generic axis specialization is bitwise-equivalent
      and lowers local 32³/64³ medians by 50.46%/51.02%.
- [x] Route native tests through the committed Nextest profile, retain
      doctests separately, and enforce RustSec/license/source policy.
- [x] Pin CI provider materialization to Atlas `05b7f5d`, whose registered
      graph contains Proteus, Asclepius, and the reconciled Hephaestus revision.
- [ ] Pass exact-head hosted CI.

## Codex — H-068 Aequitas domain units [arch] — done 2026-07-19

- [x] Store `EnergyMeV` and `VoxelSpacingMm` with Aequitas `Energy<f64>` and
      `Length<f64>` while retaining Helios validation at construction.
- [x] Preserve the MeV and millimetre public contracts through explicit unit
      conversion, compile-time layout assertions, and bounded round-trip
      properties derived from machine epsilon.
- [x] Keep `HounsfieldUnit` Helios-owned because its calibrated scale is not a
      linear SI quantity, and record the ownership decision in ADR 0001.
- [x] Verify formatting, warning-denied all-target Clippy, 17/17 configured
      Nextest tests, doctests, and warning-clean rustdoc for `helios-core`.

## Codex — H-067 Apollo lock refresh [patch] — done 2026-07-18

- [x] Replace the invalid stale one-line Apollo edit with a complete Cargo
      resolution selecting `apollo-fft` 0.25.0, Eunomia 0.4.0, Leto 0.38.2,
      and Hephaestus 0.17.0 while deleting the transitive `num-complex` package.
- [x] Run locked metadata, formatting, warning-denied workspace Clippy,
      configured workspace Nextest, doctests, and warning-clean rustdoc.
- [x] Merge the verified lock-only change and advance the Atlas Helios gitlink
      without staging fresh RITK, Themis, or root package-manager work.

**Evidence:** locked metadata and format pass; warning-denied
all-target/all-feature workspace Clippy passes; configured Nextest is 272/272;
all ten Rust library doctest targets pass with zero examples; workspace
rustdoc is warning-clean.

## Codex — H-066 direct dependency ownership [patch] — done 2026-07-17

- [x] Preserve the stale peer's removal of the unused workspace
  `num-traits` declaration; no Helios manifest or source imports it directly.
- [x] Verify the complete workspace with `cargo check --workspace --locked`.
- [x] Audit the inverse dependency tree: remaining `num-traits` paths are
  provider-owned transitive requirements through Eunomia, Gaia, Half, WGPU,
  and their dependencies, not a Helios-owned edge.

## Codex — H-011b NIST X-ray mass-attenuation data [minor] — done 2026-07-17

- [x] Pin the authoritative NIST X-ray source, supported material set, energy
  grid, and absorption-edge interpolation contract in Rustdoc and the audit.
- [x] Add one typed, allocation-free material/energy lookup API that returns
  validated `MassAttenuation` values without a consumer-side table copy.
- [x] Add source-value, boundary, interpolation, and `f32`/`f64` tests; run
  package formatting, warning-denied Clippy, nextest, doctests, and rustdoc.
  `helios-physics`: fmt and Clippy clean; nextest 37/37; doctests 0/0;
  rustdoc clean; semver checks 196/196 (baseline `origin/main`, assumed minor).

### H-003d [minor] done — oriented `VoxelGrid`

- [x] Consume Leto 0.38's checked rotation-column contract through
  `helios-math`; store one local-index-to-world `Isometry3` in `VoxelGrid`.
- [x] Preserve the axis-aligned constructor as the identity-pose case; verify
  rotated f32/f64 index/world mapping and inverse round trips.
- [x] Transform CPU projector/deposition clipping into scaled-index space;
  serialize and validate the full pose through HDF5; reject the unsupported GPU
  pose before upload rather than dropping rotation.
- [x] Regenerate the locked 0.1.0 provider graph under Rust 1.95. Package
  `fmt`, warning-denied Clippy, doctests, rustdoc, and nextest pass; nextest
  runs 104/104 domain/solver/GPU tests, including live-GPU differentials.
  Workspace example compilation and `cargo fmt --all --check` pass. Semver
  checks pass 196/196 each for `helios-domain`, `helios-math`,
  `helios-solver`, and `helios-gpu`. H-004d remains blocked on RITK's named
  `ImageOrientationPatient` attribute provider.

## Codex — H-062 DVH threshold-query audit [patch] — done 2026-07-15

- [x] Confirmed the sorted-sample invariant and established the fixed-workload
  Criterion comparison in `crates/helios-analysis/benches/dvh_queries.rs`.
- [x] Replaced repeated full scans with the standard-library binary partition
  over the canonical sorted slice; the query path allocates nothing. Samples
  containing NaN retain the prior filter semantics through an explicit fallback.
- [x] Added finite threshold boundary, NaN, and existing f32 genericity tests.
  Focused nextest passes 34/34.
- [x] Criterion evidence is recorded in
  `validation_reports/2026-07-15-dvh-query-optimization.md`; warning-denied
  Clippy, doctest, and rustdoc pass for `helios-analysis`. The workspace-wide
  DICOM graph mismatch identified during this audit is resolved by H-063.

## Codex — H-063 DICOM provider graph alignment [patch] — done 2026-07-15

- [x] Changed the direct Helios `dicom` dependency from 0.8 to 0.10, matching
  the `ritk-dicom` provider and eliminating duplicate `dicom-core` types at
  `helios-domain`'s typed attribute boundary.
- [x] Regenerated `Cargo.lock`; the resolved DICOM packages are all 0.10.0 and
  the workspace no longer carries a DICOM 0.8 package.
- [x] Locked workspace example check, all-target all-feature Clippy, workspace
  nextest, doctests, and rustdoc pass. The DICOM-focused domain suite passes
  41/41, including synthetic slice/series round trips and typed error paths.
- [x] The former four E0308 errors in `crates/helios-domain/src/dicom.rs` are
  absent after the provider graph alignment. Evidence tier: compile-time
  dependency/type verification plus value-semantic nextest coverage.

## Codex — H-064 ritk-dicom provider-boundary correction [arch] — done 2026-07-15

- [x] Removed Helios's direct `dicom` dependency and direct `dicom::core` /
  `dicom::object` imports from production and test code.
- [x] Production reads now use `ritk-dicom` canonical tags,
  `DicomAttributeRead`, parser, transfer-syntax, and frame-decoder contracts.
- [x] Synthetic slice/series tests now write a minimal Part 10 fixture from
  bytes and exercise the real `ritk-dicom` parser/decoder path without a
  consumer-side DICOM object model.
- [x] Re-ran the locked workspace gates and verified the documentation and
  lockfile: workspace examples, all-target all-feature Clippy, 261/261
  workspace nextest tests, doctests, rustdoc, and provider-backed DICOM checks
  pass.

## Codex — H-065 Moirai 0.3 provider-graph refresh [patch] — done 2026-07-15

- [x] Regenerated `Cargo.lock` after Moirai 0.3.0 retired the unowned NUMA
  iterator and benchmark. Helios does not consume that removed surface.
- [x] `cargo check --workspace --examples --all-features` resolves every
  Moirai package at 0.3.0 and compiles the complete example graph.
- [x] Evidence tier: compiler-checked dependency resolution and example
  compilation; this lock-only change adds no Helios behavior.

### H-061 done — all runnable examples and DICOM graph audited (2026-07-14)

The three existing examples (`validate_foundation_units`,
`voxel_grid_construction`, and `tomotherapy_workflow`) compile with
`cargo check --workspace --examples --all-features` and
`cargo build --workspace --examples --all-features` against the synchronized
Atlas graph. The mdBook relative-link audit passes for `docs/book/SUMMARY.md`.
Helios consumes `ritk-dicom` for parsing, typed attributes, transfer-syntax
selection, and pixel decoding; Helios has no direct DICOM dependency.
`cargo metadata --locked` passes,
and the lockfile contains no `ndarray` package entry. Evidence tier: compile-
time/build verification plus manifest/lockfile inspection; mdBook rendering is
not claimed because `mdbook` is not installed on this host.

### H-020k done — gaia-`Aabb` collimator field aperture + delivery collimation (jaw field-shaping + penumbra). Next: H-011b NIST μ/ρ / wire aperture into the dose pipeline / oriented-scatter perf — `todo`

`helios-domain::FieldAperture` (open field = gaia `Aabb`, geometric edge penumbra via box
SDF; `contains` → `Aabb::contains_point`) + `helios-simulation::collimate_frames` (scales
per-leaf fluence by aperture transmission at `(lateral_offset, couch_mm, 0)`). Oracles:
centre→1/far→0/edge→0.5/penumbra-ramp/monotone/typed-errors/f32; field-shaping consumer
(narrow aperture → edge leaves 50%, outside 0, machine state preserved, never increases
fluence). Deepens mandated gaia consumption; the last modelled delivery gap at fluence
level. (Full per-beamlet geometric occlusion in the dose ray-trace remains a follow-on.)

### (prior) H-033c done — per-structure plan evaluation demonstrated end-to-end

Integration test `per_structure_plan_evaluation_over_delivered_dose`: helical delivery →
beam-following collapsed-cone dose → central-PTV vs off-axis-OAR masked DVH → gEUD → TCP
(PTV) / NTCP (OAR). Clinical-plausibility oracles: target hotter than OAR (rotational
convergence), PTV gEUD > OAR gEUD, PTV TCP>0.5, OAR NTCP<0.5, all well-formed.

### (prior) H-033/H-033b done — radiobiology metrics + per-structure DVH outcome methods

Superseded by H-076: Asclepius now owns these laws, and Helios retains only the
zero-copy DVH receiver boundary.

New `helios-analysis::radiobiology`: `generalized_eud` (**promoted from planning** — a dose
metric, no longer gated behind the `autodiff` feature; now generic over Scalar),
`tcp_logistic` (Niemierko), `ntcp_lkb` (Lyman–Kutcher–Burman, via eunomia `erfc`). Oracles:
gEUD power-mean bounds/monotonicity, TCP [0,1]/0.5-at-TCD50/slope, **NTCP vs the normal CDF
at published Φ(±1)=0.8413/0.1587**, f32. Architectural realignment: the metric now lives in
analysis (SoC); planning's EUD objective is unchanged and its tests use an independent inline
gEUD oracle, so planning keeps its lean core+math deps (no geometry-stack pull).

### (prior) H-031d done — generalized-EUD biological objective on the coeus tape

Updated by H-076: `asclepius-coeus` now owns the gEUD tape construction.

`generalized_eud` (gEUD = (mean(D^a))^(1/a)) + `EudPenalty`/`eud_objective_gradient_autodiff`
— a one-sided gEUD hinge whose gradient w.r.t. beam weights flows by reverse-mode AD
through matmul/pow/mean (no closed form). Oracles: power-mean bounds + monotonicity +
uniform invariance; tape value == analytic gEUD; **autodiff gradient == central finite
difference** (differential oracle over the whole tape); inactive-hinge zero; typed errors.
Concurrency note: a transient peer coeus WIP (uncommitted topk.rs deletion) briefly broke
the build; detected, left untouched, self-reconciled. A test threshold I first wrote was
analytically unsound (gEUD→max is slow, (1/N)^{1/a}); corrected to rigorous power-mean
bounds, not weakened.

### (prior) H-041c done — beam-following poly-energetic dose demonstrated end-to-end

`tomotherapy_workflow` example upgraded to `accumulate_delivered_dose_anisotropic` +
`CollapsedCone::poly_forward_peaked`; ran it and **inspected the dose render**
(centrally-concentrated with smooth penumbra — correct rotation-averaged forward-peaked
helical delivery; recon μ +0.1%, self-gamma 100%). Added integration test
`beam_following_poly_energetic_dose_end_to_end`: non-negativity, energy conservation
(dose≤terma, >85% retained), positive DVH, 3%/2 mm self-gamma 100%. Examples build gate
(`cargo build --examples`) green.

### (prior) H-020j done — poly-energetic (beam-hardened) collapsed-cone kernels

**helios-solver**: `poly_forward_peaked_kernel` + `SpectralComponent` — energy-fluence-weighted
convex combination of monoenergetic `forward_peaked_kernel`s (beam hardening: harder
components reach farther downstream). **helios-simulation**: `CollapsedCone::poly_forward_peaked`
(two constructors share a private `from_beam_kernel` SSOT). Oracles: single-component reduces to
mono (kernel 1e-15, delivered dose 1e-13), weight scale-invariance + Σ=1, harder spectrum →
higher downstream/upstream ratio + downstream centroid shift, empty→identity, f32.

### (prior) H-020i done — beam-following anisotropic delivered dose (rotated cone axes)

**helios-solver**: `directional_convolve` (oriented 1-D convolution along any unit vector,
trilinear gather) + `oriented_forward_scatter` (forward-peaked collapsed cone along an
arbitrary beam direction, Gram–Schmidt lateral basis). **helios-simulation**:
`CollapsedCone` + `accumulate_delivered_dose_anisotropic` — scatters each frame's terma
along that frame's gantry direction before summing, so forward-peaking follows the
rotating beam; per-frame beamlet deposition extracted to a shared `deposit_frame_terma`
(SSOT, `accumulate_delivered_dose` behavior unchanged). Oracles: axis reduction to the
separable pipeline (1e-10), oblique downstream>upstream + lateral symmetry, conservation,
delivered-dose centroid shifts downstream under forward-peaking, linearity, f32. This
wires the H-020h capability into actual delivered dose at real gantry angles.

### (prior) H-020h done — anisotropic beam-aligned collapsed-cone scatter (grid-axis)

`forward_peaked_kernel` (asymmetric up/downstream ranges, Σ=1) +
`anisotropic_scatter_superposition` (forward-peaked along the beam axis, symmetric
lateral) on a generalized `convolve_axis_at`. The anisotropy oracle caught a real
defect — the gather was correlation, not convolution (invisible to symmetric kernels,
inverts asymmetric ones) — fixed with symmetric results bitwise-unchanged. Verified:
exact reduction to the symmetric case, downstream > upstream, lateral symmetry,
interior conservation, f32/axis-select. The last modelled dose-fidelity gap at the
axis-aligned level is closed; rotated per-gantry cones = H-020i.

### (prior) H-031c done — non-quadratic DVH-penalty objective + optimizer

`DvhPenalty` + `dvh_objective_gradient_autodiff` + `optimize_beam_weights_dvh` (feature
`autodiff`): one-sided underdose/overdose penalties with the gradient from the coeus
tape (`relu` kinks via reverse-mode AD, weights as `[1]` constant `Var`s, one backward).
Verified: tape grad == hand sub-gradient (1e-12), value cross-check, zero inside the
band, and the optimizer picks the OAR-sparing beamlet (target ≥ floor, OAR ≤ ceiling).
This is the no-closed-form capability coeus was mandated for.

### (prior) H-010b done — GPU HU→μ mapper over the Hephaestus authored-kernel seam

`helios_gpu::GpuAttenuationMapper` computes `μ = max(fma(scale, HU, offset), 0)` in
one WGSL dispatch authored as a Helios consumer kernel over
`hephaestus_core::KernelInterface` + `KernelSource<Wgsl>`. No type-specific
affine-clamp helper was required upstream: the Hephaestus seam is the support
surface. Evidence tier: differential/empirical validation. Current focused checks:
`rustup run nightly cargo nextest run -p helios-gpu attenuation` passes 5/5
(closed-form clamp, solver oracle via `helios_solver::attenuation_map`, typed error
paths) and `rustup run nightly cargo nextest run -p hephaestus-core -p
hephaestus-wgpu stream` passes 8/8. No separate throughput report is added for this
single elementwise upload/download lane; H-043/H-043b already record that isolated
elementwise GPU calls are transfer-bound and that throughput evidence belongs to the
resident on-device pipeline.

### (prior) H-031b RESOLVED — coeus autodiff consumed (last mandated component)

The peer's moirai-core refactor landed (2451715) → the parked module re-landed and
verified: `objective_gradient_autodiff` (feature `autodiff`) — coeus tape gradient ==
exact hand gradient `Aᵀ(Ax−d)` within 1e-12, zero at the least-squares optimum, typed
shape errors. **All mandated Atlas components are now consumed** (ritk, gaia,
hephaestus, moirai, coeus, consus, leto, hermes, eunomia; mnemosyne/themis + apollo
transitively). 217 `--all-features` tests pass, clippy/fmt clean. Detect-and-reconcile
note: peers refactored hephaestus ops into a dialect seam (`UnaryExpr<Wgsl>`) migrating
our `ExpNegOp`/volume kernel compatibly — helios builds unchanged.

### (prior) H-031b advanced — apollo feature-unification fix; blocked on moirai WIP

Root-caused the coeus consumption failure across three layers: (1) leto-ops E0034 ×332
under `leto/ndarray-compat` → trigger was apollo's **vestigial** workspace-level feature
request (no apollo lib/test code needs it) → **fixed upstream** (apollo f1ddf7a, whole
apollo workspace verified `--all-targets` without it); (2) with that resolved, the next
layer is **peer WIP in moirai-core** (staged `dtype/` deletion + `mod security` w/o the
file — actively churning) — their claimed scope, sequenced behind, NOT touched. The
complete autodiff module (tape gradient + differential test vs the exact hand gradient)
is parked (session scratchpad) for immediate re-landing; `DoseInfluence::rows()` +
the coeus `[patch]` block landed now. 209 `--all-features` tests pass.

### (prior) H-043b RESOLVED — resident GPU projector, 171×/371× vs CPU

Upstreamed `ray_line_integrals` to hephaestus (volume ray-integral kernel, 792ccc3 +
9354260, 4 live-GPU oracles; peer WIP untouched via disjoint files) and consumed it as
`helios_gpu::GpuProjector` (μ resident on-device, batched sinogram projection, mm→cm at
the boundary). **171× @ 90×128 / 371× @ 360×256 sinogram vs single-thread CPU** on a
128³ volume; per-ray differential vs `forward_project_ray` within 1e-3 (derived f32
bound). Closes G-18; GPU-scaling gate demonstrated on the pipeline workload. 200 default
tests pass; clippy/fmt clean. Bench + report committed.

### (prior) H-043b step 1 — fused ExpNegOp upstreamed + consumed

Upstreamed `ExpNegOp` (`exp(−x)`) to hephaestus-wgpu (commit 669a9b3; live-GPU contract
test; peer WIP in reduction.rs/device.rs untouched — disjoint-scope concurrent work, not
deferred). `beam_transmission_into` = one dispatch, no intermediate buffer: GPU +30% at
4M (373→485 Melem/s) but still PCIe-bound at 0.66–0.73× CPU (report addendum). The
remaining GPU-beats-CPU path is the resident on-device μ→projection→transmission pipeline.

### (prior) H-048 done — perf/consolidation pass (8.3× scatter kernel; MM_PER_CM SSOT)

`Volume::as_slice` zero-copy accessor (documented layout contract); `convolve_axis`
strided-slice rewrite — **8.3×/7.4× at 32³/64³, bitwise-identical** (baseline report in
validation_reports/); `save_volume_hdf5` via the slice view; `MM_PER_CM` 5 duplicates →
one helios-core SSOT constant. 207 `--all-features` tests pass, clippy/fmt clean.

### (prior) H-046 done — consus HDF5 volumetric storage (mandated consus consumed)

`helios-domain::{save_volume_hdf5, load_volume_hdf5}` (feature `storage`) archive a
`Volume` (data + grid geometry) to standard HDF5 via consus-core/hdf5/io ([patch]ed to
the local checkout; skew-free). Verified: bitwise f64 round-trip, HDF5 superblock
signature, f32 exactness, typed error paths. Adds `HeliosError::Storage`. 198 default /
**207 `--all-features`** tests pass. Consumed mandated components now: ritk, gaia,
hephaestus, moirai, consus, leto, hermes, eunomia. Remaining: coeus (H-031b),
mnemosyne/themis (indirect via leto), apollo.

### (prior) H-021b done — moirai-parallel helical-projection dispatch

`simulate_helical_sinogram` dispatches per-projection work via moirai `Adaptive`
(deterministic, order-preserving; verified at 256 projections).

### (prior) H-047 done — geometric ROI masks (per-structure DVH)

`spherical_mask`/`box_mask` predicates for `Dvh::from_volume_masked`.

### (prior) H-045 done — portal (EPID) exit dosimetry

`frame_portal_fluence` — per-leaf transmitted fluence; shared `beamlet_ray`/`gantry_basis`.

`helios-simulation::frame_portal_fluence` computes per-leaf transmitted fluence
`Ψ_leaf·exp(−τ_leaf)` for a delivery frame (delivery-verification image), sharing the
`beamlet_ray`/`gantry_basis` geometry with dose accumulation (consolidated into a
`Beamlet` struct). Verified: full transmission at μ=0, Beer–Lambert attenuation, closed
leaf 0, darkening with μ, f32. 193 default / 198 `--all-features` tests pass. The
delivery-side imaging surface (MVCT acquisition, portal dosimetry) is now complete.

### (prior) H-044b done — NCC registration (robust IGRT on low-texture images)

`register_translation_ncc` — NCC over overlap, rejecting zero-variance overlaps; recovers
a known shift on the flat-background spike phantom where SSD is ambiguous.

### (prior) H-032e + G-17 — homogeneity index; coverage link unblocked (lld)

`Dvh::homogeneity_index` (ICRU-83); `RUSTFLAGS=-Clink-arg=-fuse-ld=lld` links the
instrumented build (183 ran), region attribution pending (H-060).

### (prior) H-032c done — local-normalization gamma + low-dose cutoff

`gamma_index_3d_local` — global + local 3%/2 mm gamma; local γ=1.0 vs global γ=0.2 in
low-dose, cutoff exclusion.

### (prior) H-020g done — inverse-square fluence falloff on the divergent fan

`deposit_ray_terma_diverging` scales per-segment terma by `(SAD/r)²`; verified SAD→∞
limit + entry/exit steepening. Remaining dose fidelity: anisotropic CC kernel (H-020h).

### (prior) H-032b done — structure-masked (per-PTV/OAR) DVH

`Dvh::from_volume_masked` — the per-structure DVH clinical DVH-agreement metrics use.

`helios-analysis::Dvh::from_volume_masked` builds a per-structure DVH from a voxel-mask
predicate — how clinical DVH-agreement metrics are evaluated (target vs OAR). `from_volume`
consolidated as the unmasked case. Verified: disjoint target/OAR masks give distinct means,
single-voxel point DVH. 180 default / 185 `--all-features` tests pass. Remaining analysis:
RT-struct ROI rasterization (ritk) + local-normalization gamma (H-032c).

### (prior) H-041b done — runnable end-to-end example with inspected PNG renders

`examples/tomotherapy_workflow.rs` — full pipeline demo, renders inspected; recon μ +0.1%,
self-gamma 100%.

`helios-simulation/examples/tomotherapy_workflow.rs` runs the full pipeline and renders
`ct/mu/recon/dose.png` (Output & visual verification — **inspected**: phantom, FBP
recovery, central rotational dose falloff). Reports recon μ +0.1% of target, self-gamma
100%. `cargo build --examples` covers it in CI; output dir gitignored. 178 default / 183
`--all-features` tests pass. This completes the Sprint-5 `examples/` required artifact.

### (prior) H-041 done — end-to-end workflow validation (integration test)

Shared μ drives imaging + therapy branches with self-consistency oracles across all
layers.

`helios-simulation/tests/end_to_end.rs`: one shared μ (CT→`attenuation_map`) drives both
the imaging branch (Radon→FBP→water-ROI μ within 20%; registration recovers a known
shift) and the therapy branch (helical MLC delivery→divergent-fan dose→scatter→DVH mean>0,
3%/2 mm self-gamma 100%). Proves domain/physics/solver/imaging/analysis/simulation compose
across seams (test-only dev-deps, no production cycle). 178 default / 183 `--all-features`
tests pass. This is the integrated imaging-delivery clinical-realism workflow on
synthetic/self-consistent data; runnable example + Python on real DICOM = H-041b.

### (prior) H-044 done — IGRT rigid translation registration

`register_translation` — whole-voxel couch-shift estimate; recovers a known shift exactly.

### (prior) H-030c done — SIRT iterative reconstruction

`sirt_reconstruction` (normalized SIRT, non-negativity-projected) — a second MVCT
reconstruction method robust where FBP streaks; shared `back_project_rows` with FBP.

### (prior) H-020f done — divergent point-source fan (TomoTherapy beam geometry)

`BeamGeometry::PointSource` — beamlets diverge from a focal spot (verified parallel
limit + multi-row divergence).

`helios-simulation::BeamGeometry` seam: `accumulate_delivered_dose` now supports a
divergent point-source fan (`PointSource { source_axis_mm }`) alongside the parallel
approximation — beamlets diverge from a focal spot with depth (true TomoTherapy fan).
Verified: reduces to parallel as SAD→∞ (dose within 1e-4); off-axis beamlet sweeps ≥3
detector rows under divergence. Existing parallel oracles intact. 169 default / 174
`--all-features` tests pass. Remaining dose fidelity: anisotropic beam-aligned CC
kernel + inverse-square falloff (H-020g).

### (prior) H-004c done — multi-slice DICOM series → 3-D HU Volume

`load_ct_series` stacks a real CT/MVCT series into a 3-D HU `Volume` (sorted by z,
uniform Δz). Real-input path complete (slice + series).

`helios-domain::load_ct_series` (feature `dicom`) stacks a real DICOM CT/MVCT series
into a 3-D HU `Volume`: identical-in-plane-geometry validation, sort by
`ImagePositionPatient` z, derived uniform Δz (rejecting missing/duplicate slices),
`k`-stacking. Shares `read_slice`/`scatter_slice` with the single-slice loader. Verified
by a shuffled 3-slice synthetic round-trip + error paths. **A real CT series can now
drive the full pipeline** (μ-map → projection → recon → dose → metrics). 172 Rust tests
pass with `--all-features` (169 + 3 series). Remaining real-inputs: HU newtypes +
oriented pose (H-004d), RT-struct/registration (ritk-registration, needs burn).

### (prior) H-004b done — DICOM single-slice path (mandatory ritk consumed)

`load_ct_slice` via the `ritk-dicom` provider boundary. First consumption of the
mandatory ritk Atlas component.

### (prior) H-033b done — MVCT quantum-noise model (imaging noise/CNR sub-gate)

`add_quantum_noise` (Tyche counter-addressed normal sampling over photon statistics)
— the H-033 metrics now run on genuinely noisy reconstructions (MVCT
accuracy+noise+contrast quantified on synthetic).

### (prior) H-043 done — GPU-vs-CPU scaling study (performance instrument)

`helios-gpu/benches/transmission_throughput.rs` + `validation_reports/2026-07-01-gpu-
transmission-throughput.md` deliver the performance-gate measurement instrument with
real numbers (RTX 5080 vs Core Ultra 9 285K). **Honest finding (G-18):** the isolated
`exp(−τ)` kernel is transfer-bound — GPU reaches only ~0.5–0.72× a single-threaded CPU
loop because each call round-trips the buffer over PCIe for ~1 flop/element. GPU
throughput needs the on-device fused pipeline (H-043b); "competitive with VoLO" is not
claimed (no external reference). Benchmark is instrument-only; 160 lib tests unchanged.

### (prior) H-033 done — MVCT image-quality metrics (imaging gate, synthetic)

`helios-analysis::image_quality` adds the MVCT quality instruments — reconstruction
accuracy (`volume_rmse`/`volume_relative_l2_error`), noise (`roi_statistics`),
contrast (`michelson_contrast`), CNR (`contrast_to_noise_ratio`) — and an end-to-end
`helios-imaging` test that quantifies FBP disk-recon accuracy (interior mean within
15% of μ₀), background suppression, contrast (>0.85), CNR (>1). 160 Rust tests pass
(was 151: +8 metric oracles, +1 end-to-end). The MVCT accuracy/contrast gate is met on
synthetic phantoms; genuine quantum-noise validation (H-033b) and real-data (H-004b)
remain. **Coverage gate blocked (G-17):** `cargo llvm-cov` fails to link under the
MSYS2 GNU toolchain (profiler-runtime linker error) — coverage % unmeasurable here.

### (prior) H-020e done — collapsed-cone stage 2 (lateral scatter)

`scatter_superposition` (+ `symmetric_deposition_kernel`): separable 3-D kernel
superposition → lateral penumbra + build-up, completing the two-stage dose model.
Fidelity gaps (G-16): anisotropic CC kernel + divergent fan (H-020f).

### (prior) H-020d done — delivery→dose loop closed

`deposit_ray_terma` + `accumulate_delivered_dose`: per-frame per-leaf beamlets →
delivered-dose `Volume`, exact `w·(1−e^{−τ})` conservation oracle.

### (prior) H-040 done — 11/11 crate roster complete

`helios-python` (H-040): thin PyO3 wrappers over physics/planning, abi3-py39 wheel,
13 pytest equivalence tests green.

### (prior in-flight) H-004b `helios-domain` ritk DICOM load path — `todo`

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
`Isometry3` gains transforms), H-011d (exact Siddon), H-004b (ritk DICOM),
H-011b (NIST μ/ρ tables).

## Gate status (last run, H-020k — collimator field aperture)

| Gate | Result |
|------|--------|
| `cargo build` (whole workspace) | pass (all 11 crates) |
| `cargo nextest run` (focused H-010b) | `rustup run nightly cargo nextest run -p helios-gpu attenuation` pass: 5/5 (incl. live-GPU HU→μ solver differential) |
| `cargo nextest run --all-features` | **231 passed / 0 failed** (+7: oriented scatter + anisotropic delivery) |
| `cargo clippy -D warnings` / `cargo fmt --check` | clean (helios + hephaestus additions) |
| criterion `forward_projection_sinogram` | GPU **171×/371×** vs single-thread CPU (report committed) |
| criterion `scatter_superposition` | 8.3× @32³ / 7.4× @64³ vs recorded baseline |
| `cargo llvm-cov` (coverage %) | link unblocked via lld (183 ran instrumented); attribution empty on GNU target (G-17/H-060) |
| `pytest` (helios-python, maturin develop) | 13 passed / 0 failed |
| `cargo clippy -D warnings` | 0 code warnings |
| `cargo test --doc` | pass |
| `cargo fmt --check` | pass |
| `cargo llvm-cov` (coverage %) | **blocked** — GNU-toolchain profiler-runtime link failure (G-17) |
| GPU-vs-CPU perf (H-043) | instrument delivered; transmission kernel transfer-bound, GPU ~0.5–0.72× CPU (G-18) |

**11/11 crates — full roster delivered.** `helios-python` is a thin abi3-py39 PyO3
surface (`import helios`) over the physics/planning cores, GIL released around the
solve, verified by value-semantic pytest equivalence. The deterministic pipeline
(CT→μ→forward-projection/Radon→FBP recon; helical+MLC delivery; dose; DVH/gamma;
inverse planning) is coherent on synthetic phantoms.

Next (clinical-validation gates need real data + heavier Atlas backends): H-004b ritk
DICOM CT/MVCT load path (real inputs), H-020d delivered-dose accumulation
(delivery→dose loop), H-031b coeus-autodiff planning, H-021b moirai GPU orchestration,
H-040b numpy zero-copy Volume/attenuation/recon exposure.

Clinical-realism gate: helical synchronization ✓, binary-MLC leakage/tongue-and-
groove ✓, integrated imaging-delivery workflow (H-020c) ✓. Remaining: IGRT
registration workflows (H-041, needs ritk), delivered-dose accumulation (H-020d).
Crates 9/11; remaining helios-planning (coeus), helios-python (PyO3).

9/11 crates. Imaging round-trip works: CT→μ→Radon→FBP recovers μ. Remaining crates:
helios-planning (coeus inverse planning), helios-python (PyO3). Next: H-031 planning
or H-004b ritk DICOM (real inputs → clinical validation) or H-020b binary-MLC.

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
- **H-020h** anisotropic collapsed-cone kernel.
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
# H-076 review — Asclepius response ownership

- [x] Delete Helios's duplicate gEUD, logistic TCP, and Lyman NTCP functions.
- [x] Store the DVH response sample as Aequitas absorbed-dose quantities and
  borrow it directly into Asclepius without a conversion allocation.
- [x] Return Asclepius typed domain failures from all DVH response methods and
  migrate the end-to-end delivered-dose caller.
- [x] Delegate the Coeus gEUD tape to `asclepius-coeus` while retaining the
  Helios planning objective.
- [x] Advance Proteus and Helios to the merged Aequitas response-quantity
  identity and confirm the focused dependency graph has no duplicate Aequitas.
- [x] Keep the DVH benchmark instrument source-compatible with the historical
  scalar sample and the candidate's typed Aequitas storage.
- [x] Verify 84 focused analysis, planning, and simulation tests, including the
  analytical gEUD value, finite-difference planning gradient, typed failures,
  and masked PTV/OAR outcome workflow.
- [x] Complete exact locked all-target/all-feature check, warning-denied Clippy,
  270 workspace nextest tests, doctests, Rustdoc, examples, and supply-chain
  gates. Semver analysis classifies the intentionally removed public functions
  and fallible response signatures as the declared breaking migration.
- [x] Pin public Asclepius merge `ceb8b6d`, patch both workspace packages to
  the Atlas sibling checkout, and use Atlas merge `05b7f5d` for both checkout
  action provenance and provider-graph resolution.
- [ ] Pass exact-head hosted Rust, Python-binding, and benchmark verification.
