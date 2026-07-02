# Helios Checklist (tactical)

**Sprint target version:** `0.0.1` (Foundation → Sprint 2 physics/GPU)
**Current phase:** Phase 2 (Execution). Sprint 1 domain core complete
(`helios-core`/`math`/`domain`); Sprint 2 opened with `helios-physics`.

## Owner: claude-helios

### H-043b RESOLVED — resident GPU projector, 171×/371× vs CPU. Next: H-031b coeus-autodiff / H-020h anisotropic CC / H-010b GPU HU→μ — `todo`

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

`load_ct_slice` via `ritk-dicom` (dicom-rs, skew-free). First consumption of the
mandatory ritk Atlas component.

### (prior) H-033b done — MVCT quantum-noise model (imaging noise/CNR sub-gate)

`add_quantum_noise` (SplitMix64 photon statistics) — the H-033 metrics now run on
genuinely noisy reconstructions (MVCT accuracy+noise+contrast quantified on synthetic).

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
`Isometry3` gains transforms), H-011d (exact Siddon), H-010b (GPU HU→μ + throughput),
H-004b (ritk DICOM), H-011b (NIST μ/ρ tables).

## Gate status (last run, H-043b — resident GPU projector)

| Gate | Result |
|------|--------|
| `cargo build` (whole workspace) | pass (all 11 crates) |
| `cargo nextest run` (default) | **200 passed / 0 failed** (incl. live-GPU projector differential) |
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
