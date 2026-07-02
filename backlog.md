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
| H-004b | `helios-domain::load_ct_slice` (feature `dicom`): single-slice DICOM CT/MVCT → HU `Volume` via `ritk-dicom` (parse + rescale-calibrated decode) + geometry (Rows/Columns/PixelSpacing/ImagePositionPatient). Validated by a synthetic-DICOM round-trip through the real dicom-rs parser. **Mandatory ritk integration now consumed.** | [minor] | done | claude-helios | `crates/helios-domain/**`, `crates/helios-core/**` |
| H-004c | `helios-domain::load_ct_series`: multi-slice DICOM **series** stacking → 3-D HU `Volume` (parse+decode each slice, validate identical in-plane geometry, sort by `ImagePositionPatient` z, derive uniform Δz). Shared `read_slice`/`scatter_slice` with `load_ct_slice`. Verified by a shuffled 3-slice synthetic round-trip + empty/non-uniform error paths. | [minor] | done | claude-helios | `crates/helios-domain/**` |
| H-004d | `helios-domain`: `CtVolume`/`MvctVolume` HU-semantic newtypes + `ImageOrientationPatient` → oriented grid pose (pairs with H-003d oriented `VoxelGrid`) | [minor] | todo | — | `crates/helios-domain/**` |
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
| H-010 | `helios-gpu`: dispatch over `hephaestus_core::ComputeDevice` + hephaestus-wgpu. **Done** — `beam_transmission_into` (GPU `exp(-τ)` via NegOp+ExpOp) differentially validated vs CPU on a live adapter; `default_device`. | [minor] | done | claude-helios | `crates/helios-gpu/**` |
| H-010b | `helios-gpu`: GPU HU→μ kernel (needs a fused affine-clamp `UnaryWgslOp`, absent from hephaestus-wgpu's op set) differentially validated vs `helios-solver::attenuation_map`; GPU forward-projector; throughput benchmark vs CPU | [minor] | todo | — | `crates/helios-gpu/**` |
| H-011 | `helios-physics`: photon attenuation relations — `LinearAttenuation`/`MassAttenuation`, Beer–Lambert, HVL, HU→density calibration | [minor] | done | claude-helios | `crates/helios-physics/**` |
| H-011b | `helios-physics`: NIST XCOM μ/ρ data tables (energy-indexed, per material) loaded into `MassAttenuation` | [minor] | todo | — | `crates/helios-physics/**` |
| H-011d2 | `helios-physics`: Klein–Nishina Compton + Thomson cross-sections; `compton_mass_attenuation`/`electrons_per_gram` (μ/ρ derived from σ_KN, validated vs NIST water) | [minor] | done | claude-helios | `crates/helios-physics/**` |
| H-011d3 | `helios-physics`: KN differential cross-section + energy-transfer σ_tr (quadrature, self-validated vs closed-form total); `compton_mass_energy_transfer` kerma coefficient (validated vs NIST water μ_tr/ρ at 1 MeV) | [minor] | done | claude-helios | `crates/helios-physics/**` |
| H-055 | Decouple numeric/physics from geometry: `helios-math` `geometry` feature (default on); `helios-physics` builds without it. Keeps physics buildable during geometry-stack churn (G-14). | [arch] | done | claude-helios | `crates/helios-math/**`, `crates/helios-physics/**` |
| H-011c | `helios-solver::forward_project_ray`: ray-march optical depth ∫μ dl of a gaia `Ray` through a μ `Volume` (clip to grid `Aabb`, midpoint trilinear sampling). MVCT forward-projection / dose ray-trace core. | [minor] | done | claude-helios | `crates/helios-solver/**` |
| H-011d | `helios-solver`: exact Siddon voxel-DDA + oriented-grid (rotated pose) projection; full parallel/fan sinogram over a detector | [minor] | todo | — | `crates/helios-solver/**` |
| H-012 | `helios-solver`: GPU MVCT forward projector (Siddon/Joseph); CPU reference | [minor] | todo | — | `crates/helios-solver/**` |
| H-012b | `helios-solver`: HU→μ attenuation-map engine (CPU reference; differential oracle for the GPU kernel) | [minor] | done | claude-helios | `crates/helios-solver/**` |
| H-013a | `helios-solver::primary_fluence_parallel_x`: primary-transport stage — Beer–Lambert attenuated fluence Ψ=Ψ₀·exp(−∫μ dl), +x parallel beam. Also fixed projector optical-depth units (mm→cm). | [minor] | done | claude-helios | `crates/helios-solver/**` |
| H-013b | `helios-solver`: dose = TERMA ⊛ kernel (collapsed-cone). `dose_convolution_x` + `exponential_deposition_kernel` — exact oracles (delta identity, normalized-kernel interior conservation, physical build-up). **Verified** (G-14 resolved). | [minor] | done | claude-helios | `crates/helios-solver/**` |
| H-003c | `helios-math` re-exports adapted to the new `leto::geometry` API (Point2/Point3/Vector3/UnitVector3) + gaia Aabb/Ray; `VoxelGrid` simplified to axis-aligned. Restored full-workspace build. | [minor] | done | claude-helios | `crates/helios-math/**`, `crates/helios-domain/**` |
| H-003d | `helios-domain`: oriented `VoxelGrid` (DICOM `ImageOrientationPatient` cosines) once a rigid-transform primitive with `transform`/`inverse` exists upstream. | [minor] | todo | — | `crates/helios-domain/**` |

## Sprint 3 — Delivery

| ID | Item | Class | Status | Owner | Scope |
|----|------|-------|--------|-------|-------|
| H-020 | `helios-domain`: helical delivery kinematics (gantry rotation + couch velocity + pitch/time synchronization) — `HelicalDelivery<T>` | [minor] | done | claude-helios | `crates/helios-domain/**` |
| H-020b | `helios-domain`: binary-MLC `LeafOpenTimeSinogram` + `MlcModel` (leakage/transmission + tongue-and-groove effective fluence) | [minor] | done | claude-helios | `crates/helios-domain/**` |
| H-020c | `helios-simulation::simulate_helical_delivery`: integrated MLC + helical kinematics → time-ordered `DeliveryFrame` sequence (machine state + effective per-leaf fluence). | [minor] | done | claude-helios | `crates/helios-simulation/**` |
| H-020d | Per-leaf beamlet ray-trace: `helios-solver::deposit_ray_terma` (primary-energy terma deposition, exact `w·(1−e^{−τ})` conservation oracle) + `helios-simulation::accumulate_delivered_dose` (per-frame per-leaf beamlets → dose `Volume`). Closes the delivery→dose loop feeding DVH/gamma. | [major] | done | claude-helios | `crates/helios-solver/**`, `crates/helios-simulation/**`, `crates/helios-domain/**` |
| H-020e | Lateral scatter kernel: `helios-solver::scatter_superposition` — separable 3-D convolution of the delivered terma with normalized centred deposition kernels (`symmetric_deposition_kernel`), producing lateral penumbra + build-up. Identity-kernel differential vs the primary reference; interior energy conservation. | [major] | done | claude-helios | `crates/helios-solver/**`, `crates/helios-simulation/**` |
| H-020f | Divergent point-source fan: `BeamGeometry::{Parallel, PointSource}` seam in `accumulate_delivered_dose` — beamlets diverge from a focal spot (true TomoTherapy fan) instead of running parallel. Verified: reduces to parallel as SAD→∞; off-axis beamlet sweeps multiple detector rows (divergence). | [minor] | done | claude-helios | `crates/helios-simulation/**` |
| H-020g | Inverse-square fluence falloff along the divergent fan: `deposit_ray_terma_diverging` scales per-segment terma by `(SAD/r)²` from the focal spot; `BeamGeometry::PointSource` uses it. Verified: reduces to no-falloff as SAD→∞; steepens the entry/exit dose ratio. | [minor] | done | claude-helios | `crates/helios-solver/**`, `crates/helios-simulation/**` |
| H-020h | Anisotropic (forward-peaked, poly-energetic) beam-aligned collapsed-cone kernel replacing the separable-isotropic scatter approximation + per-leaf collimation via gaia MLC geometry; differential vs the isotropic reference | [major] | todo | — | `crates/helios-solver/**`, `crates/helios-simulation/**` |
| H-021 | `helios-simulation::simulate_helical_sinogram`: time-dependent helical MVCT acquisition — gantry rotation + couch translation drive the forward projector per projection (helix). CPU reference; moirai orchestration = H-021b. | [minor] | done | claude-helios | `crates/helios-simulation/**` |
| H-021b | `helios-simulation`: moirai-parallel projection dispatch; fan/cone-beam detector rows (full sinogram); motion modeling | [minor] | todo | — | `crates/helios-simulation/**` |
| H-022 | Binary-MLC leakage/transmission/tongue-and-groove model | [minor] | todo | — | `crates/helios-domain/**` |

## Sprint 4 — Planning & imaging

| ID | Item | Class | Status | Owner | Scope |
|----|------|-------|--------|-------|-------|
| H-030a | `helios-imaging::parallel_beam_radon` + `Sinogram`: MVCT forward projection (parallel-beam Radon). Validated vs analytical disk sinogram. | [minor] | done | claude-helios | `crates/helios-imaging/**` |
| H-030 | `helios-imaging::filtered_back_projection`: MVCT reconstruction (Ram-Lak FBP). Round-trip recovers disk μ (centre within 15%, background ~0). | [minor] | done | claude-helios | `crates/helios-imaging/**` |
| H-030c | `helios-imaging::sirt_reconstruction`: SIRT iterative reconstruction (normalized `x ← max(0, x + λ C⁻¹ Aᵀ R⁻¹(b − Ax))`) — robust to noise/sparse-angle where FBP streaks. Shared `back_project_rows` extracted (FBP + SIRT, no duplication). Verified: converges to its forward model (interior mean within 15% of μ₀), monotone error decrease, zero→zero, f32. | [minor] | done | claude-helios | `crates/helios-imaging/**` |
| H-030b | `helios-imaging`: OS-SEM/MLEM statistical reconstruction, portal dosimetry, IGRT registration workflows (via ritk). *(noise/contrast metrics = H-033/H-033b done; SIRT = H-030c done.)* | [major] | todo | — | `crates/helios-imaging/**` |
| H-031 | `helios-planning::optimize_beam_weights`: projected-gradient inverse planning (½‖Ax−d‖², x≥0) + `DoseInfluence`. Convex-convergence oracles. | [minor] | done | claude-helios | `crates/helios-planning/**` |
| H-031b | `helios-planning`: coeus-autodiff backend for general (DVH/biological, non-quadratic) objectives; replace the exact hand gradient. Needs coeus dep + patch set. | [major] | todo | — | `crates/helios-planning/**` |
| H-032 | `helios-analysis`: DVH (cumulative, Dx/Vx/mean) + 3D gamma index (Low, global norm) + pass rate | [minor] | done | claude-helios | `crates/helios-analysis/**` |
| H-032b | `helios-analysis::Dvh::from_volume_masked`: structure-masked (per-PTV/OAR) DVH via a voxel-mask predicate; `from_volume` consolidated as the unmasked case. Verified: target/OAR masks yield distinct means, single-voxel point DVH. | [minor] | done | claude-helios | `crates/helios-analysis/**` |
| H-032c | `helios-analysis::gamma_index_3d_local`: local-normalization gamma (`ΔD = criterion·D_r`) + low-dose cutoff (excludes points below threshold); shares one impl with the global variant. Verified: equals global for uniform dose, stricter in low-dose, cutoff exclusion. | [minor] | done | claude-helios | `crates/helios-analysis/**` |
| H-032d | `helios-analysis`: RT-struct ROI rasterization → DVH mask (via `ritk`), feeding `Dvh::from_volume_masked`; portal dosimetry | [minor] | todo | — | `crates/helios-analysis/**` |
| H-032e | `helios-analysis::Dvh::homogeneity_index`: ICRU-83 `HI = (D₂−D₉₈)/D₅₀` target plan-quality metric. Verified: 0 uniform, 1.92 ramp, zero-guarded. | [minor] | done | claude-helios | `crates/helios-analysis/**` |
| H-060 | Coverage % (G-17): **conclusively CI-only on this host.** `lld` unblocks the instrumented *link*; but source/region attribution is empty on `x86_64-pc-windows-gnu` for BOTH `cargo llvm-cov` (0 regions) and `grcov 0.10.5` (NaN%) from the same 145 profraw — a toolchain limitation. Obtain the % on `x86_64-pc-windows-msvc` or a Linux CI container. | [patch] | todo (CI) | — | CI/tooling |
| H-033 | `helios-analysis::image_quality`: MVCT quality metrics — reconstruction accuracy (`volume_rmse`, `volume_relative_l2_error`), noise (`roi_statistics` std), contrast (`michelson_contrast`) + CNR (`contrast_to_noise_ratio`). End-to-end recon-accuracy/contrast test in `helios-imaging` (FBP disk → metrics). | [minor] | done | claude-helios | `crates/helios-analysis/**`, `crates/helios-imaging/**` |
| H-033b | `helios-imaging::add_quantum_noise` — deterministic seeded MVCT quantum-noise model (`N=N₀e^{−τ}`, Poisson≈Gaussian draw, `τ'=−ln(N'/N₀)`) via a committed SplitMix64 PRNG + `Sinogram::from_readings`/`map_readings`. Validated: `Var(τ')≈e^{τ}/N₀` vs analytical photon statistics, noise↑ with attenuation, high-flux→clean, determinism; end-to-end noisy-recon std↑ and flux-scaling. | [minor] | done | claude-helios | `crates/helios-imaging/**` |

## Sprint 5 — End-to-end

| ID | Item | Class | Status | Owner | Scope |
|----|------|-------|--------|-------|-------|
| H-040 | `helios-python`: thin PyO3 API (`import helios`) over physics/planning; abi3-py39 wheel (maturin); 13 value-semantic pytest equivalence tests. **Completes the 11-crate roster.** | [minor] | done | claude-helios | `crates/helios-python/**` |
| H-040b | `helios-python`: expose Volume / attenuation_map / forward_project_ray / filtered_back_projection via numpy zero-copy (`numpy`/`PyArray`); needs the geometry feature + array-buffer protocol | [minor] | todo | — | `crates/helios-python/**` |
| H-041 | End-to-end workflow validation: one integration test (`helios-simulation/tests/end_to_end.rs`) where a shared μ drives both the imaging branch (Radon→FBP→registration) and the therapy branch (helical MLC delivery→divergent-fan dose→scatter→DVH/gamma), all self-consistency oracles. | [minor] | done | claude-helios | `crates/helios-simulation/tests/**` |
| H-041b | Runnable example `helios-simulation/examples/tomotherapy_workflow.rs`: CT→μ→(Radon/FBP recon) + (helical MLC delivery→divergent-fan dose→scatter)→DVH/gamma, rendering `ct/mu/recon/dose.png` (Output & visual verification — inspected: phantom, FBP recovery, central rotational dose falloff). Prints recon μ +0.1%, self-gamma 100%. | [minor] | done | claude-helios | `crates/helios-simulation/examples/**` |
| H-041c | Extend the example to consume a real DICOM series (`load_ct_series`, feature `dicom`) and expose it from Python (`helios-python`); commit a golden-image snapshot with a derived tolerance | [minor] | todo | — | `crates/helios-simulation/**`, `crates/helios-python/**` |
| H-044 | `helios-imaging::register_translation`: IGRT rigid whole-voxel translation registration (setup-error / couch-shift estimate, mean-SSD over overlap, exhaustive search). Verified: recovers a known applied shift exactly (±, zero), f32. | [minor] | done | claude-helios | `crates/helios-imaging/**` |
| H-044b | `helios-imaging::register_translation_ncc`: normalized-cross-correlation registration (robust on low-texture images — rejects zero-variance overlaps that make SSD ambiguous). Verified: recovers a known shift on the flat-background spike phantom where SSD is ambiguous, and on a textured phantom; f32. | [minor] | done | claude-helios | `crates/helios-imaging/**` |
| H-044c | `helios-imaging`: sub-voxel (interpolated-peak) translation refinement + rotation; deformable / mutual-information registration via `ritk-registration` (heavier burn build) | [major] | todo | — | `crates/helios-imaging/**` |
| H-042 | Validation report: gamma/DVH vs reference; MVCT image metrics | [minor] | todo | — | `validation_reports/**` |
| H-043 | Performance: GPU-vs-CPU scaling study — criterion benchmark of `beam_transmission_into` across sizes (`helios-gpu/benches/transmission_throughput.rs`) + quantitative report. Finding: the isolated transmission kernel is transfer-bound; GPU does not beat CPU at any tested size (RTX 5080 vs Core Ultra 9 285K). | [minor] | done | claude-helios | `crates/helios-gpu/benches/**`, `validation_reports/**` |
| H-043b | Performance: on-device fused imaging pipeline (HU→μ → forward projection → `exp(−τ)` resident on GPU, one CT upload + one sinogram download) to amortize transfer and realize GPU throughput > CPU; re-benchmark vs the H-043 baseline | [major] | todo | — | `crates/helios-gpu/**`, `crates/helios-solver/**` |
