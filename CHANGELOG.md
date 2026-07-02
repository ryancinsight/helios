# Changelog

All notable changes to Helios are documented here. Format follows
[Keep a Changelog](https://keepachangelog.com/); versioning is
[SemVer 2.0.0](https://semver.org/). Pre-1.0: minor bumps may break, documented
under a Breaking subsection.

## [0.0.1] — Unreleased (Sprint 1: Foundation)

### Added
- Cargo workspace skeleton (`edition = "2021"`, `resolver = "2"`) with
  `workspace.package`/`workspace.lints`/`workspace.dependencies` SSOT declaring the
  Atlas stack as remote git dependencies (ritk, gaia, hephaestus, moirai, coeus,
  consus, leto, hermes, mnemosyne, themis, apollo — package names verified).
- `rust-toolchain.toml` (stable + rustfmt/clippy), `.config/nextest.toml` (30 s
  slow / 60 s terminate test-time budget), `.gitignore`.
- `helios-core` crate:
  - `HeliosError` typed error enum (`thiserror`, `#[non_exhaustive]`) with
    `InvalidDomainValue { field, value, reason }`.
  - `constants` module: CODATA-2018 / ICRU-90 physical constants (speed of light,
    elementary charge, Avogadro, vacuum permittivity, electron mass/rest energy,
    classical electron radius, MeV↔J, water density, water mean excitation energy)
    with value-semantic derivation tests (mass–energy equivalence, `r_e` defining
    relation, exact SI constants).
  - Validating newtypes `EnergyMeV`, `HounsfieldUnit`, `VoxelSpacingMm`
    (`#[repr(transparent)]`, `TryFrom<f64>` boundary validation, `Display`).
- `helios-math` crate:
  - `Scalar` seam re-exported as `eunomia::RealField` (Atlas numeric SSOT), with
    `FloatElement`/`NumericElement`/`CastFrom`/`CastTo`.
  - leto linear-algebra substrate re-exported (`Vector3`, `Point3`, `Isometry3`,
    `Quaternion`, `UnitQuaternion`, `Translation3`, `Vector2`).
  - leto consumed with `default-features = false, features=["std"]` to avoid the
    leto→mnemosyne→themis version skew (see `gap_audit.md` G-10).
  - Geometry *primitives* (`Aabb`/`Ray`/intersection) are **owned by gaia**, not
    Helios; `helios-math` re-exports `gaia::{Aabb, Ray}` (H-003b) once gaia's leto
    migration is consumed via the H-050 wiring. Helios does not define its own.
- `helios-domain` crate:
  - `VoxelGrid<T: Scalar>`: anisotropic per-axis spacing + rigid leto `Isometry3`
    patient pose; `index_to_world`/`world_to_index`/`voxel_center` affine mapping;
    construction validates non-zero dims, non-overflowing voxel count, and
    finite/positive spacing (`HeliosError::InvalidDomainValue`).
  - `Volume<T: Scalar>`: dense scalar field backed by leto `Array3` (C-contiguous),
    `from_shape_fn`/`zeros`/`from_shape_vec`, `get`, and `sample_trilinear`/
    `sample_world`. Trilinear reproduces affine fields exactly (analytical oracle).
- `helios-math` (H-055): geometry vocabulary (leto substrate + gaia primitives)
  moved behind a default `geometry` feature; the `Scalar` numeric seam is always
  available. Lets the numeric/physics layers build independently of the geometry
  kernel (and of the concurrent geometry-stack churn, G-14).
- `helios-physics::compton` (H-011d2): Klein–Nishina total Compton cross-section
  and Thomson cross-section (first-principles, `r_e`/`m_e c²` from `helios-core`).
  Analytical oracles: Thomson matches the CODATA value, low-energy KN limit →
  Thomson from below, monotonic decrease with energy, σ(6 MeV) ≪ σ_T, f32
  differential vs f64 (near-α=0 cancellation documented as f64-conditioned).
  Plus `electrons_per_gram` and `compton_mass_attenuation` — the Compton μ/ρ
  *derived from first principles* (`σ_KN · N_A·Z/A`); validated against the NIST
  water value at 1 MeV (0.0707 cm²/g, Compton-dominated) to within 2% — a computed
  coefficient, not a fabricated table entry.
- `helios-physics::compton` energy transfer (H-011d3): `klein_nishina_differential`
  (dσ/dΩ), `compton_energy_transfer_cross_section` (σ_tr by quadrature of the
  differential), `compton_mean_energy_transfer_fraction`, and
  `compton_mass_energy_transfer` — the collision-kerma coefficient (μ_tr/ρ)_C,
  the dose-relevant Compton quantity. Self-validated (numeric total matches the
  closed-form σ_KN to 1e-4) and validated against NIST water μ_tr/ρ at 1 MeV
  (≈0.0310 cm²/g) to within 5%.
- `helios-physics` crate:
  - `LinearAttenuation<T>` (cm⁻¹) and `MassAttenuation<T>` (cm²/g) validated
    newtypes; `μ = (μ/ρ)·ρ` via `MassAttenuation::to_linear`.
  - Beer–Lambert `transmission(path_cm)` and `half_value_layer` (`None` for μ=0).
  - `relative_electron_density_from_hu` / `mass_density_from_hu` (first-order CT
    calibration: air→0, water→1, clamped below air).
  - Analytical tests: `T(HVL)=½`, `T(0)=1`, μ scaling with density, HU reference
    points, f32 genericity.
- `helios-gpu` crate (H-010): GPU compute over `hephaestus_core::ComputeDevice` +
  hephaestus-wgpu. `default_device` (wgpu adapter) and `beam_transmission_into` —
  MVCT detector transmission `exp(−τ)` computed on the GPU (`NegOp`+`ExpOp`),
  differentially validated against CPU `f32::exp` on a live adapter. Replicated
  hephaestus's mnemosyne/moirai/hermes `[patch]` set so the GPU dependency cluster
  resolves against the local checkout (fixes the leto→mnemosyne→themis skew, G-12).
- `helios-solver::forward_project_ray` (H-011c): MVCT forward-projection / dose
  ray-trace core — clips a gaia `Ray` to the `VoxelGrid` world `Aabb`, then
  midpoint ray-marches the trilinearly-sampled μ `Volume` to the optical depth
  `∫μ dl`. Axis-aligned grids (oriented-grid + exact Siddon tracked H-011d). 5
  analytical tests: homogeneous slab `τ=μ·L`, affine-field midpoint-exact,
  step-invariance, miss→`None`, f32. First consumer of the wired gaia geometry.
- `helios-solver::primary_fluence_parallel_x` (H-013a): dose-engine primary-
  transport stage — Beer–Lambert attenuated primary energy fluence
  `Ψ=Ψ₀·exp(−∫μ dl)` for a +x parallel beam via O(N) cumulative optical depth.
  Analytical oracles: homogeneous exponential depth curve, unattenuated entry,
  heterogeneous accumulation, f32. Kernel superposition (dose) tracked H-013b.
- **Fixed** `forward_project_ray` optical-depth units: `μ` is cm⁻¹ but the grid is
  mm, so path length is now converted mm→cm to yield a true dimensionless `τ`
  (previously 10× too large).
- `helios-physics`:
  - `projection` module: geometry-free ray line-integral reduction —
    `optical_depth(τ = Σ μᵢ·Lᵢ)` and `beam_transmission(exp(−τ))` over
    `(LinearAttenuation, length)` segments. 5 analytical tests (empty path,
    homogeneous = μ·L discretization oracle, additivity, multiplicative
    composition, f32). The geometry-coupled projector over this reduction landed
    in `helios-solver` (H-011c).
- `helios-domain::{save_volume_hdf5, load_volume_hdf5}` (H-046, feature `storage`):
  volumetric storage boundary via **consus** (the mandated Atlas storage component, now
  consumed — pure-Rust consus-core/consus-hdf5/consus-io, `[patch]`ed to the local
  checkout). Archives a dose/CT/MVCT `Volume` to a standard HDF5 file — a `volume`
  dataset (f64 LE, the Volume's own C-contiguous `(i,j,k)` order) plus a 6-element
  `geometry` dataset (spacing + origin) — and reconstructs the typed `Volume` on load.
  Adds `HeliosError::Storage` (distinct from `Dicom`). Verified: bitwise f64 round-trip
  on a distinct-per-voxel field with non-trivial spacing/origin, the file carries the
  standard HDF5 superblock signature (external-tool interoperability), f32 round-trips
  exactly through the f64 archive, and missing-file/garbage inputs are typed errors.
- `helios-simulation::simulate_helical_sinogram` moirai-parallel dispatch (H-021b): the
  independent per-projection forward projections are now dispatched through moirai's
  `Adaptive` execution policy (`map_collect_index_with` — sequential below its threshold,
  parallel above), consuming the mandated **moirai** orchestration component. The
  index-ordered collect makes the result identical to a sequential run regardless of
  scheduling (each projection is an independent read of `μ`; no reduction), verified by a
  determinism/order-preservation oracle at 256 projections. (The peer `mnemosyne-arena`
  breakage that blocked this last cycle has been reconciled; full workspace green again.)
- `helios-analysis::{spherical_mask, box_mask}` (H-047): geometric ROI mask predicates
  (sphere / axis-aligned box over a `VoxelGrid`) returning `Fn([usize;3]) -> bool`,
  directly usable as the mask for `Dvh::from_volume_masked` — per-structure DVH/statistics
  on analytically-defined ROIs (phantom inserts, simple targets) without a hand-written
  closure. Verified: radius/box voxel selection, masked-DVH mean over the ROI, f32.
  (Contour-defined ROIs via a ritk RT-struct rasterization remain H-032d.)
- `helios-simulation::frame_portal_fluence` (H-045): portal (EPID) exit dosimetry — the
  per-leaf transmitted primary fluence `Ψ_leaf · exp(−τ_leaf)` for a delivery frame, the
  image used to *verify* delivered fluence against the plan. Composes the forward
  projector with the per-leaf beamlet geometry, now factored into a shared
  `beamlet_ray`/`gantry_basis` (a `Beamlet` struct) used by both portal dosimetry and
  dose accumulation — one fan-geometry definition, no duplication. Verified: full
  transmission at μ=0, Beer–Lambert attenuation (`fluence·e^{−μ·chord}`), closed leaf
  reads 0, higher μ darkens, f32.
- `helios-imaging::register_translation_ncc` (H-044b): normalized-cross-correlation
  rigid translation registration, robust on low-texture images. It maximizes
  `NCC = Σ(m−m̄)(f−f̄)/(N·σ_m·σ_f)` over the overlap; because NCC measures correlation
  (invariant to intensity offset/scale), a shift that slides all structure out of the
  overlap leaves a near-constant, zero-variance region that is *rejected* rather than
  scored as a perfect match — curing the SSD false-minimum documented for
  `register_translation` (H-044). Verified: recovers a known shift on the flat-background
  spike phantom where plain SSD is ambiguous, on a textured phantom, and is generic over
  f32. Sub-voxel/rotation/deformable registration via ritk = H-044c.
- `helios-analysis::Dvh::homogeneity_index` (H-032e): the ICRU-83 dose homogeneity
  index `HI = (D₂ − D₉₈)/D₅₀`, a standard target plan-quality metric (lower = more
  homogeneous). Verified: `0` for uniform dose, `1.92` for a 0..99 ramp (D₂=98, D₉₈=2,
  D₅₀=50); `0`-guarded when `D₅₀` is zero.
- Tooling (G-17, coverage): the coverage-instrumentation *link* is unblocked —
  `RUSTFLAGS="-Clink-arg=-fuse-ld=lld"` (LLVM lld from the MSYS2 toolchain) links the
  instrumented binaries where the mingw bfd `ld` failed, and the full suite runs
  instrumented (183 tests, 356 profraw). A residual `cargo llvm-cov` region-attribution
  issue on the GNU target still leaves the coverage *percentage* unquantified (path
  forward: grcov / MSVC / Linux CI). Recorded, not fabricated.
- `helios-analysis::gamma_index_3d_local` (H-032c): local-normalization gamma index —
  the dose-difference criterion `ΔD = criterion·D_r` scales with the *local* reference
  dose (vs a single global value), the stricter, appropriate metric where relative
  agreement everywhere matters; reference points below `low_dose_cutoff` are excluded
  (gamma 0), the standard low-dose threshold that also avoids a vanishing `ΔD`. Shares
  one impl with the global `gamma_index_3d` (a `Norm` enum selects the mode). Verified:
  equals the global gamma for uniform dose, is strictly stricter in a low-dose region
  (local γ=1.0 vs global γ=0.2 on a two-level phantom), and excludes sub-cutoff points.
- `helios-solver::deposit_ray_terma_diverging` + divergent-fan inverse-square (H-020g):
  the divergent point-source fan now applies inverse-square fluence falloff — each
  per-segment terma is scaled by `(SAD/r)²` from the focal spot (`r` = focal-to-segment
  distance), 1 at isocentre, >1 nearer the source, <1 beyond. `BeamGeometry::PointSource`
  routes through it (parallel path unchanged; shared ray-march, no duplication).
  Verified: reduces to the energy-conserving no-falloff deposition as `SAD → ∞`, and
  steepens the entry/exit dose ratio (near-source enhancement) beyond pure attenuation.
  Anisotropic beam-aligned collapsed-cone kernel = H-020h.
- `helios-analysis::Dvh::from_volume_masked` (H-032b): structure-masked (per-PTV/OAR)
  cumulative DVH built only from voxels a mask predicate selects — the per-structure DVH
  clinical plan evaluation and DVH-agreement metrics operate on. `from_volume` is now the
  `include ≡ true` case (consolidated, one code path). Verified: disjoint target/OAR masks
  yield distinct means (2.0 vs 8.0 vs whole-volume 5.0), and a single-voxel mask is a
  point DVH. RT-struct ROI rasterization (via ritk) + local-normalization gamma = H-032c.
- Runnable end-to-end example (H-041b): `helios-simulation/examples/tomotherapy_workflow.rs`
  runs the full pipeline (CT→μ→Radon/FBP recon + helical MLC delivery→divergent-fan
  dose→collapsed-cone scatter→DVH/gamma) and renders inspectable PNGs (`ct/mu/recon/dose`)
  per the Output & visual verification discipline. On the synthetic phantom it reports
  recon water-ROI μ within +0.1% of 0.06 cm⁻¹ and a 3%/2 mm self-gamma of 100%; the
  rendered images were inspected and depict the expected structure (air/water/bone
  phantom, FBP recovery, central rotational-delivery dose falloff). Adds `image`
  (PNG-only) as a dev-dep; `cargo build --examples` covers it in CI.
- End-to-end workflow validation (H-041): `helios-simulation/tests/end_to_end.rs` — a
  single integration test where one shared attenuation volume `μ` (from a CT phantom via
  `attenuation_map`) drives BOTH the imaging branch (parallel-beam Radon → FBP → water-ROI
  μ recovery within 20%; rigid registration recovers a known couch shift) AND the therapy
  branch (helical MLC delivery → divergent-fan terma deposition → collapsed-cone scatter →
  dose, then DVH mean > 0 and 3%/2 mm self-gamma 100% pass). Proves the domain / physics /
  solver / imaging / analysis / simulation layers compose across their seams. Pulls
  `helios-imaging`/`helios-analysis` as test-only dev-deps (no production cycle). This is
  the integrated imaging-delivery clinical-realism workflow on synthetic/self-consistent
  data; a runnable `examples/` + Python program on real DICOM is H-041b.
- `helios-imaging::register_translation` (H-044): IGRT rigid setup correction — the
  whole-voxel translation `s` aligning a daily image (e.g. an MVCT reconstruction) to a
  planning reference by minimizing the mean squared difference over their overlap
  (exhaustive search over `±max_shift`). The couch-shift / setup-error estimate IGRT
  applies before delivery. Verified: recovers a known applied shift exactly (positive,
  negative, zero) on a textured phantom, f32. Assumes textured images; masked/NCC
  metric, sub-voxel refinement, and rotation/deformable registration via `ritk` = H-044b.
- `helios-imaging::sirt_reconstruction` (H-030c): SIRT iterative MVCT reconstruction —
  `x ← max(0, x + λ · C⁻¹ ⊙ Aᵀ(R⁻¹ ⊙ (b − A x)))` with `A` the Radon projector, `R⁻¹`
  per-ray chord normalization, `C⁻¹` per-voxel hit normalization, and a non-negativity
  (`μ ≥ 0`) projection. Iterative and robust to noise / sparse-angle data where FBP
  streaks. The back-projection geometry is extracted into a shared `back_project_rows`
  used by both FBP and SIRT (consolidation — FBP re-expressed as ramp-filter →
  back_project, net deletion of the duplicated loop). Verified: converges to its own
  forward model (interior mean within the 15% reconstruction tolerance; whole-image L2
  edge-Gibbs-dominated < 20%), error falls monotonically with iterations, zero sinogram
  → zero image, f32.
- `helios-simulation::BeamGeometry` + divergent-fan dose accumulation (H-020f): the
  `accumulate_delivered_dose` beam model is now a seam — `BeamGeometry::Parallel`
  (small-fan approximation, unchanged) or `BeamGeometry::PointSource { source_axis_mm }`,
  where each MLC leaf's beamlet runs from a focal spot through its isocentre-plane
  offset point so beamlets **diverge with depth** (the true TomoTherapy fan). Verified:
  the point-source fan reduces to parallel as `SAD → ∞` (matching total dose within
  1e-4), and a far off-axis beamlet that stays in one detector row when parallel sweeps
  ≥3 rows under divergence (on a 1 mm grid). Existing parallel oracles unchanged.
  Anisotropic beam-aligned CC kernel + inverse-square falloff = H-020g.
- `helios-domain::load_ct_series` (H-004c, feature `dicom`): multi-slice DICOM series
  → 3-D HU `Volume`. Parses/decodes each slice, validates an identical in-plane grid
  (Rows/Columns/PixelSpacing/in-plane origin within a 1 µm tolerance), sorts by
  `ImagePositionPatient` z, derives a uniform slice spacing (rejecting a missing/duplicate
  slice), and stacks along `k` (origin at the lowest-z slice). Refactors the shared
  per-slice parse+decode into `read_slice` + `scatter_slice` (used by both
  `load_ct_slice` and the series loader — no duplication). Verified by a *shuffled*
  3-slice synthetic round-trip (sorted stacking, Δz derived, HU per slice) plus
  single-path==single-slice equivalence and empty/non-uniform error paths. A real CT
  series can now drive the full pipeline.
- `helios-domain::load_ct_slice` (H-004b, feature `dicom`): the real-input DICOM
  boundary — parses a CT/MVCT slice with `ritk-dicom` (dicom-rs backend), decodes the
  pixel frame with its `RescaleSlope`/`RescaleIntercept` calibration to Hounsfield
  units, and maps Rows/Columns/PixelSpacing/SliceThickness/ImagePositionPatient into a
  typed HU `Volume` on an axis-aligned `VoxelGrid`. **First consumption of the mandatory
  `ritk` Atlas component.** Feature-gated so the dicom-rs parser stays out of the core
  build (a complete impl, not a stub). Adds `HeliosError::Dicom`. Verified by a
  deterministic synthetic-DICOM round-trip through the real parser (2×2 slice, raw
  [10,20,30,40] · slope 2 − 10 → HU [10,30,50,70], spacing/origin exact) and a
  missing-file error-path test. Multi-slice series stacking = H-004c.
- `helios-imaging::add_quantum_noise` (H-033b): deterministic MVCT quantum
  (photon-counting) noise model — `N = N₀·exp(−τ)`, Poisson draw (Gaussian
  approximation `N + √N·z`, exact for large counts), `τ' = −ln(N'/N₀)` — via a
  committed SplitMix64 PRNG (no external dep, reproducible from a seed). Adds
  `Sinogram::from_readings` (validated constructor) and `Sinogram::map_readings`
  (geometry-preserving combinator). Validated against analytical photon statistics:
  `Var(τ') ≈ exp(τ)/N₀` (10% ensemble tol), noise grows with attenuation, high-flux →
  clean line integrals, seed determinism, f32. An end-to-end `helios-imaging` test
  injects noise into the disk sinogram, reconstructs, and confirms interior-ROI noise
  exceeds the noiseless ripple and falls with photon flux — closing the MVCT
  *noise/CNR* sub-gate on synthetic data (the metrics from H-033 now run on genuinely
  noisy reconstructions).
- `helios-gpu/benches/transmission_throughput.rs` (H-043): GPU-vs-CPU scaling study for
  the Beer–Lambert transmission kernel (criterion, elements/s across 1 k–4 M). Delivers
  the performance-gate measurement instrument + a quantitative report
  (`validation_reports/2026-07-01-gpu-transmission-throughput.md`). Honest finding: the
  isolated `exp(−τ)` kernel is transfer-bound — on an RTX 5080 it reaches only
  ~0.5–0.72× a single-threaded CPU loop because each call round-trips the buffer over
  PCIe for ~1 flop/element (a correct roofline result). GPU throughput needs the
  on-device fused pipeline filed as H-043b; "competitive with VoLO" is not claimed (no
  external reference).
- `helios-analysis::image_quality` (H-033): quantitative MVCT image-quality metrics —
  reconstruction accuracy (`volume_rmse`, `volume_relative_l2_error` vs a ground-truth
  attenuation volume), noise (`roi_statistics` — mean + population std over a uniform
  ROI), contrast (`michelson_contrast`), and detectability (`contrast_to_noise_ratio`).
  Oracles: uniform-ROI zero noise, hand-computed mean/std, Michelson `(3,1)=0.5`, CNR
  `|10−4|/2=3`, RMSE identity + constant-offset, relative-L2 closed form (`0.25`),
  dimension-mismatch / zero-norm errors, f32. An end-to-end test in `helios-imaging`
  reconstructs the disk phantom (Radon→FBP) and quantifies interior-ROI accuracy
  (mean within 15 % of μ₀), background suppression, disk/air contrast (>0.85), and CNR
  (>1) with these metrics — the MVCT-reconstruction-accuracy/contrast gate on synthetic
  data. Stochastic quantum-noise injection (for end-to-end noise/CNR) = H-033b.
- `helios-solver::scatter_superposition` + `symmetric_deposition_kernel` (H-020e):
  stage 2 of the collapsed-cone / convolution dose model — spreads the delivered
  terma (stage 1) into dose. Separable 3-D convolution (`K = kₓ ⊗ k_y ⊗ k_z`, three
  `O(N·taps)` axis passes) with centred, `Σ=1`-normalized per-axis kernels; produces
  lateral penumbra (a beamlet's energy reaches off-line voxels) and depth build-up
  that the primary-only terma lacks. Oracles: `[1]`-kernel identity (differential vs
  the primary reference), interior point-source energy conservation, symmetric
  spread, off-axis penumbra, fluence linearity, kernel normalization/peaking, f32,
  and an end-to-end `accumulate_delivered_dose → scatter_superposition` composition
  test (zero-terma off-line voxel gains scattered dose). Separable-isotropic
  approximation; anisotropic forward-peaked CC kernel + divergent fan = H-020f.
- `helios-solver::deposit_ray_terma` + `helios-simulation::accumulate_delivered_dose`
  (H-020d): the delivery→dose loop. `deposit_ray_terma` ray-marches a gaia `Ray`
  through the μ volume depositing the primary-beam energy lost in each path segment,
  `w·(e^{−τ_before} − e^{−τ_after})`, into the nearest voxel; the per-segment losses
  telescope, so the returned total is exactly `w·(1 − e^{−τ})` (step-independent
  conservation oracle) and equals the summed voxel dose. `accumulate_delivered_dose`
  builds per-leaf beamlets from each `DeliveryFrame` (gantry angle → axial-plane
  direction, couch → z-slice, leaf index → lateral offset, effective fluence →
  weight) and sums their terma into a delivered-dose `Volume` — the input the DVH /
  gamma gates consume. Oracles: single central beamlet vs analytic `w·(1−e^{−μ·L})`,
  linearity in fluence, frame superposition, three-leaf offset fan, zero-fluence, f32.
  Adds `Volume::add_at` (bounds-checked scatter accumulation) and `Volume::sum`.
  Beamlets are parallel (small-fan approximation); divergent fan + lateral scatter =
  H-020e.
- `helios-python` crate (H-040): thin PyO3 binding surface (`import helios`) — the
  11th and final crate, completing the workspace roster. Geometry-free `f64`
  wrappers over the physics/planning cores: `thomson_cross_section`,
  `klein_nishina_cross_section`, `compton_mass_attenuation`, `mass_density_from_hu`,
  `optimize_beam_weights` (GIL released via `Python::allow_threads` around the
  iterative solve). Untrusted-input hardening at the boundary: non-finite/non-positive
  energies and shape mismatches map to Python `ValueError`. abi3-py39 cdylib
  (`maturin`); no domain logic and no other Helios crate depends on `pyo3`. Verified
  by 13 value-semantic `pytest` equivalence tests (Thomson exact, Klein–Nishina
  Thomson-limit + monotonicity, water μ/ρ vs NIST 0.0707 cm²/g, HU→density
  calibration, identity/non-negativity planning oracles, error paths) against the
  `maturin develop` module.
- `helios-planning` crate (H-031): inverse treatment planning by projected gradient
  descent — `DoseInfluence` (linear dose model `A`, `apply`/`transpose_apply`) and
  `optimize_beam_weights` minimizing `½‖A x − d‖²` over `x ≥ 0`. Convex-convergence
  oracles: identity problem → prescription, negative target → 0, monotone objective
  decrease, diagonal least-squares solution, f32. (coeus-autodiff backend for
  non-quadratic objectives = H-031b.)
- `helios-imaging::filtered_back_projection` (H-030): MVCT reconstruction by
  Ram-Lak filtered back-projection (ramp filter + linear-interpolated back-
  projection, computed in cm so it recovers μ directly). Forward→reconstruct
  round-trip on a disk phantom recovers the interior μ (centre within 15%) with
  near-zero background — the MVCT-reconstruction-accuracy capability.
- `helios-imaging` crate (H-030a): `parallel_beam_radon` + `Sinogram` — the MVCT
  forward-projection sinogram `p(θ,s) = ∫μ dl` over projection angles and signed
  detector offsets, built on the ray-march projector. Validated against the
  analytical uniform-disk sinogram `μ·2√(R²−s²)` (2% at 0.5 mm voxels),
  angle-independence, off-object zero. FBP reconstruction = H-030.
- `helios-simulation` integrated delivery (H-020c): `simulate_helical_delivery`
  ties `HelicalDelivery` kinematics to the binary-MLC `LeafOpenTimeSinogram`/
  `MlcModel` → a time-ordered `DeliveryFrame` sequence (gantry angle + couch +
  effective per-leaf fluence with leakage/T&G). `total_delivered_fluence`
  integrates the sequence. The integrated imaging-delivery-workflow layer. Oracles:
  frame count/kinematics track the sinogram, per-frame fluence matches `MlcModel`,
  all-closed → leakage-only total, all-open → full total, f32.
- `helios-simulation` crate (H-021): `simulate_helical_sinogram` — time-dependent
  helical MVCT acquisition integrating `HelicalDelivery` (gantry rotation + couch
  translation, a helix) with the forward projector: each projection rotates the
  central beam in the axial plane at the couch `z` and forward-projects through the
  μ volume → optical depth + transmission. Analytical oracles: projection count,
  axial central-ray τ = μ·chord, uniform-cube rotational symmetry (0°=90°), couch
  monotonicity, empty→full transmission, f32. CPU reference (moirai parallel
  dispatch + fan/cone detector = H-021b).
- Geometry-stack migration (H-003c): adapted Helios to the new `leto::geometry`
  API after the upstream leto rewrite settled — `helios-math` re-exports
  `Point2/Point3/Vector3/UnitVector3` (+ gaia `Aabb`/`Ray`); `VoxelGrid` simplified
  to **axis-aligned** (origin + spacing), dropping the now-reduced leto `Isometry3`
  pose (oriented grids tracked H-003d); projector pose-rotation check removed.
  Restored the full-workspace build (97 tests, all crates incl. live GPU).
- `helios-solver` dose kernel superposition (H-013b): `dose_convolution_x`
  (dose = TERMA ⊛ forward kernel) + `exponential_deposition_kernel`. Analytical
  oracles: delta-kernel identity (dose = TERMA), normalized-kernel interior energy
  conservation, physical depth-dose build-up, empty-kernel. Now verified (was
  blocked by G-14).
- Integration wiring (H-050): `[patch]` redirecting `leto`/`eunomia`/`gaia` git
  sources to the local synchronized Atlas checkout, so Helios builds against one
  consistent source and consumes gaia's **migrated leto/eunomia geometry**.
- `helios-math` (H-003b): re-exports `gaia::{Aabb, Ray}` as the Helios geometry
  vocabulary (upstream ownership); bridge test verifies a gaia `Ray` intersects a
  gaia `Aabb` through Helios. Unblocks the voxel-DDA projector (H-011c).
- `helios-domain` binary MLC (H-020b): `LeafOpenTimeSinogram` (validated per-
  projection/leaf open-time fractions) + `MlcModel` — effective transmitted
  fluence = leakage-adjusted transmission (`open + (1−open)·τ`) minus a
  tongue-and-groove edge loss where a neighbour is more closed, clamped to `[0,1]`.
  The binary-MLC leakage/tongue-and-groove clinical-realism capability. Analytical
  oracles: closed→leakage, open→1, uniform-row no T&G loss, isolated-open-leaf
  underdose, neighbour-aware sinogram application, bounds, f32.
- `helios-domain`: `HelicalDelivery<T>` — helical TomoTherapy delivery kinematics
  (gantry rotation + couch translation + pitch synchronization). Projection/time →
  gantry angle (unwrapped + wrapped) and couch position; pitch relation
  (`couch_travel_per_rotation = pitch·field_width`), couch velocity. 7 analytical
  tests (one-rotation advances angle by 2π and couch by the pitch travel,
  projection↔time agreement, half-rotation = π, monotonic couch, f32).
- `helios-analysis` crate:
  - `Dvh`: cumulative dose-volume histogram from a dose `Volume` — `min`/`max`/
    `mean`, `volume_fraction_at_dose` (Vx), `dose_at_volume_fraction` (Dx,
    nearest-rank). Verified on uniform (step) and ramp (known quantiles) fields.
  - `gamma_index_3d`: Low's 3D gamma (dose-difference / distance-to-agreement,
    global normalization) with grid + criterion validation, and `gamma_pass_rate`.
    Analytical oracles: identical→γ=0/100% pass, γ scales with dose-ratio, 2×
    criterion→fail, f32 genericity. The 3%/2 mm quality-gate machinery.
- `helios-solver` crate:
  - `attenuation_map`: deterministic per-voxel HU→μ engine mapping a CT `Volume`
    to a linear-attenuation `Volume` (cm⁻¹) via `ρ = mass_density_from_hu(HU)` and
    `μ = (μ/ρ)·ρ` (Compton-dominated MV approximation). CPU reference — the
    differential oracle for the future GPU kernel (H-010).
  - Tests: uniform-water constant μ, air→0/bone-scaling, per-voxel closed-form
    differential match over a heterogeneous field, grid preservation, f32.
- Foundation documentation: `README.md`, `ARCHITECTURE.md` (layering + Atlas
  dependency map), and PM artifacts `backlog.md`, `CHECKLIST.md`, `gap_audit.md`,
  `SPRINT_1.md`, `SPRINT_2.md`.

### Verification
- `cargo build`, `cargo clippy --all-targets --all-features -- -D warnings`,
  `cargo fmt --check`, `cargo nextest run` (13 tests) — all green.

[0.0.1]: https://github.com/ryancinsight/helios/releases/tag/v0.0.1
