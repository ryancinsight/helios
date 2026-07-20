# Changelog

All notable changes to Helios are documented here. Format follows
[Keep a Changelog](https://keepachangelog.com/); versioning is
[SemVer 2.0.0](https://semver.org/). Pre-1.0: minor bumps may break, documented
under a Breaking subsection.

## [0.1.0] — Unreleased

### Breaking

- H-069: `MassAttenuation::to_linear` now accepts Proteus `MassDensity<T>`
  instead of an unvalidated raw scalar. Proteus owns the shared material-density
  validity boundary; Helios retains mass-attenuation and CT-calibration laws.
- H-068: `EnergyMeV` and `VoxelSpacingMm` retain their Helios validation
  contracts while storing Aequitas `Energy<f64>` and `Length<f64>`. Compile-time
  layout assertions preserve their one-word transparent representation.
  Their `get` methods are no longer `const` because unit conversion is
  trait-driven. `HounsfieldUnit` remains Helios-owned because HU is a
  calibrated imaging scale rather than an SI linear unit.

### Changed

- H-074: aligned CI path-dependency checkout with Atlas commit `afd5e16`,
  whose gitlinks contain the merged Aequitas, Proteus, and Hephaestus provider
  graph. The checkout action implementation remains immutably pinned to
  `9bfb722`.
- H-073: upgraded the thin Python binding boundary to PyO3 0.29, closing
  RUSTSEC-2025-0020 and RUSTSEC-2026-0177. Compute now uses `Python::detach`
  under PyO3's corrected thread-safety contract.
- H-072: aligned the direct Aequitas and Proteus revisions with the merged
  temperature-response provider graph. Helios and Hephaestus now resolve one
  Aequitas source identity, while the radiation material contract remains
  unchanged.
- H-071: replaced the copied same-run Python benchmark check with the
  Atlas-owned phase-replicated Criterion gate pinned to merge `9bfb722`.
  Pull-request CI now compares baseline and candidate production revisions on
  one runner through ABBA and BAAB blocks, while native tests use the committed
  Nextest budget, doctests run separately, and RustSec/license/source policy is
  enforced with pinned verification tools.
- CI now materializes root and Coeus-transitive external path dependencies
  through the exact Atlas-owned provider-checkout action instead of a
  consumer-owned provider list.
- H-070: routed deterministic MVCT normal sampling through Tyche's
  counter-addressed `StandardNormal` provider, removed Helios's duplicate
  mutable SplitMix64/Box-Muller implementation, and pinned the seed-to-reading
  mapping with an exact regression vector.
- H-067: refreshed the reproducibility lock to `apollo-fft` 0.25.0, Eunomia
  0.4.0, Leto 0.38.2, and Hephaestus 0.17.0 without changing Helios source or
  manifests. The lock graph no longer contains `num-complex`.
- H-066: removed Helios's unused direct `num-traits` workspace dependency.
  Provider-owned transitive requirements remain resolved through the locked
  Eunomia, Gaia, Half, and WGPU graph.
- Updated the reproducibility pin to `apollo-fft` 0.23.0 and
  `apollo-leto-interop` 0.17.0 from the merged Apollo provider release.
- H-005: reconciled the stale binary-MLC/collimator roadmap item with the
  delivered H-020b `LeafOpenTimeSinogram`/`MlcModel` and H-020k
  gaia-backed `FieldAperture`/delivery-collimation contracts; marked duplicate
  H-022 done; restored the workspace formatter gate for the affected math and
  solver sources.
- Helios 0.1.0 declares Rust 1.95, matching its merged Mnemosyne 0.5 and
  Leto 0.38 provider graph. The lockfile is the reproducibility pin; no
  revision-qualified first-party source quarantine is introduced.
- H-065: refreshed the locked Moirai provider graph to 0.3.0 after the
  provider-owned NUMA iterator retirement; all Helios examples compile against
  the resolved release graph.
- H-064: moved Helios DICOM parsing, typed attribute access, transfer-syntax
  selection, pixel decoding, and synthetic-input verification behind the
  `ritk-dicom` public API. Removed Helios's direct `dicom` dependency; the
  dicom-rs implementation is now an internal RITK provider detail.
- H-063: aligned Helios's direct DICOM dependency with the `ritk-dicom`
  provider at version 0.10.0 and regenerated the lockfile. The workspace now
  resolves one `dicom-core` type across the DICOM boundary; the four former
  `helios-domain` E0308 errors are closed.
- H-062: `helios-analysis::Dvh::volume_fraction_at_dose` now uses a
  zero-allocation binary lower bound over its sorted sample instead of scanning
  every voxel for each threshold query. NaN-containing samples retain the
  previous filter semantics through an explicit fallback. The fixed Criterion
  comparison is recorded in
  `validation_reports/2026-07-15-dvh-query-optimization.md`.
- H-061: validated all three runnable examples against the synchronized Atlas
  provider graph and removed Helios's direct dicom-rs `ndarray` feature
  activation; `ritk-dicom` remains the pixel-decoding owner.
- Normalized the five remaining runnable-example source files with the pinned
  Rust formatter; `cargo fmt --all --check` is clean across the workspace.

### Added
- H-011b: embedded the selected NIST X-ray `μ/ρ` tables for dry air, liquid
  water, and cortical bone over 10 keV–20 MeV. `NistMaterial` returns a
  validated generic `MassAttenuation` with allocation-free, native-precision
  log-linear interpolation between edge-free table knots; it explicitly does
  not claim to reproduce XCOM's cubic-spline output.
- H-003d: `VoxelGrid` now owns a Leto `Isometry3` pose and exposes
  `VoxelGrid::oriented`; its index/world transforms apply anisotropic spacing
  in local index space and the rigid pose exactly once. Axis-aligned grids use
  the identity rotation. The CPU projector and terma deposition clip in that
  local index frame, preserving world-space millimetre path length for an
  oriented volume. HDF5 geometry now carries the three rotation columns in its
  15-value dataset and validates them on load, so storage cannot discard pose.
  The present Hephaestus field kernel has no pose metadata and rejects a
  non-identity grid before upload; it does not silently project an oriented
  volume as axis-aligned. Generic f32/f64 index/world mapping and value-semantic
  projector/storage tests cover the rigid-pose contract. DICOM
  `ImageOrientationPatient` ingestion remains blocked on RITK's named provider
  tag (H-004d), so no raw DICOM tag is duplicated in Helios.

### Breaking

- The pre-0.1 HDF5 `geometry` dataset had six values (spacing plus origin). It
  now has fifteen, appending three world-space rotation columns. Regenerate
  archived volumes with 0.1.0 so their rigid grid pose is explicit and
  validated on load.
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
- `helios-gpu::GpuAttenuationMapper` (H-010b): GPU HU→μ map as a Helios-authored
  WGSL kernel over Hephaestus `KernelInterface` + `KernelSource<Wgsl>`, computing
  `μ = max(fma(scale, HU, offset), 0)` in one dispatch. Differential oracles cover
  the closed-form clamp and `helios_solver::attenuation_map`; the upstream seam
  required no type-specific affine-clamp helper.
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
- End-to-end demonstration + coverage of the beam-following poly-energetic dose
  (H-041c): `examples/tomotherapy_workflow` now drives the therapy stage through
  `accumulate_delivered_dose_anisotropic` with a two-component (soft/hard)
  `CollapsedCone::poly_forward_peaked`, and its dose render was **inspected** (a
  centrally-concentrated distribution with smooth penumbra — the rotation-averaged
  forward-peaked helical delivery; recon μ +0.1%, self-gamma 100%). A new integration
  test `beam_following_poly_energetic_dose_end_to_end` exercises the full CT→delivery→
  anisotropic-poly-dose→DVH→gamma pipeline with analytical oracles: non-negativity,
  energy conservation (scattered dose ≤ deposited terma, >85% retained), positive DVH
  mean, and 3%/2 mm self-gamma 100%.
- Collimator field aperture + delivery collimation (H-020k). **helios-domain**:
  `FieldAperture` — the jaw-shaped open field as a gaia `Aabb` with a linear geometric
  edge penumbra. `transmission(point)` is 1 deep inside, 0 deep outside, and 0.5 on the
  geometric edge, ramping across a `±penumbra_mm` band (standard-box SDF; `contains`
  delegates to gaia `Aabb::contains_point`). **helios-simulation**: `collimate_frames`
  scales each leaf's fluence by the aperture transmission at that leaf's collimator
  coordinate `(lateral_offset, couch_mm, 0)` — jaw field-shaping + penumbra on top of the
  MLC modulation. Verified: centre→1 / far→0 / edge→0.5 / penumbra-band ramp / monotone /
  typed errors / f32; and end-to-end a narrow aperture shapes the field (edge leaves →
  50 %, outside leaves → 0, machine state preserved) while collimation never increases
  fluence. Deepens the mandated gaia-geometry consumption. Poly-energetic collapsed-cone
  kernels (H-020j). **helios-solver**:
  `poly_forward_peaked_kernel` + `SpectralComponent` — the deposition kernel is the
  energy-fluence-weighted convex combination of the monoenergetic `forward_peaked_kernel`s
  of a spectrum's components (each already Σ=1, so the sum renormalizes to Σ=1; weights
  need not be pre-normalized). Models **beam hardening**: harder (longer-range) components
  reach farther downstream. **helios-simulation**: `CollapsedCone::poly_forward_peaked`
  carries it into delivered dose (the two `CollapsedCone` constructors share a private
  `from_beam_kernel` assembly — SSOT). Verified: a single positive-weight component reduces
  to the monoenergetic kernel exactly (1e-15) and to the monoenergetic delivered dose
  (1e-13); weight scale-invariance + Σ=1; a harder-weighted spectrum raises the
  downstream/upstream kernel-mass ratio and shifts the delivered-dose beam-axis centroid
  downstream; empty spectrum → identity; f32. Per-leaf gaia collimation filed as H-020k.
- Beam-following anisotropic delivered dose (H-020i). **helios-solver**:
  `directional_convolve` — an oriented 1-D convolution along an **arbitrary** unit
  direction via trilinear resampling of the field (gather = `p − offset·d`, boundary
  samples → 0) — and `oriented_forward_scatter`, which composes it into a beam-frame
  collapsed cone (forward-peaked along the beam, symmetric across a Gram–Schmidt lateral
  basis). Reduces to `anisotropic_scatter_superposition` (H-020h) when the beam is a grid
  axis and the sample step is that axis's pitch (samples land on nodes — the differential
  oracle). **helios-simulation**: `CollapsedCone` (forward-peaked kernel config) and
  `accumulate_delivered_dose_anisotropic`, which scatters **each frame's terma along that
  frame's own gantry direction** before summing, so the forward-peaked physics follows
  the rotating helical beam (a single scatter on the pooled terma has no coherent beam
  axis). The per-frame beamlet deposition is extracted to a shared `deposit_frame_terma`
  (SSOT; `accumulate_delivered_dose` unchanged). Verified: axis-aligned reduction to the
  separable pipeline (1e-10), oblique downstream>upstream with lateral symmetry, interior
  energy conservation, forward-peaking shifts the delivered-dose beam-axis centroid
  downstream, and linearity + f32. Poly-energetic spectra + gaia per-leaf collimation are
  filed as H-020j.
- `helios-solver::{forward_peaked_kernel, anisotropic_scatter_superposition}` (H-020h):
  the beam-aligned **anisotropic** collapsed-cone scatter stage. `forward_peaked_kernel`
  builds a Σ=1 deposition kernel with *different* upstream/downstream exponential ranges
  (secondary electrons travel forward; backscatter is short-ranged), returning the
  zero-offset index; `anisotropic_scatter_superposition` applies it along the beam axis
  with the symmetric kernel laterally, on a generalized `convolve_axis_at` (explicit
  centre; the centred path delegates — one loop, no duplication). **Defect caught by the
  new oracle:** the shared gather was correlation (`src[pos+off]`), not convolution —
  invisible to every symmetric kernel, but it inverts anisotropy; fixed to
  `src[pos−off]` with symmetric results bitwise-unchanged. Verified: equal ranges reduce
  exactly to `scatter_superposition`; a point source deposits strictly more energy
  downstream than upstream while lateral symmetry holds; interior energy conserved;
  f32 + beam-axis selectability. Rotated per-gantry cone axes = H-020i.
- End-to-end per-structure plan evaluation (H-033c): a new integration test
  `per_structure_plan_evaluation_over_delivered_dose` runs the full surface — helical
  delivery → beam-following collapsed-cone dose → **masked** DVH (central target vs
  off-axis OAR spheres) → gEUD → TCP / NTCP. Oracles are clinical-plausibility +
  probability well-formedness: the central target is hotter than the off-axis OAR
  (rotational convergence), PTV gEUD > OAR gEUD, PTV TCP > 0.5 (TCD50 below the target
  gEUD), OAR NTCP < 0.5 (TD50 above the OAR gEUD) — proving the masks + DVH + radiobiology
  metrics compose over real delivered dose.
- Per-structure outcome methods on the DVH (H-033b): `Dvh::dose_sample` (zero-copy view
  of the structure's ascending-sorted doses) plus `Dvh::{generalized_eud, tcp_logistic,
  ntcp_lkb}`, which evaluate the radiobiology models on the sample the histogram already
  holds — **no dose-volume re-scan** — at the natural receiver (method form over free
  function). So a masked (PTV/OAR) DVH now reports its own gEUD and TCP/NTCP directly.
  Verified: `Dvh::generalized_eud` matches the free function on the reused sample across
  a; a uniform-dose structure's TCP/NTCP reduce to the pointwise models (0.5 at
  TCD50/TD50) and match the free functions at the structure gEUD; a masked hot-half
  structure yields higher gEUD/NTCP than the cold half.
- Radiobiology plan-evaluation metrics (H-033) in a new `helios-analysis::radiobiology`
  module: `generalized_eud` (Niemierko gEUD, **promoted here from helios-planning** — a
  dose metric belongs in analysis, not gated behind planning's `autodiff` feature; now
  generic over `Scalar` and always available), plus the outcome models built on it —
  `tcp_logistic` (Niemierko logistic TCP `1/(1+(TCD50/gEUD)^{4γ50})`) and `ntcp_lkb`
  (Lyman–Kutcher–Burman NTCP `Φ((gEUD−TD50)/(m·TD50))` via eunomia's `erfc`). Verified:
  gEUD power-mean bounds/monotonicity/uniform-invariance; TCP bounded [0,1], 0.5 at TCD50,
  monotone, slope-sharpening; NTCP matches the normal CDF at the published Φ(±1)=0.8413/
  0.1587 and Φ(0)=0.5, bounded/monotone; f32. helios-planning's EUD objective is unchanged
  and its tests now use an independent inline gEUD oracle (a differential test must not
  check code against itself), so planning keeps its lean core+math dep set.
- `helios-planning::{EudPenalty, EudKind, eud_objective_gradient_autodiff}`
  (H-031d, feature `autodiff`): a **generalized-EUD (Niemierko)** biological planning
  objective. `generalized_eud` computes `gEUD = (mean(D^a))^(1/a)` (a=1→mean, a→+∞→max
  for serial/OAR control, a→−∞→min for parallel/target coverage). The gEUD of `A·x` is
  built from differentiable `matmul`/`pow`/`mean` ops on the coeus tape, so a one-sided
  quadratic gEUD penalty (`EudPenalty` — OAR upper limit / target lower limit) has its
  gradient w.r.t. beam weights by reverse-mode AD — a gradient with **no closed form**,
  the capability the mandated coeus component exists for. Verified: gEUD power-mean bounds
  + monotonicity in a + uniform-dose invariance; the tape gEUD value matches the analytic
  `generalized_eud`; the objective gradient matches a central finite difference (the
  differential oracle over the whole gEUD-plus-penalty tape); zero gradient when the hinge
  is inactive; typed errors for a=0 / shape mismatch.
- `helios-planning::{DvhPenalty, dvh_objective_gradient_autodiff, optimize_beam_weights_dvh}`
  (H-031c, feature `autodiff`): the **non-quadratic** clinical planning objective —
  one-sided DVH-style penalties `L(x) = w_u·Σ relu(floor − A·x)² + w_o·Σ relu(A·x −
  ceiling)²` (underdose below the prescription floor and overdose above the OAR ceiling
  penalized; the band in between free) — with its gradient from the coeus tape (`relu`
  kinks handled by reverse-mode AD; weights folded as `[1]`-shaped constant `Var`s, one
  backward pass) and a projected-gradient optimizer on top. Verified: the tape gradient
  matches the hand sub-gradient `−2w_u·Aᵀrelu(floor−Ax) + 2w_o·Aᵀrelu(Ax−ceiling)`
  within 1e-12; objective value cross-checked; zero value/gradient strictly inside the
  band; and on a target/OAR toy problem the optimizer selects the OAR-sparing beamlet
  (target dose ≥ floor, OAR dose ≤ ceiling, weights ≥ 0). This is the capability the
  mandated coeus component exists for — objectives with no closed-form gradient.
- `helios-planning::objective_gradient_autodiff` (H-031b **resolved**, feature
  `autodiff`): the planning gradient `∇ₓ ½‖A·x − d‖²` computed by coeus reverse-mode
  autodiff (`Var`/`matmul`/`sub`/`mul`/`sum` tape over the MoiraiBackend) — the mandated
  **coeus** component, the last unconsumed Atlas component, is now in use. Verified: the
  tape gradient matches the exact hand gradient `Aᵀ(A·x − d)` within 1e-12 (differential
  test), is zero at the least-squares optimum, and shape mismatches are typed errors.
  Landing required two cross-repo unblocks in prior cycles: apollo's vestigial
  `leto/ndarray-compat` feature leak (apollo f1ddf7a) and the peer's moirai-core
  refactor completing (moirai 2451715). Adds `DoseInfluence::rows()`. This is the
  gradient backend that generalizes to non-quadratic (DVH/biological) objectives.
- Resident GPU forward projection (H-043b **resolved**): `helios_gpu::GpuProjector`
  uploads the attenuation volume once and forward-projects whole ray batches per
  dispatch through hephaestus's new `ray_line_integrals` volume ray-integral kernel
  (**upstreamed**, commits 792ccc3/9354260: slab-clip to the node AABB →
  `n = ceil(L/step)` midpoint trilinear samples, one thread per ray, 4 live-GPU
  analytical oracles). Measured on a 128³ volume (report
  `validation_reports/2026-07-02-gpu-projection-throughput.md`): **171× vs the
  single-thread CPU projector at a 90×128 sinogram (75.4 ms → 0.441 ms) and 371× at
  360×256 (589.6 ms → 1.591 ms)** — residency converts the GPU from the
  transfer-bound elementwise loss into a two-order-of-magnitude win on the pipeline's
  dominant workload. Per-ray differential agreement with `forward_project_ray` within
  a derived 1e-3 f32 bound (live-adapter test); misses are exactly 0 on both paths.
  Closes G-18. ("VoLO-competitive" remains unclaimed — no reference engine here.)
- Fused GPU transmission kernel (H-043b step 1): `beam_transmission_into` now dispatches
  hephaestus's fused `ExpNegOp` (`exp(−x)`) — **upstreamed to hephaestus-wgpu for this
  path** (commit 669a9b3, with a live-GPU contract test) — one dispatch and no
  intermediate device buffer, replacing the `NegOp → ExpOp` chain. Measured: +30 % GPU
  throughput at 4M elements (373→485 Melem/s); still PCIe-transfer-bound at 0.66–0.73×
  CPU (honest addendum in the H-043 validation report). Remaining H-043b scope: the
  full on-device μ→projection→transmission pipeline.
- Performance/consolidation pass (H-048): `Volume::as_slice` — a public zero-copy view
  with a documented C-contiguous `(i,j,k)` layout contract (the private alias deleted;
  one accessor). The dose engine's hottest kernel, `scatter::convolve_axis`, now iterates
  that slice with a precomputed axis stride instead of per-voxel bounds-checked
  `get().expect()` inside a `from_shape_fn` closure: **8.3× faster at 32³
  (4.31→0.52 ms) and 7.4× at 64³ (37.41→5.02 ms), bitwise-identical results** (all 35
  solver oracles unchanged; criterion baseline + roofline analysis in
  `validation_reports/2026-07-02-scatter-superposition-optimization.md`; new
  `scatter_superposition` bench committed). `save_volume_hdf5` serializes the field via
  the same slice view. `MM_PER_CM` — previously duplicated in **five** modules — is now
  a single SSOT constant in `helios-core::constants` (all sites import it).
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
  boundary — parses a CT/MVCT slice with the `ritk-dicom` provider, decodes the
  pixel frame with its `RescaleSlope`/`RescaleIntercept` calibration to Hounsfield
  units, and maps Rows/Columns/PixelSpacing/SliceThickness/ImagePositionPatient into a
  typed HU `Volume` on an axis-aligned `VoxelGrid`. **First consumption of the mandatory
  `ritk` Atlas component.** Feature-gated so the RITK DICOM provider stays out of the
  core build (a complete impl, not a stub). Adds `HeliosError::Dicom`. Verified by a
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
