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
- `helios-imaging` crate (H-030a): `parallel_beam_radon` + `Sinogram` — the MVCT
  forward-projection sinogram `p(θ,s) = ∫μ dl` over projection angles and signed
  detector offsets, built on the ray-march projector. Validated against the
  analytical uniform-disk sinogram `μ·2√(R²−s²)` (2% at 0.5 mm voxels),
  angle-independence, off-object zero. FBP reconstruction = H-030.
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
