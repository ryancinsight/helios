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
    Helios; consuming them is tracked as H-003b (blocked on gaia's leto-geometry
    migration, `gap_audit.md` G-11). Helios does not define its own.
- `helios-domain` crate:
  - `VoxelGrid<T: Scalar>`: anisotropic per-axis spacing + rigid leto `Isometry3`
    patient pose; `index_to_world`/`world_to_index`/`voxel_center` affine mapping;
    construction validates non-zero dims, non-overflowing voxel count, and
    finite/positive spacing (`HeliosError::InvalidDomainValue`).
  - `Volume<T: Scalar>`: dense scalar field backed by leto `Array3` (C-contiguous),
    `from_shape_fn`/`zeros`/`from_shape_vec`, `get`, and `sample_trilinear`/
    `sample_world`. Trilinear reproduces affine fields exactly (analytical oracle).
- `helios-physics` crate:
  - `LinearAttenuation<T>` (cm⁻¹) and `MassAttenuation<T>` (cm²/g) validated
    newtypes; `μ = (μ/ρ)·ρ` via `MassAttenuation::to_linear`.
  - Beer–Lambert `transmission(path_cm)` and `half_value_layer` (`None` for μ=0).
  - `relative_electron_density_from_hu` / `mass_density_from_hu` (first-order CT
    calibration: air→0, water→1, clamped below air).
  - Analytical tests: `T(HVL)=½`, `T(0)=1`, μ scaling with density, HU reference
    points, f32 genericity.
  - `projection` module: geometry-free ray line-integral reduction —
    `optical_depth(τ = Σ μᵢ·Lᵢ)` and `beam_transmission(exp(−τ))` over
    `(LinearAttenuation, length)` segments. 5 analytical tests (empty path,
    homogeneous = μ·L discretization oracle, additivity, multiplicative
    composition, f32). The voxel-DDA *segment generation* half awaits gaia
    geometry (G-11).
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
