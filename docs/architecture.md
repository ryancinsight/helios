# Helios Architecture

Helios (Helios-rs) is an **integrator** repository in the Atlas multi-repo stack:
a Cargo workspace that composes the Atlas foundation/infrastructure/domain crates
into a unified platform for radiation-therapy simulation/planning and radiation
imaging. It mirrors the kwavers/cfdrs model exactly — Atlas crates are consumed as
**remote git dependencies** (each carries its own `[workspace]`), and Helios's own
crates form a strictly **unidirectional, layered** dependency graph.

## Layering (unidirectional; dependencies point inward/upward only)

```
                       helios-python  (PyO3 bindings; thin, no domain logic)
                             │
        ┌────────────┬───────┴───────┬─────────────┐
   helios-planning  helios-simulation  helios-imaging   (application layers)
        │            │                 │
        └──────┬─────┴────────┬────────┘
          helios-analysis   helios-solver ── helios-gpu   (compute layers)
                     │            │              │
                helios-physics ───┘              │
                     │                           │
                helios-domain ──────────────────┘
                     │
                helios-math
                     │
                helios-core                       (foundation; depends on nothing project-local)
```

A lower layer never depends on a higher one. `helios-core` is the innermost crate;
`helios-python` is the only crate permitted to depend on `pyo3`.

## Crate responsibilities

| Crate | Responsibility | Status |
|-------|----------------|--------|
| `helios-core` | Typed errors, physical constants, validating domain newtypes, config, logging, arena hooks. | **implemented (0.0.1)** |
| `helios-math` | Numeric seam (`Scalar` = `eunomia::RealField`), leto linear-algebra substrate re-export, numerical methods. Geometry *primitives* (`Aabb`/`Ray`/mesh) are consumed from **gaia**, not defined here. | **implemented (0.0.1)** |
| `helios-domain` | Patient/imaging geometry (CT/MVCT), beam/source/sensor models, binary MLC + collimator geometry, helical delivery kinematics. Landed: `VoxelGrid` + `Volume`, including a Leto `Isometry3` oriented-grid pose; `HelicalDelivery`; binary-MLC `MlcModel`; DICOM ingest (`load_ct_slice`/`load_ct_series` → HU `Volume`, via `ritk-dicom`, feature `dicom`); and HDF5 volumetric storage (`save_volume_hdf5`/`load_volume_hdf5`, via consus, feature `storage`). HU-semantic newtypes and DICOM `ImageOrientationPatient` ingestion remain provider-sequenced; `FieldAperture` (jaw field-shaping + penumbra over a gaia `Aabb`) landed. | **partial (0.1.0)** |
| `helios-physics` | Helios-specific radiation physics that is not a shared transport law: HU→relative-density calibration and Compton cross-section/energy-transfer models. It returns Hyperion coefficient types instead of owning or re-exporting a parallel coefficient vocabulary. | **partial (0.1.0)** |
| `helios-solver` | Spatial integration and deterministic dose engines over Helios voxel grids. Landed (0.1.0): Proteus-density + Hyperion mass-to-linear HU→μ conversion, rigid-pose-aware `forward_project_ray` (∫μ dl), Hyperion-validated primary-fluence transport, transactional `deposit_ray_terma`, and collapsed-cone/scatter kernels. Hyperion owns scalar transport laws; Helios owns grid traversal, deposition, convolution, and backend differential boundaries. | **partial (0.1.0)** |
| `helios-analysis` | Dosimetric analysis (DVH, gamma), imaging quality metrics, visualization. Landed (0.0.1): cumulative + structure-masked (per-PTV/OAR) DVH, 3-D gamma index (3%/2 mm, global + local normalization with low-dose cutoff), MVCT `image_quality` metrics (RMSE / relative-L2 accuracy, ROI noise, Michelson contrast, CNR), and zero-copy evaluation of Asclepius gEUD/TCP/NTCP laws over the stored Aequitas dose sample. RT-struct rasterization + portal dosimetry pending. | **partial (0.0.1)** |
| `helios-simulation` | Time-dependent helical TomoTherapy delivery with synchronized MVCT and motion. Landed (0.0.1): Moirai-dispatched `simulate_helical_sinogram` with Hyperion transmission, `simulate_helical_delivery`, fallible per-frame beamlet dose accumulation for parallel/divergent fans, Hyperion-backed EPID exit fluence, and beam-following collapsed-cone dose. Cone-beam detector geometry and motion remain pending. | **partial (0.0.1)** |
| `helios-planning` | Inverse planning / optimization (gradient-based, multi-criteria). Landed (0.0.1): projected-gradient beam-weight optimizer (`DoseInfluence` + `optimize_beam_weights`) and the coeus-autodiff gradient backend (`objective_gradient_autodiff`, feature `autodiff`, differentially pinned to the exact hand gradient). Non-quadratic objectives on that backend: DVH-band penalty + optimizer (`DvhPenalty`/`optimize_beam_weights_dvh`) and `EudPenalty`, whose tape delegates gEUD construction to `asclepius-coeus`. TCP/NTCP objectives pending. | **partial (0.0.1)** |
| `helios-imaging` | MVCT acquisition modeling, reconstruction, and IGRT workflows. Landed (0.0.1): parallel-beam Radon forward transform (`Sinogram`), Ram-Lak FBP + SIRT iterative reconstruction, Hyperion-validated deterministic quantum noise, and rigid translation registration (`register_translation` SSD + `register_translation_ncc`). Sub-voxel/deformable registration through RITK remains pending. | **partial (0.0.1)** |
| `helios-gpu` | GPU dispatch over `hephaestus_core::ComputeDevice` + hephaestus-wgpu. `beam_transmission_into` (GPU `exp(−τ)`), `GpuAttenuationMapper` (HU→μ fused affine-clamp), and an axis-aligned `GpuProjector` landed; the latter rejects a non-identity `VoxelGrid` pose before upload until Hephaestus owns pose-bearing field geometry. Resident pipeline throughput reports live under `validation_reports/`. | **partial (0.1.0)** |
| `helios-python` | Thin PyO3 API (`import helios`): geometry-free physics/planning wrappers (Thomson/Klein–Nishina cross-sections, Compton μ/ρ, HU→density, projected-gradient beam-weight optimization). abi3-py39 wheel; GIL released around the planning solve. No domain logic. | **implemented (0.0.1)** |

Crates are created only when their layer is built (architecture_scoping growth
triggers — no speculative empty-crate scaffolding). The workspace `members` list
grows as each crate lands; `workspace.dependencies` already declares the full
Atlas set as the integration SSOT.

## Atlas dependency map

Package names below are verified against each Atlas repo's manifests. Git URLs are
the SSOT in the root `Cargo.toml` `[workspace.dependencies]`.

| Atlas component | Crates (packages) | Consumed by | Purpose in Helios |
|-----------------|-------------------|-------------|-------------------|
| **ritk** | `ritk-dicom` (**consumed**, feature `dicom`), `ritk-core`, `ritk-registration` | domain, analysis, imaging | Sole Helios-facing DICOM boundary (CT/MVCT parse, typed attributes, transfer-syntax dispatch, and rescale-calibrated decode); Helios has no direct dicom-rs dependency. RT struct/plan/dose, registration, VTK pending. |
| **gaia** | `gaia` (`Aabb`/`Ray` **consumed**) | math, domain | Geometry kernel: `Ray`/`Aabb` for projection + the `FieldAperture` collimator; binary MLC, patient surface/mesh pending. |
| **hephaestus** | `hephaestus-core`, `hephaestus-wgpu` (`-cuda`, `-metal` optional) | gpu, solver | GPU compute dispatch, WGSL pipelines, kernel caching. |
| **moirai** | `moirai`, `moirai-parallel` (`-async`, `-gpu`, `-iter`) | simulation, solver | Orchestration of time-dependent helical delivery + imaging; execution policies. |
| **coeus** | `coeus-autograd`, `coeus-tensor`, `coeus-core` (**consumed**, feature `autodiff`) | planning | Reverse-mode autodiff planning gradient (`objective_gradient_autodiff`), differentially pinned to the exact hand gradient; generalizes to non-quadratic objectives. |
| **asclepius** | `asclepius`, `asclepius-coeus` (**consumed**) | analysis, planning | Canonical gEUD, logistic TCP, Lyman NTCP, and stabilized differentiable gEUD tape construction; Helios retains DVHs and planning objectives. |
| **proteus** | `proteus` (**consumed**) | solver | Validated material density used by the HU→linear-attenuation boundary; no Helios-owned material-property newtype. |
| **hyperion** | `hyperion` (**consumed**) | physics, solver, imaging, simulation; gpu tests | Photon/optical coefficient types, NIST reference data, mass-to-linear conversion, optical depth, and Beer–Lambert transmission. Helios directly consumes these contracts and retains no compatibility facade. |
| **consus** | `consus-core`, `consus-hdf5`, `consus-io` (**consumed**, feature `storage`) | domain | Volumetric storage: `Volume` ↔ standard HDF5 archive (`save_volume_hdf5`/`load_volume_hdf5`, data plus validated rigid grid geometry). Zarr/compression pending. |
| **leto** | `leto` | math | Strided/typed array substrate. |
| **hermes** | `hermes-simd` | math | Portable SIMD for field/kernel/projection kernels. |
| **mnemosyne** | `mnemosyne-core` | core | Arena allocation and memory management for large 3D/4D datasets. |
| **themis** | `themis` | core | Optimal placement (NUMA/CPU/GPU) for large medical datasets. |
| **apollo** | `apollo` (`apollo-fft`) | solver, imaging | Spectral/transform methods for convolution kernels and reconstruction. |

## Key architectural invariants

- **Generic-first numeric seam.** All compute is parameterized through a `Scalar`
  trait (backed by `hermes`/`leto`) from first authorship; concrete numeric types
  appear only at I/O/FFI boundaries. `helios-core` constants are `f64` literals at
  their definition boundary and are converted into `T: Scalar` by callers.
- **Backend seams.** GPU/accelerator dispatch is mediated by a `ComputeBackend`
  trait (`hephaestus` wgpu/cuda/metal implementors); execution regimes (sync/async/
  parallel) by an `ExecutionPolicy` trait (`moirai`). Domain/physics code never
  imports a concrete device crate.
- **Domain purity.** Core domain entities carry no infrastructure types
  (`wgpu::Buffer`, `hdf5` handles); they stay isolation-testable. Infrastructure is
  reached only through trait seams.
- **Validating boundaries.** External input (DICOM, PyO3 args) is validated into
  typed domain newtypes at the boundary; invalid states are unrepresentable in the
  core (`EnergyMeV`, `HounsfieldUnit`, `VoxelSpacingMm`, …).

## Verification tiers

Correctness rests on the evidence hierarchy: type-level encoding (newtypes,
const-generic shapes) → property/fuzz tests → differential equivalence
(GPU vs CPU reference, epsilon-bounded per reduction order) → empirical validation
(gamma analysis, DVH, MVCT image metrics vs published TomoTherapy data). Each
constant/threshold is analytically derived and cited at its assertion site.
