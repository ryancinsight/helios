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
| `helios-domain` | Patient/imaging geometry (CT/MVCT), beam/source/sensor models, binary MLC + collimator geometry, helical delivery kinematics. Landed (0.0.1): `VoxelGrid` + `Volume`, `HelicalDelivery`, binary-MLC `MlcModel`, and DICOM ingest (`load_ct_slice` single-slice + `load_ct_series` multi-slice → 3-D HU `Volume`, via `ritk-dicom`, feature `dicom`). HU-semantic newtypes, oriented pose, and gaia MLC geometry pending. | **partial (0.0.1)** |
| `helios-physics` | Radiation interaction physics: photon/electron transport, scatter, attenuation, projection physics. Beer–Lambert attenuation + HU→density landed (0.0.1); NIST μ/ρ data, ray-marched line integral, electron transport pending. | **partial (0.0.1)** |
| `helios-solver` | GPU-accelerated deterministic dose engines (collapsed-cone / convolution-superposition) and imaging forward projectors/reconstruction. Landed (0.0.1): HU→μ engine, `forward_project_ray` (∫μ dl), primary-fluence transport, `deposit_ray_terma` (primary-energy terma deposition, exact `w·(1−e^{−τ})` conservation), `scatter_superposition` (separable 3-D kernel superposition → lateral penumbra + build-up). Anisotropic CC kernel + divergent fan + GPU backend pending. | **partial (0.0.1)** |
| `helios-analysis` | Dosimetric analysis (DVH, gamma), imaging quality metrics, visualization. Landed (0.0.1): cumulative + structure-masked (per-PTV/OAR) DVH, 3-D gamma index (3%/2 mm, global + local normalization with low-dose cutoff), and MVCT `image_quality` metrics (RMSE / relative-L2 accuracy, ROI noise, Michelson contrast, CNR). RT-struct rasterization + portal dosimetry pending. | **partial (0.0.1)** |
| `helios-simulation` | Time-dependent helical TomoTherapy delivery with synchronized MVCT and motion. Landed (0.0.1): `simulate_helical_sinogram` (MVCT acquisition), `simulate_helical_delivery` (MLC+kinematics → `DeliveryFrame`s), `accumulate_delivered_dose` (per-frame per-leaf beamlets → delivered-dose `Volume`, closing the delivery→dose loop) with a `BeamGeometry` seam (parallel or divergent point-source fan). moirai parallel dispatch + cone-beam detector + motion pending. | **partial (0.0.1)** |
| `helios-planning` | Inverse planning / optimization (gradient-based, multi-criteria). Projected-gradient beam-weight optimizer (`DoseInfluence` + `optimize_beam_weights`) landed (0.0.1); coeus-autodiff backend for non-quadratic objectives pending (H-031b). | **partial (0.0.1)** |
| `helios-imaging` | MVCT acquisition modeling, reconstruction, portal dosimetry, IGRT workflows. Landed (0.0.1): parallel-beam Radon forward transform (`Sinogram`), Ram-Lak FBP + SIRT iterative reconstruction (shared back-projector), a deterministic quantum-noise model (`add_quantum_noise`), and IGRT rigid translation registration (`register_translation`). Portal dosimetry + sub-voxel/deformable registration (ritk) pending. | **partial (0.0.1)** |
| `helios-gpu` | GPU dispatch over `hephaestus_core::ComputeDevice` + hephaestus-wgpu. `beam_transmission_into` (GPU `exp(−τ)`) landed (0.0.1), differentially validated vs CPU; more kernels + throughput benchmarks pending. | **partial (0.0.1)** |
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
| **ritk** | `ritk-dicom` (**consumed**, feature `dicom`), `ritk-core`, `ritk-registration` | domain, analysis, imaging | DICOM I/O (CT/MVCT via `ritk-dicom` — parse + rescale-calibrated decode); RT struct/plan/dose, registration, VTK pending. |
| **gaia** | `gaia` | domain | Geometry kernel: binary MLC, collimators, jaws, patient surface/mesh. |
| **hephaestus** | `hephaestus-core`, `hephaestus-wgpu` (`-cuda`, `-metal` optional) | gpu, solver | GPU compute dispatch, WGSL pipelines, kernel caching. |
| **moirai** | `moirai`, `moirai-parallel` (`-async`, `-gpu`, `-iter`) | simulation, solver | Orchestration of time-dependent helical delivery + imaging; execution policies. |
| **coeus** | `coeus-core`, `coeus-tensor`, `coeus-autograd`, `coeus-optim` | planning, solver | Tensors, autodiff for optimization/sensitivity, neural dose/imaging surrogates. |
| **consus** | `consus-core`, `consus-hdf5`, `consus-io`, `consus-compression` | analysis, imaging, simulation | Volumetric storage (HDF5/Zarr), compression of dose/CT/MVCT datasets. |
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
