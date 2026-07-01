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
| `helios-domain` | Patient/imaging geometry (CT/MVCT), beam/source/sensor models, binary MLC + collimator geometry, helical delivery kinematics. `VoxelGrid` + `Volume` landed (0.0.1); DICOM I/O and beam/MLC models pending. | **partial (0.0.1)** |
| `helios-physics` | Radiation interaction physics: photon/electron transport, scatter, attenuation, projection physics. Beer–Lambert attenuation + HU→density landed (0.0.1); NIST μ/ρ data, ray-marched line integral, electron transport pending. | **partial (0.0.1)** |
| `helios-solver` | GPU-accelerated deterministic dose engines (collapsed-cone / convolution-superposition) and imaging forward projectors/reconstruction. HU→μ material-property engine landed (0.0.1, CPU reference); dose/projector engines + GPU backend pending. | **partial (0.0.1)** |
| `helios-analysis` | Dosimetric analysis (DVH, gamma), imaging quality metrics, visualization. DVH + 3D gamma index landed (0.0.1); structure-masked DVH, local-norm gamma, imaging metrics pending. | **partial (0.0.1)** |
| `helios-simulation` | Time-dependent helical TomoTherapy delivery with synchronized MVCT and motion. `simulate_helical_sinogram` (gantry+couch → per-projection forward projection) landed (0.0.1); moirai parallel dispatch + fan/cone detector + motion pending. | **partial (0.0.1)** |
| `helios-planning` | Inverse planning / optimization (gradient-based, multi-criteria). | planned (Sprint 4) |
| `helios-imaging` | MVCT acquisition modeling, reconstruction, portal dosimetry, IGRT workflows. Parallel-beam Radon forward transform (`Sinogram`) landed (0.0.1); FBP reconstruction, portal dosimetry, IGRT pending. | **partial (0.0.1)** |
| `helios-gpu` | GPU dispatch over `hephaestus_core::ComputeDevice` + hephaestus-wgpu. `beam_transmission_into` (GPU `exp(−τ)`) landed (0.0.1), differentially validated vs CPU; more kernels + throughput benchmarks pending. | **partial (0.0.1)** |
| `helios-python` | High-level PyO3 API over simulation/planning/imaging. | planned (Sprint 5) |

Crates are created only when their layer is built (architecture_scoping growth
triggers — no speculative empty-crate scaffolding). The workspace `members` list
grows as each crate lands; `workspace.dependencies` already declares the full
Atlas set as the integration SSOT.

## Atlas dependency map

Package names below are verified against each Atlas repo's manifests. Git URLs are
the SSOT in the root `Cargo.toml` `[workspace.dependencies]`.

| Atlas component | Crates (packages) | Consumed by | Purpose in Helios |
|-----------------|-------------------|-------------|-------------------|
| **ritk** | `ritk-core`, `ritk-io`, `ritk-registration` | domain, analysis, imaging | DICOM I/O (CT, MVCT, RT struct/plan/dose), image registration, VTK visualization. |
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
