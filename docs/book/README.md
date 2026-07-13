# Helios Radiotherapy Simulation Suite

Helios is a high-performance, pure-Rust radiotherapy simulation library built on the
[Atlas physics stack](https://github.com/ryancinsight).  It covers the full treatment
pipeline — CT image processing, dose calculation, treatment planning, delivery simulation,
and clinical plan verification — with a consistent zero-copy, zero-cost-abstraction design.

## Architecture at a Glance

```
                ┌──────────────────────────────────────────────────────────┐
                │                    helios-simulation                     │
                │         (end-to-end clinical workflow orchestration)     │
                └───┬──────────────┬──────────────┬───────────────┬───────┘
                    │              │              │               │
          ┌─────────▼──┐  ┌───────▼──────┐ ┌────▼──────┐ ┌─────▼───────┐
          │helios-solver│  │helios-planning│ │helios-    │ │helios-      │
          │ (dose calc) │  │ (opt + DVH)  │ │imaging    │ │analysis     │
          └──────┬──────┘  └──────────────┘ └───────────┘ └─────────────┘
                 │
          ┌──────▼──────┐
          │helios-physics│  (attenuation, cross-sections, spectral models)
          └──────┬───────┘
                 │
          ┌──────▼──────┐
          │helios-domain│  (VoxelGrid, Volume, MLC, helical delivery)
          └──────┬───────┘
                 │
          ┌──────▼──────┐
          │helios-math  │  (Scalar seam, Point3, Vector3 via eunomia/leto)
          └──────┬───────┘
                 │
          ┌──────▼──────┐
          │helios-core  │  (EnergyMeV, HounsfieldUnit, VoxelSpacingMm, errors)
          └─────────────┘
```

## Atlas Dependencies

| Helios need | Atlas crate |
|---|---|
| Scalar field traits | `eunomia` (RealField, FloatElement) |
| Linear algebra | `leto` (Point3, Vector3, Array3) |
| SIMD kernels | `hermes-simd` |
| Parallel iteration | `moirai-parallel` |
| Memory allocation | `mnemosyne` |
| GPU tensors | `hephaestus-core` / `hephaestus-wgpu` |
| Image I/O | `ritk-image`, `ritk-dicom` |
| Storage | `consus-hdf5` |

## Getting Started

```
cargo run -p helios-core --example validate_foundation_units
cargo run -p helios-domain --example voxel_grid_construction
cargo run -p helios-simulation --example tomotherapy_workflow -- /tmp/helios_output
```

## Chapters

This book progresses from the lowest-level physics types ([Part I](foundations.md))
through imaging and dose calculation to full end-to-end workflows ([Part V](workflow_tomotherapy.md)).
Each chapter links to runnable examples that exercise the described functionality.
