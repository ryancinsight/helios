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

## Figures

The figures linked from the chapter tables of contents are deterministically
generated SVG assets committed with the chapter text. Each entry names the
example that the figure mirrors; reproducibility is enforced by the
`xtask prebook` manifest at `figures/MANIFEST.json`.

- [Photon transmission vs depth (water, 100 keV)](figures/photon_attenuation_depth.svg)
  — mirrors `helios-physics` example `photon_attenuation`.
- [HU → relative electron density calibration](figures/ct_calibration_curve.svg)
  — mirrors `helios-physics` example `photon_attenuation` (CT calibration block).
- [Single-angle sinogram profile (θ = 0)](figures/radon_sinogram_disk.svg)
  — mirrors `helios-imaging` example `radon_sinogram`.
- [Differential DVH (Gaussian phantom)](figures/dvh_curve.svg)
  — mirrors `helios-analysis` example `dvh_analysis`.
- [Central-slice dose heatmap (helical TomoTherapy)](figures/dose_slice_heatmap.svg)
  — mirrors `helios-simulation` example `tomotherapy_workflow`.
- [MLC leaf-open sinogram (21 leaves × 40 projections)](figures/helical_mlc_fluence.svg)
  — mirrors `helios-simulation` example `tomotherapy_workflow` (modulated
  central-band aperture).
- [Helios layered architecture on the Atlas stack](figures/architecture_stack.svg)
  — pure schematic, no example source.
