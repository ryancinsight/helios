# Atlas Crate Dependency Map

Helios sits at the top of the Atlas dependency graph.

## Helios Layer Map

`
helios-simulation   ← top-level integrator
    ├── helios-imaging       (Radon, FBP)
    ├── helios-physics       (attenuation, spectra)
    ├── helios-solver        (attenuation maps)
    ├── helios-analysis      (DVH, gamma)
    ├── helios-domain        (VoxelGrid, Volume)
    ├── helios-math          (Point3, Scalar)
    └── helios-gpu           (hephaestus wgpu kernels)

Atlas Foundation
    ├── coeus-core / coeus-tensor  (tensor ops, autodiff)
    ├── leto                       (CPU arrays)
    ├── hephaestus-wgpu            (GPU arrays)
    ├── moirai                     (async runtime + parallelism)
    ├── mnemosyne                  (memory allocator)
    ├── hermes-simd                (SIMD kernels)
    ├── eunomia                    (numeric traits)
    ├── apollo-fft                 (FFT)
    └── themis                     (topology / placement)
`

## Third-Party Dependencies

| Crate | Purpose | Will be replaced |
|---|---|---|
| dicom (dicom-rs) | DICOM I/O | itk-dicom native path |
| image | PNG output in examples | Retained (output only) |
| nyhow | Error propagation | Retained |

## Further Reading

- [API Reference Index](appendix_api.md)