# Chapter 7 — Filtered Back Projection

<!-- generated-figure-start -->
![Figure 7.1 — Filtered Back Projection](figures/ch07/fig01_7_filtered_back_projection.svg)
*Figure 7.1 — Filtered Back Projection*
<!-- generated-figure-end -->

Filtered back projection (FBP) reconstructs a 2-D image from a sinogram by:

1. **Ramp filtering** each projection in the Fourier domain: |ω| · P̂(θ, ω)
2. **Back-projecting** filtered projections across all angles

```text
use helios_imaging::filtered_back_projection;

// `sinogram` is a `Sinogram<T>`; `recon` is the target `VoxelGrid<T>`.
let recon = filtered_back_projection(&sinogram, &recon_grid);
```

## Frequency Filter

The ramp filter amplifies high frequencies, sharpening edges.
A Hanning window can be applied to suppress noise:

```text
H(ω) = |ω| · W(ω),   W(ω) = 0.5 + 0.5 cos(π|ω|/ω_max)
```

## Ramp Filtering

Ramp filtering is the reconstruction's only super-linear-cost stage: the Ram-Lak
kernel spans the **full** detector width, so applying it as a spatial convolution
costs `O(n_ang · n_off²)`.

`helios_imaging::ramp` performs it as a zero-padded linear convolution through
the Atlas transform provider (`apollo-fft`) instead, at
`O(n_ang · n_off log n_off)`. The pad length is `N ≥ 3·n_off − 2` — exactly the
length of the linear convolution of an `n_off`-sample row with the
`2·n_off − 1`-tap kernel — so no output sample is touched by circular
wrap-around, and the transform path **reproduces** the spatial convolution
rather than approximating it. The spatial form is retained as
`RampMethod::Direct`, the reference the transform path is differentially tested
against.

`filtered_back_projection` selects the transform path automatically above the
crossover width. `ramp_filter_rows` is the stage on its own, for callers that
need to choose the method explicitly.

## Accuracy

FBP is an approximate inversion for finite angular sampling:
- **Reconstruction RMSE** depends on 
_angles and detector pitch
- For the 61-voxel phantom in the tomotherapy example: RMSE < 0.005 cm⁻¹

## Further Reading

- [Example: Tomotherapy Workflow](examples/tomotherapy_workflow.md)
- [Parallel-Beam Radon Transform](imaging_radon.md)
