# Chapter 8 — MVCT and Correction Workflows

Megavoltage CT (MVCT) is acquired on a TomoTherapy unit using the 3.5 MV treatment
beam before each fraction for patient setup verification.

## MVCT vs kVCT

| Property | kVCT (diagnostic) | MVCT (treatment) |
|---|---|---|
| Energy | 80–140 kVp | 3.5 MV |
| HU range | Full clinical | Reduced contrast |
| Spatial resolution | 0.5 mm | ~1 mm |
| Use | Treatment planning | Daily alignment |

## Workflow in Helios

MVCT reconstruction uses the same FBP pipeline as kVCT, with a different
beam-hardening correction appropriate for MV photons:

`
ust
let mvct_sinogram = parallel_beam_radon(&mvct_projection, n_angles);
let mvct_recon = filtered_back_projection(&mvct_sinogram, n_angles, nx);
`

## Further Reading

- [Filtered Back Projection](imaging_fbp.md)
- [Adaptive Radiotherapy](workflow_adaptive.md)