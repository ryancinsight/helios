# Appendix C — API Reference Index

Generated documentation is available via cargo doc:

`ash
# Helios full workspace docs
cd path/to/helios
cargo doc --workspace --no-deps --open
`

## Key Public APIs by Crate

### helios-core
- EnergyMeV, HounsfieldUnit, VoxelSpacingMm — validated domain types
- HeliosError — unified error type

### helios-domain  
- VoxelGrid — axis-aligned 3D grid with coordinate mapping
- Volume<T> — dense scalar field over a VoxelGrid

### helios-imaging
- parallel_beam_radon — forward Radon transform
- iltered_back_projection — FBP reconstruction

### helios-simulation
- simulate_helical_delivery — TomoTherapy terma simulation
- ccumulate_delivered_dose_anisotropic — Collapsed-cone dose
- CollapsedCone — poly-energetic dose kernel configuration

### helios-analysis
- Dvh — Dose-Volume Histogram
- gamma_index_3d — 3D gamma metric
- gamma_pass_rate — fraction of passing voxels
- 
oi_statistics — ROI mean/std/min/max

## Further Reading

- [Changelog](appendix_changelog.md)