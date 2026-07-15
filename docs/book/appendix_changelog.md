# Changelog

## Unreleased

### Added
- Full TomoTherapy end-to-end simulation pipeline
- Parallel-beam Radon transform and FBP reconstruction
- Collapsed-cone poly-energetic dose calculation
- DVH, gamma index, and ROI statistics analysis
- Three validated examples: foundation units, voxel grid, tomotherapy workflow
- 7-part mdBook structure with complete chapter scaffolding
- Atlas-native dependency graph (zero third-party math libraries)
- GPU sparse kernel support via hephaestus-wgpu

### Dependencies
- Migrated from nalgebra/ndarray to leto (CPU arrays)
- Migrated from tokio/rayon to moirai (async + parallel)
- unomia for numeric trait unification
- coeus for tensor operations and autodiff
- hephaestus for GPU acceleration

## v0.0.1

- Initial project scaffold
- helios-core domain types
- helios-math geometry primitives