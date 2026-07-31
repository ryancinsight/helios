# Summary

[Introduction](README.md)

# Part I — Foundations

- [1. Physics Domain Types and Safety Boundaries](foundations.md)
  - [Example: Validating Foundation Units](examples/validate_foundation_units.md)
- [2. Voxel Grids and Volumetric Data](domain_geometry.md)
  - [Example: VoxelGrid and Volume Construction](examples/voxel_grid_construction.md)
- [3. Scalar Fields and Numeric Abstractions](numerics.md)
- [4. Memory and Allocation: Mnemosyne Integration](memory.md)

# Part II — CT Imaging and Attenuation

- [5. Hounsfield Units and Attenuation Maps](imaging_ct.md)
  - [Example: Photon Attenuation Physics](examples/photon_attenuation.md)
- [6. Parallel-Beam Radon Transform](imaging_radon.md)
  - [Example: Radon Sinogram](examples/radon_sinogram.md)
- [7. Filtered Back Projection](imaging_fbp.md)
  - [Example: FBP Reconstruction](examples/fbp_reconstruction.md)
- [8. MVCT and Correction Workflows](imaging_mvct.md)
  - [Example: SIRT Iterative Reconstruction](examples/sirt_reconstruction.md)
  - [Example: IGRT Setup Correction via Registration](examples/mvct_registration.md)

# Part III — Dose Calculation

- [9. Mass Attenuation and Photon Cross Sections](dose_attenuation.md)
  - [Example: Compton Scattering Physics](examples/compton_physics.md)
- [10. Terma and Energy Deposition](dose_terma.md)
- [11. Collapsed-Cone Convolution](dose_collapsed_cone.md)
  - [Example: Collapsed-Cone 3-D Dose Engine](examples/collapsed_cone_3d.md)
- [12. Beam Hardening and Poly-Energetic Spectra](dose_spectra.md)

# Part IV — Treatment Delivery and Planning

- [13. MLC Models and Leaf Sequencing](planning_mlc.md)
- [14. Helical Delivery and Sinogram](planning_helical.md)
- [15. Dose-Volume Histograms](planning_dvh.md)
  - [Example: DVH Analysis](examples/dvh_analysis.md)
  - [Example: DVH-Constrained Beam-Weight Optimization](examples/dvh_optimization.md)
- [16. Gamma Index and Plan Verification](planning_gamma.md)
  - [Example: Gamma Index Comparison](examples/gamma_index.md)

# Part V — End-to-End Clinical Workflows

- [17. TomoTherapy End-to-End Workflow](workflow_tomotherapy.md)
  - [Example: Tomotherapy Workflow](examples/tomotherapy_workflow.md)
- [18. LINAC-Based Step-and-Shoot Delivery](workflow_linac.md)
  - [Example: LINAC Dose Accumulation](examples/linac_dose_accumulation.md)
- [19. Adaptive Radiotherapy with MVCT](workflow_adaptive.md)
  - [Example: Adaptive RT Workflow](examples/adaptive_rt_workflow.md)

# Part VI — GPU Acceleration

- [20. GPU Backend Overview: Hephaestus Integration](gpu_overview.md)
- [21. GPU-Accelerated Dose Kernels](gpu_dose.md)
  - [Example: GPU Attenuation Map and Forward Projection](examples/gpu_attenuation_projection.md)
- [22. Coeus Tensor Operations for Dose Grids](gpu_coeus.md)

# Part VII — Atlas Stack Migration

- [23. Migration Overview](migration_overview.md)
- [24. Eunomia: Numeric Traits](migration_eunomia.md)
- [25. Leto: Arrays and Linalg](migration_arrays.md)
- [26. Leto: Geometry](migration_geometry.md)
- [27. Hermes: SIMD Lanes](migration_simd.md)
- [28. Mnemosyne and Themis: Memory](migration_memory.md)
- [29. Moirai: Concurrency](migration_concurrency.md)
- [30. Apollo: FFT](migration_fft.md)
- [31. Leto: GAT Tiling](migration_gat_tiles.md)
- [32. Coeus: Tensors and Autodiff](migration_coeus.md)
- [33. Ritk: Image I/O](migration_image_io.md)
- [34. Migration Validation](migration_validation.md)

# Part VIII — Validation and Benchmarking

- [35. Reference Phantoms and Ground Truth](validation_phantoms.md)
- [36. Analytical Solutions and Regression Tests](validation_regression.md)
  - [Example: Regression and Analytical Validation](examples/validation_regression.md)
- [37. Clinical Protocol Compliance](validation_clinical.md)
  - [Example: Clinical Protocol Validation](examples/validation_clinical.md)

# Appendix

- [A. Atlas Crate Dependency Map](appendix_dependencies.md)
- [B. Atlas Glossary](appendix_glossary.md)
- [C. API Reference Index](appendix_api.md)
- [D. Changelog](appendix_changelog.md)
- [E. Book Organization Forward Roadmap](BOOK_ORGANIZATION.md)
