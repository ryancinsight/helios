# Summary

[Introduction](README.md)

## Deterministic figures

- [Photon transmission vs depth](figures/photon_attenuation_depth.svg)
- [HU to relative electron density calibration](figures/ct_calibration_curve.svg)
- [Single-angle Radon sinogram profile](figures/radon_sinogram_disk.svg)
- [Differential DVH](figures/dvh_curve.svg)
- [Central-slice dose heatmap](figures/dose_slice_heatmap.svg)
- [MLC leaf-open sinogram](figures/helical_mlc_fluence.svg)
- [Helios layered architecture](figures/architecture_stack.svg)

---

# Part I — Foundations

- [1. Physics Domain Types and Safety Boundaries](foundations.md)
  - [Example: Validating Foundation Units](examples/validate_foundation_units.md)
- [2. Voxel Grids and Volumetric Data](domain_geometry.md)
  - [Example: VoxelGrid and Volume Construction](examples/voxel_grid_construction.md)
- [3. Scalar Fields and Numeric Abstractions](numerics.md)
- [4. Memory and Allocation: Mnemosyne Integration](memory.md)

---

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

---

# Part III — Dose Calculation

- [9. Mass Attenuation and Photon Cross Sections](dose_attenuation.md)
  - [Example: Compton Scattering Physics](examples/compton_physics.md)
- [10. Terma and Energy Deposition](dose_terma.md)
- [11. Collapsed-Cone Convolution](dose_collapsed_cone.md)
  - [Example: Collapsed-Cone 3-D Dose Engine](examples/collapsed_cone_3d.md)
- [12. Beam Hardening and Poly-Energetic Spectra](dose_spectra.md)

---

# Part IV — Treatment Delivery and Planning

- [13. MLC Models and Leaf Sequencing](planning_mlc.md)
- [14. Helical Delivery and Sinogram](planning_helical.md)
- [15. Dose-Volume Histograms](planning_dvh.md)
  - [Example: DVH Analysis](examples/dvh_analysis.md)
  - [Example: DVH-Constrained Beam-Weight Optimization](examples/dvh_optimization.md)
- [16. Gamma Index and Plan Verification](planning_gamma.md)
  - [Example: Gamma Index Comparison](examples/gamma_index.md)

---

# Part V — End-to-End Clinical Workflows

- [17. TomoTherapy End-to-End Workflow](workflow_tomotherapy.md)
  - [Example: Tomotherapy Workflow](examples/tomotherapy_workflow.md)
- [18. LINAC-Based Step-and-Shoot Delivery](workflow_linac.md)
  - [Example: LINAC Dose Accumulation](examples/linac_dose_accumulation.md)
- [19. Adaptive Radiotherapy with MVCT](workflow_adaptive.md)
  - [Example: Adaptive RT Workflow](examples/adaptive_rt_workflow.md)

---

# Part VI — GPU Acceleration

- [20. GPU Backend Overview: Hephaestus Integration](gpu_overview.md)
- [21. GPU-Accelerated Dose Kernels](gpu_dose.md)
  - [Example: GPU Attenuation Map and Forward Projection](examples/gpu_attenuation_projection.md)
- [22. Coeus Tensor Operations for Dose Grids](gpu_coeus.md)

---

# Part VII — Validation and Benchmarking

- [23. Reference Phantoms and Ground Truth](validation_phantoms.md)
- [24. Analytical Solutions and Regression Tests](validation_regression.md)
  - [Example: Regression and Analytical Validation](examples/validation_regression.md)
- [25. Clinical Protocol Compliance](validation_clinical.md)
  - [Example: Clinical Protocol Validation](examples/validation_clinical.md)

---

# Part VIII — Atlas Stack Integration (Migration Reference)

Helios is one primary consumer of the Atlas stack (the other is kwavers).
This part documents the Atlas-first design and migration surface area.

- [26. Migration Overview: ndarray/nalgebra/burn → Atlas](migration_overview.md)
- [27. Eunomia: Numeric Trait Unification](migration_eunomia.md)
- [28. Leto: Arrays and Linear Algebra](migration_arrays.md)
- [29. Leto: Geometry — VoxelGrid, MLC, Beam Isometries](migration_geometry.md)
- [30. Hermes: SIMD Lanes and Vectorized Kernels](migration_simd.md)
- [31. Mnemosyne and Themis: Memory](migration_memory.md)
- [32. Moirai: Concurrency](migration_concurrency.md)
- [33. Apollo: FFT and Spectral Methods](migration_fft.md)
- [34. Leto: GAT-Based Tiling and Lending Iterators](migration_gat_tiles.md)
- [35. Coeus: Tensors and Autodiff](migration_coeus.md)
- [36. Ritk: Image I/O — DICOM, NIfTI, PNG](migration_image_io.md)
- [37. Migration Validation: TG-119 and Atlas Parity](migration_validation.md)

---

# Appendices

- [A. Atlas Crate Dependency Map](appendix_dependencies.md)
- [B. Atlas Glossary](appendix_glossary.md)
- [C. API Reference Index](appendix_api.md)
- [D. Changelog](appendix_changelog.md)
- [E. Book Organization Forward Roadmap](BOOK_ORGANIZATION.md)

- [Stray test figure (drift fixture)](figures/atlas_drift_fixture_test_only_no_such_figure.svg)
