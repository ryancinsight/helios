# Summary

[Introduction](README.md)

---

# Part I — Foundations

- [Physics Domain Types and Safety Boundaries](foundations.md)
  - [Example: Validating Foundation Units](examples/validate_foundation_units.md)
- [Voxel Grids and Volumetric Data](domain_geometry.md)
  - [Example: VoxelGrid and Volume Construction](examples/voxel_grid_construction.md)
- [Scalar Fields and Numeric Abstractions](numerics.md)
- [Memory and Allocation: Mnemosyne Integration](memory.md)

---

# Part II — CT Imaging and Attenuation

- [Hounsfield Units and Attenuation Maps](imaging_ct.md)
  - [Example: Photon Attenuation Physics](examples/photon_attenuation.md)
- [Parallel-Beam Radon Transform](imaging_radon.md)
  - [Example: Radon Sinogram](examples/radon_sinogram.md)
- [Filtered Back Projection](imaging_fbp.md)
  - [Example: FBP Reconstruction](examples/fbp_reconstruction.md)
- [MVCT and Correction Workflows](imaging_mvct.md)
  - [Example: SIRT Iterative Reconstruction](examples/sirt_reconstruction.md)
  - [Example: IGRT Setup Correction via Registration](examples/mvct_registration.md)

---

# Part III — Dose Calculation

- [Mass Attenuation and Photon Cross Sections](dose_attenuation.md)
  - [Example: Compton Scattering Physics](examples/compton_physics.md)
- [Terma and Energy Deposition](dose_terma.md)
- [Collapsed-Cone Convolution](dose_collapsed_cone.md)
  - [Example: Collapsed-Cone 3-D Dose Engine](examples/collapsed_cone_3d.md)
- [Beam Hardening and Poly-Energetic Spectra](dose_spectra.md)

---

# Part IV — Treatment Delivery and Planning

- [MLC Models and Leaf Sequencing](planning_mlc.md)
- [Helical Delivery and Sinogram](planning_helical.md)
- [Dose-Volume Histograms](planning_dvh.md)
  - [Example: DVH Analysis](examples/dvh_analysis.md)
  - [Example: DVH-Constrained Beam-Weight Optimization](examples/dvh_optimization.md)
- [Gamma Index and Plan Verification](planning_gamma.md)
  - [Example: Gamma Index Comparison](examples/gamma_index.md)

---

# Part V — End-to-End Clinical Workflows

- [TomoTherapy End-to-End Workflow](workflow_tomotherapy.md)
  - [Example: Tomotherapy Workflow](examples/tomotherapy_workflow.md)
- [LINAC-Based Step-and-Shoot Delivery](workflow_linac.md)
  - [Example: LINAC Dose Accumulation](examples/linac_dose_accumulation.md)
- [Adaptive Radiotherapy with MVCT](workflow_adaptive.md)
  - [Example: Adaptive RT Workflow](examples/adaptive_rt_workflow.md)

---

# Part VI — GPU Acceleration

- [GPU Backend Overview: Hephaestus Integration](gpu_overview.md)
- [GPU-Accelerated Dose Kernels](gpu_dose.md)
  - [Example: GPU Attenuation Map and Forward Projection](examples/gpu_attenuation_projection.md)
- [Coeus Tensor Operations for Dose Grids](gpu_coeus.md)

---

# Part VII — Validation and Benchmarking

- [Reference Phantoms and Ground Truth](validation_phantoms.md)
- [Analytical Solutions and Regression Tests](validation_regression.md)
  - [Example: Regression and Analytical Validation](examples/validation_regression.md)
- [Clinical Protocol Compliance](validation_clinical.md)
  - [Example: Clinical Protocol Validation](examples/validation_clinical.md)

---

# Part VIII — Atlas Stack Integration (Migration Reference)

Helios is one primary consumer of the Atlas stack (the other is kwavers).
This part documents the Atlas-first design and migration surface area.

- [Migration Overview: ndarray/nalgebra/burn → Atlas](migration_overview.md)
- [Eunomia: Numeric Trait Unification](migration_eunomia.md)
- [Leto: Arrays and Linear Algebra](migration_arrays.md)
- [Leto: Geometry — VoxelGrid, MLC, Beam Isometries](migration_geometry.md)
- [Hermes: SIMD Lanes and Vectorized Kernels](migration_simd.md)
- [Mnemosyne and Themis: Memory](migration_memory.md)
- [Moirai: Concurrency](migration_concurrency.md)
- [Apollo: FFT and Spectral Methods](migration_fft.md)
- [Leto: GAT-Based Tiling and Lending Iterators](migration_gat_tiles.md)
- [Coeus: Tensors and Autodiff](migration_coeus.md)
- [Ritk: Image I/O — DICOM, NIfTI, PNG](migration_image_io.md)
- [Migration Validation: TG-119 and Atlas Parity](migration_validation.md)

---

# Appendices

- [Atlas Crate Dependency Map](appendix_dependencies.md)
- [Atlas Glossary](appendix_glossary.md)
- [API Reference Index](appendix_api.md)
- [Changelog](appendix_changelog.md)
- [Book Organization Forward Roadmap](BOOK_ORGANIZATION.md)
