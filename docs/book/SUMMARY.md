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

---

# Part III — Dose Calculation

- [Mass Attenuation and Photon Cross Sections](dose_attenuation.md)
- [Terma and Energy Deposition](dose_terma.md)
- [Collapsed-Cone Convolution](dose_collapsed_cone.md)
- [Beam Hardening and Poly-Energetic Spectra](dose_spectra.md)

---

# Part IV — Treatment Delivery and Planning

- [MLC Models and Leaf Sequencing](planning_mlc.md)
- [Helical Delivery and Sinogram](planning_helical.md)
- [Dose-Volume Histograms](planning_dvh.md)
  - [Example: DVH Analysis](examples/dvh_analysis.md)
- [Gamma Index and Plan Verification](planning_gamma.md)
  - [Example: Gamma Index Comparison](examples/gamma_index.md)

---

# Part V — End-to-End Clinical Workflows

- [TomoTherapy End-to-End Workflow](workflow_tomotherapy.md)
  - [Example: Tomotherapy Workflow](examples/tomotherapy_workflow.md)
- [LINAC-Based Step-and-Shoot Delivery](workflow_linac.md)
- [Adaptive Radiotherapy with MVCT](workflow_adaptive.md)

---

# Part VI — GPU Acceleration

- [GPU Backend Overview: Hephaestus Integration](gpu_overview.md)
- [GPU-Accelerated Dose Kernels](gpu_dose.md)
- [Coeus Tensor Operations for Dose Grids](gpu_coeus.md)

---

# Part VII — Validation and Benchmarking

- [Reference Phantoms and Ground Truth](validation_phantoms.md)
- [Analytical Solutions and Regression Tests](validation_regression.md)
- [Clinical Protocol Compliance](validation_clinical.md)

---

# Appendices

- [Atlas Crate Dependency Map](appendix_dependencies.md)
- [API Reference Index](appendix_api.md)
- [Changelog](appendix_changelog.md)
