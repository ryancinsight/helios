# Appendix E — Book Organization Forward Roadmap

## Overview

This document outlines the planned book structure for the Helios radiation therapy 
simulation suite, following the kwavers model. The book will document the physics, 
mathematics, and implementation of Helios using the Atlas stack.

> **Authoritative TOC:** [`SUMMARY.md`](SUMMARY.md) is the single source of truth for
> the *current* book (the file mdBook renders). This document is the *forward roadmap*
> — its Parts/Chapters describe intended future expansion and may run ahead of the
> chapters and examples that exist today.

---

## Book Title

**Helios: Radiation Therapy Simulation, Planning, and Imaging**

---

## Parts and Chapters

### Part I: Foundations

#### Chapter 1: Introduction to Helios
- Overview of radiation therapy simulation
- Atlas stack integration
- Architecture principles (SRP, SoC, SSOT, DIP, DRY)

#### Chapter 2: Numerical Methods
- Finite difference methods (FDM)
- Monte Carlo methods
- Linear system solvers
- Time integration schemes

#### Chapter 3: Radiation Physics
- Photon interactions
- Electron transport
- Dose calculation algorithms
- Monte Carlo for radiation transport

#### Chapter 4: Patient Modeling
- DICOM RT structure sets
- Voxelized phantoms
- Geometry representation (Gaia)

### Part II: Core Simulation

#### Chapter 5: Dose Calculation
- Pencil beam algorithms
- Monte Carlo dose calculation
- Dose convolution methods
- GPU acceleration with Hephaestus

#### Chapter 6: Treatment Planning
- Beam optimization
- IMRT/VMAT planning
- Dose-volume constraints
- Optimization algorithms

#### Chapter 7: Monte Carlo Methods
- Photon and electron transport
- Variance reduction techniques
- Parallelization with Moirai
- Validation and benchmarking

### Part III: Imaging and Verification

#### Chapter 8: Cone Beam CT
- Image reconstruction
- Filtered back projection
- Iterative reconstruction
- Image quality metrics

#### Chapter 9: Dose Verification
- DVH analysis
- Gamma comparison
- 3D dose visualization
- QA phantoms

#### Chapter 10: Image Guidance
- Real-time tracking
- Motion management
- Image registration
- fiducial-based alignment

### Part IV: Advanced Applications

#### Chapter 11: Brachytherapy
- High-dose-rate (HDR) planning
- Dwell time optimization
- Source modeling
- Inverse planning

#### Chapter 12: Proton and Ion Therapy
- Bragg peak modeling
- Range calculation
- LET and biological dose
- Optimization with uncertainty

#### Chapter 13: Adaptive Therapy
- Daily adaptation
- Deformation modeling
- Fast replanning
- Real-time adaptation

### Part V: Atlas Stack Integration

#### Chapter 14: Memory Management with Mnemosyne
- Memory allocation strategies
- NUMA locality
- Placement-aware data structures

#### Chapter 15: Geometry with Gaia and Leto
- Patient geometry representation
- Beam modeling
- Dose grid operations

#### Chapter 16: Parallel Computing with Moirai
- Monte Carlo parallelization
- Dose calculation parallelization
- GPU-accelerated parallelization

#### Chapter 17: GPU Computing with Hephaestus
- GPU dose calculation kernels
- Memory management
- Kernel caching

#### Chapter 18: Python Integration
- PyO3 bindings
- API design
- Example workflows

### Part VI: Validation and Benchmarking

#### Chapter 19: Validation Framework
- TG-119 compliance
- AAPM test cases
- Comparison with clinical TPS

#### Chapter 20: Performance Benchmarking
- Monte Carlo benchmarks
- Dose calculation benchmarks
- Parallel scaling studies

---

## Example Structure

Following the kwavers pattern, each chapter will have corresponding examples:

```
crates/helios/examples/
├── chapter_5_dose_calculation/
│   ├── basic_dose.rs
│   └── gpu_dose.rs
├── chapter_6_treatment_planning/
│   ├── basic_optimization.rs
│   └── imrt_planning.rs
└── ...
```

---

## Implementation Status

### Currently Implemented

- `helios-core` - Core abstractions and data structures
- `helios-math` - Mathematical primitives (uses leto for arrays)
- `helios-domain` - Domain-specific types and geometry
- `helios-physics` - Physics models and equations
- `helios-solver` - Forward and inverse solvers
- `helios-analysis` - Analysis and validation tools
- `helios-gpu` - GPU acceleration backend
- `helios-simulation` - Simulation orchestration
- `helios-imaging` - Imaging and reconstruction
- `helios-planning` - Treatment planning algorithms
- `helios-python` - Python bindings

### Migration Status

- **Leto**: COMPLETE - All array operations use leto
- **Gaia**: COMPLETE - Geometry uses Gaia primitives
- **Moirai**: COMPLETE - Concurrency uses Moirai
- **Hephaestus**: COMPLETE - GPU uses Hephaestus
- **Apollo**: COMPLETE - FFT uses Apollo
- **Coeus**: COMPLETE - Tensors and autodiff use Coeus
- **Mnemosyne**: COMPLETE - Memory uses Mnemosyne
- **Themis**: COMPLETE - Placement uses Themis
- **Hermes**: COMPLETE - SIMD uses Hermes

---

## Build Instructions

```bash
# Build the book
mdbook build docs/book

# Serve the book (with hot reload)
mdbook serve docs/book
```

## Organization Principles

1. **Executable Examples**: Each chapter has corresponding examples in `crates/helios/examples/`
2. **Deep Hierarchy**: Code organized by SRP with redundancy-free structure
3. **Single Source of Truth**: Atlas crates provide SSOT for all foundational operations
4. **Zero-Cost Abstractions**: Heavy use of const generics, ZSTs, and phantom types
5. **DRY**: No duplicated logic across layers
6. **DIP**: Depend on abstractions, not concrete implementations
7. **Cow**: Copy-on-write where appropriate for safe mutation

---

## Timeline

### Phase 1: Core Implementation
- [x] Implement all helios crates
- [x] Migrate to Atlas stack (leto, gaia, moirai, hephaestus)
- [x] Python bindings (pyo3)
- [x] Unit tests and benchmarks

### Phase 2: Documentation
- [x] Create book structure (authoritative TOC in `SUMMARY.md`; this document is the forward roadmap)
- [~] Write chapter content (Parts I–VII drafted; Advanced Applications chapters pending)
- [x] Create executable examples (8 chapter-mapped examples in `crates/*/examples/`)
- [x] Validate all examples compile *and run* (each passes its analytical oracle)

### Phase 3: Validation
- [ ] Clinical validation studies
- [ ] TG-119 compliance testing
- [ ] Performance benchmarking
- [ ] Peer review

---

## References

- TG-119: AAPM protocol for IMRT QA
- AAPM TG-185: Monte Carlo for proton therapy
- AAPM TG-215: Dose calculation verification
- AAPM TG-218: VMAT QA

---

## Related Documentation

- [kwavers Book](https://github.com/ryancinsight/kwavers/docs/book) - Computational Ultrasound
- [Atlas Stack Documentation](https://github.com/ryancinsight/atlas) - Core infrastructure
