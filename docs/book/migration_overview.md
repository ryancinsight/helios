# Chapter 26 — Migration Overview: ndarray/nalgebra/burn → Atlas

Helios is one of two primary consumers of the Atlas stack (the other is
kwavers).  This Part documents the destination crates, the principles they
encode, and the per-crate migration surface that each module touches.

## Why Atlas

The Atlas stack replaces six categories of third-party dependency with one
first-party monolith:

| Legacy dependency | Atlas replacement | Reason |
|---|---|---|
| `ndarray`, `nalgebra::Matrix` | `leto::NdArray<T,D>` (CPU) | one dense type across all crates |
| `nalgebra::Point3/Vector3/Isometry3` | `leto::Point3/Vector3/Isometry3` | same geometry on CPU/GPU |
| `tokio`, `rayon` | `moirai` executor + task graph | one runtime, async + parallel |
| `packed_simd` | `hermes-simd` lanes | portable SSE/AVX/NEON/WASM |
| `rustfft` / `realfft` | `apollo` forward FFT + autodiff | spectral w/ gradient pipeline |
| `num-traits`, `num-complex` | `eunomia::RealField/ComplexField` | one trait frontier |
| `burn::Tensor`, `tch` | `coeus::Tensor<T>` + autograd | shared autodiff backends |
| `image`, `dicom` | `ritk-image`, `ritk-dicom` | DICOM + NIfTI + PNG unified |
| `alloc` / `jemalloc` / `mimalloc` | `mnemosyne::Arena` + `themis` | NUMA-aware, zero-fragmentation |

Atlas enforces **SRP** (one well-defined surface per crate), **SSOT** (no
overlap), **DIP** (depend on [`eunomia`] traits, not impls), **DRY** (no
duplicated kernels), and zero-cost abstractions (ZSTs, `PhantomData`,
const generics, GATs).

## Migration Status (Helios)

| Module | Status | Notes |
|---|---|---|
| `helios-core` | COMPLETE | `EnergyMeV`/`HounsfieldUnit`/`VoxelSpacingMm` validating newtypes |
| `helios-math` | COMPLETE | scalar seam via `eunomia::RealField` |
| `helios-domain` | COMPLETE | VoxelGrid + Volume over `leto::NdArray` |
| `helios-physics` | COMPLETE | attenuation tables ported to `eunomia` traits |
| `helios-solver` | COMPLETE | collapsed-cone convolution via `apollo` (when spectral) |
| `helios-analysis` | COMPLETE | DVH / Gamma index over `leto::NdArray<F, Ix3>` |
| `helios-gpu` | COMPLETE | `hephaestus::Backend` selection (`Wgpu`, `Cuda`) |
| `helios-imaging` | COMPLETE | FBP, SIRT, Radon transform over `leto::NdArray<F,Ix2/Ix3>` |
| `helios-planning` | COMPLETE | DVH-constrained optimization via `coeus::Tensor` |
| `helios-simulation` | COMPLETE | orchestration via `moirai::Executor` |
| `helios-python` | COMPLETE | PyO3 boundary, vectors re-exported as NumPy arrays |

Helios is in the **cleanup** phase — the only remaining work is parity
benchmarks and legacy-crate pruning (see
[`BOOK_ORGANIZATION.md`](BOOK_ORGANIZATION.md)).

## How To Read This Part

The sub-chapters are organized **from the trait layer outward**.

1. [Eunomia: Numeric Traits](migration_eunomia.md) — the trait frontier every other crate depends on.
2. [Leto: Arrays and Linalg](migration_arrays.md) — CPU dense and sparse, plus geometry.
3. [Leto: Geometry](migration_geometry.md) — points, vectors, isometries, VoxelGrid, MLC.
4. [Hermes: SIMD Lanes](migration_simd.md) — vectorized Radon inverse, attenuation map sweeps.
5. [Mnemosyne and Themis: Memory](migration_memory.md) — dose-volume storage, NUMA-aware pools.
6. [Moirai: Concurrency](migration_concurrency.md) — accelerated tomographic reconstruction parallelization.
7. [Apollo: FFT](migration_fft.md) — collapsed-cone spectral path.
8. [Leto: GAT Tiling](migration_gat_tiles.md) — sliding-window tile streaming over dose grids.
9. [Coeus: Tensors and Autodiff](migration_coeus.md) — DVH-constrained optimization and adjoint reconstruction.
10. [Ritk: Image I/O](migration_image_io.md) — DICOM, NIfTI, PNG ingestion.
11. [Migration Validation](migration_validation.md) — TG-119 compliance + dose-vs-baseline parity.

## Performance Contract

Atlas promises the same per-flop throughput as the legacy crates but
**strictly better constants**:

- Zero heap allocations on hot paths (typed arenas via `mnemosyne`).
- Zero virtual dispatch on numeric kernels (specialization via `eunomia`).
- Zero copy for read-only dose-volume views (`Cow` over `leto::NdArray`).
- Zero abstraction tax (GATs for lifetime variance, ZSTs for typing).

Where Atlas does not yet match legacy throughput, the
[Migration Validation](migration_validation.md) chapter records the gap
and governs the cleanup pass that closes it.

## Cross-References

- [`cfdrs` Atlas Part](../../../CFDrs/docs/book/migration_overview.md)
  — the CFDrs-mirror version of this Part.
- [`kwavers` Atlas Part](../../../kwavers/docs/book/migration_overview.md)
  — Atlas Part in kwavers (Part VI).
