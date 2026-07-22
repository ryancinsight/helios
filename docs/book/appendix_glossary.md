# Atlas Glossary

A consolidated glossary spanning every Atlas crate that Helios, CFDrs, and
Kwavers consume.  Where a term has an Atlas-crate-specific meaning, the
relevant crate is named in parentheses.  See also
[`CFDrs` glossary](../../CFDrs/docs/book/appendix_glossary.md) and
[`Kwavers` glossary](../../kwavers/docs/book/appendix_glossary.md) for
the cross-book shared vocabulary.

## A

**Alloc, Arena Allocation** (`mnemosyne`).  A type-erased memory region
that owns its chunks and recycles them in O(1) on `reset()`.  Helios uses
arenas for dose-volume storage and per-patient transient scratch.

**Apollo** (`apollo`).  Atlas forward-only FFT crate.  Apollo exposes
`FftPlan::forward_real` and `forward_complex`.  Inverse is done via
forward + complex-conjugate, with autodiff-tape preservation.

**Atlas Stack**.  Unified first-party mathematics, memory, SIMD, and
concurrency stack — `eunomia` / `leto` / `hephaestus` / `coeus` /
`apollo` / `hermes-simd` / `mnemosyne` / `themis` / `moirai` / `ritk`
/ `consus`.

**Autograd / Autodiff** (`coeus`).  Reverse-mode automatic
differentiation.  Coeus tracks an autodiff tape so that any
`Tensor<T>` participates in a backward graph.

## B

**Backend (Compute)** (`hephaestus`).  A `hephaestus::Backend` is one
of `Cpu`, `Wgpu`, or `Cuda`.  Atlas consumers carry a backend as a type
parameter so the same source compiles for any backend.

**BeamIsometry** (`leto`).  Helios-specific typed composition of
gantry/collimator/couch `Isometry3`.  Replaces legacy matrix-product
chains.

## C

**Coeus** (`coeus`).  Atlas tensor + autodiff crate.  Equivalent to
PyTorch/JAX/Burn from the Atlas perspective.

**CowArray** (`leto`).  A `Cow<'a, NdArray<T, D>>` used for
read-then-write slicing without cloning the underlying storage.

**Compile-time dispatch**.  Atlas uses `#[cfg(target_arch)]` +
`PhantomData` so the same source compiles to specialized code per
(`Scalar`, `Backend`) pair, no virtual dispatch.

**ComputeBackend**.  A `eunomia::Backend` of `Cpu`/`Wgpu`/`Cuda`.

## D

**DICOM** (`ritk`).  Digital Imaging and Communications in Medicine.
Ritk parses DICOM into typed `DicomElement`s.

**DIP / DIF**.  Dependency-Inversion Principle: depend on traits
(`eunomia::RealField`, `moirai::Executor`), not impls.

**DRY**.  Don't Repeat Yourself.  Atlas enforces DRY at the trait
frontier.

## E

**Eunomia** (`eunomia`).  Atlas numeric-trait crate that unifies
`num-traits` + `num-complex` + custom precision types.

**Executor** (`moirai`).  Unified async + parallel worker pool
replacing `tokio::Runtime` and `rayon::ThreadPool`.

## F

**Foundations**.  The lowest chapter: type definitions, validation
contracts, dimension checking.

**FftPlan** (`apollo`).  Pre-computed FFT plan for a fixed shape.  Used
by every Atlas consumer running spectral methods.

**FloatElement** (`eunomia`).  Atlas trait bound for floating-point
scalars that want Atlas-wide SIMD / autodiff / MPI integration.

## H

**Hephaestus** (`hephaestus`).  Atlas GPU crate — wgpu (cross-platform)
and CUDA backends.  Exposes `WgpuBackend`, `CudaBackend`.

**Hermes** (`hermes-simd`).  Atlas SIMD crate with portable SSE2 /
AVX2 / AVX-512 / NEON / WASM lane dispatch.

**Host–Device Sync**.  Boundary between `leto` (CPU) and `hephaestus`
(GPU).  Mediated by `hephaestus::sync_*`.

## I

**Isometry3** (`leto`).  3-D rigid transformation (rotation +
translation).  Equivalent to `nalgebra::Isometry3`, but uses Leto's
quaternion representation.

## K

**KernelCache** (`hephaestus`).  Pre-compiled shader / CUDA kernels,
keyed by `(shape, dtype)`.  Avoids recompilation per launch.

## L

**Lending Iterator**.  Atlas iterator that yields a value whose
lifetime is **shorter** than `&self`.  Encoded via `type Item<'a>` (RFC 1598).

**Leto** (`leto`).  Atlas CPU dense/sparse storage crate —
`NdArray<T, D>`, `CowArray`, `CsrMatrix<T>`, `DMatrix<T>`, plus the
geometry types `Point3<T>`, `Vector3<T>`, `Isometry3<T>`.

## M

**Mnemosyne** (`mnemosyne`).  Atlas arena allocator crate —
`Arena`, `ScratchArena`, with chunked growth.

**Moirai** (`moirai`).  Atlas async + parallel crate.  Replaces
`tokio::Runtime` + `rayon::ThreadPool`.

**MLC** (`helios-domain`).  Multi-Leaf Collimator typed geometry,
replacing legacy `Vec<(f64, f64)>` extents and runtime Varian/Elekta enums.

## N

**NdArray** (`leto`).  Atlas dense n-dimensional array.  Replaces
`ndarray::Array*` and `nalgebra::DMatrix/Tensor`.

**NUMA Placement** (`themis`).  Per-`Arena` binding to a physical
core or NUMA node.

## P

**PhantomData**.  Zero-sized marker type forcing the compiler to
honour capability requirements.

**PyO3 boundary**.  Helios Rust ↔ Python boundary.

## R

**RealField** (`eunomia`).  Atlas trait bound for floating-point-
real-type.

**Ritk** (`ritk`).  Atlas image-toolkit crate — DICOM, NIfTI, PNG.

## S

**SIMD Lane** (`hermes-simd`).  A `SimdLane`-implementing type
(`Avx2Lane<F>`, `Avx512Lane<F>`, …).

**Spectral Tile** (`apollo`, `leto`).  Frequency-domain partition
struct used by tiled spectral solvers.

**SRP**.  Single Responsibility Principle.

**SSOT**.  Single Source of Truth.

## T

**TaskGraph** (`moirai`).  DAG of dependent compute tasks, executed by
an `Executor`.

**Tensor** (`coeus`).  Autodiff-aware n-dimensional array.

**TG-119**.  AAPM Task Group 119 IMRT commissioning phantom.  Helios's
canonical clinical validation phantom.

**Themis** (`themis`).  Atlas NUMA / physical-core placement crate.

**TileStreaming** (`leto`).  Trait with a GAT-based iterator that
yields `&Tile<'a>` over an Atlas volume.

**Typed Boundary**.  Boundary condition carried as a compile-time enum
(`BoundaryKind::Lid`, `::Wall`, …) rather than runtime string tags.

## V

**VoxelGrid** (`helios-domain`).  Helios typed grid — `shape [u32;3]`,
`spacing_mm [F;3]`, `origin_mm Point3<F>` — replaces ad-hoc `Vec<f64>`
triples.

## Z

**Zero-Copy**.  Atlas pattern — every read-only consumer takes
`Cow<'_, T>` instead of `T`.

**Zero-Cost Abstraction**.  Atlas pattern — every abstraction (ZSTs,
GATs, `PhantomData`, const generics) compiles away to zero cycles.

**ZST**.  Zero-Sized Type.  Atlas uses `PhantomData<F>` and similar
ZSTs to encode capability at compile time.
