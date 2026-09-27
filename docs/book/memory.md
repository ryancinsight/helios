# Chapter 4 — Memory and Allocation: Mnemosyne Integration

<!-- generated-figure-start -->
![Figure 4.1 — Memory and Allocation: Mnemosyne Integration](figures/ch04/fig01_4_memory_and_allocation_mnemosyne_integration.svg)
*Figure 4.1 — Memory and Allocation: Mnemosyne Integration*
<!-- generated-figure-end -->

Helios allocates every dense physics array through the leto array substrate:
`Volume<T>` owns a C-contiguous `leto::Array3<T>` (`crates/helios-domain/src/volume.rs`),
and sinograms, dose grids, and terma volumes are built once and then read
through borrowed slices.

## The Workspace Allocator — `Mnemosyne`

Helios adopts the Atlas memory subsystem as its allocation package. The
allocator is named in exactly one place, `helios_core::memory`, which re-exports
`Mnemosyne` (a `GlobalAlloc` implementation) together with the thread-local
scratch surface (`ScratchPool`, `ScratchBank`, `AlignedVec`, `ScratchElement`).
Every other crate depends on that seam rather than on `mnemosyne` directly, so
the choice can be changed by editing one file. The seam is gated on the
default-on `mnemosyne-memory` feature.

A **program** opts in with one line at its crate root:

```rust,ignore
helios_core::install_global_allocator!();
```

That covers the `xtask` binary, all examples, and the criterion benches. A
**library** never installs it: choosing the allocator is the application's
decision, not a dependency's, and a library that installs one cannot coexist
with a consumer that has its own.

`helios-python` is deliberately excluded. It is a `cdylib` extension module
loaded *into* the CPython process, so a global allocator there would commandeer
the host interpreter's allocator for every Python object allocation — a
side effect on a process Helios does not own. The module links the memory
subsystem but leaves the allocator alone.

### Why a global allocator rather than an allocator-backed container

`leto` exposes a `mnemosyne-alloc` feature that backs `Array` with
`MnemosyneStorage`. Helios does not use it. `MnemosyneStorage` aligns to
`align_of::<T>()` — exactly what `Vec` does — so it buys no cache-line property
that `AlignedVec` does not already provide; and it implements neither `Clone`
nor `Debug`, while `Volume` derives both, so adopting it would remove a public
capability to obtain a property the global allocator already delivers. Installing
the global allocator routes *every* allocation in the program — `Vec`, `Box`, and
the `VecStorage` inside `leto::Array` alike — through Mnemosyne, which is a
strictly larger claim than swapping one container's storage.

## Aligned Storage — `mnemosyne-arena` in the planning kernel

`helios-planning` consumes `mnemosyne-arena`'s `AlignedVec<T>` as the backing
store of the dense dose-influence matrix (`DoseInfluence::data` in
`crates/helios-planning/src/optimize.rs`). The buffer *start* is 64-byte
cache-line aligned; rows are packed without padding, so row `i` begins
`i · beamlets · size_of::<T>()` bytes past it and is itself line-aligned only
when `beamlets · size_of::<T>()` is a multiple of 64. Rows are not padded: the
optimizer's measured gain at 57 beamlets per row (misaligned rows) is at least
its gain at 64 (aligned rows), so row alignment is not its source. It is the
only structure that opts into explicit aligned storage; every other dense array
goes through the leto array substrate described above.

Arena-backed *placement* for the large intermediate buffers is still tracked
work, not current behaviour: `mnemosyne-core` is declared in the workspace
dependency SSOT (`Cargo.toml`) but no member consumes it. The placement
contract Helios does consume is Themis
(`PlacementHint`/`MemoryTier` in `crates/helios-gpu/src/{attenuation,projection,transmission}.rs`),
which hints GPU-resident buffers.

## Layout Policy

| Data | Layout | Rationale |
|---|---|---|
| Volumetric arrays | C-contiguous (row-major) | Cache-friendly 3D iteration |
| Sinogram | Row per angle | Independent-angle parallelism |
| Dose grid | C-contiguous | Same as CT for subtraction |
| Dose-influence matrix | Row-major, unpadded rows, 64-byte aligned buffer start | Row-wise kernels over one contiguous buffer |

## Zero-Copy Slicing

Volume::as_slice() returns a `&[T]` borrow from the underlying leto::Array3
without allocation. Kernels operate on borrowed slices, enabling zero-copy
pipelines from CT → μ → terma → dose.

## Measuring

The allocator is part of the measurement configuration, not a neutral backdrop.
Call `helios_core::memory::warm_current_thread()` before opening a timing window:
it flushes thread-local initialisation traffic — options parsing, arena segment
acquisition, per-thread allocator setup — so the window contains the code under
test rather than the allocator's first-touch cost.

## Further Reading

- [Scalar Fields and Numeric Abstractions](numerics.md)
- [mnemosyne crate](https://github.com/ryancinsight/Mnemosyne)
