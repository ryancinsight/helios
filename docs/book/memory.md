# Chapter 4 — Memory and Allocation: Mnemosyne Integration

Helios uses mnemosyne for arena-based, zero-fragmentation allocation
of large physics arrays. Large intermediate buffers (sinograms, dose grids,
terma volumes) are stack-allocated from thread-local arenas.

## Arena Allocation

`
ust
use mnemosyne::Arena;

let arena = Arena::with_capacity(256 * 1024 * 1024); // 256 MiB
let buffer: &mut [f64] = arena.alloc_slice(64 * 64 * 64)?;
`

## Layout Policy

| Data | Layout | Rationale |
|---|---|---|
| Volumetric arrays | C-contiguous (row-major) | Cache-friendly 3D iteration |
| Sinogram | Row per angle | Independent-angle parallelism |
| Dose grid | C-contiguous | Same as CT for subtraction |

## Zero-Copy Slicing

Volume::as_slice() returns a &[T] borrow from the underlying leto::Array3
without allocation. Kernels operate on borrowed slices, enabling zero-copy
pipelines from CT → μ → terma → dose.

## Further Reading

- [Scalar Fields and Numeric Abstractions](numerics.md)
- [mnemosyne crate](https://github.com/ryancinsight/Mnemosyne)
