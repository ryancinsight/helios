# Chapter 28 — Leto: Arrays and Linear Algebra

The bulk of helios's tensor arithmetic — CT volumes, attenuation maps,
dose grids, sparse system matrices — migrates from `ndarray::Array*` and
`nalgebra::DMatrix` to **Leto's** unified `NdArray` and `CsrMatrix` types.
Leto is the **single source of truth** for dense and sparse CPU storage.

## The NdArray Type

```rust
pub struct NdArray<T: RealField, D: Dimension> {
    storage: Vec<T>,
    shape:   D,
}
```

`NdArray<T, D>` is the only dense type Atlas consumers see.  It replaces:

| Legacy | Atlas |
|---|---|
| `ndarray::Array1<f64>` | `leto::NdArray<f64, Ix1>` |
| `ndarray::Array2<f64>` | `leto::NdArray<f64, Ix2>` |
| `ndarray::Array3<f64>` | `leto::NdArray<f64, Ix3>` |
| `ndarray::Array4<f64>` | `leto::NdArray<f64, Ix4>` |
| `nalgebra::DMatrix<f64>` | `leto::DMatrix<f64>` (alias for `NdArray<f64, Ix2>`) |

The `IxN` dimension types are const-generic-shaped where possible.  Storage
is allocated via [`mnemosyne::Arena`].

## CowArray: Zero-Copy Read-Only Views

```rust
pub struct CowArray<'a, T: RealField, D: Dimension> {
    inner: Cow<'a, NdArray<T, D>>,
}
```

`CowArray` returns **borrowed** data when the underlying storage can be
shared and **owned** data when a write requires a copy.  Helios adjoint
reconstruction passes `CowArray` everywhere — reads dominate the hot path,
so the borrow case dominates and there is **zero copy** on the hot path.

## CSR Sparse Storage

Spectral and FE stiffness matrices port to:

```rust
pub struct CsrMatrix<T: RealField> {
    rows:   usize,
    cols:   usize,
    indptr: Vec<usize>,
    indices: Vec<u32>,
    data:   Vec<T>,
}
```

`CsrMatrix<T>` is monomorphised per `T`; `f32` and `f64` builds produce
distinct kernels.  `SpMV` and `SpMM` kernels route through
[`hermes-simd`] where SIMD-width permits.

## Migration Procedure

1. **Identify the boundary.** Every module has a `from_data` API that
   ingests the from-`ndarray` type.  Replace with `&NdArray<F, Ix3>`.
2. **Replace imports.** `use ndarray::Array3;` becomes
   `use leto::NdArray; use eunomia::RealField;`.
3. **Update the access pattern.** `arr[[i,j,k]]` becomes `arr.get([i,j,k])?`
   (returns `Option<&F>`); bulk passes use `arr.as_slice()`.
4. **GPU passthrough.** Anywhere `hephaestus::GpuArray<F,D>` is the desired
   backend, the same code compiles by setting `B = Wgpu`.

## Validation Examples

- [`fbp_reconstruction`](examples/fbp_reconstruction.md) — FBP exchanges
  `NdArray<f64, Ix2>` for sinogram and recon grid.
- [`radon_sinogram`](examples/radon_sinogram.md) — same.
- [`gpu_attenuation_projection`](examples/gpu_attenuation_projection.md) —
  `GpuArray<f64, Ix3>` exercises the same code path.
- [`compton_physics`](examples/compton_physics.md) — mass attenuation
  coefficients over `NdArray`.

## Further Reading

- [`leto` source](../../../leto/crates/)
- [Coeus: Tensors and Autodiff](migration_coeus.md) — autodiff-aware
  tensor siblings.
- [Leto: Geometry](migration_geometry.md)
