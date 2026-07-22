# Chapter 30 — Hermes: SIMD Lanes and Vectorized Kernels

Helios's hot paths — Radon inverse, attenuation map sweeps, FBP back-
projection, sparse kernels — migrate from hand-written SIMD intrinsics
and `packed_simd` to **Hermes**' portable SIMD lanes.  Hermes detects the
host CPU once at runtime (via `is_x86_feature_detected!`) and routes every
kernel through the best lane available.

## Lane Hierarchy

```rust
pub trait SimdLane: Copy + Send + Sync + 'static {
    const LANES: usize;
    type Scalar: RealField;
    fn splat(v: Self::Scalar) -> Self;
    fn add(self, rhs: Self) -> Self;
    fn mul(self, rhs: Self) -> Self;
    ...
}

pub struct SseLane<F>(PhantomData<F>);
pub struct Avx2Lane<F>(PhantomData<F>);
pub struct Avx512Lane<F>(PhantomData<F>);
pub struct NeonLane<F>(PhantomData<F>);
pub struct WasmLane<F>(PhantomData<F>);
```

Helios's CFG-time dispatch (in `helios-physics`'s build script):

```rust
#[cfg(target_arch = "x86_64")]
use hermes_simd::Avx2Lane as ActiveLane;
#[cfg(target_arch = "aarch64")]
use hermes_simd::NeonLane as ActiveLane;
```

The same source compiles to SSE on a Haswell server, AVX-512 on a Sapphire
Rapids node, and NEON on an Apple-Silicon workstation.  No `#ifdef`
blocks inside solver code.

## Migration From packed_simd

| Legacy | Atlas |
|---|---|
| `packed_simd::f64x4` | `hermes_simd::Avx2Lane<f64>` |
| `packed_simd_2::f32x8` | `hermes_simd::Avx2Lane<f32>` |
| manual `#[target_feature(enable = "avx2")]` | `Avx2Lane::add(...)` |
| runtime `is_x86_feature_detected!` | compile-time `#[cfg]` |

The shape of helios SIMD kernels usually shrinks from "instruction-led"
to "shape-led":

```rust
use hermes_simd::SimdLane;

#[inline]
pub fn dose_kernel_sum<F: FloatElement>(terma: &[F], scatter: &[F]) -> F {
    let a = F::Lane::load(terma);
    let s = F::Lane::load(scatter);
    a.mul(s).reduce_add()
}
```

## Kernel Catalogue (Helios)

| Kernel | Hermes spelling |
|---|---|
| FBP ramp filter | `F::Lane::convolve(line, ram_lak_kernel)` |
| Radon sinogram scatter | `F::Lane::stencil(grid, weights)` |
| Sparse SpMV | `F::Lane::scatter_mul(indices, vals, accum)` |
| DVH histogram | `F::Lane::conditional(cond, then_, else_)` |
| Dose-volume scatter-add | `F::Lane::scatter_add(indices, vals, accum)` |

## Validation Examples

- [`fbp_reconstruction`](examples/fbp_reconstruction.md) — FBP ramps
  through `Avx2Lane` kernel.
- [`gpu_detection` (Carrier)](../../CFDrs/docs/book/examples/gpu_detection.md) —
  when GPU is enabled, this routes through `hephaestus` lanes; otherwise
  Hermes handles the same problem on CPU.
- [`parallel_beam_radon`](../../../helios/crates/helios-imaging/examples/radon_sinogram.rs) —
  hermes-vectorized Radon forward.

## Further Reading

- [`hermes-simd` source](../../../hermes/crates/hermes-simd/)
- [Leto: Arrays](migration_arrays.md)
- [Leto: GAT Tiling](migration_gat_tiles.md)
