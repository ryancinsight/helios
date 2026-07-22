# Chapter 35 — Coeus: Tensors and Autodiff

Helios uses **Coeus** as the autodiff-aware tensor crate.  Coeus replaces
`burn::Tensor`, `tch`, and `candle-core` with a single first-party
trait frontier (built on `eunomia::RealField`) and an autodiff tape that
integrates with the rest of the Atlas stack.

## Coeus Surface

```rust
pub struct Tensor<T: RealField, B: Backend = Cpu> {
    storage:    Storage<T, B>,
    grad:       Option<TapeArc<T>>,
    requires_grad: bool,
}

pub trait CoeusScalar: RealField + ... { /* f32 / f64 / f16 / bf16 rules */ }

impl<T: CoeusScalar, B: Backend> Tensor<T, B> {
    pub fn from_storage(storage: NdArray<T, IxN>) -> Self;
    pub fn backward(&self, gradient: Tensor<T, B>);
    pub fn logits(&self) -> &NdArray<T, IxN>;
}
```

Helios currently exercises `coeus` on three workflows:

1. **DVH-constrained beam-weight optimization** ([dvh_optimization](examples/dvh_optimization.md))
   — gradients through `Tensor<f64, Cpu>::logits` flow back to beam
   weights.
2. **MVCT registration** ([mvct_registration](examples/mvct_registration.md))
   — mutual-information loss with `Tensor` gradient tape.
3. **`tomotherapy_workflow`** ([tomotherapy_workflow](examples/tomotherapy_workflow.md))
   — full treatment pipeline finished with a `Tensor`-tracked loss
   backprop through dose accumulation.

## Migration Procedure

| Legacy | Atlas |
|---|---|
| `burn::Tensor::from_data(...)` | `coeus::Tensor::from_storage(NdArray)` |
| `burn::Tensor::backward()` | `coeus::Tensor::backward(grad)` |
| `tch::Tensor::new()` | `coeus::Tensor::from_storage(...)` |
| `candle_core::Tensor` | `coeus::Tensor` |
| `Tensor::requires_grad_(true)` | `Tensor::with_grad_tape()` |
| ad-hoc double-floating forward | `forward(&[CoeusScalar])` |

## Trait Frontier

Coeus depends on `eunomia::RealField`, just like every Atlas crate.  This
is the **single trait frontier**: a `Tensor<f64, Cpu>` and a `Tensor<f32,
Wgpu>` both implement `CoeusScalar::from_f64`, with monomorphized kernels
and no virtual dispatch.

The autodiff tape is **monomorphized** per `(Scalar, Backend)` pair — there
is no separate dtype-specific kernel, and the backward graph is a
`TapeArc<T>` carrying the **same** scalar type as the forward pass.

## GPU Pipeline

```rust
let dose_gpu: Tensor<f64, Wgpu> = Tensor::from_storage(dose_ndarr).to_wgpu();
let loss = dose_gpu.sum();
loss.backward();  // gradients flow back from GPU to host
```

Coeus's GPU kernel suite routes through `hephaestus::GpuCsrMatrix` for
sparse operations and `hephaestus::GpuKernel::x_predefined_convolution`
for the dense path.

## Validation Examples

- [`dvh_optimization`](examples/dvh_optimization.md) — DVH-constrained
  optimization with `Tensor` gradients.
- [`mvct_registration`](examples/mvct_registration.md) — mutual info
  loss with autodiff tape.
- [`tomotherapy_workflow`](examples/tomotherapy_workflow.md) — end-to-end
  autodiff-aware pipeline.
- [`adaptive_rt_workflow`](examples/adaptive_rt_workflow.md) — adaptive
  replanning with per-stop gradient steps.

## Further Reading

- [`coeus` source](../../../coeus/)
- [Leto: Arrays](migration_arrays.md)
- [Moirai: Concurrency](migration_concurrency.md)
