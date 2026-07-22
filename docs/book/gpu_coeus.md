# Chapter 22 — Coeus Tensor Operations for Dose Grids

coeus-core provides the tensor substrate for GPU-resident dose grids.
All operations are expressed through the ComputeBackend trait, making
them portable across CPU (MoiraiBackend) and GPU backends.

## Tensor-Based Dose Operations

`
ust
use coeus_core::ComputeBackend;
use coeus_tensor::Tensor;

fn scale_dose<B: ComputeBackend>(
    dose: &Tensor<f32, B>,
    prescription_gy: f32,
) -> Tensor<f32, B> {
    dose.clone() * prescription_gy
}
`

## autodiff Integration

coeus-autograd enables gradient computation through dose calculation,
supporting inverse planning and machine-learning surrogate training.

## Further Reading

- [GPU Backend Overview](gpu_overview.md)
- [Coeus crate](https://github.com/ryancinsight/Coeus)