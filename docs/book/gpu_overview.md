# GPU Backend Overview: Hephaestus Integration

Helios GPU acceleration is provided by hephaestus-wgpu (cross-platform
WebGPU/Vulkan/Metal/DX12) or hephaestus-cuda (NVIDIA CUDA).

## Backend Selection

`ust
#[cfg(feature = \"gpu-wgpu\")]
use hephaestus_wgpu::WgpuBackend as GpuBackend;

#[cfg(not(feature = \"gpu-wgpu\"))]
use coeus_core::MoiraiBackend as GpuBackend; // CPU fallback
`

## GPU Tensor Operations

All GPU arrays are hephaestus::Array<T, B, D>, sharing the same API as
leto::Array<T, _, D> for seamless CPU↔GPU switching.

## Sparse GPU Support

Both hephaestus-wgpu and hephaestus-cuda include:
- GpuCsrMatrix<T> — compressed sparse row storage
- spmv — sparse matrix × dense vector
- spmm — sparse matrix × dense matrix

## Further Reading

- [GPU-Accelerated Dose Kernels](gpu_dose.md)
- [Coeus Tensor Operations](gpu_coeus.md)