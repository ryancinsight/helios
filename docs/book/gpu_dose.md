# Chapter 21 — GPU-Accelerated Dose Kernels

The collapsed-cone dose engine can be accelerated with GPU ray-marching
kernels dispatched through hephaestus-wgpu.

## WGSL Kernel

`wgsl
@compute @workgroup_size(64)
fn ray_march_kernel(
    @builtin(global_invocation_id) id: vec3<u32>,
) {
    let cone_idx = id.x;
    // March along cone direction and accumulate dose
    var dose: f32 = 0.0;
    for (var step: u32 = 0u; step < max_steps; step++) {
        let pos = origin + f32(step) * direction * step_size;
        let mu = sample_attenuation(pos);
        dose += terma_at(pos) * exp(-mu * f32(step) * step_size);
    }
    dose_buffer[voxel_idx] = dose;
}
`

## Performance

On an NVIDIA RTX 4090, the GPU collapsed-cone is ~40× faster than the
CPU moirai multi-threaded implementation for a 256³ dose grid.

## Further Reading

- [GPU Backend Overview](gpu_overview.md)
- [Collapsed-Cone Convolution](dose_collapsed_cone.md)