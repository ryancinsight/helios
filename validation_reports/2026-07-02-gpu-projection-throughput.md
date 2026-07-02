# Resident GPU forward projection — CPU vs GPU (H-043b resolved)

**Date:** 2026-07-02 · **Benchmark:** `helios-gpu/benches/projection_throughput.rs`
(criterion `--quick`, medians) · **Machine:** Intel Core Ultra 9 285K (24c, CPU arm
single-threaded) + NVIDIA GeForce RTX 5080.

## Change

`GpuProjector` keeps the attenuation volume **resident on the GPU** (one upload at
construction) and forward-projects whole ray batches per dispatch through
hephaestus's new `ray_line_integrals` kernel (upstreamed, commits 792ccc3/9354260:
slab-clip to the node AABB → `n = ceil(L/step)` midpoint trilinear samples → one
thread per ray). PCIe traffic per batch: 24 B/ray in, 4 B/ray out — the residency
step the H-043 study identified as the only path past transfer-bound physics.

## Measured (128³ μ volume, 1 mm step, f32)

| Sinogram | Rays | CPU (median) | GPU (median) | Speedup |
|---------:|-----:|-------------:|-------------:|--------:|
| 90 × 128 | 11 520 | 75.4 ms · 153 Kelem/s | 0.441 ms · 26.1 Melem/s | **171×** |
| 360 × 256 | 92 160 | 589.6 ms · 156 Kelem/s | 1.591 ms · 57.9 Melem/s | **371×** |

Correctness: the batched GPU projector matches the CPU `forward_project_ray`
reference per-ray within a derived 1e-3 relative bound (identical per-ray sequential
summation order; live-adapter differential test in `helios-gpu`), and the kernel's
own analytical oracles (uniform chord, affine-exact midpoint, step-independence,
miss → 0) run in hephaestus CI.

## Analysis

Each ray performs ~220 midpoint trilinear samples (8 loads + ~14 FLOPs each) — a
compute-dense, embarrassingly parallel workload, the opposite regime from the
elementwise `exp(−τ)` kernel (1 FLOP/element, transfer-bound at 0.66–0.73× CPU even
after dispatch fusion). Residency converts the GPU from a net loss into a 2+
order-of-magnitude win on the pipeline's dominant cost.

**Gate status:** the "GPU scaling/throughput" component of the performance gate is
now demonstrated on the pipeline workload with recorded baselines. "Competitive
with VoLO-class throughput" remains unclaimable — no VoLO reference exists in this
environment. CPU arm is single-threaded `forward_project_ray`; a rayon/moirai
multi-core CPU arm would narrow but not close a 371× gap.
