# GPU vs CPU throughput — Beer–Lambert transmission kernel (H-043)

**Date:** 2026-07-01 · **Benchmark:** `helios-gpu/benches/transmission_throughput.rs`
(criterion, `--quick`) · **Kernel:** `beam_transmission_into` — `out[i] = exp(−τ[i])`.

## Machine class

| Component | Detail |
|-----------|--------|
| CPU | Intel Core Ultra 9 285K — 24 cores / 24 threads |
| GPU | NVIDIA GeForce RTX 5080 (wgpu default adapter) |
| CPU arm | single-threaded serial `f32::exp` loop (the differential reference) |
| GPU arm | `helios_gpu::beam_transmission_into`: host→device upload → `NegOp`+`ExpOp` dispatch → device→host download, **per call** |

## Measured throughput (median, higher = better)

| Elements | CPU (Melem/s) | GPU (Melem/s) | GPU/CPU |
|---------:|--------------:|--------------:|--------:|
| 1 024 | 686 | 4.3 | 0.006× |
| 16 384 | 714 | 54.7 | 0.077× |
| 262 144 | 662 | 323 | 0.49× |
| 1 048 576 | 732 | 524 | 0.72× |
| 4 194 304 | 633 | 373 | 0.59× |

(`--quick` single-estimate medians; indicative, not CI-grade statistics. Full
criterion baselines with confidence intervals are a CI step.)

## Analysis (roofline)

`exp(−τ)` is ~1 arithmetic op per element and reads/writes one `f32` each — it is
**memory-bandwidth-bound**, arithmetic intensity ≈ 0.25 flop/byte. The GPU arm pays a
fixed per-call cost (buffer alloc + PCIe H2D of `4N` bytes + dispatch launch + D2H of
`4N` bytes). For a kernel this cheap, that transfer/launch overhead dominates:

- Small `N` (≤16 k): GPU is 13–160× **slower** — pure launch/transfer overhead.
- Large `N` (1 M): GPU throughput peaks (~524 Melem/s) as the fixed cost amortizes, but
  still **0.72×** the CPU because every call round-trips 8 MB over PCIe for one flop
  each; at 4 M the 16 MB H2D+D2H download cost pulls it back to 0.59×.
- A single-threaded, non-SIMD CPU loop already beats the GPU at every size here.

**Conclusion:** offloading an isolated, transfer-bound elementwise kernel to the GPU is
a net loss on this hardware — the result is correct and expected, not a defect. GPU
throughput for Helios requires **keeping data resident on-device and fusing** the
imaging pipeline (HU→μ map → ray-march/forward projection → `exp(−τ)`) so the τ buffer
is produced and consumed on the GPU without per-op round-trips, amortizing one upload
(CT) and one download (sinogram) over many kernels. Filed as **H-043b** (on-device
fused pipeline); it is the change that would let the GPU arm exceed the CPU.

## Gate status

- Performance gate (GPU scaling/timing): the **scaling study instrument is delivered**
  and produces real numbers. "Competitive with VoLO-class throughput" is **not**
  claimed — no VoLO reference is available in this environment, and the current isolated
  kernel is transfer-bound (H-043b is the path to GPU-favourable throughput). Evidence
  tier: empirical (criterion, this machine); no external reference engine.

The benchmark is a measurement instrument: optimization work changes the *kernel /
pipeline*, never the benchmark body or its timed region.

## Addendum (2026-07-02): fused `ExpNegOp` kernel

`beam_transmission_into` now runs hephaestus's fused `ExpNegOp` (`exp(-x)`, upstreamed
for this path) — **one** dispatch and no intermediate device buffer, replacing the
`NegOp → ExpOp` chain. Re-measured (criterion `--quick`, same machine/instrument):

| Elements | CPU (Melem/s) | GPU chained (07-01) | GPU fused (07-02) |
|---------:|--------------:|--------------------:|------------------:|
| 262 144 | 738 | 323 | 392 |
| 1 048 576 | 732 | 524 | 450 |
| 4 194 304 | 667 | 373 | **485 (+30 %)** |

Conclusion sharpened, not changed: fusing the dispatch chain helps at large sizes but
the isolated kernel remains **PCIe-transfer-bound** — the GPU arm still peaks at
~0.66–0.73× the single-threaded CPU. Removing a dispatch cannot beat the physics of
round-tripping 4 bytes/element each way for ~1 flop. The only remaining path to
GPU > CPU on this workload is the full on-device pipeline (μ-map → projection →
transmission resident on the accelerator; one CT upload, one sinogram download),
which is the remaining H-043b scope.
