# Scatter-superposition kernel optimization — baseline comparison (H-048)

**Date:** 2026-07-02 · **Benchmark:** `helios-solver/benches/scatter_superposition.rs`
(criterion `--quick`, median of the reported interval) · **Machine:** Intel Core
Ultra 9 285K (24c), single-threaded kernel.

## Change

`scatter::convolve_axis` — the dose engine's per-voxel hot path (`O(N·taps)` per
axis, three passes) — was rewritten from per-voxel `Volume::get(i,j,k).expect(…)`
access inside a `from_shape_fn` closure (three-axis bounds check + flat-index
recomputation + closure indirection per tap) to direct iteration of the volume's
zero-copy `as_slice()` view with a precomputed axis stride (one strided slice read
per tap). Tap summation order is unchanged, so results are **bitwise identical**
(all 35 solver oracles pass unchanged — the differential guarantee).

## Measured (5-tap kernel per axis, f64)

| Volume | Before (median) | After (median) | Speedup |
|-------:|----------------:|---------------:|--------:|
| 32³ | 4.31 ms · 7.61 Melem/s | 0.52 ms · 62.9 Melem/s | **8.3×** |
| 64³ | 37.41 ms · 7.01 Melem/s | 5.02 ms · 52.2 Melem/s | **7.4×** |

Criterion's own change estimate at 32³: −87.9 % time / +725 % throughput.

## Analysis

The kernel is bound by per-element index arithmetic and bounds checking, not
floating-point arithmetic (5 mul-adds/element vs ~20 checked-index operations in
the old form). Removing the redundant checks and the per-element closure lets the
inner `k` pass run over contiguous memory. The remaining single slice bounds check
per tap is retained (`#![forbid(unsafe_code)]`); further gains would come from
windowed iterators on the interior region, not `unsafe`.

Benchmark body unchanged before/after (measurement-instrument discipline); the
baseline was captured on the prior kernel in the same session, same machine.
