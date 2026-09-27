//! Throughput benchmark for projected-gradient beam-weight optimization
//! (`optimize_beam_weights`): each iteration is one `A · x` and one `Aᵀ · r`
//! over the dense dose-influence matrix.
//!
//! Shapes are `voxels × beamlets`. `256 × 16` is the small-problem regime where
//! per-iteration allocation is a visible share of the step; the two `4096`-voxel
//! shapes differ only in row stride: 64 `f64` beamlets give 512-byte rows that
//! stay cache-line aligned, 57 give 456-byte rows that do not, so the pair
//! separates a row-alignment effect from everything else.
//!
//! Measurement instrument only: optimization changes the kernel, never this
//! body. Baselines are recorded in the delivering PR.

use criterion::{BenchmarkId, Criterion, Throughput};
use helios_planning::{optimize_beam_weights, DoseInfluence};
use std::hint::black_box;

/// `(voxels, beamlets)` problem shapes.
const SHAPES: &[(usize, usize)] = &[(256, 16), (4096, 57), (4096, 64)];

/// Projected-gradient iterations per measured call.
const ITERATIONS: usize = 10;

fn bench_optimizer(c: &mut Criterion) {
    let mut group = c.benchmark_group("optimize_beam_weights");
    for &(voxels, beamlets) in SHAPES {
        // Non-uniform non-negative entries so no arithmetic short-circuits.
        let rows: Vec<f64> = (0..voxels * beamlets)
            .map(|i| ((i * 37 % 101) as f64 + 1.0) / 101.0)
            .collect();
        let influence = DoseInfluence::from_rows(voxels, beamlets, rows)
            .expect("invariant: rows.len() == voxels * beamlets");
        let prescription = vec![1.0; voxels];
        // ‖AᵀA‖ ≤ ‖A‖_F² < voxels · beamlets for entries in (0, 1], so this step
        // satisfies the convergence bound `step < 2/‖AᵀA‖`.
        let step = 1.0 / (voxels * beamlets) as f64;

        group.throughput(Throughput::Elements(
            (ITERATIONS * voxels * beamlets) as u64,
        ));
        group.bench_function(
            BenchmarkId::from_parameter(format!("{voxels}x{beamlets}")),
            |b| {
                b.iter(|| {
                    optimize_beam_weights(
                        black_box(&influence),
                        black_box(&prescription),
                        black_box(ITERATIONS),
                        black_box(step),
                    )
                });
            },
        );
    }
    group.finish();
}

fn main() {
    let mut criterion = Criterion::default().configure_from_args();
    bench_optimizer(&mut criterion);
    criterion.final_summary();
}
