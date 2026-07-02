//! Throughput benchmark for the collapsed-cone scatter stage
//! (`scatter_superposition`), the dose engine's per-voxel hot path
//! (`O(N · taps)` per axis pass, three passes).
//!
//! Measurement instrument only: optimization changes the kernel, never this
//! body. Baselines are recorded in the corresponding CHANGELOG/commit entry.
#![allow(missing_docs)] // criterion_group! generates an undocumented harness item.

use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use helios_domain::{Volume, VoxelGrid};
use helios_math::Point3;
use helios_solver::{scatter_superposition, symmetric_deposition_kernel};

fn bench_scatter(c: &mut Criterion) {
    let mut group = c.benchmark_group("scatter_superposition");
    for &n in &[32usize, 64] {
        let grid = VoxelGrid::axis_aligned([n, n, n], [2.0, 2.0, 2.0], Point3::new(0.0, 0.0, 0.0))
            .expect("grid");
        // Non-uniform terma so no arithmetic short-circuits.
        let terma =
            Volume::from_shape_fn(grid, |idx| (idx[0] * 7 + idx[1] * 3 + idx[2]) as f64 * 1e-4);
        let kernel = symmetric_deposition_kernel(0.6_f64, 0.2, 2); // 5 taps/axis

        group.throughput(Throughput::Elements((n * n * n) as u64));
        group.bench_with_input(BenchmarkId::from_parameter(n), &n, |b, _| {
            b.iter(|| scatter_superposition(black_box(&terma), &kernel, &kernel, &kernel));
        });
    }
    group.finish();
}

criterion_group!(benches, bench_scatter);
criterion_main!(benches);
