//! Repeated cumulative-DVH threshold-query benchmark.
//!
//! The scan arm is the pre-optimization reference: it counts qualifying
//! values directly from the sorted sample. The production arm measures the
//! public query used by plan metrics. Both arms execute the same fixed query
//! workload and retain the accumulated value so Criterion cannot elide it.
#![allow(missing_docs)]

use aequitas::systems::si::quantities::AbsorbedDose;
use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use helios_analysis::Dvh;
use helios_domain::{Volume, VoxelGrid};
use helios_math::Point3;

const SAMPLE_EDGE: usize = 64;
const QUERY_COUNT: usize = 1_024;

fn sample_value([i, j, k]: [usize; 3]) -> f64 {
    (i * 17 + j * 5 + k) as f64 * 0.001
}

fn scan_volume_fraction(sample: &[f64], dose: f64) -> f64 {
    let at_least = sample.iter().filter(|&&value| value >= dose).count();
    at_least as f64 / sample.len() as f64
}

fn workload() -> (Dvh<f64>, Box<[f64]>, Vec<f64>) {
    let grid = VoxelGrid::axis_aligned(
        [SAMPLE_EDGE; 3],
        [1.0, 1.0, 1.0],
        Point3::new(0.0, 0.0, 0.0),
    )
    .expect("benchmark grid is valid");
    let sample = (0..SAMPLE_EDGE)
        .flat_map(|i| {
            (0..SAMPLE_EDGE)
                .flat_map(move |j| (0..SAMPLE_EDGE).map(move |k| sample_value([i, j, k])))
        })
        .collect();
    let volume = Volume::from_shape_fn(grid, sample_value);
    let dvh = Dvh::from_volume(&volume);
    let queries = (0..QUERY_COUNT)
        .map(|index| (index % (SAMPLE_EDGE * 17)) as f64 * 0.001)
        .collect();
    (dvh, sample, queries)
}

fn bench_dvh_queries(c: &mut Criterion) {
    let (dvh, sample, queries) = workload();
    let mut group = c.benchmark_group("dvh_volume_fraction_at_dose");
    group.throughput(Throughput::Elements(QUERY_COUNT as u64));

    group.bench_function(BenchmarkId::new("scan_reference", QUERY_COUNT), |b| {
        b.iter(|| {
            let total = queries.iter().fold(0.0, |acc, &dose| {
                acc + scan_volume_fraction(black_box(&sample), black_box(dose))
            });
            black_box(total)
        });
    });
    group.bench_function(BenchmarkId::new("production", QUERY_COUNT), |b| {
        b.iter(|| {
            let total = queries.iter().fold(0.0, |acc, &dose| {
                acc + black_box(
                    dvh.volume_fraction_at_dose(AbsorbedDose::from_base(black_box(dose))),
                )
            });
            black_box(total)
        });
    });
    group.finish();
}

criterion_group!(benches, bench_dvh_queries);
criterion_main!(benches);
