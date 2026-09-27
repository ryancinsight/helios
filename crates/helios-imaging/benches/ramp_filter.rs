//! Throughput benchmark for the FBP ramp-filter stage (`ramp_filter_rows`) — the
//! reconstruction's only super-linear-cost step.
//!
//! Each detector width is measured twice: the frequency-domain path (the adopted
//! default) and the spatial convolution it replaced. The pair is the measured
//! criterion for adopting the Atlas transform provider (`apollo-fft`) at this
//! member (backlog H-113).
//!
//! Measurement instrument only: optimization changes the kernel, never this
//! body. Baselines are recorded in the corresponding CHANGELOG/commit entry.
#![allow(missing_docs)] // criterion_group! generates an undocumented harness item.

use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use helios_imaging::{ram_lak_kernel, ramp_filter_rows, RampMethod};

/// Detector widths from a small MVCT panel to a large flat panel. `64` is the
/// transform crossover itself, so the constant is measured rather than asserted;
/// the rest are all above it, where `Auto` selects the frequency-domain path.
const WIDTHS: &[usize] = &[64, 181, 512, 1024];

/// Projection count held fixed so the detector width is the only variable.
const ANGLES: usize = 64;

/// Detector spacing in cm, matching the FBP stage's own units.
const DS_CM: f64 = 0.2;

fn bench_ramp_filter(c: &mut Criterion) {
    let mut group = c.benchmark_group("ramp_filter_rows");
    for &n_off in WIDTHS {
        let kernel = ram_lak_kernel(n_off, DS_CM);
        // Non-uniform rows so no arithmetic short-circuits.
        let rows: Vec<f64> = (0..ANGLES)
            .flat_map(|a| (0..n_off).map(move |k| ((a * n_off + k) as f64 * 0.017).sin() + 0.25))
            .collect();

        group.throughput(Throughput::Elements((ANGLES * n_off) as u64));
        for (label, method) in [("fft", RampMethod::Auto), ("direct", RampMethod::Direct)] {
            group.bench_function(BenchmarkId::new(label, n_off), |b| {
                b.iter(|| {
                    let mut buffer = rows.clone();
                    ramp_filter_rows(&mut buffer, ANGLES, n_off, &kernel, DS_CM, method);
                    black_box(buffer)
                });
            });
        }
    }
    group.finish();
}

criterion_group!(benches, bench_ramp_filter);
criterion_main!(benches);
