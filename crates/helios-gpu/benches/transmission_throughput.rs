//! GPU-vs-CPU throughput / scaling study for the Beer–Lambert transmission
//! kernel (`exp(−τ)`), the imaging-forward-projection hot path.
//!
//! Measures the same computation on the CPU (`f32::exp` per element) and on the
//! GPU (`helios_gpu::beam_transmission_into`, upload → dispatch → download) across
//! input sizes, reporting elements/second. The crossover — where GPU throughput
//! overtakes CPU once the per-dispatch overhead is amortized — is the scaling
//! result the performance gate tracks. Absolute numbers are machine-specific;
//! the recorded baseline (`validation_reports/`) names the adapter/CPU.
//!
//! This is a measurement instrument: optimization changes the *kernel*, never the
//! benchmark body. If no GPU adapter is present the GPU arm is skipped (CPU-only).
#![allow(missing_docs)] // criterion_group! generates an undocumented harness item.

use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use helios_gpu::{beam_transmission_into, default_device};

/// CPU reference: element-wise `exp(−τ)` (the differential oracle for the kernel).
fn cpu_transmission(tau: &[f32], out: &mut [f32]) {
    for (o, &t) in out.iter_mut().zip(tau) {
        *o = (-t).exp();
    }
}

fn bench_transmission(c: &mut Criterion) {
    // Representative sinogram/projection-buffer sizes (elements).
    const SIZES: [usize; 5] = [1_024, 16_384, 262_144, 1_048_576, 4_194_304];

    let device = default_device().ok();
    if device.is_none() {
        eprintln!("helios-gpu bench: no GPU adapter available — running CPU arm only");
    }

    let mut group = c.benchmark_group("beam_transmission");
    for &n in &SIZES {
        // Optical depths spanning a realistic [0, 1) range.
        let tau: Vec<f32> = (0..n).map(|i| (i % 100) as f32 * 0.01).collect();
        let mut out = vec![0.0f32; n];

        group.throughput(Throughput::Elements(n as u64));
        group.bench_with_input(BenchmarkId::new("cpu", n), &n, |b, _| {
            b.iter(|| cpu_transmission(black_box(&tau), black_box(&mut out)));
        });

        if let Some(device) = device.as_ref() {
            group.bench_with_input(BenchmarkId::new("gpu", n), &n, |b, _| {
                b.iter(|| {
                    beam_transmission_into(device, black_box(&tau), black_box(&mut out))
                        .expect("gpu transmission dispatch");
                });
            });
        }
    }
    group.finish();
}

criterion_group!(benches, bench_transmission);
criterion_main!(benches);
