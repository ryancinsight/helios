//! Full-sinogram forward-projection throughput: CPU ray loop vs the resident
//! GPU projector (H-043b). The workload the residency argument is about — the
//! μ volume uploads once; each batch round-trips only the rays in and one
//! scalar per ray out.
#![allow(missing_docs)] // criterion_group! generates an undocumented harness item.

use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use helios_domain::{Volume, VoxelGrid};
use helios_gpu::GpuProjector;
use helios_math::{Point3, Ray, Vector3};
use helios_solver::forward_project_ray;

const N: usize = 128; // 128³ voxels, 1 mm spacing
const STEP_MM: f32 = 1.0;

fn mu_volume() -> Volume<f32> {
    let grid = VoxelGrid::<f32>::axis_aligned([N, N, N], [1.0; 3], Point3::new(0.0, 0.0, 0.0))
        .expect("grid");
    Volume::from_shape_fn(grid, |idx| {
        0.0002 * idx[0] as f32 + 0.00015 * idx[1] as f32 + 0.0001 * idx[2] as f32 + 0.02
    })
}

/// Parallel-beam sinogram geometry: `n_angles × n_offsets` rays in the axial
/// plane through the volume centre, packed as [ox,oy,oz,dx,dy,dz] per ray.
fn sinogram_rays(n_angles: usize, n_offsets: usize) -> Vec<f32> {
    let centre = (N as f32 - 1.0) / 2.0;
    let mut rays = Vec::with_capacity(n_angles * n_offsets * 6);
    for a in 0..n_angles {
        let theta = a as f32 * std::f32::consts::PI / n_angles as f32;
        let (dx, dy) = (theta.cos(), theta.sin());
        let (px, py) = (-dy, dx); // detector axis
        for o in 0..n_offsets {
            let s = (o as f32 - (n_offsets as f32 - 1.0) / 2.0) * (N as f32 / n_offsets as f32);
            let (cx, cy) = (centre + s * px, centre + s * py);
            rays.extend_from_slice(&[cx - 300.0 * dx, cy - 300.0 * dy, centre, dx, dy, 0.0]);
        }
    }
    rays
}

fn bench_projection(c: &mut Criterion) {
    let mu = mu_volume();
    let device = helios_gpu::default_device().ok();
    if device.is_none() {
        eprintln!("projection bench: no GPU adapter — CPU arm only");
    }
    let projector = device
        .as_ref()
        .map(|d| GpuProjector::new(d, &mu).expect("upload"));

    let mut group = c.benchmark_group("forward_projection_sinogram");
    for &(n_angles, n_offsets) in &[(90usize, 128usize), (360, 256)] {
        let rays = sinogram_rays(n_angles, n_offsets);
        let n_rays = n_angles * n_offsets;
        group.throughput(Throughput::Elements(n_rays as u64));

        group.bench_with_input(
            BenchmarkId::new("cpu", format!("{n_angles}x{n_offsets}")),
            &n_rays,
            |b, _| {
                b.iter(|| {
                    let mut acc = 0.0f32;
                    for r in rays.chunks_exact(6) {
                        let tau = Ray::try_from_direction(
                            Point3::new(r[0], r[1], r[2]),
                            Vector3::new(r[3], r[4], r[5]),
                        )
                        .and_then(|ray| forward_project_ray(black_box(&mu), &ray, STEP_MM))
                        .unwrap_or(0.0);
                        acc += tau;
                    }
                    acc
                });
            },
        );

        if let (Some(device), Some(projector)) = (device.as_ref(), projector.as_ref()) {
            let mut out = vec![0.0f32; n_rays];
            group.bench_with_input(
                BenchmarkId::new("gpu", format!("{n_angles}x{n_offsets}")),
                &n_rays,
                |b, _| {
                    b.iter(|| {
                        projector
                            .project_into(device, black_box(&rays), STEP_MM, &mut out)
                            .expect("gpu projection");
                        out[0]
                    });
                },
            );
        }
    }
    group.finish();
}

criterion_group!(benches, bench_projection);
criterion_main!(benches);
