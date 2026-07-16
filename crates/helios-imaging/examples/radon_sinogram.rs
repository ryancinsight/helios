//! Helios parallel-beam Radon transform example.
//!
//! Constructs a uniform-density circular disk phantom, projects it over
//! a set of angles with the parallel-beam Radon transform, and validates
//! the sinogram against the analytical line-integral through a disk.
//!
//! Run with: cargo run --example radon_sinogram -p helios-imaging

use helios_domain::{Volume, VoxelGrid};
use helios_imaging::parallel_beam_radon;
use helios_math::Point3;

// ── Phantom ──────────────────────────────────────────────────────────────────

/// Disk phantom: uniform μ₀ inside `radius_mm`, zero outside.
/// Grid is `n×n×1` with voxel spacing `spacing_mm`.
fn disk_phantom(mu0: f64, radius_mm: f64, n: usize, spacing_mm: f64) -> Volume<f64> {
    let origin = Point3::new(0.0, 0.0, 0.0);
    let grid =
        VoxelGrid::axis_aligned([n, n, 1], [spacing_mm; 3], origin).expect("valid phantom grid");
    let centre_mm = (n as f64 - 1.0) * spacing_mm / 2.0;
    Volume::from_shape_fn(grid, move |[i, j, _k]| {
        let x = i as f64 * spacing_mm - centre_mm;
        let y = j as f64 * spacing_mm - centre_mm;
        if x * x + y * y <= radius_mm * radius_mm {
            mu0
        } else {
            0.0_f64
        }
    })
}

// ── Analytical oracle ─────────────────────────────────────────────────────────

/// Chord length through a disk of `radius` at signed detector offset `s` (mm).
fn chord_mm(radius: f64, s: f64) -> f64 {
    let s2 = s * s;
    let r2 = radius * radius;
    if s2 >= r2 {
        0.0
    } else {
        2.0 * (r2 - s2).sqrt()
    }
}

fn main() {
    const MU0: f64 = 0.2; // cm⁻¹ (water at ~100 keV)
    const RADIUS_MM: f64 = 30.0; // 3 cm radius disk
    const N_VOXELS: usize = 161;
    const SPACING: f64 = 0.5; // 0.5 mm voxels

    let phantom = disk_phantom(MU0, RADIUS_MM, N_VOXELS, SPACING);

    // Projection angles: 4 cardinal angles.
    let angles: Vec<f64> = (0..4)
        .map(|i| i as f64 * std::f64::consts::PI / 4.0)
        .collect();

    // Detector offsets centred on the rotation axis.
    let n_det = 81;
    let det_extent = 50.0_f64; // ±50 mm
    let offsets: Vec<f64> = (0..n_det)
        .map(|i| -det_extent + i as f64 * (2.0 * det_extent) / (n_det - 1) as f64)
        .collect();

    let sinogram = parallel_beam_radon(&phantom, &angles, &offsets, 200.0, SPACING);
    println!(
        "Sinogram dimensions: {} angles × {} offsets",
        angles.len(),
        offsets.len()
    );

    // Validate the θ=0 projection against the analytical chord * μ₀.
    // (Line integrals are in mm·cm⁻¹ because μ is per-cm and spacings are mm.)
    let mut max_err_pct = 0.0_f64;
    for (d, &s) in offsets.iter().enumerate() {
        let analytical_mm = chord_mm(RADIUS_MM, s);
        let analytical = MU0 * analytical_mm / 10.0; // mm → cm
        let measured = sinogram.get(0, d) as f64;
        if analytical > 0.01 {
            let err_pct = ((measured - analytical) / analytical).abs() * 100.0;
            max_err_pct = max_err_pct.max(err_pct);
        }
    }
    println!("θ=0 max relative error vs analytical chord: {max_err_pct:.2}%");
    assert!(
        max_err_pct < 5.0,
        "Radon projection deviates >5% from analytical: {max_err_pct:.2}%"
    );

    println!("✓  Radon sinogram matches analytical chord integral within 5%");
}
