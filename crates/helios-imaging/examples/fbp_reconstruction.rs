//! Helios FBP (filtered back-projection) reconstruction example.
//!
//! Builds a disk phantom, generates a sinogram via the parallel-beam Radon
//! transform, and reconstructs with the Ram-Lak ramp-filtered back-projection.
//! Validates that the reconstructed attenuation at the disk centre matches μ₀
//! to within 15% — the typical FBP discretization tolerance.
//!
//! Run with: cargo run --example fbp_reconstruction -p helios-imaging

use helios_domain::{Volume, VoxelGrid};
use helios_imaging::{filtered_back_projection, parallel_beam_radon};
use helios_math::Point3;

// ── Phantom helpers ───────────────────────────────────────────────────────────

fn disk_phantom(mu0: f64, radius_mm: f64) -> Volume<f64> {
    let n = 161;
    let spacing = 0.5;
    let origin = Point3::new(0.0, 0.0, 0.0);
    let grid =
        VoxelGrid::axis_aligned([n, n, 1], [spacing; 3], origin).expect("valid phantom grid");
    let centre = (n as f64 - 1.0) * spacing / 2.0;
    Volume::from_shape_fn(grid, move |[i, j, _k]| {
        let dx = i as f64 * spacing - centre;
        let dy = j as f64 * spacing - centre;
        if (dx * dx + dy * dy).sqrt() <= radius_mm {
            mu0
        } else {
            0.0_f64
        }
    })
}

fn uniform_angles(n: usize) -> Vec<f64> {
    (0..n)
        .map(|a| a as f64 * std::f64::consts::PI / n as f64)
        .collect()
}

fn uniform_offsets(half_mm: f64, n: usize) -> Vec<f64> {
    let ds = 2.0 * half_mm / (n - 1) as f64;
    (0..n).map(|j| -half_mm + j as f64 * ds).collect()
}

// ── Main ──────────────────────────────────────────────────────────────────────

fn main() {
    const MU0: f64 = 0.04; // cm⁻¹ (soft tissue at 6 MV)
    const RADIUS_MM: f64 = 25.0;
    const N_ANGLES: usize = 180;
    const N_OFFSETS: usize = 181;

    println!("Building disk phantom  μ₀={MU0} cm⁻¹  R={RADIUS_MM} mm");
    let phantom = disk_phantom(MU0, RADIUS_MM);

    let angles = uniform_angles(N_ANGLES);
    let offsets = uniform_offsets(45.0, N_OFFSETS);

    println!("Running parallel-beam Radon transform  ({N_ANGLES} angles × {N_OFFSETS} offsets)");
    let sino = parallel_beam_radon(&phantom, &angles, &offsets, 400.0, 0.25);

    // Reconstruction grid: 41×41×1 at 2 mm.
    let recon_grid = VoxelGrid::axis_aligned([41, 41, 1], [2.0; 3], Point3::new(0.0, 0.0, 0.0))
        .expect("valid reconstruction grid");

    println!("Running FBP reconstruction  (41×41×1 at 2 mm)");
    let recon = filtered_back_projection(&sino, &recon_grid);

    // --- centre voxel ---
    let centre_mu = recon.get(20, 20, 0).expect("centre voxel in range");
    let err_pct = ((centre_mu - MU0) / MU0).abs() * 100.0;
    println!("Centre voxel  μ_recon={centre_mu:.4}  μ_truth={MU0}  error={err_pct:.1}%");
    assert!(err_pct < 15.0, "FBP centre error {err_pct:.1}% > 15%");

    // --- background corner ---
    let bg = recon.get(2, 2, 0).expect("corner voxel in range");
    println!("Background corner  μ={bg:.4}  (expect ≈0)");
    assert!(
        bg.abs() < 0.1 * MU0,
        "FBP background {bg:.4} should be near zero"
    );

    println!("✓  FBP recovered disk attenuation within {err_pct:.1}% of truth");
}
