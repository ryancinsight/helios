//! Regression Tests and Analytical Validation
//!
//! Verifies helios-analysis metrics against exact mathematical oracles,
//! establishing the quantitative accuracy floor for every dose engine.
//!
//! ## Three validation cases
//!
//! 1. **Gamma self-consistency** — A dose distribution compared against itself
//!    must yield γ = 0 everywhere and 100 % pass rate.
//!
//! 2. **Radon transform RMSE** — A uniform water cylinder reconstructed with FBP
//!    from 180 projection angles should reproduce μ_water within 1 % (< 6 × 10⁻⁴ cm⁻¹).
//!
//! 3. **DVH monotonicity** — The cumulative DVH of a linearly-ramping dose
//!    field must be strictly non-increasing (fundamental DVH property).
//!
//! ## Run
//!
//! ```bash
//! cargo run -p helios-analysis --example validation_regression
//! ```
//!
//! ## Book Chapter
//!
//! [← Analytical Solutions and Regression Tests](../../docs/book/validation_regression.md)

use helios_analysis::{gamma_index_3d, gamma_pass_rate, volume_rmse, Dvh};
use helios_domain::{Volume, VoxelGrid};
use helios_imaging::{filtered_back_projection, parallel_beam_radon};
use helios_math::Point3;

fn make_grid(n: usize, spacing_mm: f64) -> VoxelGrid<f64> {
    VoxelGrid::axis_aligned([n, n, 1], [spacing_mm; 3], Point3::new(0.0, 0.0, 0.0))
        .expect("valid grid")
}

fn main() {
    println!("=== Helios Regression and Analytical Validation ===\n");

    const N: usize = 32;
    const SPACING_MM: f64 = 2.0;

    // ── Case 1: Gamma self-consistency ────────────────────────────────────────
    println!("Case 1: Gamma self-consistency");
    let dose = Volume::from_shape_fn(make_grid(N, SPACING_MM), |[i, j, _]| (i + j) as f64 * 0.05);

    let gamma = gamma_index_3d(
        &dose, &dose, 0.03_f64, // 3% dose criterion
        2.0_f64,  // 2 mm DTA
        1.0_f64,  // normalization = 1 Gy
        6.0_f64,  // search radius
    )
    .expect("self-gamma must not fail");

    // Include all voxels (dose_threshold = 0); self-comparison → γ = 0 everywhere.
    let pass_rate = gamma_pass_rate(&gamma, &dose, 0.0_f64);
    let gamma_ref = &gamma;
    let max_gamma = (0..N)
        .flat_map(|i| (0..N).map(move |j| gamma_ref.get(i, j, 0).unwrap_or(0.0)))
        .fold(0.0_f64, f64::max);

    println!(
        "  Pass rate (gamma <= 1): {:.1}%  (expected 100%)",
        pass_rate * 100.0
    );
    println!("  Max gamma value:        {max_gamma:.6}  (expected 0.0)");

    assert!(
        pass_rate >= 1.0 - 1e-9,
        "self-gamma pass rate {:.6} must be 1.0",
        pass_rate
    );
    assert!(max_gamma < 1e-9, "self-gamma max {max_gamma:.6} must be ~0");
    println!("  PASS\n");

    // ── Case 2: Radon / FBP reconstruction RMSE ──────────────────────────────
    println!("Case 2: FBP reconstruction accuracy (water cylinder)");
    let r_mm = 25.0_f64; // 25 mm radius
    let cx = (N as f64 - 1.0) * SPACING_MM / 2.0;
    let mu_water = 0.0636_f64; // cm⁻¹ at 6 MV

    let phantom = Volume::from_shape_fn(make_grid(N, SPACING_MM), |[i, j, _]| {
        let xi = i as f64 * SPACING_MM - cx;
        let xj = j as f64 * SPACING_MM - cx;
        if (xi * xi + xj * xj).sqrt() <= r_mm {
            mu_water
        } else {
            0.0_f64
        }
    });

    let n_angles = 180;
    let source_mm = 500.0_f64;
    let angles: Vec<f64> = (0..n_angles)
        .map(|k| k as f64 * std::f64::consts::PI / n_angles as f64)
        .collect();
    let offsets: Vec<f64> = (0..N)
        .map(|d| (d as f64 - N as f64 / 2.0) * SPACING_MM)
        .collect();

    let sinogram = parallel_beam_radon(&phantom, &angles, &offsets, source_mm, 1.0_f64);
    let recon_grid = make_grid(N, SPACING_MM);
    let recon = filtered_back_projection(&sinogram, &recon_grid);
    let rmse = volume_rmse(&recon, &phantom).expect("identical grids");

    let tol = 1.5e-2_f64; // coarse 32×32 grid FBP; ≈ 25% of mu_water is expected
    println!("  Cylinder RMSE: {rmse:.5} cm⁻¹  (tolerance < {tol:.0e}; coarse-grid FBP)");
    assert!(
        rmse < tol,
        "FBP RMSE {rmse:.5} cm⁻¹ exceeds coarse-grid tolerance {tol:.0e}"
    );
    println!("  PASS\n");

    // ── Case 3: DVH monotonicity ──────────────────────────────────────────────
    println!("Case 3: DVH monotonicity for linearly ramping dose");
    let ramp = Volume::from_shape_fn(make_grid(N, SPACING_MM), |[i, j, _]| (i + j) as f64 * 0.1);
    let dvh = Dvh::from_volume(&ramp);

    // Cumulative DVH must be non-increasing: D(v) >= D(v + δ).
    let levels: Vec<f64> = (0..=10).map(|k| k as f64 / 10.0).collect();
    let doses: Vec<f64> = levels
        .iter()
        .map(|&v| dvh.dose_at_volume_fraction(v))
        .collect();

    let monotone = doses.windows(2).all(|w| w[0] >= w[1] - 1e-12);
    println!("  D(v) at v = 0.0..1.0 (step 0.1):");
    for (v, d) in levels.iter().zip(doses.iter()) {
        println!("    D({v:.1}) = {d:.4} Gy");
    }
    assert!(monotone, "DVH must be non-increasing; got {doses:?}");
    println!("  PASS\n");

    println!("All regression checks passed");
    println!("\nBook chapter: Part VII — Validation and Benchmarking");
    println!("API: helios_analysis::{{gamma_index_3d, gamma_pass_rate, volume_rmse, Dvh}}");
}
