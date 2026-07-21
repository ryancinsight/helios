//! SIRT Iterative CT Reconstruction
//!
//! Demonstrates the **Simultaneous Iterative Reconstruction Technique** (SIRT)
//! on a synthetic 64×64 attenuation phantom with a bone-like insert.
//!
//! Where filtered back-projection ([`fbp_reconstruction`]) is a one-shot
//! analytical inverse, SIRT is an algebraic reconstruction that iteratively
//! tightens the solution toward the least-squares fit of the sinogram —
//! particularly useful for noisy, limited-angle, or undersampled data.
//!
//! ## Algorithm: Normalized SIRT
//!
//! ```text
//! x ← max(0, x + λ · C⁻¹ ⊙ Aᵀ( R⁻¹ ⊙ (b − Ax) ))
//! ```
//!
//! - `A` — parallel-beam Radon forward projector
//! - `Aᵀ` — back-projector (transpose)
//! - `R⁻¹` — per-ray chord-length normalization
//! - `C⁻¹` — per-voxel hit-weight normalization
//! - `λ = 1.0` — standard relaxation (convergence-stable for 0 < λ < 2)
//! - Non-negativity projection encodes `μ ≥ 0`
//!
//! For consistent (noiseless) data the residual decreases monotonically toward
//! the least-squares solution.
//!
//! ## Stages
//!
//! 1. Build a synthetic water/bone phantom with `Volume::from_shape_fn`
//! 2. Acquire a 32-angle sinogram with `parallel_beam_radon`
//! 3. Reconstruct with `sirt_reconstruction` (10 iterations)
//! 4. Compare RMSE vs the FBP baseline from `filtered_back_projection`
//!
//! ## Run
//!
//! ```bash
//! cargo run -p helios-imaging --example sirt_reconstruction
//! ```
//!
//! ## Book Chapter
//!
//! [← MVCT and Correction Workflows](../../docs/book/imaging_mvct.md)

use helios_domain::{Volume, VoxelGrid};
use helios_imaging::{
    filtered_back_projection, parallel_beam_radon, sirt_reconstruction,
};
use helios_math::Point3;

fn main() {
    println!("=== SIRT Iterative CT Reconstruction ===\n");

    // ── 1. Synthetic phantom ──────────────────────────────────────────────────
    //
    // 64×64×1 slab (single axial slice for speed); 2 mm isotropic voxels.
    let nx = 64usize;
    let ny = 64usize;
    let nz = 1usize;
    let voxel_mm = 2.0_f64;

    let grid = VoxelGrid::axis_aligned(
        [nx, ny, nz],
        [voxel_mm; 3],
        Point3::new(0.0, 0.0, 0.0),
    )
    .expect("valid phantom grid");

    // Water: μ ≈ 0.02 cm⁻¹ (6 MV Compton regime, density ≈ 1.0 g/cm³)
    // Bone insert (20×20 centre block): μ ≈ 0.04 cm⁻¹
    let centre_x = nx / 2;
    let centre_y = ny / 2;
    let bone_half = 10usize;

    let phantom = Volume::from_shape_fn(grid, |[i, j, _k]| {
        let in_bone = i >= centre_x.saturating_sub(bone_half)
            && i < centre_x + bone_half
            && j >= centre_y.saturating_sub(bone_half)
            && j < centre_y + bone_half;
        if in_bone { 0.04_f64 } else { 0.02_f64 }
    });

    let phantom_max = 0.04_f64;
    println!("Phantom: {nx}×{ny}×{nz} voxels, {voxel_mm:.0} mm isotropic");
    println!("  Water μ = 0.02 cm⁻¹, Bone μ = 0.04 cm⁻¹\n");

    // ── 2. Forward projection (sinogram acquisition) ───────────────────────────
    //
    // 32 projection angles uniformly in [0, π), 64 detector bins, fan source at
    // 500 mm SAD (simulates TomoTherapy-scale geometry).
    let n_angles = 32usize;
    let n_detectors = 64usize;
    let source_mm = 500.0_f64;
    let step_mm = 1.0_f64;

    let angles: Vec<f64> = (0..n_angles)
        .map(|k| k as f64 * std::f64::consts::PI / n_angles as f64)
        .collect();
    let offsets: Vec<f64> = (0..n_detectors)
        .map(|d| (d as f64 - n_detectors as f64 / 2.0) * voxel_mm)
        .collect();

    let sinogram = parallel_beam_radon(&phantom, &angles, &offsets, source_mm, step_mm);
    println!(
        "Sinogram acquired: {} angles × {} detectors",
        n_angles, n_detectors
    );

    // ── 3. FBP baseline reconstruction ────────────────────────────────────────
    let fbp = filtered_back_projection(&sinogram, &grid);
    let fbp_rmse = rmse(&phantom, &fbp);
    println!("\nFBP (1 pass, analytic):");
    println!("  RMSE vs phantom = {fbp_rmse:.5} cm⁻¹  ({:.1}% of peak)",
        fbp_rmse / phantom_max * 100.0);

    // ── 4. SIRT iterative reconstruction ──────────────────────────────────────
    let iterations = 10;
    let relaxation = 1.0_f64;

    let sirt = sirt_reconstruction(&sinogram, &grid, source_mm, step_mm, iterations, relaxation);
    let sirt_rmse = rmse(&phantom, &sirt);
    println!("\nSIRT ({iterations} iterations, λ={relaxation}):");
    println!("  RMSE vs phantom = {sirt_rmse:.5} cm⁻¹  ({:.1}% of peak)",
        sirt_rmse / phantom_max * 100.0);

    // ── 5. Convergence comparison ─────────────────────────────────────────────
    println!("\nRMSE comparison:");
    println!("  FBP  : {fbp_rmse:.5} cm⁻¹");
    println!("  SIRT : {sirt_rmse:.5} cm⁻¹");
    if sirt_rmse < fbp_rmse {
        println!("  → SIRT improves on FBP by {:.1}% for this undersampled geometry",
            (fbp_rmse - sirt_rmse) / fbp_rmse * 100.0);
    } else {
        println!("  → FBP is competitive for uniform-angle full sampling");
    }

    // ── 6. Bone contrast check ────────────────────────────────────────────────
    let bone_val = sirt.get(centre_x, centre_y, 0).expect("centre within grid");
    let water_val = sirt.get(0, 0, 0).expect("corner within grid");
    println!("\nSIRT bone/water contrast:");
    println!("  Bone  (centre): {bone_val:.5} cm⁻¹  (phantom = 0.04)");
    println!("  Water (corner): {water_val:.5} cm⁻¹  (phantom = 0.02)");

    // Physics plausibility: reconstructed values must be non-negative and
    // contrast must be positive.
    assert!(
        bone_val >= 0.0,
        "SIRT must produce non-negative attenuation; got {bone_val:.5}"
    );
    assert!(
        bone_val > water_val,
        "Bone must have higher attenuation than water; bone={bone_val:.5} water={water_val:.5}"
    );
    assert!(
        sirt_rmse.is_finite(),
        "SIRT RMSE must be finite"
    );

    println!("\nAll physics checks passed ✓");
    println!("\nBook chapter: Part II — CT Imaging and Attenuation (MVCT Workflows)");
    println!("API: helios_imaging::{{parallel_beam_radon, filtered_back_projection, sirt_reconstruction}}");
}

/// Root-mean-square error between two volumes over all voxels.
fn rmse(reference: &Volume<f64>, reconstruction: &Volume<f64>) -> f64 {
    let [nx, ny, nz] = reference.grid().dims();
    let mut sum_sq = 0.0;
    let mut n = 0usize;
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                if let (Some(r), Some(rec)) = (reference.get(i, j, k), reconstruction.get(i, j, k))
                {
                    let d = rec - r;
                    sum_sq += d * d;
                    n += 1;
                }
            }
        }
    }
    if n == 0 { return 0.0; }
    (sum_sq / n as f64).sqrt()
}
