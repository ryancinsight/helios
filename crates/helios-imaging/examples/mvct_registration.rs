//! IGRT Rigid Setup Correction via Translation Registration
//!
//! Demonstrates **integer-voxel translation registration** for image-guided
//! radiation therapy (IGRT) using the MVCT sinogram and reconstruction pipeline.
//!
//! In clinical IGRT a daily MVCT is acquired before delivery, aligned to the
//! planning CT, and the resulting couch shift corrects the patient setup error.
//! This example simulates that workflow on a synthetic 32×32×32 phantom with a
//! known applied shift.
//!
//! ## Stages
//!
//! 1. Build a reference (planning) phantom with bone inserts
//! 2. Create a "daily" image by shifting it 3 voxels along each axis
//! 3. Run `register_translation` to recover the shift
//! 4. Assert the recovered shift matches exactly
//! 5. Apply the correction and measure residual alignment error
//!
//! ## Algorithm
//!
//! Exhaustive whole-voxel search minimizing:
//! ```text
//! cost(s) = mean_v ( moving(v) − fixed(v − s) )²
//! ```
//! over shifts `s ∈ [−max_shift, max_shift]³`.
//!
//! ## Run
//!
//! ```bash
//! cargo run -p helios-imaging --example mvct_registration
//! ```
//!
//! ## Book Chapter
//!
//! [← MVCT and Correction Workflows](../../docs/book/imaging_mvct.md)

use helios_domain::{Volume, VoxelGrid};
use helios_imaging::register_translation;
use helios_math::Point3;

fn main() {
    println!("=== IGRT Rigid Setup Correction: Translation Registration ===\n");

    // ── 1. Planning reference phantom ─────────────────────────────────────────
    //
    // 32×32×32 voxels, 3 mm isotropic. A 10×10×10 block of bone (μ = 0.04 cm⁻¹)
    // centred in water background (μ = 0.02 cm⁻¹).
    let n = 32usize;
    let voxel_mm = 3.0_f64;
    let grid = VoxelGrid::axis_aligned([n, n, n], [voxel_mm; 3], Point3::new(0.0, 0.0, 0.0))
        .expect("valid planning grid");

    let bone_half = 5usize;
    let cx = n / 2;
    let cy = n / 2;
    let cz = n / 2;

    let reference = Volume::from_shape_fn(grid, |[i, j, k]| {
        let in_bone = i >= cx.saturating_sub(bone_half)
            && i < cx + bone_half
            && j >= cy.saturating_sub(bone_half)
            && j < cy + bone_half
            && k >= cz.saturating_sub(bone_half)
            && k < cz + bone_half;
        if in_bone { 0.04_f64 } else { 0.02_f64 }
    });

    println!("Planning CT: {n}×{n}×{n} voxels, {voxel_mm:.0} mm isotropic");
    println!("  Bone μ = 0.04 cm⁻¹ in {bone_half}×{bone_half}×{bone_half} central block\n");

    // ── 2. Simulate daily MVCT with applied setup error ────────────────────────
    //
    // Shift: [+3, −2, +1] voxels — the "patient was mispositioned".
    let true_shift: [isize; 3] = [3, -2, 1];

    // Build the "daily" volume by applying the inverse shift to sampling:
    // daily(v) = reference(v − shift)  ⟹  shifted by +shift in voxel space.
    let daily = Volume::from_shape_fn(grid, |[i, j, k]| {
        let si = i as isize - true_shift[0];
        let sj = j as isize - true_shift[1];
        let sk = k as isize - true_shift[2];
        if si < 0 || sj < 0 || sk < 0 {
            return 0.0_f64;
        }
        reference
            .get(si as usize, sj as usize, sk as usize)
            .unwrap_or(0.0)
    });

    println!("Daily MVCT simulated with setup error: [{}, {}, {}] voxels",
        true_shift[0], true_shift[1], true_shift[2]);

    // ── 3. Registration ────────────────────────────────────────────────────────
    //
    // Search up to ±5 voxels per axis (covers clinical ≤15 mm at 3 mm resolution).
    let max_shift = [5usize; 3];
    let detected_shift = register_translation(&reference, &daily, max_shift);

    println!("\nRegistration result:");
    println!("  True shift:     [{}, {}, {}] voxels",
        true_shift[0], true_shift[1], true_shift[2]);
    println!("  Detected shift: [{}, {}, {}] voxels",
        detected_shift[0], detected_shift[1], detected_shift[2]);

    // ── 4. Validate exact recovery ────────────────────────────────────────────
    for axis in 0..3 {
        assert_eq!(
            detected_shift[axis], true_shift[axis],
            "Axis {axis}: detected {} ≠ true {}",
            detected_shift[axis], true_shift[axis]
        );
    }
    println!("  Exact recovery: ✓");

    // ── 5. Couch correction and residual alignment ─────────────────────────────
    //
    // Correction: corrected(v) = daily(v + detected_shift).
    // Derivation: register_translation finds s s.t. daily(v) ≈ reference(v − s).
    // So reference(v) = daily(v + s). Applying +detected_shift to the lookup
    // index in daily recovers the reference at each voxel.
    let corrected = Volume::from_shape_fn(grid, |[i, j, k]| {
        let ci = i as isize + detected_shift[0];
        let cj = j as isize + detected_shift[1];
        let ck = k as isize + detected_shift[2];
        if ci < 0 || cj < 0 || ck < 0 {
            return 0.0_f64;
        }
        daily
            .get(ci as usize, cj as usize, ck as usize)
            .unwrap_or(0.0)
    });

    // Residual: corrected vs reference (should be 0 in the interior)
    let residual = alignment_rmse(&reference, &corrected, n);
    println!("\nPost-correction alignment:");
    println!("  Residual RMSE (interior) = {residual:.2e} cm⁻¹");
    println!("  (Expected ~0 for exact integer shift in noise-free phantom)");

    let couch_shift_mm: Vec<f64> = detected_shift
        .iter()
        .map(|&s| s as f64 * voxel_mm)
        .collect();
    println!("\nCouch correction (mm): [{:.1}, {:.1}, {:.1}]",
        couch_shift_mm[0], couch_shift_mm[1], couch_shift_mm[2]);

    assert!(
        residual < 1e-6,
        "Interior residual {residual:.2e} should be ~0 after exact correction"
    );

    println!("\nAll IGRT alignment checks passed ✓");
    println!("\nBook chapter: Part II — CT Imaging and Attenuation (MVCT Workflows)");
    println!("API: helios_imaging::register_translation");
}

/// RMSE over the interior voxels (avoiding border affected by the shift padding).
fn alignment_rmse(reference: &Volume<f64>, corrected: &Volume<f64>, n: usize) -> f64 {
    let border = 5usize; // skip voxels near the boundary where shift causes zero-padding
    let mut sum_sq = 0.0;
    let mut count = 0usize;
    for i in border..(n - border) {
        for j in border..(n - border) {
            for k in border..(n - border) {
                let r = reference.get(i, j, k).unwrap_or(0.0);
                let c = corrected.get(i, j, k).unwrap_or(0.0);
                let d = c - r;
                sum_sq += d * d;
                count += 1;
            }
        }
    }
    if count == 0 { return 0.0; }
    (sum_sq / count as f64).sqrt()
}
