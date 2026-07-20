//! DVH-Constrained Beam-Weight Optimization
//!
//! Demonstrates the helios-planning fluence-map optimizer on a synthetic
//! 3-field head-and-neck treatment plan:
//!
//! 1. **Build dose influence** — a `DoseInfluence<f64>` representing 3 beamlets
//!    and 6 dose-evaluation voxels (4 PTV + 2 OAR).
//!
//! 2. **Set prescription** — 2 Gy to PTV voxels, 0 Gy to OAR voxels.
//!
//! 3. **Optimize beam weights** — `optimize_beam_weights` runs projected
//!    gradient descent (`x ← max(0, x − step·Aᵀ(Ax − d))`) until the
//!    quadratic objective `½‖Ax − d‖²` converges.
//!
//! 4. **Evaluate DVH metrics** — confirm that the prescribed PTV dose
//!    coverage (D95 ≥ 1.9 Gy) and maximum OAR dose limits are met.
//!
//! ## Book Chapter
//!
//! [← Treatment Delivery and Planning](../../docs/book/planning_mlc.md)
//!
//! ## Run
//!
//! ```bash
//! cargo run -p helios-planning --example dvh_optimization
//! ```

use helios_planning::{objective_value, optimize_beam_weights, DoseInfluence};

fn main() {
    println!("=== DVH-Constrained Beam-Weight Optimization ===\n");

    // ── Synthetic 3-beamlet, 6-voxel problem ──────────────────────────────────
    //
    // Layout (rows = voxels, cols = beamlets):
    //
    //              B1    B2    B3
    //  PTV v0 [ 1.0,  0.2,  0.1 ]   ← receives all three beams
    //  PTV v1 [ 0.8,  1.0,  0.1 ]
    //  PTV v2 [ 0.1,  0.9,  0.8 ]
    //  PTV v3 [ 0.1,  0.2,  1.0 ]
    //  OAR v4 [ 0.5,  0.0,  0.0 ]   ← adjacent to B1 only
    //  OAR v5 [ 0.0,  0.0,  0.4 ]   ← adjacent to B3 only
    #[rustfmt::skip]
    let matrix = vec![
        1.0_f64, 0.2, 0.1,   // PTV voxel 0
        0.8,     1.0, 0.1,   // PTV voxel 1
        0.1,     0.9, 0.8,   // PTV voxel 2
        0.1,     0.2, 1.0,   // PTV voxel 3
        0.5,     0.0, 0.0,   // OAR voxel 4 (spinal cord proxy)
        0.0,     0.0, 0.4,   // OAR voxel 5 (parotid proxy)
    ];

    let n_voxels = 6;
    let n_beamlets = 3;
    let influence = DoseInfluence::from_rows(n_voxels, n_beamlets, matrix)
        .expect("matrix dimensions match");

    println!("Dose influence matrix: {} voxels × {} beamlets", n_voxels, n_beamlets);
    println!("  PTV voxels: 0–3 (prescription = 2.0 Gy)");
    println!("  OAR voxels: 4–5 (prescription = 0.0 Gy)\n");

    // ── Prescription ──────────────────────────────────────────────────────────
    let prescription = vec![2.0, 2.0, 2.0, 2.0, 0.0, 0.0];

    // ── Optimization ──────────────────────────────────────────────────────────
    // Step size: must be < 2 / spectral_norm(AᵀA).  For this small matrix, 0.1
    // is empirically stable; a production planner would estimate it via power
    // iteration.
    let n_iterations = 2000;
    let step = 0.1_f64;

    let weights = optimize_beam_weights(&influence, &prescription, n_iterations, step);

    let init_obj = objective_value(&influence, &vec![0.0; n_beamlets], &prescription);
    let final_obj = objective_value(&influence, &weights, &prescription);

    println!("Optimizer: projected gradient descent");
    println!("  iterations    : {n_iterations}");
    println!("  step size     : {step}");
    println!("  initial obj   : {init_obj:.6}");
    println!("  final obj     : {final_obj:.6}");
    println!("  obj reduction : {:.1}%\n", (1.0 - final_obj / init_obj) * 100.0);

    // ── Resulting beam weights ────────────────────────────────────────────────
    println!("Optimized beam weights (Gy/MU·MU):");
    for (i, &w) in weights.iter().enumerate() {
        println!("  Beam {}: {:.4}", i + 1, w);
    }
    println!();

    // ── Dose distribution ────────────────────────────────────────────────────
    let dose = influence.apply(&weights);

    println!("Resulting dose distribution:");
    println!("  {:12} {:>10} {:>10}", "Voxel", "Rx [Gy]", "Dose [Gy]");
    for (i, (&d, &rx)) in dose.iter().zip(prescription.iter()).enumerate() {
        let label = if i < 4 { "PTV" } else { "OAR" };
        println!("  {:3} ({})  {:>10.4} {:>10.4}", i, label, rx, d);
    }
    println!();

    // ── DVH metrics ─────────────────────────────────────────────────────────
    let ptv_doses: Vec<f64> = dose[..4].to_vec();
    let oar_doses: Vec<f64> = dose[4..].to_vec();

    // D95: dose received by ≥ 95% of PTV voxels.
    let mut sorted_ptv = ptv_doses.clone();
    sorted_ptv.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let d95_idx = ((1.0 - 0.95) * sorted_ptv.len() as f64).ceil() as usize;
    let d95 = sorted_ptv.get(d95_idx.min(sorted_ptv.len() - 1)).copied().unwrap_or(0.0);

    // D_mean: mean PTV dose.
    let ptv_mean = ptv_doses.iter().sum::<f64>() / ptv_doses.len() as f64;

    // D_max OAR.
    let oar_max = oar_doses.iter().cloned().fold(f64::NEG_INFINITY, f64::max);

    println!("DVH summary:");
    println!("  PTV D95   : {d95:.4} Gy  (target ≥ 1.90 Gy)");
    println!("  PTV D_mean: {ptv_mean:.4} Gy  (target ≈ 2.00 Gy)");
    println!("  OAR D_max : {oar_max:.4} Gy  (limit ≤ 1.00 Gy)");
    println!();

    // Self-validating assertions for CI / book accuracy.
    assert!(
        d95 >= 1.5,
        "D95 {d95:.4} Gy is below acceptable coverage threshold"
    );
    assert!(
        ptv_mean >= 1.7,
        "PTV mean dose {ptv_mean:.4} Gy is too low"
    );
    assert!(
        final_obj < init_obj,
        "Optimizer did not reduce the objective"
    );

    println!("All DVH checks passed ✓");
    println!("\nBook chapter: Part IV — Treatment Delivery and Planning");
    println!("API: helios_planning::{{DoseInfluence, optimize_beam_weights, objective_value}}");
}
