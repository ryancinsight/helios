//! Adaptive Radiotherapy Workflow with MVCT-Guided Setup Correction
//!
//! Demonstrates the daily IGRT–ART decision loop on a synthetic phantom:
//!
//! 1. **Planning phase** — build a reference CT phantom and compute the
//!    planned 4-field box dose from `helios-simulation`.
//! 2. **Daily MVCT** — simulate a 3-voxel patient setup error, reconstruct
//!    the daily image, and use `register_translation` to recover the shift.
//! 3. **Corrected delivery** — apply the couch-shift correction and recompute
//!    dose on the shifted anatomy.
//! 4. **Adaptive decision gate** — compare the planned vs corrected dose via
//!    DVH metrics and the 3 %/2 mm gamma index; decide whether to adapt.
//!
//! ## Clinical Context
//!
//! In adaptive radiotherapy (ART) the planning dose must be periodically
//! re-evaluated against the daily patient geometry. The gamma index and DVH
//! comparison are the quantitative gates: if gamma pass-rate ≥ 95 % and DVH
//! deviation is small the fraction proceeds; otherwise an online replan is
//! triggered.
//!
//! ## Run
//!
//! ```bash
//! cargo run -p helios-simulation --example adaptive_rt_workflow
//! ```
//!
//! ## Book Chapter
//!
//! [← Adaptive Radiotherapy with MVCT](../../docs/book/workflow_adaptive.md)

use helios_analysis::{gamma_index_3d, gamma_pass_rate, roi_statistics, Dvh};
use helios_domain::{Volume, VoxelGrid};
use helios_imaging::register_translation;
use helios_math::Point3;
use helios_physics::MassAttenuation;
use helios_simulation::{accumulate_delivered_dose, BeamGeometry, DeliveryFrame};
use helios_solver::attenuation_map;

// ── Phantom parameters ────────────────────────────────────────────────────────
const N: usize = 32;
const SPACING_MM: f64 = 2.0;
const MU_WATER_CM: f64 = 0.0636;
const MU_BONE_CM: f64 = 0.24;
const N_LEAVES: usize = 16;
const LEAF_WIDTH_MM: f64 = 2.0;

fn make_grid() -> VoxelGrid<f64> {
    VoxelGrid::axis_aligned([N, N, 1], [SPACING_MM; 3], Point3::new(0.0, 0.0, 0.0))
        .expect("valid grid")
}

/// Build the synthetic CT phantom (HU): water background with a 10×10×1 bone insert.
fn planning_ct() -> Volume<f64> {
    Volume::from_shape_fn(make_grid(), |[i, j, _]| {
        let bone = i >= 11 && i < 21 && j >= 11 && j < 21;
        if bone { 500.0_f64 } else { 0.0_f64 } // bone ≈ 500 HU, water = 0 HU
    })
}

/// Shift the CT phantom by `[si, sj, 0]` voxels (setup error simulation).
fn shift_phantom(ct: &Volume<f64>, si: isize, sj: isize) -> Volume<f64> {
    let [nx, ny, _nz] = ct.grid().dims();
    Volume::from_shape_fn(make_grid(), |[i, j, k]| {
        let src_i = i as isize - si;
        let src_j = j as isize - sj;
        if src_i < 0 || src_j < 0 || src_i >= nx as isize || src_j >= ny as isize {
            return 0.0_f64;
        }
        ct.get(src_i as usize, src_j as usize, k).unwrap_or(0.0)
    })
}

fn build_attenuation_map(ct: &Volume<f64>) -> Volume<f64> {
    attenuation_map(
        ct,
        MassAttenuation::new(MU_WATER_CM).expect("valid μ/ρ"),
        1.0_f64,
    )
}

fn four_field_box(fluence: f64) -> Vec<DeliveryFrame<f64>> {
    [0.0_f64, 90.0, 180.0, 270.0]
        .iter()
        .enumerate()
        .map(|(idx, &deg)| DeliveryFrame {
            projection: idx,
            gantry_angle_rad: deg.to_radians(),
            couch_mm: 0.0,
            leaf_fluence: vec![fluence; N_LEAVES],
        })
        .collect()
}

fn compute_dose(mu: &Volume<f64>) -> Volume<f64> {
    accumulate_delivered_dose(
        &four_field_box(1.0),
        mu,
        BeamGeometry::Parallel { standoff_mm: 500.0 },
        LEAF_WIDTH_MM,
        0.5_f64,
    )
}

fn main() {
    println!("=== Adaptive Radiotherapy Workflow ===\n");

    // ── 1. Planning phase ─────────────────────────────────────────────────────
    println!("Phase 1: Planning");
    let plan_ct = planning_ct();
    let plan_mu = build_attenuation_map(&plan_ct);
    let plan_dose = compute_dose(&plan_mu);
    let plan_dvh = Dvh::from_volume(&plan_dose);
    let d_mean_plan = plan_dvh.mean();
    let d_max_plan = plan_dvh.max();
    println!("  Planned dose:  mean = {d_mean_plan:.4e}  max = {d_max_plan:.4e}");

    // ── 2. Daily MVCT — simulate 3-voxel setup error ──────────────────────────
    println!("\nPhase 2: Daily MVCT (simulated setup error)");
    let applied_shift = [2isize, 1]; // 2 voxels in x, 1 in y
    let daily_ct = shift_phantom(&plan_ct, applied_shift[0], applied_shift[1]);

    // Registration: recover the shift from daily vs reference images
    let recovered = register_translation(
        &plan_ct,
        &daily_ct,
        [4, 4, 0], // search ±4 voxels in x/y, 0 in z (2D phantom)
    );
    println!(
        "  Applied shift:    [{}, {}]",
        applied_shift[0], applied_shift[1]
    );
    println!("  Recovered shift:  [{}, {}]", recovered[0], recovered[1]);

    let registration_ok = recovered[0] == applied_shift[0] && recovered[1] == applied_shift[1];
    if registration_ok {
        println!("  ✓ Registration exact");
    } else {
        println!("  ✗ Registration mismatch — couch shift may be approximate");
    }

    // ── 3. Corrected delivery ─────────────────────────────────────────────────
    println!("\nPhase 3: Corrected delivery (couch shift applied)");
    // Apply correction: shift anatomy back by the recovered displacement
    let corrected_ct = shift_phantom(&daily_ct, -recovered[0], -recovered[1]);
    let corrected_mu = build_attenuation_map(&corrected_ct);
    let corrected_dose = compute_dose(&corrected_mu);
    let corrected_dvh = Dvh::from_volume(&corrected_dose);
    let d_mean_corr = corrected_dvh.mean();
    let d_max_corr = corrected_dvh.max();
    println!("  Corrected dose: mean = {d_mean_corr:.4e}  max = {d_max_corr:.4e}");

    // ── 4. Adaptive decision gate ─────────────────────────────────────────────
    println!("\nPhase 4: Adaptive decision gate");

    // Gamma comparison: corrected vs planned
    let gamma = gamma_index_3d(
        &plan_dose,
        &corrected_dose,
        0.03_f64, // 3% dose criterion
        2.0_f64,  // 2 mm DTA
        d_max_plan,
        6.0_f64,  // search radius
    )
    .expect("identical grids");
    let pass_rate = gamma_pass_rate(&gamma, &plan_dose, 0.0_f64);

    // DVH mean deviation
    let mean_dev = if d_mean_plan > 1e-12 {
        (d_mean_corr - d_mean_plan).abs() / d_mean_plan
    } else {
        0.0
    };

    // Dose uniformity in planning target volume (bone insert region)
    let ptv_stats = roi_statistics(&plan_dose, [11, 11, 0], [21, 21, 1]);
    let ptv_stats_corr = roi_statistics(&corrected_dose, [11, 11, 0], [21, 21, 1]);

    println!("  Gamma pass-rate (3%/2mm):    {:.1}%", pass_rate * 100.0);
    println!("  Mean dose deviation:         {:.2}%", mean_dev * 100.0);
    println!("  PTV dose (planned):  mean = {:.4e}  noise = {:.4e}", ptv_stats.mean, ptv_stats.std);
    println!("  PTV dose (corrected): mean = {:.4e}  noise = {:.4e}", ptv_stats_corr.mean, ptv_stats_corr.std);

    // Clinical acceptance criteria
    let gamma_ok = pass_rate >= 0.95;
    let dvh_ok = mean_dev < 0.05; // < 5% mean dose deviation

    println!("\n  Gamma ≥ 95%:     {}", if gamma_ok { "✓ PASS" } else { "✗ FAIL" });
    println!("  DVH deviation < 5%: {}", if dvh_ok { "✓ PASS" } else { "✗ FAIL" });

    let decision = if gamma_ok && dvh_ok {
        "PROCEED — couch-corrected delivery within clinical tolerance"
    } else {
        "REPLAN — dose deviations exceed ART thresholds"
    };
    println!("\n  Decision: {decision}");

    // Verify core properties
    assert!(plan_dvh.max() > 1e-12, "Planned dose must be positive");
    assert!(corrected_dvh.max() > 1e-12, "Corrected dose must be positive");

    println!("\nAll adaptive RT workflow steps completed successfully ✓");
    println!("\nBook chapter: Part V — Adaptive Radiotherapy with MVCT");
    println!(
        "API: helios_imaging::register_translation + helios_simulation::accumulate_delivered_dose\n     + helios_analysis::{{gamma_index_3d, gamma_pass_rate, roi_statistics, Dvh}}"
    );
}
