//! Clinical plan validation: multi-structure evaluation with DVH, biological
//! outcome, and image quality metrics.
//!
//! Demonstrates a head-and-neck treatment plan scenario with three structures:
//!
//! - **PTV** (planning target volume) at ~60 Gy prescription
//! - **Parotid L** (organ at risk) receiving a lower dose
//! - **Spinal cord** (serial organ, strict dose limit)
//!
//! For each structure the example computes:
//!
//! 1. DVH coverage metrics: D₉₅, D_mean, homogeneity index
//! 2. Biological outcome: gEUD, TCP (logistic), NTCP (Lyman-Kutcher-Burman)
//! 3. Image quality assessment: ROI statistics and contrast-to-noise ratio
//!    between PTV and parotid regions
//!
//! Run with: `cargo run --example validation_clinical -p helios-analysis`

use helios_analysis::{
    contrast_to_noise_ratio, gamma_index_3d, gamma_pass_rate, michelson_contrast, roi_statistics,
    Dvh,
};
use helios_domain::{Volume, VoxelGrid};
use helios_math::Point3;

// ── Phantom construction ─────────────────────────────────────────────────────

/// Synthetic dose volume for a head-and-neck plan: a central PTV hotspot with
/// lateral parotid-like regions receiving lower dose and a posterior cord strip.
/// Synthetic dose volume for a head-and-neck plan: a flat-top PTV core at
/// prescription dose with Gaussian penumbral falloff, low-dose parotid and
/// spinal cord regions.
fn hn_phantom(dims: [usize; 3], spacing_mm: f64) -> Volume<f64> {
    let origin = Point3::new(0.0, 0.0, 0.0);
    let grid = VoxelGrid::axis_aligned(dims, [spacing_mm; 3], origin).expect("valid phantom grid");
    let cx = (dims[0] as f64 - 1.0) * spacing_mm / 2.0;
    let cy = (dims[1] as f64 - 1.0) * spacing_mm / 2.0;
    let ptv_radius = 15.0; // mm
    let penumbra = 4.0; // mm penumbral width

    Volume::from_shape_fn(grid, move |[i, j, k]| {
        let _z = k as f64 * spacing_mm;
        let dx = i as f64 * spacing_mm - cx;
        let dy = j as f64 * spacing_mm - cy;
        let r = (dx * dx + dy * dy).sqrt();

        // PTV: flat 60 Gy core inside r ≤ 15 mm, Gaussian falloff in penumbra
        let ptv = if r <= ptv_radius {
            60.0
        } else if r <= ptv_radius + penumbra {
            60.0 * (-(r - ptv_radius).powi(2) / (2.0 * 2.0 * 2.0)).exp()
        } else {
            0.0
        };

        // Parotid L: lateral region at y ≈ 60 mm receiving scattered dose ~12 Gy
        let parotid = if (dy - 20.0).abs() < 6.0 {
            12.0 * (-(dy - 20.0).powi(2) / (2.0 * 4.0 * 4.0)).exp()
        } else {
            0.0
        };

        // Spinal cord: posterior strip at ~8 Gy
        let cord = if (dx - 22.0).abs() < 4.0 {
            8.0 * (-(dx - 22.0).powi(2) / (2.0 * 2.0 * 2.0)).exp()
        } else {
            0.0
        };

        ptv + parotid + cord
    })
}

// ── Main ──────────────────────────────────────────────────────────────────────

fn main() {
    const PRESCRIPTION_GY: f64 = 60.0;
    const DIMS: [usize; 3] = [41, 41, 41];
    const SPACING_MM: f64 = 2.0;

    println!("=== Clinical Plan Validation ===\n");

    let dose = hn_phantom(DIMS, SPACING_MM);
    let cx = (DIMS[0] as f64 - 1.0) * SPACING_MM / 2.0;
    let cy = (DIMS[1] as f64 - 1.0) * SPACING_MM / 2.0;

    // ── Structure masks ───────────────────────────────────────────────────────
    // PTV: central sphere r ≤ 15 mm
    let ptv_mask = |idx: [usize; 3]| -> bool {
        let x = idx[0] as f64 * SPACING_MM - cx;
        let y = idx[1] as f64 * SPACING_MM - cy;
        x * x + y * y <= 15.0 * 15.0
    };
    // Parotid L: lateral region centred at y ≈ 60 mm (20 mm from PTV centre)
    let parotid_mask = |idx: [usize; 3]| -> bool {
        let x = idx[0] as f64 * SPACING_MM;
        let y = idx[1] as f64 * SPACING_MM;
        (56.0..=64.0).contains(&y) && (cx - 5.0..=cx + 5.0).contains(&x)
    };
    // Spinal cord: posterior strip centred at x ≈ 62 mm (22 mm from PTV centre)
    let cord_mask = |idx: [usize; 3]| -> bool {
        let x = idx[0] as f64 * SPACING_MM;
        let y = idx[1] as f64 * SPACING_MM;
        (58.0..=66.0).contains(&x) && (cy - 3.0..=cy + 3.0).contains(&y)
    };

    // ── DVH construction ──────────────────────────────────────────────────────
    let ptv_dvh = Dvh::from_volume_masked(&dose, ptv_mask);
    let parotid_dvh = Dvh::from_volume_masked(&dose, parotid_mask);
    let cord_dvh = Dvh::from_volume_masked(&dose, cord_mask);

    println!("Structure         Voxels   D_mean (Gy)   D₉₅ (Gy)    HI");
    println!("──────────────────────────────────────────────────────────");
    report_dvh("PTV", &ptv_dvh, PRESCRIPTION_GY);
    report_dvh("Parotid L", &parotid_dvh, PRESCRIPTION_GY);
    report_dvh("Spinal cord", &cord_dvh, PRESCRIPTION_GY);

    // ── Clinical coverage assertions ──────────────────────────────────────────
    let ptv_d95 = ptv_dvh.dose_at_volume_fraction(0.95);
    // Flat-top PTV: D₉₅ should be near prescription (60 Gy)
    assert!(
        ptv_d95 >= PRESCRIPTION_GY * 0.90,
        "PTV D₉₅ {ptv_d95:.2} Gy is below 90% of prescription"
    );

    let parotid_mean = parotid_dvh.mean();
    // Parotid mean should be well below prescription
    assert!(
        parotid_mean < PRESCRIPTION_GY * 0.50,
        "Parotid mean dose {parotid_mean:.2} Gy exceeds 50% of prescription"
    );

    let cord_max = cord_dvh.max();
    // Cord dose should be well below the 50 Gy hard limit
    assert!(
        cord_max < 50.0,
        "Spinal cord max dose {cord_max:.2} Gy exceeds 50 Gy hard limit"
    );

    println!("\n✓  DVH coverage passes clinical criteria\n");

    // ── Biological outcome ────────────────────────────────────────────────────
    println!("Biological outcome models (Asclepius laws via DVH dose sample):");

    // gEUD: volume-effect parameter a > 0 emphasizes hot voxels (parallel organs),
    // a < 0 emphasizes cold voxels (serial organs). Positive a avoids 0^negative
    // when masks include zero-dose voxels.
    let ptv_geud = ptv_dvh.generalized_eud(5.0).expect("valid gEUD");
    let parotid_geud = parotid_dvh.generalized_eud(1.0).expect("valid gEUD");
    let cord_geud = cord_dvh.generalized_eud(3.0).expect("valid gEUD");
    println!("  PTV gEUD (a=5):     {ptv_geud:.2} Gy");
    println!("  Parotid gEUD (a=1): {parotid_geud:.2} Gy");
    println!("  Cord gEUD (a=3):    {cord_geud:.2} Gy");

    // TCP: tumour control probability (logistic model)
    let tcp = ptv_dvh.tcp_logistic(5.0, 55.0, 3.0).expect("valid TCP");
    println!("  PTV TCP:            {tcp:.4}  ({:.1}%)", tcp * 100.0);
    assert!(tcp > 0.5, "PTV TCP {tcp:.4} below 50%");

    // NTCP: normal-tissue complication probability (Lyman-Kutcher-Burman)
    let parotid_ntcp = parotid_dvh.ntcp_lkb(1.0, 30.0, 0.15).expect("valid NTCP");
    let cord_ntcp = cord_dvh.ntcp_lkb(3.0, 45.0, 0.12).expect("valid NTCP");
    println!(
        "  Parotid NTCP:       {parotid_ntcp:.4}  ({:.1}%)",
        parotid_ntcp * 100.0
    );
    println!(
        "  Cord NTCP:          {cord_ntcp:.4}  ({:.1}%)",
        cord_ntcp * 100.0
    );
    assert!(
        parotid_ntcp < 0.30,
        "Parotid NTCP {parotid_ntcp:.4} exceeds 30%"
    );
    assert!(cord_ntcp < 0.05, "Cord NTCP {cord_ntcp:.4} exceeds 5%");

    println!("\n✓  Biological outcome within clinical tolerance\n");

    // ── Image quality assessment ──────────────────────────────────────────────
    println!("Image quality metrics (ROI statistics over structure regions):");

    // PTV ROI: central block where dose is flat at ~60 Gy
    let ptv_stats = roi_statistics(&dose, [16, 16, 16], [25, 25, 25]);
    // Parotid ROI: lateral region centred at y ≈ 60 mm
    let parotid_stats = roi_statistics(&dose, [18, 28, 18], [23, 32, 23]);
    // Cord ROI: posterior strip centred at x ≈ 62 mm
    let cord_stats = roi_statistics(&dose, [29, 19, 19], [33, 22, 22]);

    println!(
        "  PTV ROI:        mean={:.2} Gy, σ={:.4} Gy",
        ptv_stats.mean, ptv_stats.std
    );
    println!(
        "  Parotid ROI:    mean={:.2} Gy, σ={:.4} Gy",
        parotid_stats.mean, parotid_stats.std
    );
    println!(
        "  Cord ROI:       mean={:.2} Gy, σ={:.4} Gy",
        cord_stats.mean, cord_stats.std
    );

    // Michelson contrast between PTV and parotid signal levels
    let contrast = michelson_contrast(ptv_stats.mean, parotid_stats.mean);
    println!("  PTV/Parotid Michelson contrast: {contrast:.4}");

    // CNR: PTV vs background (parotid as background, PTV std as noise proxy)
    let cnr = if ptv_stats.std > 0.0 {
        contrast_to_noise_ratio(ptv_stats.mean, parotid_stats.mean, ptv_stats.std)
    } else {
        f64::INFINITY
    };
    println!("  PTV/Parotid CNR:  {cnr:.2}");

    println!("\n✓  Image quality metrics computed\n");

    // ── Gamma index plan verification ─────────────────────────────────────────
    println!("Gamma index verification (3%/2 mm global vs identical plan):");
    let gamma = gamma_index_3d(
        &dose,
        &dose, // self-comparison (identical → 100% pass)
        0.03,
        2.0,
        PRESCRIPTION_GY,
        6.0,
    )
    .expect("self-gamma must succeed");
    let pass = gamma_pass_rate(&gamma, &dose, PRESCRIPTION_GY * 0.10);
    println!("  Self-comparison pass rate: {pass:.1}%");
    assert!(pass >= 0.999, "Self-gamma pass rate {pass:.4} < 99.9%");
    println!("  ✓  Gamma self-consistency verified\n");

    // ── Summary ───────────────────────────────────────────────────────────────
    println!("All clinical validation checks passed");
    println!("API: helios_analysis::{{Dvh, gamma_index_3d, gamma_pass_rate, roi_statistics,");
    println!("       michelson_contrast, contrast_to_noise_ratio}}");
}

fn report_dvh(name: &str, dvh: &Dvh<f64>, _rx: f64) {
    let d95 = dvh.dose_at_volume_fraction(0.95);
    let mean = dvh.mean();
    let hi = dvh.homogeneity_index();
    println!(
        "  {name:<16} {:>6}   {mean:>8.2}     {d95:>7.2}    {hi:.4}",
        dvh.count()
    );
}
