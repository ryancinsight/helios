//! Helios gamma index plan-comparison example.
//!
//! Demonstrates the 3%/2 mm global-normalization gamma index — the standard
//! clinical criterion for comparing a calculated dose distribution against a
//! reference (measurement or Monte Carlo). A pass rate ≥ 95% is the typical
//! clinical acceptance threshold.
//!
//! Constructs a reference Gaussian dose volume and two evaluated plans:
//! - "Good plan": identical to reference (expect 100% pass rate, γ ≈ 0).
//! - "Degraded plan": same Gaussian shifted by 1.5 mm (expect <100% pass rate).
//!
//! Run with: cargo run --example gamma_index -p helios-analysis

use aequitas::systems::si::{
    quantities::{AbsorbedDose, Length},
    units::{Gray, Millimeter},
};
use helios_analysis::{gamma_index_3d, gamma_pass_rate};
use helios_domain::{Volume, VoxelGrid};
use helios_math::Point3;

// ── Helpers ───────────────────────────────────────────────────────────────────

fn gaussian_dose(
    dims: [usize; 3],
    spacing_mm: f64,
    peak: f64,
    sigma_mm: f64,
    shift_mm: [f64; 3],
) -> Volume<f64> {
    let origin = Point3::new(0.0, 0.0, 0.0);
    let grid = VoxelGrid::axis_aligned(dims, [spacing_mm; 3], origin).expect("valid dose grid");
    let centre = [
        (dims[0] as f64 - 1.0) * spacing_mm / 2.0 + shift_mm[0],
        (dims[1] as f64 - 1.0) * spacing_mm / 2.0 + shift_mm[1],
        (dims[2] as f64 - 1.0) * spacing_mm / 2.0 + shift_mm[2],
    ];
    let two_sig_sq = 2.0 * sigma_mm * sigma_mm;
    Volume::from_shape_fn(grid, move |[i, j, k]| {
        let dx = i as f64 * spacing_mm - centre[0];
        let dy = j as f64 * spacing_mm - centre[1];
        let dz = k as f64 * spacing_mm - centre[2];
        peak * (-(dx * dx + dy * dy + dz * dz) / two_sig_sq).exp()
    })
}

// ── Main ──────────────────────────────────────────────────────────────────────

fn main() {
    const DIMS: [usize; 3] = [21, 21, 21];
    const SPACING_MM: f64 = 2.0; // 2 mm voxels
    const PEAK_GY: f64 = 60.0;
    const SIGMA_MM: f64 = 10.0;
    const NORM_DOSE: f64 = PEAK_GY; // global normalization = max dose
    const DOSE_THRESHOLD: f64 = PEAK_GY * 0.10; // 10% low-dose cutoff
                                                // 3%/2 mm clinical criterion
    const DOSE_DIFF: f64 = 0.03;
    const DTA_MM: f64 = 2.0;
    const SEARCH_MM: f64 = 6.0; // 3× DTA
    let normalization_dose = AbsorbedDose::from_unit::<Gray>(NORM_DOSE);
    let dose_threshold = AbsorbedDose::from_unit::<Gray>(DOSE_THRESHOLD);
    let dta = Length::from_unit::<Millimeter>(DTA_MM);
    let search_radius = Length::from_unit::<Millimeter>(SEARCH_MM);

    // Reference plan: centred Gaussian
    let reference = gaussian_dose(DIMS, SPACING_MM, PEAK_GY, SIGMA_MM, [0.0; 3]);

    // --- Case 1: identical plan (γ should be ≈ 0, pass rate = 100%) ---
    let identical = gaussian_dose(DIMS, SPACING_MM, PEAK_GY, SIGMA_MM, [0.0; 3]);
    let gamma_identical = gamma_index_3d(
        &reference,
        &identical,
        DOSE_DIFF,
        dta,
        normalization_dose,
        search_radius,
    )
    .expect("gamma computation succeeded");
    let pass_identical = gamma_pass_rate(&gamma_identical, &reference, dose_threshold);
    println!(
        "Identical plan:  pass rate = {:.1}%  (expect 100%)",
        pass_identical * 100.0
    );
    assert!(
        pass_identical >= 0.999,
        "Identical plans must achieve 100% pass rate"
    );

    // --- Case 2: shifted plan (1.5 mm shift, within 2 mm DTA — expect high pass) ---
    let shifted_small = gaussian_dose(DIMS, SPACING_MM, PEAK_GY, SIGMA_MM, [1.5, 0.0, 0.0]);
    let gamma_shifted_small = gamma_index_3d(
        &reference,
        &shifted_small,
        DOSE_DIFF,
        dta,
        normalization_dose,
        search_radius,
    )
    .expect("gamma computation succeeded");
    let pass_shifted_small = gamma_pass_rate(&gamma_shifted_small, &reference, dose_threshold);
    println!(
        "1.5 mm shift:    pass rate = {:.1}%  (expect high)",
        pass_shifted_small * 100.0
    );

    // --- Case 3: large shift (4 mm — exceeds DTA, expect lower pass rate) ---
    let shifted_large = gaussian_dose(DIMS, SPACING_MM, PEAK_GY, SIGMA_MM, [4.0, 0.0, 0.0]);
    let gamma_shifted_large = gamma_index_3d(
        &reference,
        &shifted_large,
        DOSE_DIFF,
        dta,
        normalization_dose,
        search_radius,
    )
    .expect("gamma computation succeeded");
    let pass_shifted_large = gamma_pass_rate(&gamma_shifted_large, &reference, dose_threshold);
    println!(
        "4.0 mm shift:    pass rate = {:.1}%  (expect <100%)",
        pass_shifted_large * 100.0
    );
    assert!(
        pass_shifted_large < pass_shifted_small,
        "Larger shift must produce worse pass rate"
    );

    println!(
        "\nGamma criterion: {}%/{} mm (global normalization = {NORM_DOSE} Gy)",
        (DOSE_DIFF * 100.0) as usize,
        DTA_MM
    );
    println!("✓  Gamma index correctly ranks plan quality: identical > small shift > large shift");
}
