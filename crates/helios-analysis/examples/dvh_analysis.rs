//! Helios dose-volume histogram (DVH) analysis example.
//!
//! Constructs a synthetic dose distribution representing a tumour target volume
//! (PTV) receiving ~60 Gy and surrounding critical structures, then builds DVHs
//! and reports clinical evaluation metrics:
//!
//! - `D₉₅` (dose to 95% of the PTV — the coverage metric)
//! - `D_mean` (mean dose)
//! - Homogeneity Index (ICRU-83: `(D₂ − D₉₈) / D₅₀`)
//! - Volume fraction above prescription
//!
//! Run with: cargo run --example dvh_analysis -p helios-analysis

use helios_analysis::Dvh;
use helios_domain::{Volume, VoxelGrid};
use helios_math::Point3;

// ── Phantom dose construction ─────────────────────────────────────────────────

/// Synthetic dose volume: a Gaussian beam depositing a prescribed dose at
/// the grid centre, dropping off with 3-D Gaussian rolloff.
fn gaussian_dose_phantom(
    dims: [usize; 3],
    spacing_mm: f64,
    peak_gy: f64,
    sigma_mm: f64,
) -> Volume<f64> {
    let origin = Point3::new(0.0, 0.0, 0.0);
    let grid =
        VoxelGrid::axis_aligned(dims, [spacing_mm; 3], origin).expect("valid phantom grid");
    let centre = [
        (dims[0] as f64 - 1.0) * spacing_mm / 2.0,
        (dims[1] as f64 - 1.0) * spacing_mm / 2.0,
        (dims[2] as f64 - 1.0) * spacing_mm / 2.0,
    ];
    let two_sig_sq = 2.0 * sigma_mm * sigma_mm;
    Volume::from_shape_fn(grid, move |[i, j, k]| {
        let dx = i as f64 * spacing_mm - centre[0];
        let dy = j as f64 * spacing_mm - centre[1];
        let dz = k as f64 * spacing_mm - centre[2];
        peak_gy * (-(dx * dx + dy * dy + dz * dz) / two_sig_sq).exp()
    })
}

// ── Main ──────────────────────────────────────────────────────────────────────

fn main() {
    const PRESCRIPTION_GY: f64 = 60.0; // 60 Gy prescription dose
    const PEAK_GY: f64 = 63.0; // slight hot spot at centre
    const SIGMA_MM: f64 = 25.0; // beam sigma (Gaussian width)

    println!("Building synthetic dose phantom  peak={PEAK_GY} Gy  σ={SIGMA_MM} mm");
    let dose = gaussian_dose_phantom([31, 31, 31], 2.0, PEAK_GY, SIGMA_MM);

    // --- Full-volume DVH ---
    let dvh = Dvh::from_volume(&dose);
    println!("\nFull-volume DVH  ({} voxels)", dvh.count());
    println!("  D_min  = {:.2} Gy", dvh.min());
    println!("  D_mean = {:.2} Gy", dvh.mean());
    println!("  D_max  = {:.2} Gy", dvh.max());

    // Clinical coverage: V_95 and D_95 relative to prescription
    let d95 = dvh.dose_at_volume_fraction(0.95);
    let v95_pct = dvh.volume_fraction_at_dose(PRESCRIPTION_GY * 0.95) * 100.0;
    println!("\n  D₉₅  = {d95:.2} Gy  ({:.1}% of Rx)", d95 / PRESCRIPTION_GY * 100.0);
    println!("  V₉₅%  = {v95_pct:.1}%  (volume receiving ≥ 95% Rx = {:.1} Gy)", PRESCRIPTION_GY * 0.95);

    // ICRU-83 homogeneity index (lower = better)
    let hi = dvh.homogeneity_index();
    println!("  HI    = {hi:.4}  (0 = perfectly uniform)");

    // PTV-masked DVH: inner 15×15×15 sphere of ~40 mm radius
    let centre_idx: [usize; 3] = [15, 15, 15];
    let ptv_mask_radius: usize = 7; // voxels (~14 mm)
    let ptv_dvh = Dvh::from_volume_masked(&dose, |idx| {
        let di = (idx[0] as isize - centre_idx[0] as isize).pow(2);
        let dj = (idx[1] as isize - centre_idx[1] as isize).pow(2);
        let dk = (idx[2] as isize - centre_idx[2] as isize).pow(2);
        (di + dj + dk) as usize <= ptv_mask_radius * ptv_mask_radius
    });
    let ptv_d95 = ptv_dvh.dose_at_volume_fraction(0.95);
    let ptv_mean = ptv_dvh.mean();
    println!("\nPTV-masked DVH  ({} voxels inside r≤{} voxel radius)", ptv_dvh.count(), ptv_mask_radius);
    println!("  PTV D_mean = {ptv_mean:.2} Gy");
    println!("  PTV D₉₅   = {ptv_d95:.2} Gy");

    assert!(
        ptv_d95 >= PRESCRIPTION_GY * 0.90,
        "PTV D₉₅ {ptv_d95:.2} Gy is below 90% of prescription"
    );
    assert!(
        ptv_mean >= PRESCRIPTION_GY * 0.95,
        "PTV mean dose {ptv_mean:.2} Gy is below 95% of prescription"
    );

    println!("\n✓  DVH coverage metrics pass clinical criteria");
}
