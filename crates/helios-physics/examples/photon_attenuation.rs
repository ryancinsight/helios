//! Helios photon attenuation physics example.
//!
//! Demonstrates the Beer–Lambert transmission model, mass attenuation, and
//! CT-number-to-density calibration from `helios-physics`. Covers:
//!
//! 1. **LinearAttenuation** — validated μ (cm⁻¹) type with transmission law
//! 2. **Half-value layer** — HVL = ln(2)/μ for a monoenergetic beam
//! 3. **MassAttenuation** → LinearAttenuation — scaling μ/ρ by mass density
//! 4. **CT calibration** — HU → relative electron density → mass density
//!
//! Run with: cargo run --example photon_attenuation -p helios-physics

use helios_physics::{
    mass_density_from_hu, relative_electron_density_from_hu, LinearAttenuation, MassAttenuation,
};

fn main() {
    // ── 1. Beer–Lambert transmission ─────────────────────────────────────────

    // Water linear attenuation at ~100 keV: μ ≈ 0.171 cm⁻¹
    let mu_water = LinearAttenuation::new(0.171_f64).expect("valid μ");
    println!("Water at 100 keV:  μ = {:.4} cm⁻¹", mu_water.get());

    let hvl_cm = mu_water.half_value_layer().expect("non-zero μ");
    println!("Half-value layer:  HVL = {hvl_cm:.3} cm  (should be ≈ 4.05 cm)");
    assert!(
        (hvl_cm - 4.05).abs() < 0.1,
        "HVL mismatch: got {hvl_cm:.3} cm"
    );

    // Transmission through 10 cm of water
    let path_cm = 10.0_f64;
    let t = mu_water.transmission(path_cm);
    let t_analytical = (-mu_water.get() * path_cm).exp();
    println!("Transmission(10 cm water) = {t:.4}  (analytical = {t_analytical:.4})");
    assert!((t - t_analytical).abs() < 1e-10, "Transmission mismatch");

    // ── 2. Exponential attenuation vs depth ──────────────────────────────────

    println!("\nAttenuation depth profile (water, 100 keV):");
    println!("  Depth (cm)  |  Transmitted fraction");
    for &depth in &[0.0, 1.0, 2.0, 5.0, 10.0, 20.0] {
        let frac = mu_water.transmission(depth);
        let bar: String = "#".repeat((frac * 40.0) as usize);
        println!("  {:10.1}  |  {:.4}  {bar}", depth, frac);
    }

    // ── 3. MassAttenuation → LinearAttenuation ───────────────────────────────

    // Soft tissue at 100 keV: μ/ρ ≈ 0.169 cm²/g, ρ ≈ 1.06 g/cm³
    let mu_rho_tissue = MassAttenuation::new(0.169_f64).expect("valid μ/ρ");
    let density_tissue = 1.06_f64; // g/cm³
    let mu_tissue = mu_rho_tissue
        .to_linear(density_tissue)
        .expect("valid density");
    println!("\nSoft tissue  μ/ρ = {:.4} cm²/g  ρ = {density_tissue} g/cm³  → μ = {:.4} cm⁻¹",
        mu_rho_tissue.get(), mu_tissue.get());
    assert!(
        (mu_tissue.get() - 0.169 * 1.06).abs() < 1e-6,
        "μ = (μ/ρ)ρ identity failed"
    );

    // ── 4. CT number → density calibration ───────────────────────────────────

    println!("\nCT calibration (HU → relative electron density → mass density):");
    let water_density = 1.0_f64; // g/cm³
    let test_materials: &[(&str, f64)] = &[
        ("Air",        -1000.0),
        ("Lung",        -800.0),
        ("Adipose",     -100.0),
        ("Water",          0.0),
        ("Soft tissue",   50.0),
        ("Compact bone", 700.0),
    ];
    for &(name, hu) in test_materials {
        let rho_rel = relative_electron_density_from_hu(hu);
        let rho = mass_density_from_hu(hu, water_density);
        println!("  {:<16}  HU={:>6.0}  ρ_rel={:.3}  ρ={:.3} g/cm³", name, hu, rho_rel, rho);
    }

    // Water at 0 HU → ρ_rel = 1.0
    let water_rho = relative_electron_density_from_hu(0.0_f64);
    assert!((water_rho - 1.0).abs() < 1e-10, "Water HU→ρ should be exactly 1");

    // Air at -1000 HU → ρ_rel = 0.0
    let air_rho = relative_electron_density_from_hu(-1000.0_f64);
    assert!(air_rho == 0.0, "Air HU→ρ should be 0 (clamped)");

    println!("\n✓  Beer–Lambert, HVL, mass attenuation, and CT calibration all verified");
}
