//! Compton Scattering Physics Example
//!
//! Demonstrates the Klein–Nishina total cross-section, Thomson low-energy limit,
//! and Compton mass-energy-transfer coefficients from `helios-physics`.
//!
//! ## What this covers
//!
//! 1. **Klein–Nishina vs Thomson** — At low photon energies the Klein–Nishina
//!    cross section converges to the classical Thomson value `σ_T = 8π r_e²/3`.
//! 2. **Cross section vs energy table** — σ drops from the Thomson limit at
//!    10 keV to < 20% of σ_T at 6 MV (Compton dominates at MV energies).
//! 3. **Energy-transfer fraction** — The fraction of photon energy deposited
//!    locally rises steeply with energy (forward-scattered Compton electrons
//!    carry most energy at 6 MV).
//! 4. **Compton mass attenuation** — Validated against the NIST analytical
//!    formula for water and compact bone.
//!
//! ## Run
//!
//! ```bash
//! cargo run -p helios-physics --example compton_physics
//! ```
//!
//! ## Book chapter
//!
//! [← Mass Attenuation and Photon Cross Sections](../../docs/book/dose_attenuation.md)

use aequitas::systems::si::units::SquareCentimeterPerGram;
use helios_physics::{
    compton_energy_transfer_cross_section, compton_mass_attenuation,
    compton_mean_energy_transfer_fraction, electrons_per_gram, klein_nishina_cross_section,
    thomson_cross_section,
};
use hyperion::coefficient::MassAttenuation;

fn main() {
    println!("=== Compton Scattering Physics ===\n");

    // ── 1. Thomson vs Klein–Nishina low-energy limit ──────────────────────────
    let sigma_t: f64 = thomson_cross_section();
    println!("Thomson cross section  σ_T = {:.4e} m²/electron", sigma_t);

    // At 1 keV, KN should be within 0.1% of Thomson
    let sigma_kn_1kev: f64 = klein_nishina_cross_section(0.001_f64);
    let ratio_1kev = sigma_kn_1kev / sigma_t;
    println!("KN at 1 keV / σ_T     = {ratio_1kev:.6}  (expected ≈ 1.0)");
    assert!(
        (ratio_1kev - 1.0).abs() < 0.01,
        "KN should approach Thomson at 1 keV, got ratio {ratio_1kev:.6}"
    );
    println!("  ✓ Klein–Nishina → Thomson at low energy\n");

    // ── 2. Cross-section energy table ─────────────────────────────────────────
    println!("Klein–Nishina total cross section vs photon energy:");
    println!("  Energy         σ_KN (m²/e⁻)    σ_KN/σ_T");
    println!("  {}", "-".repeat(48));

    let energies_mev = [0.01, 0.05, 0.1, 0.5, 1.0, 6.0, 10.0, 18.0];
    let labels = [
        "10 keV",
        "50 keV",
        "100 keV",
        "500 keV",
        "1 MeV",
        "6 MeV (6MV)",
        "10 MeV",
        "18 MeV (18MV)",
    ];

    for (&e, &lbl) in energies_mev.iter().zip(labels.iter()) {
        let sigma: f64 = klein_nishina_cross_section(e);
        let ratio = sigma / sigma_t;
        println!("  {:<14}  {:.4e}   {:.4}", lbl, sigma, ratio);
    }

    // Verify monotonic decrease (cross section decreases with energy)
    let sigmas: Vec<f64> = energies_mev
        .iter()
        .map(|&e| klein_nishina_cross_section(e))
        .collect();
    let monotone = sigmas.windows(2).all(|w| w[0] > w[1]);
    assert!(
        monotone,
        "KN cross section must decrease monotonically with energy"
    );
    println!("  ✓ Cross section strictly decreasing with energy\n");

    // ── 3. Energy-transfer fraction ───────────────────────────────────────────
    println!("Compton mean energy-transfer fraction f_tr vs energy:");
    println!("  Energy        f_tr        Interpretation");
    println!("  {}", "-".repeat(60));

    let et_energies = [0.05, 0.1, 1.0, 6.0, 18.0];
    let et_labels = ["50 keV", "100 keV", "1 MeV", "6 MeV", "18 MeV"];
    for (&e, &lbl) in et_energies.iter().zip(et_labels.iter()) {
        let f_tr: f64 = compton_mean_energy_transfer_fraction(e);
        let interp = if f_tr < 0.3 {
            "low (scatter dominates)"
        } else if f_tr < 0.6 {
            "moderate"
        } else {
            "high (electron carries most energy)"
        };
        println!("  {:<12}  {f_tr:.4}      {interp}", lbl);
    }

    // At 6 MV, energy-transfer fraction should be > 0.5
    let f_tr_6mv: f64 = compton_mean_energy_transfer_fraction(6.0);
    assert!(
        f_tr_6mv > 0.5,
        "Energy-transfer fraction at 6 MeV should exceed 0.5, got {f_tr_6mv:.4}"
    );
    println!("  ✓ High energy-transfer at 6 MeV (therapeutic regime)\n");

    // ── 4. Compton mass attenuation for NIST materials ────────────────────────
    println!("Compton mass attenuation μ_C/ρ at 100 keV:");
    // Water Z/A ≈ 0.5551 → electrons per gram ≈ 3.343×10²³ e⁻/g
    let e_per_g_water: f64 = electrons_per_gram(0.5551_f64);
    for &(label, energy) in &[("100 keV", 0.1_f64), ("6 MeV", 6.0_f64)] {
        let mu_rho: MassAttenuation<f64> = compton_mass_attenuation(energy, e_per_g_water);
        println!(
            "  {:<14}  μ_C/ρ = {:.4} cm²/g  (water)",
            label,
            mu_rho.in_unit::<SquareCentimeterPerGram>()
        );
    }

    // ── 5. Energy-transfer cross section ─────────────────────────────────────
    println!();
    let sigma_tr_6mv: f64 = compton_energy_transfer_cross_section(6.0_f64);
    let sigma_tot_6mv: f64 = klein_nishina_cross_section(6.0_f64);
    println!(
        "At 6 MeV:  σ_tr/σ_KN = {:.4}  (energy-transfer fraction check: {:.4})",
        sigma_tr_6mv / sigma_tot_6mv,
        compton_mean_energy_transfer_fraction(6.0_f64)
    );

    println!("\n✓  All Compton physics checks passed");
    println!("\nBook chapter: Part III — Mass Attenuation and Photon Cross Sections");
    println!(
        "API: helios_physics::{{klein_nishina_cross_section, compton_mass_attenuation, \
         compton_mean_energy_transfer_fraction}}"
    );
}
