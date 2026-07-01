//! Compton scattering — the Klein–Nishina total cross-section.
//!
//! Compton (incoherent) scattering dominates photon interactions at MV energies
//! (the TomoTherapy 6 MV regime), so its cross-section drives attenuation and
//! energy transfer there. The Klein–Nishina total cross-section **per electron**
//! is, with `α = E / m_e c²` the photon energy in electron-rest-mass units,
//!
//! ```text
//! σ_KN = 2π r_e² { (1+α)/α² [ 2(1+α)/(1+2α) − ln(1+2α)/α ]
//!                  + ln(1+2α)/(2α) − (1+3α)/(1+2α)² }
//! ```
//!
//! As `α → 0` this reduces to the classical Thomson cross-section
//! `σ_T = (8/3)π r_e²` (the low-energy analytical oracle used in the tests).
//! Result units are m² per electron.

use crate::attenuation::MassAttenuation;
use helios_core::constants::{
    AVOGADRO_PER_MOL, CLASSICAL_ELECTRON_RADIUS_M, ELECTRON_REST_ENERGY_MEV,
};
use helios_math::{NumericElement, Scalar};

/// Square centimetres per square metre (`1 m² = 10⁴ cm²`).
const CM2_PER_M2: f64 = 1.0e4;

/// Classical Thomson total cross-section per electron, `σ_T = (8/3)π r_e²` (m²).
///
/// The zero-energy limit of [`klein_nishina_cross_section`].
#[must_use]
pub fn thomson_cross_section<T: Scalar>() -> T {
    let r_e = T::from_f64(CLASSICAL_ELECTRON_RADIUS_M);
    T::from_f64(8.0 / 3.0) * T::PI * r_e * r_e
}

/// Klein–Nishina total Compton cross-section per electron (m²) for a photon of
/// energy `photon_energy_mev` (MeV).
///
/// # Panics
/// Does not panic for finite positive energies. Very small energies are handled
/// by the closed form (which is numerically stable away from `α = 0`); for the
/// exact `α → 0` limit use [`thomson_cross_section`].
#[must_use]
pub fn klein_nishina_cross_section<T: Scalar>(photon_energy_mev: T) -> T {
    let r_e = T::from_f64(CLASSICAL_ELECTRON_RADIUS_M);
    let one = <T as NumericElement>::ONE;
    let two = T::from_f64(2.0);
    let three = T::from_f64(3.0);

    let alpha = photon_energy_mev * T::from_f64(ELECTRON_REST_ENERGY_MEV).recip();
    let one_plus_two_alpha = one + two * alpha;
    let ln_term = one_plus_two_alpha.ln();

    let term1 = (one + alpha)
        * (alpha * alpha).recip()
        * (two * (one + alpha) * one_plus_two_alpha.recip() - ln_term * alpha.recip());
    let term2 = ln_term * (two * alpha).recip();
    let term3 = -(one + three * alpha) * (one_plus_two_alpha * one_plus_two_alpha).recip();

    two * T::PI * r_e * r_e * (term1 + term2 + term3)
}

/// Electrons per gram for a material with effective `z_over_a` = ⟨Z/A⟩
/// (`N_A · Z/A`). Water is ≈0.5551 → ≈3.343×10²³ e⁻/g.
#[must_use]
pub fn electrons_per_gram<T: Scalar>(z_over_a: T) -> T {
    T::from_f64(AVOGADRO_PER_MOL) * z_over_a
}

/// Compton contribution to the mass attenuation coefficient (cm²/g), derived
/// from first principles as `(μ/ρ)_C = σ_KN(E) · (electrons per gram)`.
///
/// This connects [`klein_nishina_cross_section`] (m²/electron, converted to cm²)
/// to [`MassAttenuation`]. In the MV regime Compton dominates, so for water this
/// reproduces the tabulated total `μ/ρ` closely (the test validates against the
/// NIST value at 1 MeV) — a derived coefficient, not a fabricated table entry.
#[must_use]
pub fn compton_mass_attenuation<T: Scalar>(
    photon_energy_mev: T,
    electrons_per_gram: T,
) -> MassAttenuation<T> {
    let sigma_cm2 = klein_nishina_cross_section(photon_energy_mev) * T::from_f64(CM2_PER_M2);
    let mu_over_rho = sigma_cm2 * electrons_per_gram;
    MassAttenuation::new(mu_over_rho)
        .expect("invariant: Compton cross-section and electron density are non-negative")
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    // CODATA Thomson cross-section: 6.6524587321e-29 m². Independent numeric check
    // of the (8/3)π r_e² formula against the published value.
    #[test]
    fn thomson_matches_published_value() {
        assert_relative_eq!(
            thomson_cross_section::<f64>(),
            6.652_458_732_1e-29,
            epsilon = 1e-38,
            max_relative = 1e-6
        );
    }

    #[test]
    fn low_energy_limit_approaches_thomson_from_below() {
        // At α ≈ 0.002 (E = 1 keV) Klein–Nishina must be within ~1% of Thomson and
        // strictly below it (σ_KN = σ_T(1 − 2α + …)).
        //
        // Evaluated in f64: the closed form has a (1+α)/α²·[big−big] cancellation
        // that is ill-conditioned as α→0, so this near-limit check requires f64
        // precision (f32 genericity is checked at a well-conditioned energy below).
        let sigma_t = thomson_cross_section::<f64>();
        let sigma_kn = klein_nishina_cross_section(0.001_f64);
        assert!(sigma_kn < sigma_t, "KN must be below Thomson");
        assert_relative_eq!(sigma_kn / sigma_t, 1.0, max_relative = 1e-2);
    }

    #[test]
    fn cross_section_decreases_with_energy() {
        // Compton cross-section falls monotonically across the diagnostic→MV range.
        let energies = [0.03_f64, 0.1, 0.5, 1.25, 6.0, 18.0];
        let sigmas: Vec<f64> = energies
            .iter()
            .map(|&e| klein_nishina_cross_section(e))
            .collect();
        for pair in sigmas.windows(2) {
            assert!(pair[1] < pair[0], "σ must decrease with energy: {pair:?}");
            assert!(pair[1] > 0.0, "σ must stay positive");
        }
    }

    #[test]
    fn six_mv_is_well_below_thomson() {
        // At 6 MeV Compton is far into the relativistic regime; σ ≪ σ_T.
        let ratio = klein_nishina_cross_section(6.0_f64) / thomson_cross_section::<f64>();
        assert!(ratio < 0.2, "σ_KN(6 MeV)/σ_T = {ratio} should be ≪ 1");
        assert!(ratio > 0.0);
    }

    #[test]
    fn electrons_per_gram_matches_water() {
        // Water ⟨Z/A⟩ = 0.5551 → 3.343×10²³ e⁻/g.
        assert_relative_eq!(
            electrons_per_gram(0.5551_f64),
            3.343e23,
            max_relative = 1e-3
        );
    }

    #[test]
    fn derived_compton_mu_over_rho_matches_water_nist_at_1mev() {
        // First-principles Compton μ/ρ for water at 1 MeV vs the NIST total
        // (0.0707 cm²/g), which is Compton-dominated (~99.8%) at this energy.
        // The derived Compton-only value must match to within ~2%.
        let epg = electrons_per_gram(0.5551_f64);
        let mu_over_rho = compton_mass_attenuation(1.0_f64, epg).get();
        assert_relative_eq!(mu_over_rho, 0.0707, max_relative = 2e-2);
    }

    #[test]
    fn derived_mu_over_rho_decreases_with_energy() {
        let epg = electrons_per_gram(0.5551_f64);
        let low = compton_mass_attenuation(0.5_f64, epg).get();
        let high = compton_mass_attenuation(6.0_f64, epg).get();
        assert!(
            high < low && high > 0.0,
            "μ/ρ Compton must fall with energy"
        );
    }

    #[test]
    fn cross_section_is_generic_over_scalar_f32() {
        // Evaluate at a well-conditioned energy (E = 0.1 MeV, α ≈ 0.2 — no near-α=0
        // cancellation) and check the f32 path reproduces the f64 physics and stays
        // below Thomson. Differential f32-vs-f64 within f32 precision.
        let sigma_kn_f32 = klein_nishina_cross_section(0.1_f32);
        let sigma_kn_f64 = klein_nishina_cross_section(0.1_f64);
        assert!(sigma_kn_f32 < thomson_cross_section::<f32>());
        assert_relative_eq!(
            f64::from(sigma_kn_f32) / sigma_kn_f64,
            1.0,
            max_relative = 1e-4
        );
    }
}
