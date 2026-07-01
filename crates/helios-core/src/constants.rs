//! Physical constants for radiation transport and dosimetry.
//!
//! Values are SI/CODATA-2018 recommended values or ICRU-recommended reference
//! data, quoted at their published precision. They are `f64` at the definition
//! boundary; higher layers convert into `T: Scalar` at call sites. Relationships
//! between constants are checked by the value-semantic tests at the bottom of
//! this module, which serve as the derivation record for the quoted values.
//!
//! # References
//! - CODATA 2018: <https://physics.nist.gov/cuu/Constants/> (fundamental constants).
//! - ICRU Report 90 (2016): key data for measurement standards — liquid-water
//!   mean excitation energy.

/// Speed of light in vacuum, `c` (m·s⁻¹). Exact by SI definition (2019).
pub const SPEED_OF_LIGHT_M_PER_S: f64 = 299_792_458.0;

/// Elementary charge, `e` (C). Exact by SI definition (2019).
pub const ELEMENTARY_CHARGE_C: f64 = 1.602_176_634e-19;

/// Avogadro constant, `N_A` (mol⁻¹). Exact by SI definition (2019).
pub const AVOGADRO_PER_MOL: f64 = 6.022_140_76e23;

/// Vacuum electric permittivity, `ε₀` (F·m⁻¹). CODATA 2018.
pub const VACUUM_PERMITTIVITY_F_PER_M: f64 = 8.854_187_812_8e-12;

/// Electron rest mass, `m_e` (kg). CODATA 2018.
pub const ELECTRON_MASS_KG: f64 = 9.109_383_701_5e-31;

/// Electron rest energy, `m_e c²` (MeV). CODATA 2018: 0.510 998 950 00 MeV.
///
/// Cross-checked against `ELECTRON_MASS_KG`, `SPEED_OF_LIGHT_M_PER_S`, and
/// [`MEV_TO_JOULE`] in the module tests.
pub const ELECTRON_REST_ENERGY_MEV: f64 = 0.510_998_950_00;

/// Classical electron radius, `r_e` (m). CODATA 2018.
///
/// Cross-checked against the defining relation
/// `r_e = e² / (4π ε₀ m_e c²)` in the module tests.
pub const CLASSICAL_ELECTRON_RADIUS_M: f64 = 2.817_940_326_2e-15;

/// Energy conversion factor, joules per megaelectronvolt (J·MeV⁻¹).
///
/// Equal to `ELEMENTARY_CHARGE_C × 10⁶` by definition of the electronvolt; the
/// module tests assert this identity exactly.
pub const MEV_TO_JOULE: f64 = ELEMENTARY_CHARGE_C * 1.0e6;

/// Dosimetric reference density of liquid water (g·cm⁻³).
///
/// Unit-density water is the ICRU dosimetric reference medium; kV/MV imaging and
/// dose engines calibrate CT numbers and stopping powers against it.
pub const WATER_DENSITY_G_PER_CM3: f64 = 1.0;

/// Mean excitation energy of liquid water, `I` (eV). ICRU Report 90 (2016).
///
/// Supersedes the ICRU-37 value of 75.0 eV; used in Bethe stopping-power
/// evaluation for electron/proton transport.
pub const WATER_MEAN_EXCITATION_ENERGY_EV: f64 = 78.0;

#[cfg(test)]
mod tests {
    use super::*;

    /// The electronvolt conversion is defined as `e` scaled to mega-, so the
    /// identity is exact in IEEE-754 (single multiply by an exactly
    /// representable power of ten is the only rounding, shared by both sides).
    #[test]
    fn mev_to_joule_is_elementary_charge_scaled() {
        assert_eq!(MEV_TO_JOULE, ELEMENTARY_CHARGE_C * 1.0e6);
    }

    /// Electron rest energy derived from mass–energy equivalence must agree with
    /// the quoted `ELECTRON_REST_ENERGY_MEV`.
    ///
    /// Tolerance derivation: the three inputs (`m_e`, `c`, `e`) each carry
    /// CODATA relative uncertainty ≤ 3×10⁻¹⁰; the product/quotient propagates to
    /// ≲ 1×10⁻⁹, and IEEE-754 rounding over ~4 operations adds ≲ 4·ε ≈ 9×10⁻¹⁶.
    /// A relative bound of 1×10⁻⁸ covers both with margin.
    #[test]
    fn electron_rest_energy_matches_mass_energy_equivalence() {
        let rest_energy_joule = ELECTRON_MASS_KG * SPEED_OF_LIGHT_M_PER_S * SPEED_OF_LIGHT_M_PER_S;
        let rest_energy_mev = rest_energy_joule / MEV_TO_JOULE;
        let rel_err = (rest_energy_mev - ELECTRON_REST_ENERGY_MEV).abs() / ELECTRON_REST_ENERGY_MEV;
        assert!(
            rel_err < 1.0e-8,
            "m_e c² = {rest_energy_mev} MeV vs quoted {ELECTRON_REST_ENERGY_MEV} MeV (rel_err {rel_err:e})"
        );
    }

    /// Classical electron radius must satisfy `r_e = e² / (4π ε₀ m_e c²)`.
    ///
    /// Tolerance derivation: inputs carry CODATA relative uncertainty ≤ 1.5×10⁻¹⁰
    /// (ε₀) with the rest bounded tighter; error propagates to ≲ 3×10⁻¹⁰, and the
    /// ~6-operation evaluation adds ≲ 6·ε. A relative bound of 1×10⁻⁸ is
    /// comfortably conservative.
    #[test]
    fn classical_electron_radius_matches_defining_relation() {
        let m_e_c2 = ELECTRON_MASS_KG * SPEED_OF_LIGHT_M_PER_S * SPEED_OF_LIGHT_M_PER_S;
        let r_e = ELEMENTARY_CHARGE_C * ELEMENTARY_CHARGE_C
            / (4.0 * core::f64::consts::PI * VACUUM_PERMITTIVITY_F_PER_M * m_e_c2);
        let rel_err = (r_e - CLASSICAL_ELECTRON_RADIUS_M).abs() / CLASSICAL_ELECTRON_RADIUS_M;
        assert!(
            rel_err < 1.0e-8,
            "r_e = {r_e} m vs quoted {CLASSICAL_ELECTRON_RADIUS_M} m (rel_err {rel_err:e})"
        );
    }

    /// Guard the exactly-defined SI constants against accidental edits.
    #[test]
    fn si_defined_constants_hold_exact_values() {
        assert_eq!(SPEED_OF_LIGHT_M_PER_S, 299_792_458.0);
        assert_eq!(ELEMENTARY_CHARGE_C, 1.602_176_634e-19);
        assert_eq!(AVOGADRO_PER_MOL, 6.022_140_76e23);
    }
}
