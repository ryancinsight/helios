//! Photon attenuation: Beer–Lambert transmission and CT-density calibration.
//!
//! The narrow-beam Beer–Lambert law gives the transmitted fraction of a
//! monoenergetic photon beam through a homogeneous slab of thickness `x`:
//!
//! ```text
//! I(x) / I₀ = exp(−μ·x)
//! ```
//!
//! where `μ` is the **linear attenuation coefficient** (cm⁻¹). `μ = (μ/ρ)·ρ`
//! relates it to the material-specific **mass attenuation coefficient** `μ/ρ`
//! (cm²/g) and mass density `ρ` (g/cm³). Tabulated material data is exposed
//! through the [`tables`] module and constructs these validated types directly.

use helios_core::HeliosError;
use helios_math::{NumericElement, Scalar};

/// Tabulated material mass-attenuation data.
pub mod tables;

pub use tables::NistMaterial;

/// Linear attenuation coefficient `μ` in cm⁻¹.
///
/// Validated finite and non-negative at construction (a negative coefficient
/// would model gain, which is unphysical for passive attenuation).
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
#[repr(transparent)]
pub struct LinearAttenuation<T: Scalar>(T);

impl<T: Scalar> LinearAttenuation<T> {
    /// Construct from a coefficient in cm⁻¹.
    ///
    /// # Errors
    /// Returns [`HeliosError::InvalidDomainValue`] if `mu` is non-finite or negative.
    pub fn new(mu: T) -> Result<Self, HeliosError> {
        if !mu.is_finite() {
            return Err(HeliosError::InvalidDomainValue {
                field: "LinearAttenuation",
                value: mu.to_f64(),
                reason: "attenuation coefficient must be finite",
            });
        }
        if mu < <T as NumericElement>::ZERO {
            return Err(HeliosError::InvalidDomainValue {
                field: "LinearAttenuation",
                value: mu.to_f64(),
                reason: "attenuation coefficient must be non-negative",
            });
        }
        Ok(Self(mu))
    }

    /// The coefficient in cm⁻¹.
    #[must_use]
    pub fn get(self) -> T {
        self.0
    }

    /// Narrow-beam transmitted fraction `exp(−μ·path)` through `path_length_cm`.
    ///
    /// `path_length_cm` is the geometric path through the medium (cm) and is
    /// expected to be non-negative; the pure Beer–Lambert value is returned
    /// without clamping so callers can compose path segments.
    #[must_use]
    pub fn transmission(self, path_length_cm: T) -> T {
        (-(self.0 * path_length_cm)).exp()
    }

    /// Half-value layer `ln(2)/μ` (cm): the thickness that halves intensity.
    ///
    /// Returns `None` when `μ = 0` (a non-attenuating medium has no finite HVL).
    #[must_use]
    pub fn half_value_layer(self) -> Option<T> {
        if self.0 == <T as NumericElement>::ZERO {
            return None;
        }
        Some(T::LN_2 * self.0.recip())
    }
}

/// Mass attenuation coefficient `μ/ρ` in cm²/g.
///
/// Validated finite and non-negative at construction. Multiply by a mass density
/// to obtain a [`LinearAttenuation`].
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
#[repr(transparent)]
pub struct MassAttenuation<T: Scalar>(T);

impl<T: Scalar> MassAttenuation<T> {
    /// Construct from a coefficient in cm²/g.
    ///
    /// # Errors
    /// Returns [`HeliosError::InvalidDomainValue`] if the value is non-finite or negative.
    pub fn new(mu_over_rho: T) -> Result<Self, HeliosError> {
        if !mu_over_rho.is_finite() {
            return Err(HeliosError::InvalidDomainValue {
                field: "MassAttenuation",
                value: mu_over_rho.to_f64(),
                reason: "mass attenuation coefficient must be finite",
            });
        }
        if mu_over_rho < <T as NumericElement>::ZERO {
            return Err(HeliosError::InvalidDomainValue {
                field: "MassAttenuation",
                value: mu_over_rho.to_f64(),
                reason: "mass attenuation coefficient must be non-negative",
            });
        }
        Ok(Self(mu_over_rho))
    }

    /// The coefficient in cm²/g.
    #[must_use]
    pub fn get(self) -> T {
        self.0
    }

    /// Linear attenuation `μ = (μ/ρ)·ρ` for mass density `density_g_cm3` (g/cm³).
    ///
    /// # Errors
    /// Returns [`HeliosError::InvalidDomainValue`] if `density_g_cm3` is non-finite
    /// or negative.
    pub fn to_linear(self, density_g_cm3: T) -> Result<LinearAttenuation<T>, HeliosError> {
        if !density_g_cm3.is_finite() || density_g_cm3 < <T as NumericElement>::ZERO {
            return Err(HeliosError::InvalidDomainValue {
                field: "MassAttenuation::to_linear::density",
                value: density_g_cm3.to_f64(),
                reason: "mass density must be finite and non-negative",
            });
        }
        LinearAttenuation::new(self.0 * density_g_cm3)
    }
}

/// Relative electron/mass density from a CT number, `max(0, 1 + HU/1000)`.
///
/// First-order linear CT calibration: air (`−1000 HU`) → `0`, water (`0 HU`) → `1`,
/// and it scales linearly above. Clamped at zero below air. A scanner-specific
/// bilinear (stoichiometric) calibration is a later refinement; this is the
/// vendor-independent baseline.
#[must_use]
pub fn relative_electron_density_from_hu<T: Scalar>(hu: T) -> T {
    let ratio = T::ONE + hu * T::from_f64(1.0e-3);
    ratio.max_scalar(<T as NumericElement>::ZERO)
}

/// Mass density (g/cm³) from a CT number via [`relative_electron_density_from_hu`]
/// scaled by the reference water density `water_density_g_cm3`.
#[must_use]
pub fn mass_density_from_hu<T: Scalar>(hu: T, water_density_g_cm3: T) -> T {
    relative_electron_density_from_hu(hu) * water_density_g_cm3
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn linear_attenuation_rejects_negative_and_nonfinite() {
        assert!(LinearAttenuation::new(-0.1_f64).is_err());
        assert!(LinearAttenuation::new(f64::NAN).is_err());
        assert!(LinearAttenuation::new(f64::INFINITY).is_err());
        assert_eq!(LinearAttenuation::new(0.0_f64).unwrap().get(), 0.0);
        assert_eq!(LinearAttenuation::new(0.5_f64).unwrap().get(), 0.5);
    }

    #[test]
    fn transmission_through_zero_thickness_is_unity() {
        let mu = LinearAttenuation::new(0.2_f64).unwrap();
        assert_eq!(mu.transmission(0.0), 1.0);
    }

    #[test]
    fn transmission_at_half_value_layer_is_one_half() {
        // Analytical oracle: by definition T(HVL) = 1/2 for any μ > 0.
        let mu = LinearAttenuation::new(0.35_f64).unwrap();
        let hvl = mu.half_value_layer().expect("mu > 0 has a finite HVL");
        assert_relative_eq!(mu.transmission(hvl), 0.5, epsilon = 1e-12);
    }

    #[test]
    fn half_value_layer_is_none_for_non_attenuating_medium() {
        assert_eq!(
            LinearAttenuation::new(0.0_f64).unwrap().half_value_layer(),
            None
        );
    }

    #[test]
    fn mass_to_linear_applies_density_and_beer_lambert() {
        // μ = (μ/ρ)·ρ; value-semantic on the defining relation (not a memorized
        // NIST digit): 0.06 cm²/g at ρ = 1 g/cm³ → μ = 0.06 cm⁻¹.
        let mu_over_rho = MassAttenuation::new(0.06_f64).unwrap();
        let mu = mu_over_rho.to_linear(1.0).unwrap();
        assert_relative_eq!(mu.get(), 0.06, epsilon = 1e-15);
        // Doubling density doubles μ.
        assert_relative_eq!(
            mu_over_rho.to_linear(2.0).unwrap().get(),
            0.12,
            epsilon = 1e-15
        );
        // Transmission through 10 cm matches exp(−0.6).
        assert_relative_eq!(mu.transmission(10.0), (-0.6_f64).exp(), epsilon = 1e-12);
    }

    #[test]
    fn mass_to_linear_rejects_bad_density() {
        let m = MassAttenuation::new(0.06_f64).unwrap();
        assert!(m.to_linear(-1.0).is_err());
        assert!(m.to_linear(f64::NAN).is_err());
    }

    #[test]
    fn hu_calibration_hits_reference_points() {
        // air, water, and a linear point above water.
        assert_relative_eq!(
            relative_electron_density_from_hu(-1000.0_f64),
            0.0,
            epsilon = 1e-15
        );
        assert_relative_eq!(
            relative_electron_density_from_hu(0.0_f64),
            1.0,
            epsilon = 1e-15
        );
        assert_relative_eq!(
            relative_electron_density_from_hu(1000.0_f64),
            2.0,
            epsilon = 1e-15
        );
        // Below air is clamped to zero.
        assert_relative_eq!(
            relative_electron_density_from_hu(-1200.0_f64),
            0.0,
            epsilon = 1e-15
        );
    }

    #[test]
    fn mass_density_scales_by_water_density() {
        // With unit-density water, HU=0 → 1 g/cm³; HU=1000 → 2 g/cm³.
        assert_relative_eq!(mass_density_from_hu(0.0_f64, 1.0), 1.0, epsilon = 1e-15);
        assert_relative_eq!(mass_density_from_hu(1000.0_f64, 1.0), 2.0, epsilon = 1e-15);
    }

    #[test]
    fn attenuation_is_generic_over_scalar_f32() {
        let mu = LinearAttenuation::new(0.5_f32).unwrap();
        let hvl = mu.half_value_layer().unwrap();
        assert_relative_eq!(mu.transmission(hvl), 0.5_f32, epsilon = 1e-6);
    }

    proptest::proptest! {
        /// Beer–Lambert transmission is bounded in `[0, 1]` for any non-negative
        /// coefficient and path length. (At extreme optical depth `exp(−μx)`
        /// underflows to exactly 0 — physical extinction — so 0 is included.)
        #[test]
        fn transmission_is_bounded_unit_interval(
            mu_val in 0.0_f64..20.0,
            path in 0.0_f64..100.0,
        ) {
            let mu = LinearAttenuation::new(mu_val).unwrap();
            let t = mu.transmission(path);
            proptest::prop_assert!((0.0..=1.0).contains(&t), "T={t}");
        }

        /// Transmission decreases monotonically with path length (more material →
        /// less transmitted).
        #[test]
        fn transmission_decreases_with_path(
            mu_val in 1e-3_f64..20.0,
            p1 in 0.0_f64..50.0,
            extra in 0.0_f64..50.0,
        ) {
            let mu = LinearAttenuation::new(mu_val).unwrap();
            proptest::prop_assert!(mu.transmission(p1 + extra) <= mu.transmission(p1));
        }

        /// Beer–Lambert composes multiplicatively over concatenated path segments:
        /// `T(x₁+x₂) = T(x₁)·T(x₂)`.
        #[test]
        fn transmission_composes_over_segments(
            mu_val in 0.0_f64..10.0,
            x1 in 0.0_f64..30.0,
            x2 in 0.0_f64..30.0,
        ) {
            let mu = LinearAttenuation::new(mu_val).unwrap();
            let whole = mu.transmission(x1 + x2);
            let parts = mu.transmission(x1) * mu.transmission(x2);
            proptest::prop_assert!((whole - parts).abs() <= 1e-12 * (1.0 + whole));
        }

        /// CT-density calibration is monotonic non-decreasing in HU (denser tissue
        /// → higher density).
        #[test]
        fn density_is_monotonic_in_hu(a in -1000.0_f64..3000.0, b in -1000.0_f64..3000.0) {
            let (da, db) = (mass_density_from_hu(a, 1.0), mass_density_from_hu(b, 1.0));
            proptest::prop_assert_eq!(a <= b, da <= db);
        }
    }
}
