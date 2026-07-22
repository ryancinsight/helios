//! CT-number to relative-density calibration.

use helios_math::{NumericElement, Scalar};

/// Relative electron and mass density from a CT number, `max(0, 1 + HU/1000)`.
///
/// This first-order calibration maps air (`-1000 HU`) to zero and water
/// (`0 HU`) to one. Scanner-specific stoichiometric calibration remains a
/// Helios imaging concern rather than a photon-transport law.
#[must_use]
pub fn relative_electron_density_from_hu<T: Scalar>(hu: T) -> T {
    let ratio = T::ONE + hu * T::from_f64(1.0e-3);
    ratio.max_scalar(<T as NumericElement>::ZERO)
}

/// Mass density in `g/cm^3` from CT number and reference water density.
#[must_use]
pub fn mass_density_from_hu<T: Scalar>(hu: T, water_density_g_cm3: T) -> T {
    relative_electron_density_from_hu(hu) * water_density_g_cm3
}

#[cfg(test)]
mod tests {
    use super::*;
    use eunomia::assert_relative_eq;

    #[test]
    fn calibration_hits_air_water_and_dense_reference_points() {
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
        assert_relative_eq!(
            relative_electron_density_from_hu(-1200.0_f64),
            0.0,
            epsilon = 1e-15
        );
    }

    #[test]
    fn mass_density_scales_by_reference_water_density() {
        assert_relative_eq!(mass_density_from_hu(0.0_f64, 1.0), 1.0, epsilon = 1e-15);
        assert_relative_eq!(mass_density_from_hu(1000.0_f64, 1.0), 2.0, epsilon = 1e-15);
    }

    proptest::proptest! {
        /// The first-order calibration is monotone in its valid clinical range.
        #[test]
        fn density_is_monotonic_in_hu(a in -1000.0_f64..3000.0, b in -1000.0_f64..3000.0) {
            let (da, db) = (mass_density_from_hu(a, 1.0), mass_density_from_hu(b, 1.0));
            proptest::prop_assert_eq!(a <= b, da <= db);
        }
    }
}
