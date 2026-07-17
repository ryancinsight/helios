//! NIST photon mass-attenuation tables for selected reference materials.
//!
//! The embedded knots are the `μ/ρ` column (cm²/g) from NIST's **X-Ray Mass
//! Attenuation Coefficients** tables for [dry air], [liquid water], and
//! [cortical bone]. The common 10 keV–20 MeV range intentionally excludes the
//! absorption-edge rows in the air and bone tables. The lookup is therefore
//! defined at every stored knot and uses log-linear interpolation only between
//! adjacent edge-free knots.
//!
//! This is not a reproduction of XCOM output: XCOM uses log-log cubic-spline
//! fits for total attenuation coefficients and separately handles absorption
//! edges. The bounded, allocation-free table API is suitable for the imaging
//! and therapy photon-energy range while keeping that distinction explicit.
//!
//! [dry air]: https://physics.nist.gov/PhysRefData/XrayMassCoef/ComTab/air.html
//! [liquid water]: https://physics.nist.gov/PhysRefData/XrayMassCoef/ComTab/water.html
//! [cortical bone]: https://physics.nist.gov/PhysRefData/XrayMassCoef/ComTab/bone.html

use helios_core::HeliosError;
use helios_math::Scalar;

use super::MassAttenuation;

const KNOT_COUNT: usize = 28;
const MINIMUM_ENERGY_MEV: f64 = 0.01;
const MAXIMUM_ENERGY_MEV: f64 = 20.0;

/// Photon-energy knots in MeV shared by every embedded material table.
const PHOTON_ENERGY_MEV: [f64; KNOT_COUNT] = [
    0.01, 0.015, 0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0,
    1.25, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 15.0, 20.0,
];

const DRY_AIR_MASS_ATTENUATION: [f64; KNOT_COUNT] = [
    5.120, 1.614, 0.7779, 0.3538, 0.2485, 0.2080, 0.1875, 0.1662, 0.1541, 0.1356, 0.1233, 0.1067,
    0.09549, 0.08712, 0.08055, 0.07074, 0.06358, 0.05687, 0.05175, 0.04447, 0.03581, 0.03079,
    0.02751, 0.02522, 0.02225, 0.02045, 0.01810, 0.01705,
];

const LIQUID_WATER_MASS_ATTENUATION: [f64; KNOT_COUNT] = [
    5.329, 1.673, 0.8096, 0.3756, 0.2683, 0.2269, 0.2059, 0.1837, 0.1707, 0.1505, 0.1370, 0.1186,
    0.1061, 0.09687, 0.08956, 0.07865, 0.07072, 0.06323, 0.05754, 0.04942, 0.03969, 0.03403,
    0.03031, 0.02770, 0.02429, 0.02219, 0.01941, 0.01813,
];

const CORTICAL_BONE_MASS_ATTENUATION: [f64; KNOT_COUNT] = [
    28.51, 9.032, 4.001, 1.331, 0.6655, 0.4242, 0.3148, 0.2229, 0.1855, 0.1480, 0.1309, 0.1113,
    0.09908, 0.09022, 0.08332, 0.07308, 0.06566, 0.05871, 0.05346, 0.04607, 0.03745, 0.03257,
    0.02946, 0.02734, 0.02467, 0.02314, 0.02132, 0.02068,
];

/// Material whose NIST mass-attenuation table is embedded in Helios.
///
/// The three tables share one energy grid, so selecting a material changes
/// only the source coefficient slice; lookup remains allocation-free and
/// executes interpolation in the caller's native [`Scalar`] precision.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum NistMaterial {
    /// Dry air near sea level.
    DryAir,
    /// Liquid water.
    LiquidWater,
    /// Cortical bone from ICRU-44.
    CorticalBone,
}

impl NistMaterial {
    /// Return the material mass attenuation coefficient at `photon_energy_mev`.
    ///
    /// Exact table knots return the stored NIST value. Between adjacent knots,
    /// this evaluates a power law in native `T` arithmetic: `ln(μ/ρ)` is linear
    /// in `ln(E)`. The selected 10 keV–20 MeV range contains no absorption edge
    /// for these embedded material tables, so interpolation never crosses a
    /// represented discontinuity.
    ///
    /// # Errors
    ///
    /// Returns [`HeliosError::InvalidDomainValue`] when `photon_energy_mev` is
    /// non-finite or lies outside the embedded 0.01–20 MeV interval.
    pub fn mass_attenuation<T: Scalar>(
        self,
        photon_energy_mev: T,
    ) -> Result<MassAttenuation<T>, HeliosError> {
        if !photon_energy_mev.is_finite() {
            return Err(invalid_energy(
                photon_energy_mev,
                "photon energy must be finite",
            ));
        }

        let minimum_energy = T::from_f64(MINIMUM_ENERGY_MEV);
        let maximum_energy = T::from_f64(MAXIMUM_ENERGY_MEV);
        if photon_energy_mev < minimum_energy || photon_energy_mev > maximum_energy {
            return Err(invalid_energy(
                photon_energy_mev,
                "photon energy must be within the embedded 0.01–20 MeV table range",
            ));
        }

        let coefficients = self.coefficients();
        let upper_knot = upper_knot(photon_energy_mev);
        if photon_energy_mev == T::from_f64(PHOTON_ENERGY_MEV[upper_knot]) {
            return MassAttenuation::new(T::from_f64(coefficients[upper_knot]));
        }

        let lower_knot = upper_knot - 1;
        let lower_energy = T::from_f64(PHOTON_ENERGY_MEV[lower_knot]);
        let upper_energy = T::from_f64(PHOTON_ENERGY_MEV[upper_knot]);
        let lower_coefficient = T::from_f64(coefficients[lower_knot]);
        let upper_coefficient = T::from_f64(coefficients[upper_knot]);

        let log_fraction = (photon_energy_mev.ln() - lower_energy.ln())
            * (upper_energy.ln() - lower_energy.ln()).recip();
        let log_coefficient = lower_coefficient.ln()
            + (upper_coefficient.ln() - lower_coefficient.ln()) * log_fraction;
        MassAttenuation::new(log_coefficient.exp())
    }

    fn coefficients(self) -> &'static [f64; KNOT_COUNT] {
        match self {
            Self::DryAir => &DRY_AIR_MASS_ATTENUATION,
            Self::LiquidWater => &LIQUID_WATER_MASS_ATTENUATION,
            Self::CorticalBone => &CORTICAL_BONE_MASS_ATTENUATION,
        }
    }
}

fn upper_knot<T: Scalar>(photon_energy_mev: T) -> usize {
    let mut upper_knot = 0;
    while upper_knot + 1 < KNOT_COUNT
        && photon_energy_mev > T::from_f64(PHOTON_ENERGY_MEV[upper_knot])
    {
        upper_knot += 1;
    }
    upper_knot
}

fn invalid_energy<T: Scalar>(photon_energy_mev: T, reason: &'static str) -> HeliosError {
    HeliosError::InvalidDomainValue {
        field: "NistMaterial::mass_attenuation::photon_energy_mev",
        value: photon_energy_mev.to_f64(),
        reason,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn assert_exact_knot<T: Scalar>() {
        assert_eq!(
            NistMaterial::DryAir
                .mass_attenuation(T::from_f64(1.0))
                .unwrap()
                .get(),
            T::from_f64(0.06358)
        );
        assert_eq!(
            NistMaterial::LiquidWater
                .mass_attenuation(T::from_f64(1.0))
                .unwrap()
                .get(),
            T::from_f64(0.07072)
        );
        assert_eq!(
            NistMaterial::CorticalBone
                .mass_attenuation(T::from_f64(1.0))
                .unwrap()
                .get(),
            T::from_f64(0.06566)
        );
    }

    #[test]
    fn returns_official_material_knot_values_at_one_mev() {
        assert_exact_knot::<f32>();
        assert_exact_knot::<f64>();
    }

    #[test]
    fn returns_the_first_and_last_embedded_knots() {
        assert_eq!(
            NistMaterial::LiquidWater
                .mass_attenuation(0.01_f64)
                .unwrap()
                .get(),
            5.329
        );
        assert_eq!(
            NistMaterial::LiquidWater
                .mass_attenuation(20.0_f64)
                .unwrap()
                .get(),
            0.01813
        );
    }

    #[test]
    fn rejects_energies_outside_the_embedded_interval() {
        assert_eq!(
            NistMaterial::LiquidWater.mass_attenuation(0.009_f64),
            Err(HeliosError::InvalidDomainValue {
                field: "NistMaterial::mass_attenuation::photon_energy_mev",
                value: 0.009,
                reason: "photon energy must be within the embedded 0.01–20 MeV table range",
            })
        );
        assert_eq!(
            NistMaterial::LiquidWater.mass_attenuation(20.1_f64),
            Err(HeliosError::InvalidDomainValue {
                field: "NistMaterial::mass_attenuation::photon_energy_mev",
                value: 20.1,
                reason: "photon energy must be within the embedded 0.01–20 MeV table range",
            })
        );
    }

    #[test]
    fn rejects_nonfinite_energy_with_the_failing_contract() {
        match NistMaterial::DryAir.mass_attenuation(f64::NAN) {
            Err(HeliosError::InvalidDomainValue {
                field,
                value,
                reason,
            }) => {
                assert_eq!(field, "NistMaterial::mass_attenuation::photon_energy_mev");
                assert!(value.is_nan(), "the reported value must preserve NaN");
                assert_eq!(reason, "photon energy must be finite");
            }
            result => panic!("expected a non-finite energy error, got {result:?}"),
        }
    }

    fn assert_log_interpolation<T: Scalar>() {
        let lower_energy = T::from_f64(0.1);
        let upper_energy = T::from_f64(0.15);
        let geometric_midpoint = (lower_energy * upper_energy).sqrt();
        let expected = (T::from_f64(0.1707) * T::from_f64(0.1505)).sqrt();
        let actual = NistMaterial::LiquidWater
            .mass_attenuation(geometric_midpoint)
            .unwrap()
            .get();
        let roundoff_bound = T::from_f64(32.0) * T::EPSILON * expected;
        assert!(
            (actual - expected).abs() <= roundoff_bound,
            "actual={actual:?}, expected={expected:?}, bound={roundoff_bound:?}"
        );
    }

    #[test]
    fn interpolates_log_linearly_in_each_native_scalar_precision() {
        // At the geometric mean of two energies, log-linear interpolation must
        // equal the geometric mean of their coefficients. The bound covers the
        // finite sequence of log, subtraction, division, multiplication, and exp
        // operations (32 ulps is conservative for those operations in f32/f64).
        assert_log_interpolation::<f32>();
        assert_log_interpolation::<f64>();
    }
}
