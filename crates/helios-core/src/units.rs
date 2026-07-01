//! Validating domain newtypes.
//!
//! Each type is the single entry point through which a raw `f64` becomes a
//! domain value: construction goes through `TryFrom<f64>`, which enforces the
//! invariant so an invalid value cannot exist in memory. The types are
//! `#[repr(transparent)]` — zero-cost over the underlying `f64` and ABI-stable
//! for later FFI/PyO3 exposure — and implement `Display` for diagnostics.

use crate::error::HeliosError;
use core::fmt;

/// Beam energy in megaelectronvolts (MeV).
///
/// Valid range: finite and strictly positive. Photon and electron beam
/// nominal energies (e.g. TomoTherapy's 6 MV MV beam) are positive by
/// definition; zero and non-finite values are rejected at construction.
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
#[repr(transparent)]
pub struct EnergyMeV(f64);

impl EnergyMeV {
    /// Returns the underlying value in MeV.
    #[must_use]
    pub const fn get(self) -> f64 {
        self.0
    }
}

impl TryFrom<f64> for EnergyMeV {
    type Error = HeliosError;

    fn try_from(value: f64) -> Result<Self, Self::Error> {
        if !value.is_finite() {
            return Err(HeliosError::InvalidDomainValue {
                field: "EnergyMeV",
                value,
                reason: "energy must be finite",
            });
        }
        if value <= 0.0 {
            return Err(HeliosError::InvalidDomainValue {
                field: "EnergyMeV",
                value,
                reason: "energy must be strictly positive",
            });
        }
        Ok(Self(value))
    }
}

impl fmt::Display for EnergyMeV {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} MeV", self.0)
    }
}

/// A CT number in Hounsfield units (HU).
///
/// Valid range: finite and within `[-1024, 31743]`. The lower bound is the
/// conventional data floor (air/vacuum ≈ −1000 HU, −1024 the 12-bit floor); the
/// upper bound spans extended-range CT (16-bit signed offset). MVCT values fall
/// well inside this window.
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
#[repr(transparent)]
pub struct HounsfieldUnit(f64);

impl HounsfieldUnit {
    /// Inclusive lower bound of the representable HU range.
    pub const MIN: f64 = -1024.0;
    /// Inclusive upper bound of the representable HU range (extended CT).
    pub const MAX: f64 = 31_743.0;

    /// Returns the underlying CT number in HU.
    #[must_use]
    pub const fn get(self) -> f64 {
        self.0
    }
}

impl TryFrom<f64> for HounsfieldUnit {
    type Error = HeliosError;

    fn try_from(value: f64) -> Result<Self, Self::Error> {
        if !value.is_finite() {
            return Err(HeliosError::InvalidDomainValue {
                field: "HounsfieldUnit",
                value,
                reason: "CT number must be finite",
            });
        }
        if !(Self::MIN..=Self::MAX).contains(&value) {
            return Err(HeliosError::InvalidDomainValue {
                field: "HounsfieldUnit",
                value,
                reason: "CT number outside representable HU range [-1024, 31743]",
            });
        }
        Ok(Self(value))
    }
}

impl fmt::Display for HounsfieldUnit {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} HU", self.0)
    }
}

/// Voxel spacing along one axis, in millimetres (mm).
///
/// Valid range: finite and strictly positive. A zero or negative spacing would
/// make voxel volume ill-defined and break ray/voxel intersection in the
/// projectors, so it is rejected at construction.
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
#[repr(transparent)]
pub struct VoxelSpacingMm(f64);

impl VoxelSpacingMm {
    /// Returns the underlying spacing in millimetres.
    #[must_use]
    pub const fn get(self) -> f64 {
        self.0
    }
}

impl TryFrom<f64> for VoxelSpacingMm {
    type Error = HeliosError;

    fn try_from(value: f64) -> Result<Self, Self::Error> {
        if !value.is_finite() {
            return Err(HeliosError::InvalidDomainValue {
                field: "VoxelSpacingMm",
                value,
                reason: "spacing must be finite",
            });
        }
        if value <= 0.0 {
            return Err(HeliosError::InvalidDomainValue {
                field: "VoxelSpacingMm",
                value,
                reason: "spacing must be strictly positive",
            });
        }
        Ok(Self(value))
    }
}

impl fmt::Display for VoxelSpacingMm {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} mm", self.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn energy_accepts_positive_and_preserves_value() {
        let e = EnergyMeV::try_from(6.0).expect("6 MeV is valid");
        assert_eq!(e.get(), 6.0);
    }

    #[test]
    fn energy_rejects_zero_negative_and_nonfinite() {
        assert!(EnergyMeV::try_from(0.0).is_err());
        assert!(EnergyMeV::try_from(-1.0).is_err());
        assert!(EnergyMeV::try_from(f64::NAN).is_err());
        assert!(EnergyMeV::try_from(f64::INFINITY).is_err());
    }

    #[test]
    fn energy_error_names_field_and_reason() {
        let err = EnergyMeV::try_from(-2.5).unwrap_err();
        assert_eq!(
            err,
            HeliosError::InvalidDomainValue {
                field: "EnergyMeV",
                value: -2.5,
                reason: "energy must be strictly positive",
            }
        );
    }

    #[test]
    fn hounsfield_accepts_boundaries_and_rejects_outside() {
        assert_eq!(
            HounsfieldUnit::try_from(HounsfieldUnit::MIN)
                .expect("min is inclusive")
                .get(),
            HounsfieldUnit::MIN
        );
        assert_eq!(
            HounsfieldUnit::try_from(HounsfieldUnit::MAX)
                .expect("max is inclusive")
                .get(),
            HounsfieldUnit::MAX
        );
        assert!(HounsfieldUnit::try_from(HounsfieldUnit::MIN - 1.0).is_err());
        assert!(HounsfieldUnit::try_from(HounsfieldUnit::MAX + 1.0).is_err());
        assert!(HounsfieldUnit::try_from(f64::NAN).is_err());
    }

    #[test]
    fn hounsfield_accepts_typical_water_value() {
        // Liquid water is ~0 HU by CT calibration convention.
        assert_eq!(HounsfieldUnit::try_from(0.0).expect("water").get(), 0.0);
    }

    #[test]
    fn spacing_accepts_positive_and_rejects_nonpositive() {
        assert_eq!(
            VoxelSpacingMm::try_from(0.976_562_5)
                .expect("typical MVCT in-plane spacing")
                .get(),
            0.976_562_5
        );
        assert!(VoxelSpacingMm::try_from(0.0).is_err());
        assert!(VoxelSpacingMm::try_from(-0.5).is_err());
        assert!(VoxelSpacingMm::try_from(f64::NEG_INFINITY).is_err());
    }

    #[test]
    fn display_carries_unit_suffix() {
        assert_eq!(EnergyMeV::try_from(6.0).unwrap().to_string(), "6 MeV");
        assert_eq!(VoxelSpacingMm::try_from(1.5).unwrap().to_string(), "1.5 mm");
        assert_eq!(
            HounsfieldUnit::try_from(-1000.0).unwrap().to_string(),
            "-1000 HU"
        );
    }

    proptest::proptest! {
        /// Any finite positive value round-trips through `EnergyMeV` unchanged.
        #[test]
        fn energy_roundtrips_finite_positive(v in 1e-6_f64..1e4_f64) {
            let e = EnergyMeV::try_from(v).expect("finite positive is valid");
            proptest::prop_assert_eq!(e.get(), v);
        }

        /// Every in-range HU value is accepted and preserved; the newtype never
        /// perturbs the stored number.
        #[test]
        fn hounsfield_preserves_in_range(v in HounsfieldUnit::MIN..=HounsfieldUnit::MAX) {
            let hu = HounsfieldUnit::try_from(v).expect("in range");
            proptest::prop_assert_eq!(hu.get(), v);
        }

        /// `EnergyMeV` rejects every non-positive value (the strict-positivity
        /// invariant holds over the whole ≤0 range, not just sampled points).
        #[test]
        fn energy_rejects_all_non_positive(v in -1e6_f64..=0.0_f64) {
            proptest::prop_assert!(EnergyMeV::try_from(v).is_err());
        }

        /// `HounsfieldUnit` rejects any value strictly above the representable max.
        #[test]
        fn hounsfield_rejects_above_max(v in (HounsfieldUnit::MAX + 1.0)..1e9_f64) {
            proptest::prop_assert!(HounsfieldUnit::try_from(v).is_err());
        }

        /// `HounsfieldUnit`'s `PartialOrd` mirrors the underlying numeric order:
        /// `a ≤ b` ⇒ `HU(a) ≤ HU(b)`. Preserving order matters for windowing and
        /// thresholding on CT numbers.
        #[test]
        fn hounsfield_preserves_ordering(
            a in HounsfieldUnit::MIN..=HounsfieldUnit::MAX,
            b in HounsfieldUnit::MIN..=HounsfieldUnit::MAX,
        ) {
            let (ha, hb) = (
                HounsfieldUnit::try_from(a).unwrap(),
                HounsfieldUnit::try_from(b).unwrap(),
            );
            proptest::prop_assert_eq!(a <= b, ha <= hb);
        }

        /// Any finite positive spacing round-trips through `VoxelSpacingMm`.
        #[test]
        fn spacing_roundtrips_finite_positive(v in 1e-6_f64..1e3_f64) {
            let s = VoxelSpacingMm::try_from(v).expect("finite positive is valid");
            proptest::prop_assert_eq!(s.get(), v);
        }
    }
}
