//! Helical TomoTherapy delivery kinematics.
//!
//! In helical delivery the gantry rotates continuously while the couch
//! translates the patient through the bore, so the source traces a helix in the
//! patient frame. The rotation is discretized into a fixed number of
//! **projections** per rotation (51 on TomoTherapy), each with its own binary-MLC
//! leaf pattern.
//!
//! The couch advance per gantry rotation is set by the **pitch** — the couch
//! travel per rotation expressed in units of the field width (the jaw opening at
//! isocentre):
//!
//! ```text
//! pitch = couch_travel_per_rotation / field_width
//! ```
//!
//! This module provides the deterministic mapping from projection index (or time)
//! to gantry angle and couch position — the "helical synchronization" that the
//! delivery simulation and MVCT acquisition are driven by.

use aequitas::systems::si::{
    quantities::{Angle, Dimensionless, Length, Time, Velocity},
    units::{Millimeter, Radian, Second},
};
use helios_core::HeliosError;
use helios_math::{NumericElement, Scalar};

/// Helical delivery geometry and timing.
///
/// Construct with [`HelicalDelivery::new`]; the mapping methods are pure
/// functions of the projection index or elapsed time.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HelicalDelivery<T: Scalar> {
    projections_per_rotation: usize,
    field_width_mm: Length<T>,
    pitch: Dimensionless<T>,
    gantry_period_s: Time<T>,
    start_gantry_angle_rad: Angle<T>,
    start_couch_mm: Length<T>,
}

impl<T: Scalar> HelicalDelivery<T> {
    /// Construct a helical delivery.
    ///
    /// - `projections_per_rotation`: gantry projections per full rotation (51 on
    ///   TomoTherapy).
    /// - `field_width_mm`: jaw opening at isocentre (e.g. 25 mm).
    /// - `pitch`: couch travel per rotation ÷ field width (typically 0.2–0.5).
    /// - `gantry_period_s`: time for one full gantry rotation.
    /// - `start_gantry_angle_rad` / `start_couch_mm`: pose at projection 0.
    ///
    /// # Errors
    /// Returns [`HeliosError::InvalidDomainValue`] if `projections_per_rotation`
    /// is zero, or any of `field_width_mm`, `pitch`, `gantry_period_s` is
    /// non-finite or not strictly positive.
    pub fn new(
        projections_per_rotation: usize,
        field_width_mm: Length<T>,
        pitch: Dimensionless<T>,
        gantry_period_s: Time<T>,
        start_gantry_angle_rad: Angle<T>,
        start_couch_mm: Length<T>,
    ) -> Result<Self, HeliosError> {
        if projections_per_rotation == 0 {
            return Err(HeliosError::InvalidDomainValue {
                field: "HelicalDelivery::projections_per_rotation",
                value: 0.0,
                reason: "must be non-zero",
            });
        }
        for (value, field) in [
            (
                field_width_mm.in_unit::<Millimeter>(),
                "HelicalDelivery::field_width_mm",
            ),
            (pitch.into_base(), "HelicalDelivery::pitch"),
            (
                gantry_period_s.in_unit::<Second>(),
                "HelicalDelivery::gantry_period_s",
            ),
        ] {
            if !value.is_finite() || value <= <T as NumericElement>::ZERO {
                return Err(HeliosError::InvalidDomainValue {
                    field,
                    value: value.to_f64(),
                    reason: "must be finite and strictly positive",
                });
            }
        }
        if !start_gantry_angle_rad.in_unit::<Radian>().is_finite()
            || !start_couch_mm.in_unit::<Millimeter>().is_finite()
        {
            return Err(HeliosError::InvalidDomainValue {
                field: "HelicalDelivery::start_pose",
                value: start_gantry_angle_rad.in_unit::<Radian>().to_f64(),
                reason: "start angle and couch position must be finite",
            });
        }
        Ok(Self {
            projections_per_rotation,
            field_width_mm,
            pitch,
            gantry_period_s,
            start_gantry_angle_rad,
            start_couch_mm,
        })
    }

    /// Projections per full gantry rotation.
    #[must_use]
    pub const fn projections_per_rotation(&self) -> usize {
        self.projections_per_rotation
    }

    /// Pitch (couch travel per rotation ÷ field width).
    #[must_use]
    pub fn pitch(&self) -> Dimensionless<T> {
        self.pitch
    }

    /// Field width at isocentre (mm).
    #[must_use]
    pub fn field_width_mm(&self) -> Length<T> {
        self.field_width_mm
    }

    /// Couch travel per full gantry rotation (mm): `pitch · field_width`.
    #[must_use]
    pub fn couch_travel_per_rotation_mm(&self) -> Length<T> {
        Length::from_base(self.pitch.into_base() * self.field_width_mm.into_base())
    }

    /// Couch advance per projection (mm).
    #[must_use]
    pub fn couch_advance_per_projection_mm(&self) -> Length<T> {
        Length::from_base(
            self.couch_travel_per_rotation_mm().into_base()
                * self.projections_per_rotation_recip().into_base(),
        )
    }

    /// Constant couch velocity (mm/s): couch travel per rotation ÷ gantry period.
    #[must_use]
    pub fn couch_velocity_mm_per_s(&self) -> Velocity<T> {
        Velocity::from_base(
            self.couch_travel_per_rotation_mm().into_base() / self.gantry_period_s.into_base(),
        )
    }

    /// Gantry angle (rad) at a projection index, unwrapped (monotonically
    /// increasing across rotations).
    #[must_use]
    pub fn gantry_angle_rad(&self, projection: usize) -> Angle<T> {
        Angle::from_base(
            self.start_gantry_angle_rad.in_unit::<Radian>()
                + T::TAU
                    * T::from_f64(projection as f64)
                    * self.projections_per_rotation_recip().into_base(),
        )
    }

    /// Gantry angle wrapped into `[0, 2π)`.
    #[must_use]
    pub fn gantry_angle_wrapped_rad(&self, projection: usize) -> Angle<T> {
        let angle = self.gantry_angle_rad(projection).in_unit::<Radian>();
        let turns = (angle * T::TAU.recip()).floor();
        Angle::from_base(angle - T::TAU * turns)
    }

    /// Couch position (mm) at a projection index.
    #[must_use]
    pub fn couch_position_mm(&self, projection: usize) -> Length<T> {
        Length::from_base(
            self.start_couch_mm.into_base()
                + T::from_f64(projection as f64)
                    * self.couch_advance_per_projection_mm().into_base(),
        )
    }

    /// Elapsed time (s) at a projection index.
    #[must_use]
    pub fn time_s(&self, projection: usize) -> Time<T> {
        Time::from_base(
            T::from_f64(projection as f64)
                * self.gantry_period_s.in_unit::<Second>()
                * self.projections_per_rotation_recip().into_base(),
        )
    }

    /// Gantry angle (rad) at continuous time `t` (s), unwrapped.
    #[must_use]
    pub fn gantry_angle_at_time_rad(&self, t: Time<T>) -> Angle<T> {
        Angle::from_base(
            self.start_gantry_angle_rad.in_unit::<Radian>()
                + T::TAU * t.in_unit::<Second>() / self.gantry_period_s.in_unit::<Second>(),
        )
    }

    /// Couch position (mm) at continuous time `t` (s).
    #[must_use]
    pub fn couch_position_at_time_mm(&self, t: Time<T>) -> Length<T> {
        Length::from_base(
            self.start_couch_mm.into_base()
                + self.couch_velocity_mm_per_s().into_base() * t.into_base(),
        )
    }

    #[inline]
    fn projections_per_rotation_recip(&self) -> Dimensionless<T> {
        Dimensionless::from_base(T::from_f64(self.projections_per_rotation as f64).recip())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use eunomia::assert_relative_eq;

    fn angle(value: f64) -> Angle<f64> {
        Angle::from_unit::<Radian>(value)
    }

    fn delivery() -> HelicalDelivery<f64> {
        // 51 projections/rotation, 25 mm field, pitch 0.4, 10 s/rotation.
        HelicalDelivery::new(
            51,
            Length::from_unit::<Millimeter>(25.0),
            Dimensionless::from_base(0.4),
            Time::from_unit::<Second>(10.0),
            angle(0.0),
            Length::from_unit::<Millimeter>(0.0),
        )
        .expect("valid delivery")
    }

    #[test]
    fn rejects_invalid_parameters() {
        assert!(HelicalDelivery::new(
            0,
            Length::from_unit::<Millimeter>(25.0),
            Dimensionless::from_base(0.4),
            Time::from_unit::<Second>(10.0),
            angle(0.0),
            Length::from_unit::<Millimeter>(0.0),
        )
        .is_err());
        assert!(HelicalDelivery::new(
            51,
            Length::from_unit::<Millimeter>(0.0),
            Dimensionless::from_base(0.4),
            Time::from_unit::<Second>(10.0),
            angle(0.0),
            Length::from_unit::<Millimeter>(0.0),
        )
        .is_err());
        assert!(HelicalDelivery::new(
            51,
            Length::from_unit::<Millimeter>(25.0),
            Dimensionless::from_base(-0.4),
            Time::from_unit::<Second>(10.0),
            angle(0.0),
            Length::from_unit::<Millimeter>(0.0),
        )
        .is_err());
        assert!(HelicalDelivery::new(
            51,
            Length::from_unit::<Millimeter>(25.0),
            Dimensionless::from_base(0.4),
            Time::from_unit::<Second>(f64::NAN),
            angle(0.0),
            Length::from_unit::<Millimeter>(0.0),
        )
        .is_err());
    }

    #[test]
    fn pitch_relation_holds() {
        let d = delivery();
        // couch travel per rotation ÷ field width == pitch.
        assert_relative_eq!(
            d.couch_travel_per_rotation_mm().in_unit::<Millimeter>()
                / d.field_width_mm().in_unit::<Millimeter>(),
            d.pitch().into_base(),
            epsilon = 1e-15
        );
        // pitch 0.4 × 25 mm = 10 mm per rotation.
        assert_relative_eq!(
            d.couch_travel_per_rotation_mm().in_unit::<Millimeter>(),
            10.0,
            epsilon = 1e-13
        );
    }

    #[test]
    fn one_full_rotation_advances_angle_by_tau_and_couch_by_travel() {
        let d = delivery();
        let ppr = d.projections_per_rotation();
        // After exactly one rotation (projection = ppr): angle += 2π.
        assert_relative_eq!(
            d.gantry_angle_rad(ppr).in_unit::<Radian>() - d.gantry_angle_rad(0).in_unit::<Radian>(),
            core::f64::consts::TAU,
            epsilon = 1e-13
        );
        // Wrapped angle returns to the start.
        assert_relative_eq!(
            d.gantry_angle_wrapped_rad(ppr).in_unit::<Radian>(),
            0.0,
            epsilon = 1e-12
        );
        // Couch advanced by exactly one rotation's travel.
        assert_relative_eq!(
            d.couch_position_mm(ppr).in_unit::<Millimeter>()
                - d.couch_position_mm(0).in_unit::<Millimeter>(),
            d.couch_travel_per_rotation_mm().in_unit::<Millimeter>(),
            epsilon = 1e-12
        );
        assert_relative_eq!(
            d.couch_position_mm(ppr).in_unit::<Millimeter>(),
            10.0,
            epsilon = 1e-12
        );
    }

    #[test]
    fn half_rotation_is_pi() {
        // A start angle plus half of 51 projections... use an even count for an
        // exact half. 50 projections/rotation → projection 25 is a half-turn.
        let d = HelicalDelivery::new(
            50,
            Length::from_unit::<Millimeter>(25.0),
            Dimensionless::from_base(0.4),
            Time::from_unit::<Second>(10.0),
            angle(0.0),
            Length::from_unit::<Millimeter>(0.0),
        )
        .unwrap();
        assert_relative_eq!(
            d.gantry_angle_rad(25).in_unit::<Radian>(),
            core::f64::consts::PI,
            epsilon = 1e-13
        );
    }

    #[test]
    fn projection_and_time_parameterizations_agree() {
        let d = delivery();
        for p in [0usize, 1, 17, 51, 102] {
            let t = d.time_s(p);
            assert_relative_eq!(
                d.gantry_angle_rad(p).in_unit::<Radian>(),
                d.gantry_angle_at_time_rad(t).in_unit::<Radian>(),
                epsilon = 1e-12
            );
            assert_relative_eq!(
                d.couch_position_mm(p).in_unit::<Millimeter>(),
                d.couch_position_at_time_mm(t).in_unit::<Millimeter>(),
                epsilon = 1e-12
            );
        }
        // One rotation takes exactly the gantry period.
        assert_relative_eq!(d.time_s(51).in_unit::<Second>(), 10.0, epsilon = 1e-13);
    }

    #[test]
    fn couch_advance_is_monotonic() {
        let d = delivery();
        let mut prev = d.couch_position_mm(0).in_unit::<Millimeter>();
        for p in 1..200 {
            let z = d.couch_position_mm(p).in_unit::<Millimeter>();
            assert!(z > prev, "couch must advance monotonically");
            prev = z;
        }
    }

    #[test]
    fn kinematics_are_generic_over_scalar_f32() {
        let d = HelicalDelivery::<f32>::new(
            51,
            Length::from_unit::<Millimeter>(25.0),
            Dimensionless::from_base(0.4),
            Time::from_unit::<Second>(10.0),
            Angle::from_unit::<Radian>(0.0),
            Length::from_unit::<Millimeter>(0.0),
        )
        .unwrap();
        assert_relative_eq!(
            d.couch_travel_per_rotation_mm().in_unit::<Millimeter>(),
            10.0_f32,
            epsilon = 1e-4
        );
        assert_relative_eq!(
            d.gantry_angle_rad(51).in_unit::<Radian>(),
            core::f32::consts::TAU,
            epsilon = 1e-4
        );
    }
}
