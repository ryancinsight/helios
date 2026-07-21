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

use helios_core::HeliosError;
use helios_math::{NumericElement, Scalar};

/// Helical delivery geometry and timing.
///
/// Construct with [`HelicalDelivery::new`]; the mapping methods are pure
/// functions of the projection index or elapsed time.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HelicalDelivery<T: Scalar> {
    projections_per_rotation: usize,
    field_width_mm: T,
    pitch: T,
    gantry_period_s: T,
    start_gantry_angle_rad: T,
    start_couch_mm: T,
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
        field_width_mm: T,
        pitch: T,
        gantry_period_s: T,
        start_gantry_angle_rad: T,
        start_couch_mm: T,
    ) -> Result<Self, HeliosError> {
        if projections_per_rotation == 0 {
            return Err(HeliosError::InvalidDomainValue {
                field: "HelicalDelivery::projections_per_rotation",
                value: 0.0,
                reason: "must be non-zero",
            });
        }
        for (value, field) in [
            (field_width_mm, "HelicalDelivery::field_width_mm"),
            (pitch, "HelicalDelivery::pitch"),
            (gantry_period_s, "HelicalDelivery::gantry_period_s"),
        ] {
            if !value.is_finite() || value <= <T as NumericElement>::ZERO {
                return Err(HeliosError::InvalidDomainValue {
                    field,
                    value: value.to_f64(),
                    reason: "must be finite and strictly positive",
                });
            }
        }
        if !start_gantry_angle_rad.is_finite() || !start_couch_mm.is_finite() {
            return Err(HeliosError::InvalidDomainValue {
                field: "HelicalDelivery::start_pose",
                value: start_gantry_angle_rad.to_f64(),
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
    pub fn pitch(&self) -> T {
        self.pitch
    }

    /// Field width at isocentre (mm).
    #[must_use]
    pub fn field_width_mm(&self) -> T {
        self.field_width_mm
    }

    /// Couch travel per full gantry rotation (mm): `pitch · field_width`.
    #[must_use]
    pub fn couch_travel_per_rotation_mm(&self) -> T {
        self.pitch * self.field_width_mm
    }

    /// Couch advance per projection (mm).
    #[must_use]
    pub fn couch_advance_per_projection_mm(&self) -> T {
        self.couch_travel_per_rotation_mm() * self.projections_per_rotation_recip()
    }

    /// Constant couch velocity (mm/s): couch travel per rotation ÷ gantry period.
    #[must_use]
    pub fn couch_velocity_mm_per_s(&self) -> T {
        self.couch_travel_per_rotation_mm() * self.gantry_period_s.recip()
    }

    /// Gantry angle (rad) at a projection index, unwrapped (monotonically
    /// increasing across rotations).
    #[must_use]
    pub fn gantry_angle_rad(&self, projection: usize) -> T {
        self.start_gantry_angle_rad
            + T::TAU * T::from_f64(projection as f64) * self.projections_per_rotation_recip()
    }

    /// Gantry angle wrapped into `[0, 2π)`.
    #[must_use]
    pub fn gantry_angle_wrapped_rad(&self, projection: usize) -> T {
        let angle = self.gantry_angle_rad(projection);
        let turns = (angle * T::TAU.recip()).floor();
        angle - T::TAU * turns
    }

    /// Couch position (mm) at a projection index.
    #[must_use]
    pub fn couch_position_mm(&self, projection: usize) -> T {
        self.start_couch_mm
            + T::from_f64(projection as f64) * self.couch_advance_per_projection_mm()
    }

    /// Elapsed time (s) at a projection index.
    #[must_use]
    pub fn time_s(&self, projection: usize) -> T {
        T::from_f64(projection as f64)
            * self.gantry_period_s
            * self.projections_per_rotation_recip()
    }

    /// Gantry angle (rad) at continuous time `t` (s), unwrapped.
    #[must_use]
    pub fn gantry_angle_at_time_rad(&self, t: T) -> T {
        self.start_gantry_angle_rad + T::TAU * t * self.gantry_period_s.recip()
    }

    /// Couch position (mm) at continuous time `t` (s).
    #[must_use]
    pub fn couch_position_at_time_mm(&self, t: T) -> T {
        self.start_couch_mm + self.couch_velocity_mm_per_s() * t
    }

    #[inline]
    fn projections_per_rotation_recip(&self) -> T {
        T::from_f64(self.projections_per_rotation as f64).recip()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use eunomia::assert_relative_eq;

    fn delivery() -> HelicalDelivery<f64> {
        // 51 projections/rotation, 25 mm field, pitch 0.4, 10 s/rotation.
        HelicalDelivery::new(51, 25.0, 0.4, 10.0, 0.0, 0.0).expect("valid delivery")
    }

    #[test]
    fn rejects_invalid_parameters() {
        assert!(HelicalDelivery::new(0, 25.0, 0.4, 10.0, 0.0, 0.0).is_err());
        assert!(HelicalDelivery::new(51, 0.0, 0.4, 10.0, 0.0, 0.0).is_err());
        assert!(HelicalDelivery::new(51, 25.0, -0.4, 10.0, 0.0, 0.0).is_err());
        assert!(HelicalDelivery::new(51, 25.0, 0.4, f64::NAN, 0.0, 0.0).is_err());
    }

    #[test]
    fn pitch_relation_holds() {
        let d = delivery();
        // couch travel per rotation ÷ field width == pitch.
        assert_relative_eq!(
            d.couch_travel_per_rotation_mm() / d.field_width_mm(),
            d.pitch(),
            epsilon = 1e-15
        );
        // pitch 0.4 × 25 mm = 10 mm per rotation.
        assert_relative_eq!(d.couch_travel_per_rotation_mm(), 10.0, epsilon = 1e-13);
    }

    #[test]
    fn one_full_rotation_advances_angle_by_tau_and_couch_by_travel() {
        let d = delivery();
        let ppr = d.projections_per_rotation();
        // After exactly one rotation (projection = ppr): angle += 2π.
        assert_relative_eq!(
            d.gantry_angle_rad(ppr) - d.gantry_angle_rad(0),
            core::f64::consts::TAU,
            epsilon = 1e-13
        );
        // Wrapped angle returns to the start.
        assert_relative_eq!(d.gantry_angle_wrapped_rad(ppr), 0.0, epsilon = 1e-12);
        // Couch advanced by exactly one rotation's travel.
        assert_relative_eq!(
            d.couch_position_mm(ppr) - d.couch_position_mm(0),
            d.couch_travel_per_rotation_mm(),
            epsilon = 1e-12
        );
    }

    #[test]
    fn half_rotation_is_pi() {
        // A start angle plus half of 51 projections... use an even count for an
        // exact half. 50 projections/rotation → projection 25 is a half-turn.
        let d = HelicalDelivery::new(50, 25.0, 0.4, 10.0, 0.0, 0.0).unwrap();
        assert_relative_eq!(
            d.gantry_angle_rad(25),
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
                d.gantry_angle_rad(p),
                d.gantry_angle_at_time_rad(t),
                epsilon = 1e-12
            );
            assert_relative_eq!(
                d.couch_position_mm(p),
                d.couch_position_at_time_mm(t),
                epsilon = 1e-12
            );
        }
        // One rotation takes exactly the gantry period.
        assert_relative_eq!(d.time_s(51), 10.0, epsilon = 1e-13);
    }

    #[test]
    fn couch_advance_is_monotonic() {
        let d = delivery();
        let mut prev = d.couch_position_mm(0);
        for p in 1..200 {
            let z = d.couch_position_mm(p);
            assert!(z > prev, "couch must advance monotonically");
            prev = z;
        }
    }

    #[test]
    fn kinematics_are_generic_over_scalar_f32() {
        let d = HelicalDelivery::<f32>::new(51, 25.0, 0.4, 10.0, 0.0, 0.0).unwrap();
        assert_relative_eq!(d.couch_travel_per_rotation_mm(), 10.0_f32, epsilon = 1e-4);
        assert_relative_eq!(
            d.gantry_angle_rad(51),
            core::f32::consts::TAU,
            epsilon = 1e-4
        );
    }
}
