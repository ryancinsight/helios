//! Collimator field aperture: the jaw-defined open field with geometric edge
//! penumbra, over a gaia [`Aabb`].
//!
//! A [`FieldAperture`] models the secondary collimator (jaws) that shape the
//! open field in **collimator coordinates** — lateral leaf offset and couch axis
//! — independent of gantry angle (the collimator rotates with the gantry). Its
//! open region is a gaia `Aabb`; [`transmission`](FieldAperture::transmission)
//! is `1` deep inside, `0` deep outside, and `0.5` on the geometric edge,
//! ramping linearly across a penumbra band (the finite-source field-edge
//! blur). This is the field-shaping factor a delivery applies to each beamlet's
//! fluence.

use aequitas::systems::si::quantities::Length;
use helios_core::HeliosError;
use helios_math::{Aabb, GeometryScalar, NumericElement, Point3};

/// A rectangular collimator field aperture (a gaia `Aabb` open region) with a
/// linear geometric penumbra at its edges.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FieldAperture<T: GeometryScalar> {
    open: Aabb<T>,
    penumbra: Length<T>,
}

/// `|x|` for a `GeometryScalar`, via the ordered field (exact, no `f64` round-trip).
#[inline]
fn abs<T: GeometryScalar>(x: T) -> T {
    if x < <T as NumericElement>::ZERO {
        -x
    } else {
        x
    }
}

/// `max(a, b)` via the ordered field.
#[inline]
fn max2<T: GeometryScalar>(a: T, b: T) -> T {
    if a > b {
        a
    } else {
        b
    }
}

impl<T: GeometryScalar> FieldAperture<T> {
    /// Build an aperture from an open-field gaia `Aabb` and a positive penumbra
    /// half-width (mm).
    ///
    /// # Errors
    /// [`HeliosError::InvalidDomainValue`] if `penumbra` is not finite and
    /// positive, or the box is degenerate (`min > max` on any axis).
    pub fn new(open: Aabb<T>, penumbra: Length<T>) -> Result<Self, HeliosError> {
        let penumbra_m = penumbra.into_base();
        if !(penumbra_m.to_f64().is_finite() && penumbra_m > <T as NumericElement>::ZERO) {
            return Err(HeliosError::InvalidDomainValue {
                field: "FieldAperture::penumbra",
                value: penumbra_m.to_f64(),
                reason: "penumbra half-width must be finite and positive",
            });
        }
        if open.min.x > open.max.x || open.min.y > open.max.y || open.min.z > open.max.z {
            return Err(HeliosError::InvalidDomainValue {
                field: "FieldAperture::open",
                value: 0.0,
                reason: "aperture box has min greater than max on some axis",
            });
        }
        Ok(Self { open, penumbra })
    }

    /// Build a rectangular aperture centred at `centre` with per-axis half-widths
    /// `half` (mm) and typed penumbra width.
    ///
    /// # Errors
    /// As [`new`](Self::new); also if any half-width is negative.
    pub fn rectangular(
        centre: Point3<T>,
        half: [T; 3],
        penumbra: Length<T>,
    ) -> Result<Self, HeliosError> {
        let zero = <T as NumericElement>::ZERO;
        if half.iter().any(|&h| h < zero) {
            return Err(HeliosError::InvalidDomainValue {
                field: "FieldAperture::rectangular::half",
                value: 0.0,
                reason: "aperture half-widths must be non-negative",
            });
        }
        let min = Point3::new(centre.x - half[0], centre.y - half[1], centre.z - half[2]);
        let max = Point3::new(centre.x + half[0], centre.y + half[1], centre.z + half[2]);
        Self::new(Aabb::new(min, max), penumbra)
    }

    /// The open-field gaia `Aabb`.
    #[must_use]
    pub fn open(&self) -> &Aabb<T> {
        &self.open
    }

    /// Whether `p` is strictly inside the fully-open core (delegates to gaia
    /// [`Aabb::contains_point`]). Points in the penumbra band or outside are not
    /// "contained" even though they may transmit partially.
    #[must_use]
    pub fn contains(&self, p: &Point3<T>) -> bool {
        self.open.contains_point(p)
    }

    /// Signed distance from `p` to the open box (the standard AABB SDF): negative
    /// inside, `0` on the boundary, positive outside.
    fn signed_distance(&self, p: &Point3<T>) -> T {
        let zero = <T as NumericElement>::ZERO;
        let c = self.open.center();
        let hx = (self.open.max.x - self.open.min.x) * <T as GeometryScalar>::from_f64(0.5);
        let hy = (self.open.max.y - self.open.min.y) * <T as GeometryScalar>::from_f64(0.5);
        let hz = (self.open.max.z - self.open.min.z) * <T as GeometryScalar>::from_f64(0.5);
        // q_i = |p_i − c_i| − h_i : >0 outside that axis's slab, <0 inside.
        let q = [
            abs(p.x - c.x) - hx,
            abs(p.y - c.y) - hy,
            abs(p.z - c.z) - hz,
        ];
        let pos = |v: T| if v > zero { v } else { zero };
        let outside_sq = pos(q[0]) * pos(q[0]) + pos(q[1]) * pos(q[1]) + pos(q[2]) * pos(q[2]);
        let outside = outside_sq.sqrt();
        // Inside distance: the least-negative q (nearest face), clamped ≤ 0.
        let max_q = max2(max2(q[0], q[1]), q[2]);
        let inside = if max_q < zero { max_q } else { zero };
        outside + inside
    }

    /// Beamlet transmission at `p` (collimator coordinates): `1` deep inside the
    /// field, `0` deep outside, `0.5` on the geometric edge, ramping linearly
    /// across the `±penumbra` band. Always in `[0, 1]`.
    #[must_use]
    pub fn transmission(&self, p: &Point3<T>) -> T {
        let half = <T as GeometryScalar>::from_f64(0.5);
        let two = <T as GeometryScalar>::from_f64(2.0);
        let sdf = self.signed_distance(p);
        // Aequitas stores SI base metres; collimator coordinates in this
        // geometry module remain millimetres.
        let penumbra_mm = self.penumbra.into_base() * <T as GeometryScalar>::from_f64(1.0e3);
        // 0.5 at sdf=0; +0.5 per penumbra inside; −0.5 per penumbra outside.
        let t = half - sdf * (two * penumbra_mm).recip();
        let zero = <T as NumericElement>::ZERO;
        let one = <T as NumericElement>::ONE;
        if t < zero {
            zero
        } else if t > one {
            one
        } else {
            t
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use helios_math::ShippedScalar;
    use eunomia::assert_relative_eq;

    // 20 mm × 20 mm × 20 mm field centred at the origin, 2 mm penumbra.
    fn aperture() -> FieldAperture<f64> {
        FieldAperture::rectangular(
            Point3::new(0.0, 0.0, 0.0),
            [10.0, 10.0, 10.0],
            Length::from_base(2.0e-3),
        )
        .unwrap()
    }

    #[test]
    fn centre_transmits_fully_and_far_outside_blocks() {
        let a = aperture();
        assert_relative_eq!(
            a.transmission(&Point3::new(0.0, 0.0, 0.0)),
            1.0,
            epsilon = 1e-12
        );
        assert_relative_eq!(
            a.transmission(&Point3::new(100.0, 0.0, 0.0)),
            0.0,
            epsilon = 1e-12
        );
    }

    #[test]
    fn geometric_edge_is_half_and_penumbra_band_ramps() {
        let a = aperture();
        // On the +x face (x = 10): sdf = 0 → 50 %.
        assert_relative_eq!(
            a.transmission(&Point3::new(10.0, 0.0, 0.0)),
            0.5,
            epsilon = 1e-12
        );
        // One penumbra (2 mm) inside the face → ~100 %; one outside → ~0 %.
        assert_relative_eq!(
            a.transmission(&Point3::new(8.0, 0.0, 0.0)),
            1.0,
            epsilon = 1e-12
        );
        assert_relative_eq!(
            a.transmission(&Point3::new(12.0, 0.0, 0.0)),
            0.0,
            epsilon = 1e-12
        );
        // Halfway through the band → 75 % / 25 %.
        assert_relative_eq!(
            a.transmission(&Point3::new(9.0, 0.0, 0.0)),
            0.75,
            epsilon = 1e-12
        );
        assert_relative_eq!(
            a.transmission(&Point3::new(11.0, 0.0, 0.0)),
            0.25,
            epsilon = 1e-12
        );
    }

    #[test]
    fn transmission_is_monotone_across_the_edge_and_bounded() {
        let a = aperture();
        let mut prev = 1.0;
        for step in 0..=40 {
            let x = 6.0 + step as f64 * 0.2; // 6 → 14 mm, crossing the 10 mm edge
            let t = a.transmission(&Point3::new(x, 0.0, 0.0));
            assert!((0.0..=1.0).contains(&t), "transmission {t} out of [0,1]");
            assert!(t <= prev + 1e-12, "transmission must not increase outward");
            prev = t;
        }
    }

    #[test]
    fn contains_matches_the_gaia_core_and_rejects_penumbra() {
        let a = aperture();
        assert!(a.contains(&Point3::new(0.0, 0.0, 0.0)));
        assert!(!a.contains(&Point3::new(11.0, 0.0, 0.0))); // in penumbra, not core
    }

    #[test]
    fn invalid_penumbra_and_box_are_typed_errors() {
        let unit = Point3::new(0.0, 0.0, 0.0);
        assert!(
            FieldAperture::rectangular(unit, [10.0, 10.0, 10.0], Length::from_base(0.0)).is_err()
        );
        assert!(
            FieldAperture::rectangular(unit, [10.0, 10.0, 10.0], Length::from_base(-1.0e-3),)
                .is_err()
        );
        assert!(FieldAperture::new(
            Aabb::new(Point3::new(5.0, 0.0, 0.0), Point3::new(-5.0, 1.0, 1.0)),
            Length::from_base(2.0e-3)
        )
        .is_err());
    }

    /// Aperture transmission evaluates a penumbra profile. The three sampled
    /// points are its plateau, its half-value edge, and its cut-off; only the
    /// edge carries evaluation error, and sixteen ulps of `T` bounds it. The
    /// cut-off is compared against zero, where a relative bound is vacuous, so
    /// that assertion keeps an absolute epsilon derived from the same budget.
    const APERTURE_ULPS: f64 = 16.0;

    /// Asserts aperture plateau, half-edge, and cut-off in one scalar width.
    fn aperture_transmission_falls_from_plateau_to_cutoff<T>()
    where
        T: ShippedScalar + helios_math::GeometryScalar,
    {
        // GeometryScalar and FloatElement both define `from_f64`; name the
        // one that performs the literal conversion so the call is unambiguous.
        let cast = <T as helios_math::FloatElement>::from_f64;
        let zero = cast(0.0);
        let half_width = cast(5.0);
        let aperture = FieldAperture::rectangular(
            Point3::new(zero, zero, zero),
            [half_width; 3],
            Length::from_base(cast(1.0e-3)),
        )
        .expect("positive aperture extent and penumbra");
        let tolerance = T::EPSILON * cast(APERTURE_ULPS);

        // Centre: fully open.
        assert_relative_eq!(
            aperture.transmission(&Point3::new(zero, zero, zero)),
            cast(1.0),
            max_relative = tolerance
        );
        // On the edge: half transmission by symmetry of the penumbra.
        assert_relative_eq!(
            aperture.transmission(&Point3::new(half_width, zero, zero)),
            cast(0.5),
            max_relative = tolerance
        );
        // Far outside: fully blocked. Relative tolerance cannot express nearness
        // to zero, so this one is absolute.
        assert_relative_eq!(
            aperture.transmission(&Point3::new(cast(20.0), zero, zero)),
            zero,
            epsilon = tolerance
        );
    }

    #[test]
    fn aperture_transmission_falls_from_plateau_to_cutoff_in_single_precision() {
        aperture_transmission_falls_from_plateau_to_cutoff::<f32>();
    }

    #[test]
    fn aperture_transmission_falls_from_plateau_to_cutoff_in_double_precision() {
        aperture_transmission_falls_from_plateau_to_cutoff::<f64>();
    }
}
