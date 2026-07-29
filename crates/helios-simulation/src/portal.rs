//! Portal (EPID) exit dosimetry: the transmitted primary fluence per MLC leaf.
//!
//! An electronic portal imaging device behind the patient measures the fluence
//! that survives attenuation, `Ψ_leaf · exp(−τ_leaf)`, for each open leaf of a
//! delivery frame — the signal used to *verify* that the delivered fluence matches
//! the plan. This composes the per-leaf beamlet geometry (shared with dose
//! accumulation) and the [`forward_project_ray`](helios_solver::forward_project_ray)
//! optical depth, so a closed leaf reads 0, an unattenuated leaf reads its full
//! fluence, and attenuation multiplies each reading by `exp(−τ)`.

use crate::delivery::DeliveryFrame;
use crate::dose_accumulation::{BeamGeometry, beamlet_ray, gantry_basis};
use aequitas::systems::si::{
    quantities::{Dimensionless, EnergyPerArea, Length},
    units::Millimeter,
};
use eunomia::UnitScalar;
use helios_domain::Volume;
use helios_math::{GeometryScalar, NumericElement};
use helios_solver::forward_project_ray;
use hyperion::{TransportError, quantity::OpticalDepth};

#[cfg(test)]
use aequitas::systems::si::{quantities::Angle, units::Radian};

/// Portal exit fluence per MLC leaf for one delivery `frame` through `mu`.
///
/// Returns a vector aligned with `frame.leaf_fluence`: entry `l` is the delivered
/// leaf fluence attenuated by the beamlet's optical depth,
/// `leaf_fluence[l] · exp(−τ_l)`. `geometry`/`leaf_width`/`step` are as in
/// [`accumulate_delivered_dose`](crate::accumulate_delivered_dose). A leaf whose
/// beamlet misses the volume reads its full (unattenuated) fluence.
///
/// # Errors
///
/// Returns [`TransportError`] if a beamlet produces a negative or non-finite
/// optical depth.
pub fn frame_portal_fluence<T: GeometryScalar + UnitScalar>(
    frame: &DeliveryFrame<T>,
    mu: &Volume<T>,
    geometry: BeamGeometry<T>,
    leaf_width: Length<T>,
    step: Length<T>,
) -> Result<Vec<EnergyPerArea<T>>, TransportError<T>> {
    let zero = <T as NumericElement>::ZERO;
    let step_mm = step.in_unit::<Millimeter>();
    let (centre, dir, perp) = gantry_basis(mu.grid(), frame.gantry_angle_rad);
    frame
        .leaf_fluence
        .iter()
        .enumerate()
        .map(|(leaf, fluence)| {
            let fluence_base = *fluence.as_base();
            if fluence_base <= zero {
                return Ok(EnergyPerArea::from_base(zero)); // closed leaf: no exit signal.
            }
            let tau = beamlet_ray(centre, dir, perp, frame, leaf, leaf_width, geometry)
                .and_then(|beamlet| forward_project_ray(mu, &beamlet.ray, step_mm))
                .unwrap_or(zero);
            let transmission: Dimensionless<T> = OpticalDepth::new(Dimensionless::from_base(tau))?
                .transmission()
                .into_quantity();
            let delivered_fluence = EnergyPerArea::from_base(fluence_base);
            let exit_fluence: EnergyPerArea<T> = delivered_fluence * transmission;
            Ok(exit_fluence)
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use helios_math::ShippedScalar;
    use eunomia::assert_relative_eq;
    use helios_domain::VoxelGrid;
    use helios_math::Point3;

    // Uniform-μ cube: 9³ voxels, 2 mm spacing → central axial chord 16 mm = 1.6 cm.
    fn uniform_cube(mu_val: f64) -> Volume<f64> {
        let grid = VoxelGrid::axis_aligned([9, 9, 9], [2.0, 2.0, 2.0], Point3::new(0.0, 0.0, 0.0))
            .expect("grid");
        Volume::from_shape_fn(grid, move |_| mu_val)
    }

    fn single_leaf_frame(fluence: f64) -> DeliveryFrame<f64> {
        DeliveryFrame {
            projection: 0,
            gantry_angle_rad: Angle::from_unit::<Radian>(0.0),
            couch: Length::from_unit::<Millimeter>(8.0),
            leaf_fluence: vec![EnergyPerArea::from_base(fluence)],
        }
    }

    fn length(value: f64) -> Length<f64> {
        Length::from_unit::<Millimeter>(value)
    }

    #[test]
    fn no_attenuation_transmits_full_fluence() {
        let portal = frame_portal_fluence(
            &single_leaf_frame(2.0),
            &uniform_cube(0.0),
            BeamGeometry::Parallel {
                standoff: length(500.0),
            },
            length(2.0),
            length(0.25),
        )
        .expect("valid attenuation volume");
        assert_relative_eq!(*portal[0].as_base(), 2.0, epsilon = 1e-9);
    }

    #[test]
    fn uniform_medium_attenuates_by_beer_lambert() {
        // Central +x leaf through μ = 0.05 cm⁻¹, chord 1.6 cm → τ = 0.08.
        // Portal = fluence·exp(−0.08).
        let portal = frame_portal_fluence(
            &single_leaf_frame(3.0),
            &uniform_cube(0.05),
            BeamGeometry::Parallel {
                standoff: length(500.0),
            },
            length(2.0),
            length(0.25),
        )
        .expect("valid attenuation volume");
        assert_relative_eq!(
            *portal[0].as_base(),
            3.0 * (-0.05 * 1.6_f64).exp(),
            epsilon = 1e-9
        );
    }

    #[test]
    fn closed_leaf_reads_zero_and_more_attenuation_darkens() {
        let mu_lo = uniform_cube(0.05);
        let mu_hi = uniform_cube(0.20);
        let geom = BeamGeometry::Parallel {
            standoff: length(500.0),
        };
        // A closed leaf (0 fluence) among open ones reads exactly 0.
        let frame = DeliveryFrame {
            projection: 0,
            gantry_angle_rad: Angle::from_unit::<Radian>(0.0),
            couch: length(8.0),
            leaf_fluence: vec![1.0, 0.0, 1.0]
                .into_iter()
                .map(EnergyPerArea::from_base)
                .collect(),
        };
        let portal = frame_portal_fluence(&frame, &mu_lo, geom, length(2.0), length(0.25))
            .expect("valid attenuation volume");
        assert_relative_eq!(*portal[1].as_base(), 0.0, epsilon = 1e-15);
        // Higher μ darkens the transmitted signal (central leaf).
        let lo = *frame_portal_fluence(
            &single_leaf_frame(1.0),
            &mu_lo,
            geom,
            length(2.0),
            length(0.25),
        )
        .expect("valid attenuation volume")[0]
            .as_base();
        let hi = *frame_portal_fluence(
            &single_leaf_frame(1.0),
            &mu_hi,
            geom,
            length(2.0),
            length(0.25),
        )
        .expect("valid attenuation volume")[0]
            .as_base();
        assert!(
            hi < lo && hi > 0.0,
            "more attenuation must darken: {hi} !< {lo}"
        );
    }

    /// Portal fluence is Beer-Lambert transmission through the traversed path,
    /// scaled by the incident leaf fluence: an accumulation over the ray, one
    /// `exp`, and one product. `exp` is not required by IEEE 754 to be correctly
    /// rounded, so 64 ulps of `T` bounds the chain including the accumulation.
    const PORTAL_TRANSMISSION_ULPS: f64 = 64.0;

    /// Asserts portal transmission through a uniform slab in one scalar width.
    fn portal_fluence_follows_beer_lambert<T>()
    where
        T: ShippedScalar + helios_math::GeometryScalar,
    {
        // GeometryScalar and FloatElement both define `from_f64`; name the
        // one that performs the literal conversion so the call is unambiguous.
        let cast = <T as helios_math::FloatElement>::from_f64;
        let zero = cast(0.0);
        let spacing = cast(2.0);
        let grid = VoxelGrid::<T>::axis_aligned([9, 9, 9], [spacing; 3], Point3::new(zero, zero, zero))
            .expect("valid axis-aligned grid");

        let mu = cast(0.05);
        let attenuation = Volume::from_shape_fn(grid, |_| mu);
        let incident = cast(2.0);
        let frame = DeliveryFrame {
            projection: 0,
            gantry_angle_rad: Angle::from_unit::<Radian>(zero),
            couch: Length::from_unit::<Millimeter>(cast(8.0)),
            leaf_fluence: vec![EnergyPerArea::from_base(incident)],
        };

        let portal = frame_portal_fluence(
            &frame,
            &attenuation,
            BeamGeometry::Parallel {
                standoff: Length::from_unit::<Millimeter>(cast(500.0)),
            },
            Length::from_unit::<Millimeter>(spacing),
            Length::from_unit::<Millimeter>(cast(0.25)),
        )
        .expect("valid attenuation volume");

        // 1.6 cm of uniform mu between the source and the detector.
        let expected = incident * (-(mu * cast(1.6))).exp();
        assert_relative_eq!(
            *portal[0].as_base(),
            expected,
            max_relative = T::EPSILON * cast(PORTAL_TRANSMISSION_ULPS)
        );
    }

    #[test]
    fn portal_fluence_follows_beer_lambert_in_single_precision() {
        portal_fluence_follows_beer_lambert::<f32>();
    }

    #[test]
    fn portal_fluence_follows_beer_lambert_in_double_precision() {
        portal_fluence_follows_beer_lambert::<f64>();
    }

    #[test]
    fn negative_projected_optical_depth_is_rejected() {
        let error = frame_portal_fluence(
            &single_leaf_frame(1.0),
            &uniform_cube(-0.05),
            BeamGeometry::Parallel {
                standoff: length(500.0),
            },
            length(2.0),
            length(0.25),
        )
        .expect_err("negative optical depth must fail");
        assert!(matches!(
            error,
            TransportError::InvalidValue {
                field: hyperion::ValueKind::OpticalDepth,
                ..
            }
        ));
    }
}
