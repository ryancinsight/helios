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
use crate::dose_accumulation::{beamlet_ray, gantry_basis, BeamGeometry};
use aequitas::systems::si::quantities::{Dimensionless, EnergyPerArea};
use helios_domain::Volume;
use helios_math::{GeometryScalar, NumericElement};
use helios_solver::forward_project_ray;
use hyperion::{quantity::OpticalDepth, TransportError};

/// Portal exit fluence per MLC leaf for one delivery `frame` through `mu`.
///
/// Returns a vector aligned with `frame.leaf_fluence`: entry `l` is the delivered
/// leaf fluence attenuated by the beamlet's optical depth,
/// `leaf_fluence[l] · exp(−τ_l)`. `geometry`/`leaf_width_mm`/`step_mm` are as in
/// [`accumulate_delivered_dose`](crate::accumulate_delivered_dose). A leaf whose
/// beamlet misses the volume reads its full (unattenuated) fluence.
///
/// # Errors
///
/// Returns [`TransportError`] if a beamlet produces a negative or non-finite
/// optical depth.
pub fn frame_portal_fluence<T: GeometryScalar>(
    frame: &DeliveryFrame<T>,
    mu: &Volume<T>,
    geometry: BeamGeometry<T>,
    leaf_width_mm: T,
    step_mm: T,
) -> Result<Vec<T>, TransportError<T>> {
    let zero = <T as NumericElement>::ZERO;
    let (centre, dir, perp) = gantry_basis(mu.grid(), frame.gantry_angle_rad);
    frame
        .leaf_fluence
        .iter()
        .enumerate()
        .map(|(leaf, &fluence)| {
            if fluence <= zero {
                return Ok(zero); // closed leaf: no exit signal.
            }
            let tau = beamlet_ray(centre, dir, perp, frame, leaf, leaf_width_mm, geometry)
                .and_then(|beamlet| forward_project_ray(mu, &beamlet.ray, step_mm))
                .unwrap_or(zero);
            let transmission: Dimensionless<T> = OpticalDepth::new(Dimensionless::from_base(tau))?
                .transmission()
                .into_quantity();
            let delivered_fluence = EnergyPerArea::from_base(fluence);
            let exit_fluence: EnergyPerArea<T> = delivered_fluence * transmission;
            Ok(exit_fluence.into_base())
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
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
            gantry_angle_rad: 0.0,
            couch_mm: 8.0,
            leaf_fluence: vec![fluence],
        }
    }

    #[test]
    fn no_attenuation_transmits_full_fluence() {
        let portal = frame_portal_fluence(
            &single_leaf_frame(2.0),
            &uniform_cube(0.0),
            BeamGeometry::Parallel { standoff_mm: 500.0 },
            2.0,
            0.25,
        )
        .expect("valid attenuation volume");
        assert_relative_eq!(portal[0], 2.0, epsilon = 1e-9);
    }

    #[test]
    fn uniform_medium_attenuates_by_beer_lambert() {
        // Central +x leaf through μ = 0.05 cm⁻¹, chord 1.6 cm → τ = 0.08.
        // Portal = fluence·exp(−0.08).
        let portal = frame_portal_fluence(
            &single_leaf_frame(3.0),
            &uniform_cube(0.05),
            BeamGeometry::Parallel { standoff_mm: 500.0 },
            2.0,
            0.25,
        )
        .expect("valid attenuation volume");
        assert_relative_eq!(portal[0], 3.0 * (-0.05 * 1.6_f64).exp(), epsilon = 1e-9);
    }

    #[test]
    fn closed_leaf_reads_zero_and_more_attenuation_darkens() {
        let mu_lo = uniform_cube(0.05);
        let mu_hi = uniform_cube(0.20);
        let geom = BeamGeometry::Parallel { standoff_mm: 500.0 };
        // A closed leaf (0 fluence) among open ones reads exactly 0.
        let frame = DeliveryFrame {
            projection: 0,
            gantry_angle_rad: 0.0,
            couch_mm: 8.0,
            leaf_fluence: vec![1.0, 0.0, 1.0],
        };
        let portal = frame_portal_fluence(&frame, &mu_lo, geom, 2.0, 0.25)
            .expect("valid attenuation volume");
        assert_relative_eq!(portal[1], 0.0, epsilon = 1e-15);
        // Higher μ darkens the transmitted signal (central leaf).
        let lo = frame_portal_fluence(&single_leaf_frame(1.0), &mu_lo, geom, 2.0, 0.25)
            .expect("valid attenuation volume")[0];
        let hi = frame_portal_fluence(&single_leaf_frame(1.0), &mu_hi, geom, 2.0, 0.25)
            .expect("valid attenuation volume")[0];
        assert!(
            hi < lo && hi > 0.0,
            "more attenuation must darken: {hi} !< {lo}"
        );
    }

    #[test]
    fn portal_is_generic_over_scalar_f32() {
        let grid =
            VoxelGrid::<f32>::axis_aligned([9, 9, 9], [2.0, 2.0, 2.0], Point3::new(0.0, 0.0, 0.0))
                .unwrap();
        let mu = Volume::from_shape_fn(grid, |_| 0.05_f32);
        let frame = DeliveryFrame {
            projection: 0,
            gantry_angle_rad: 0.0_f32,
            couch_mm: 8.0,
            leaf_fluence: vec![2.0_f32],
        };
        let portal = frame_portal_fluence(
            &frame,
            &mu,
            BeamGeometry::Parallel { standoff_mm: 500.0 },
            2.0,
            0.25,
        )
        .expect("valid attenuation volume");
        assert_relative_eq!(portal[0], 2.0_f32 * (-0.05_f32 * 1.6).exp(), epsilon = 1e-5);
    }

    #[test]
    fn negative_projected_optical_depth_is_rejected() {
        let error = frame_portal_fluence(
            &single_leaf_frame(1.0),
            &uniform_cube(-0.05),
            BeamGeometry::Parallel { standoff_mm: 500.0 },
            2.0,
            0.25,
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
