//! Delivered-dose accumulation: ray-trace every [`DeliveryFrame`] beamlet into a
//! dose [`Volume`].
//!
//! Closes the delivery→dose loop. Each frame carries the machine state (gantry
//! angle, couch position) and the effective per-leaf fluence actually delivered
//! (leakage + tongue-and-groove already applied by the MLC model). This kernel
//! turns that time-ordered fluence into a spatial dose distribution by depositing
//! each leaf's beamlet terma through the attenuation volume, summed over all
//! frames — the input the DVH / gamma analysis consumes.
//!
//! # Beam geometry
//! A first-order helical TomoTherapy fan: at gantry angle `θ` the beam travels
//! along the axial-plane direction `d = (cosθ, sinθ, 0)`; each binary-MLC leaf is
//! a beamlet parallel to `d`, laterally offset along the in-plane perpendicular
//! `p = (−sinθ, cosθ, 0)` by `(leaf − centre)·leaf_width`, at the couch `z` slice.
//! Beamlets are parallel (small-fan approximation); a divergent point-source fan
//! and per-leaf collimation via gaia are a later increment.

use crate::delivery::DeliveryFrame;
use helios_domain::Volume;
use helios_math::{GeometryScalar, NumericElement, Point3, Ray, Vector3};
use helios_solver::deposit_ray_terma;

/// Accumulate the delivered dose from a helical-delivery `frames` sequence into a
/// dose [`Volume`] over the same grid as the attenuation volume `mu`.
///
/// `source_distance_mm` places each beamlet origin behind isocentre; `leaf_width_mm`
/// is the inter-leaf lateral pitch; `step_mm` is the ray-march sampling step. Dose
/// is linear in the per-leaf fluence, so scaling all fluence scales the dose and
/// the contributions of independent frames/leaves superpose (the test oracles).
#[must_use]
pub fn accumulate_delivered_dose<T: GeometryScalar>(
    frames: &[DeliveryFrame<T>],
    mu: &Volume<T>,
    source_distance_mm: T,
    leaf_width_mm: T,
    step_mm: T,
) -> Volume<T> {
    let zero = <T as NumericElement>::ZERO;
    let grid = *mu.grid();
    let [nx, ny, nz] = grid.dims();
    let centre = grid.voxel_center((nx - 1) / 2, (ny - 1) / 2, (nz - 1) / 2);
    let mut dose = Volume::zeros(grid);

    for frame in frames {
        let theta = frame.gantry_angle_rad;
        let (cos, sin) = (theta.cos(), theta.sin());
        // Central-axis direction and in-plane lateral (leaf-offset) direction.
        let dir = Vector3::new(cos, sin, zero);
        let perp = Vector3::new(-sin, cos, zero);

        let leaves = frame.leaf_fluence.len();
        let centre_leaf = <T as GeometryScalar>::from_f64((leaves as f64 - 1.0) * 0.5);
        for (leaf, &weight) in frame.leaf_fluence.iter().enumerate() {
            if weight <= zero {
                continue; // closed/leak-free leaf deposits nothing.
            }
            let offset =
                (<T as GeometryScalar>::from_f64(leaf as f64) - centre_leaf) * leaf_width_mm;
            // Beamlet parallel to `dir`, offset laterally along `perp`, at the
            // couch z (dir.z = perp.z = 0), origin standing off behind isocentre.
            let origin = Point3::new(
                centre.x + perp.x * offset - dir.x * source_distance_mm,
                centre.y + perp.y * offset - dir.y * source_distance_mm,
                frame.couch_mm,
            );
            if let Some(ray) = Ray::try_from_direction(origin, dir) {
                let _deposited = deposit_ray_terma(&mut dose, mu, &ray, weight, step_mm);
            }
        }
    }
    dose
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use helios_domain::VoxelGrid;

    // Uniform-μ cube: 9³ voxels, 2 mm spacing → node box [0,16] mm, centre 8 mm,
    // axial chord 16 mm = 1.6 cm.
    fn uniform_cube(mu_val: f64) -> Volume<f64> {
        let grid = VoxelGrid::axis_aligned([9, 9, 9], [2.0, 2.0, 2.0], Point3::new(0.0, 0.0, 0.0))
            .expect("grid");
        Volume::from_shape_fn(grid, move |_| mu_val)
    }

    // A single-leaf frame at gantry angle θ, couch z, fluence w. One leaf → the
    // beamlet is on the central axis (zero lateral offset).
    fn frame(gantry_angle_rad: f64, couch_mm: f64, w: f64) -> DeliveryFrame<f64> {
        DeliveryFrame {
            projection: 0,
            gantry_angle_rad,
            couch_mm,
            leaf_fluence: vec![w],
        }
    }

    // Analytic primary energy removed by one central axial beamlet of weight w.
    fn expected_axial_energy(mu_val: f64, w: f64) -> f64 {
        w * (1.0 - (-mu_val * 1.6_f64).exp())
    }

    #[test]
    fn single_central_beamlet_matches_analytic_energy() {
        // θ=0 → +x beamlet through the centre (couch z = 8 mm). Total delivered
        // dose = w·(1 − e^{−μ·1.6}).
        let mu = uniform_cube(0.05);
        let dose = accumulate_delivered_dose(&[frame(0.0, 8.0, 2.0)], &mu, 500.0, 2.0, 0.25);
        assert_relative_eq!(dose.sum(), expected_axial_energy(0.05, 2.0), epsilon = 1e-9);
    }

    #[test]
    fn zero_fluence_delivers_zero_dose() {
        let mu = uniform_cube(0.05);
        let dose = accumulate_delivered_dose(&[frame(0.0, 8.0, 0.0)], &mu, 500.0, 2.0, 0.5);
        assert_relative_eq!(dose.sum(), 0.0, epsilon = 1e-15);
    }

    #[test]
    fn dose_is_linear_in_fluence() {
        // Doubling every leaf fluence doubles the dose voxelwise.
        let mu = uniform_cube(0.05);
        let d1 = accumulate_delivered_dose(&[frame(0.0, 8.0, 1.0)], &mu, 500.0, 2.0, 0.25);
        let d2 = accumulate_delivered_dose(&[frame(0.0, 8.0, 2.0)], &mu, 500.0, 2.0, 0.25);
        assert_relative_eq!(d2.sum(), 2.0 * d1.sum(), epsilon = 1e-12);
        assert_relative_eq!(
            d2.get(0, 4, 4).unwrap(),
            2.0 * d1.get(0, 4, 4).unwrap(),
            epsilon = 1e-12
        );
    }

    #[test]
    fn independent_frames_superpose() {
        // Dose from [A, B] equals dose from [A] plus dose from [B]. Use two
        // different gantry angles so the beams occupy different voxels.
        let mu = uniform_cube(0.05);
        let a = frame(0.0, 8.0, 1.0);
        let b = frame(std::f64::consts::FRAC_PI_2, 8.0, 1.0);
        let together = accumulate_delivered_dose(&[a.clone(), b.clone()], &mu, 500.0, 2.0, 0.25);
        let da = accumulate_delivered_dose(&[a], &mu, 500.0, 2.0, 0.25);
        let db = accumulate_delivered_dose(&[b], &mu, 500.0, 2.0, 0.25);
        assert_relative_eq!(together.sum(), da.sum() + db.sum(), epsilon = 1e-12);
    }

    #[test]
    fn multi_leaf_fan_sums_offset_beamlets() {
        // Three open leaves at θ=0 place parallel +x beamlets at y = 6, 8, 10 mm
        // (perp = +y, 2 mm pitch). Each crosses the full 1.6 cm chord, so the
        // total energy is 3 × the single-beamlet value.
        let mu = uniform_cube(0.05);
        let f = DeliveryFrame {
            projection: 0,
            gantry_angle_rad: 0.0,
            couch_mm: 8.0,
            leaf_fluence: vec![1.0, 1.0, 1.0],
        };
        let dose = accumulate_delivered_dose(&[f], &mu, 500.0, 2.0, 0.25);
        assert_relative_eq!(
            dose.sum(),
            3.0 * expected_axial_energy(0.05, 1.0),
            epsilon = 1e-9
        );
        // The three beamlets land in distinct leaf rows (y = 3, 4, 5).
        for j in [3usize, 4, 5] {
            assert!(
                dose.get(0, j, 4).unwrap() > 0.0,
                "row {j} should be irradiated"
            );
        }
        assert_relative_eq!(dose.get(0, 1, 4).unwrap(), 0.0, epsilon = 1e-15);
    }

    #[test]
    fn accumulation_is_generic_over_scalar_f32() {
        let grid =
            VoxelGrid::<f32>::axis_aligned([9, 9, 9], [2.0, 2.0, 2.0], Point3::new(0.0, 0.0, 0.0))
                .unwrap();
        let mu = Volume::from_shape_fn(grid, |_| 0.05_f32);
        let f = DeliveryFrame {
            projection: 0,
            gantry_angle_rad: 0.0_f32,
            couch_mm: 8.0,
            leaf_fluence: vec![2.0_f32],
        };
        let dose = accumulate_delivered_dose(&[f], &mu, 500.0, 2.0, 0.25);
        let expected = 2.0_f32 * (1.0 - (-0.05_f32 * 1.6).exp());
        assert_relative_eq!(dose.sum(), expected, epsilon = 1e-5);
    }
}
