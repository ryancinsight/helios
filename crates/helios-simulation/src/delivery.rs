//! Integrated helical delivery: MLC modulation synchronized with the gantry/couch.
//!
//! Ties the [`HelicalDelivery`] kinematics (gantry rotation + couch translation)
//! to the binary-MLC [`LeafOpenTimeSinogram`]/[`MlcModel`], producing a
//! time-ordered sequence of [`DeliveryFrame`]s — each the delivery machine state
//! (gantry angle, couch position) together with the effective per-leaf fluence
//! actually delivered at that projection (leakage + tongue-and-groove applied).
//! This is the integrated imaging/delivery-workflow layer that a dose
//! accumulation (per-leaf beamlet ray-trace) builds on.

use helios_domain::{HelicalDelivery, LeafOpenTimeSinogram, MlcModel};
use helios_math::{NumericElement, Scalar};

/// The machine state and delivered fluence at one projection of a helical
/// delivery.
#[derive(Debug, Clone, PartialEq)]
pub struct DeliveryFrame<T: Scalar> {
    /// Projection index.
    pub projection: usize,
    /// Gantry angle (rad).
    pub gantry_angle_rad: T,
    /// Couch position (mm).
    pub couch_mm: T,
    /// Effective transmitted fluence per leaf (leakage + tongue-and-groove
    /// applied to the leaf-open-time pattern).
    pub leaf_fluence: Vec<T>,
}

/// Build the time-ordered delivery sequence: one [`DeliveryFrame`] per projection
/// of `leaf_open_times`, with the gantry/couch state from `delivery` and the
/// per-leaf effective fluence from `mlc`.
#[must_use]
pub fn simulate_helical_delivery<T: Scalar>(
    delivery: &HelicalDelivery<T>,
    leaf_open_times: &LeafOpenTimeSinogram<T>,
    mlc: &MlcModel<T>,
) -> Vec<DeliveryFrame<T>> {
    let (projections, leaves) = leaf_open_times.dims();
    let effective = mlc.effective_fluence_sinogram(leaf_open_times);
    (0..projections)
        .map(|p| DeliveryFrame {
            projection: p,
            gantry_angle_rad: delivery.gantry_angle_rad(p),
            couch_mm: delivery.couch_position_mm(p),
            leaf_fluence: effective[p * leaves..(p + 1) * leaves].to_vec(),
        })
        .collect()
}

/// Total delivered fluence integrated over all frames and leaves (a proxy for
/// total monitor units / beam-on).
#[must_use]
pub fn total_delivered_fluence<T: Scalar>(frames: &[DeliveryFrame<T>]) -> T {
    frames.iter().fold(<T as NumericElement>::ZERO, |acc, f| {
        acc + f
            .leaf_fluence
            .iter()
            .copied()
            .fold(<T as NumericElement>::ZERO, |a, x| a + x)
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    fn delivery() -> HelicalDelivery<f64> {
        HelicalDelivery::new(51, 25.0, 0.4, 10.0, 0.0, 0.0).expect("delivery")
    }

    #[test]
    fn frame_count_and_kinematics_track_the_sinogram() {
        let lot = LeafOpenTimeSinogram::from_fractions(3, 4, vec![0.5; 12]).unwrap();
        let mlc = MlcModel::new(0.01_f64, 0.0).unwrap();
        let del = delivery();
        let frames = simulate_helical_delivery(&del, &lot, &mlc);
        assert_eq!(frames.len(), 3);
        for (p, frame) in frames.iter().enumerate() {
            assert_eq!(frame.projection, p);
            assert_relative_eq!(
                frame.gantry_angle_rad,
                del.gantry_angle_rad(p),
                epsilon = 1e-12
            );
            assert_relative_eq!(frame.couch_mm, del.couch_position_mm(p), epsilon = 1e-12);
            assert_eq!(frame.leaf_fluence.len(), 4);
        }
    }

    #[test]
    fn frame_fluence_matches_mlc_model() {
        // Row [1,0,1,0] with T&G 0.1, no leakage — verify the frame reproduces the
        // MlcModel's neighbour-aware effective fluence.
        let lot = LeafOpenTimeSinogram::from_fractions(1, 4, vec![1.0, 0.0, 1.0, 0.0]).unwrap();
        let mlc = MlcModel::new(0.0_f64, 0.1).unwrap();
        let frames = simulate_helical_delivery(&delivery(), &lot, &mlc);
        let expected = mlc.effective_fluence_sinogram(&lot);
        assert_eq!(frames[0].leaf_fluence, expected);
        // leaf0: open=1, left=self(1), right=0 → 1 − 0.1·½·1 = 0.95.
        assert_relative_eq!(frames[0].leaf_fluence[0], 0.95, epsilon = 1e-15);
    }

    #[test]
    fn all_closed_delivers_only_leakage() {
        let lot = LeafOpenTimeSinogram::from_fractions(2, 5, vec![0.0; 10]).unwrap();
        let mlc = MlcModel::new(0.01_f64, 0.1).unwrap();
        let frames = simulate_helical_delivery(&delivery(), &lot, &mlc);
        for frame in &frames {
            for &f in &frame.leaf_fluence {
                assert_relative_eq!(f, 0.01, epsilon = 1e-15);
            }
        }
        // Total = 2 projections · 5 leaves · 0.01 leakage.
        assert_relative_eq!(
            total_delivered_fluence(&frames),
            2.0 * 5.0 * 0.01,
            epsilon = 1e-13
        );
    }

    #[test]
    fn all_open_delivers_full_fluence() {
        let lot = LeafOpenTimeSinogram::from_fractions(2, 5, vec![1.0; 10]).unwrap();
        let mlc = MlcModel::new(0.01_f64, 0.1).unwrap();
        let frames = simulate_helical_delivery(&delivery(), &lot, &mlc);
        // Uniformly open → no T&G loss, full transmission everywhere.
        assert_relative_eq!(total_delivered_fluence(&frames), 10.0, epsilon = 1e-13);
    }

    #[test]
    fn delivery_is_generic_over_scalar_f32() {
        let lot = LeafOpenTimeSinogram::from_fractions(1, 3, vec![1.0_f32, 0.0, 1.0]).unwrap();
        let mlc = MlcModel::new(0.0_f32, 0.1).unwrap();
        let del = HelicalDelivery::<f32>::new(51, 25.0, 0.4, 10.0, 0.0, 0.0).unwrap();
        let frames = simulate_helical_delivery(&del, &lot, &mlc);
        assert_relative_eq!(frames[0].leaf_fluence[0], 0.95_f32, epsilon = 1e-6);
    }
}
