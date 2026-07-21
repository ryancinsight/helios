//! Integrated helical delivery: MLC modulation synchronized with the gantry/couch.
//!
//! Ties the [`HelicalDelivery`] kinematics (gantry rotation + couch translation)
//! to the binary-MLC [`LeafOpenTimeSinogram`]/[`MlcModel`], producing a
//! time-ordered sequence of [`DeliveryFrame`]s — each the delivery machine state
//! (gantry angle, couch position) together with the effective per-leaf fluence
//! actually delivered at that projection (leakage + tongue-and-groove applied).
//! This is the integrated imaging/delivery-workflow layer that a dose
//! accumulation (per-leaf beamlet ray-trace) builds on.

use helios_domain::{FieldAperture, HelicalDelivery, LeafOpenTimeSinogram, MlcModel};
use helios_math::{GeometryScalar, NumericElement, Point3, Scalar};

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

/// Apply a secondary-collimator [`FieldAperture`] to a delivery: scale each leaf's
/// fluence by the aperture transmission at that leaf's collimator coordinate
/// `(lateral_offset, couch_mm, 0)`, where
/// `lateral_offset = (leaf − centre_leaf)·leaf_width_mm`.
///
/// This is the jaw field-shaping + geometric edge penumbra, applied on top of the
/// per-leaf MLC modulation already in `frames`: leaves outside the open field are
/// blocked (fluence → 0), those inside pass unchanged (× 1), and those straddling
/// the field edge are partially transmitted. Returns new frames; the machine state
/// (gantry, couch) is preserved.
#[must_use]
pub fn collimate_frames<T: GeometryScalar>(
    frames: &[DeliveryFrame<T>],
    aperture: &FieldAperture<T>,
    leaf_width_mm: T,
) -> Vec<DeliveryFrame<T>> {
    let zero = <T as NumericElement>::ZERO;
    frames
        .iter()
        .map(|frame| {
            let leaves = frame.leaf_fluence.len();
            let centre = <T as GeometryScalar>::from_f64((leaves as f64 - 1.0) * 0.5);
            let leaf_fluence = frame
                .leaf_fluence
                .iter()
                .enumerate()
                .map(|(leaf, &fluence)| {
                    let offset =
                        (<T as GeometryScalar>::from_f64(leaf as f64) - centre) * leaf_width_mm;
                    let point = Point3::new(offset, frame.couch_mm, zero);
                    fluence * aperture.transmission(&point)
                })
                .collect();
            DeliveryFrame {
                projection: frame.projection,
                gantry_angle_rad: frame.gantry_angle_rad,
                couch_mm: frame.couch_mm,
                leaf_fluence,
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use eunomia::assert_relative_eq;

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

    fn open_frame() -> DeliveryFrame<f64> {
        DeliveryFrame {
            projection: 3,
            gantry_angle_rad: 1.2,
            couch_mm: 4.0,
            leaf_fluence: vec![1.0; 9],
        }
    }

    #[test]
    fn wide_aperture_leaves_fluence_unchanged() {
        // Aperture far larger than the leaf bank → every leaf fully open (× 1).
        let ap = FieldAperture::rectangular(Point3::new(0.0, 4.0, 0.0), [500.0, 500.0, 500.0], 2.0)
            .unwrap();
        let out = collimate_frames(&[open_frame()], &ap, 5.0);
        assert_eq!(out[0].leaf_fluence, vec![1.0; 9]);
    }

    #[test]
    fn narrow_aperture_shapes_the_field_with_edge_penumbra() {
        // 9 leaves, 5 mm pitch → lateral offsets −20..+20 mm. Field half-width
        // 10 mm, 2 mm penumbra: the ±10 mm edge lands on leaves 2 and 6 (→ 50 %),
        // interior leaves 3–5 fully open, and leaves 0–1 / 7–8 fully blocked.
        let ap = FieldAperture::rectangular(Point3::new(0.0, 4.0, 0.0), [10.0, 500.0, 500.0], 2.0)
            .unwrap();
        let out = collimate_frames(&[open_frame()], &ap, 5.0);
        let f = &out[0].leaf_fluence;
        assert_relative_eq!(f[4], 1.0, epsilon = 1e-12); // centre, open
        assert_relative_eq!(f[3], 1.0, epsilon = 1e-12);
        assert_relative_eq!(f[2], 0.5, epsilon = 1e-12); // on the field edge
        assert_relative_eq!(f[6], 0.5, epsilon = 1e-12);
        assert_relative_eq!(f[0], 0.0, epsilon = 1e-12); // outside the field
        assert_relative_eq!(f[8], 0.0, epsilon = 1e-12);
        assert!(
            f.iter().sum::<f64>() < 9.0,
            "collimation must remove fluence"
        );
        // Machine state is preserved.
        assert_eq!(out[0].projection, 3);
        assert_relative_eq!(out[0].gantry_angle_rad, 1.2, epsilon = 1e-15);
        assert_relative_eq!(out[0].couch_mm, 4.0, epsilon = 1e-15);
    }

    #[test]
    fn collimation_never_increases_fluence() {
        // Transmission ∈ [0,1], so each leaf's fluence can only be reduced.
        let ap = FieldAperture::rectangular(Point3::new(3.0, 4.0, 0.0), [6.0, 500.0, 500.0], 1.5)
            .unwrap();
        let frame = DeliveryFrame {
            leaf_fluence: vec![0.9, 0.3, 1.0, 0.7, 0.5, 0.8, 0.2, 1.0, 0.6],
            ..open_frame()
        };
        let out = collimate_frames(std::slice::from_ref(&frame), &ap, 5.0);
        for (before, after) in frame.leaf_fluence.iter().zip(&out[0].leaf_fluence) {
            assert!(*after <= *before + 1e-12 && *after >= 0.0);
        }
    }
}
