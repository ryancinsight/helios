    #![expect(
        clippy::unwrap_used,
        reason = "ratchet HELIOS-UNWRAP-1: pre-existing debt"
    )]
    use super::*;
    use eunomia::assert_relative_eq;
    use helios_domain::VoxelGrid;
    use helios_math::ShippedScalar;

    // Uniform-μ cube: 9³ voxels, 2 mm spacing → node box [0,16] mm, centre 8 mm,
    // axial chord 16 mm = 1.6 cm.
    fn uniform_cube(mu_val: f64) -> Volume<f64> {
        let grid = VoxelGrid::axis_aligned([9, 9, 9], [2.0, 2.0, 2.0], Point3::new(0.0, 0.0, 0.0))
            .expect("grid");
        Volume::from_shape_fn(grid, move |_| mu_val)
    }

    // A single-leaf frame at gantry angle θ, couch z, fluence w. One leaf → the
    // beamlet is on the central axis (zero lateral offset).
    fn length(value: f64) -> Length<f64> {
        Length::from_unit::<Millimeter>(value)
    }

    fn length_cm(value: f64) -> Length<f64> {
        Length::from_base(value * 0.01)
    }

    fn relative_weight(value: f64) -> aequitas::systems::si::quantities::Dimensionless<f64> {
        aequitas::systems::si::quantities::Dimensionless::from_base(value)
    }

    fn angle(value: f64) -> Angle<f64> {
        Angle::from_unit::<Radian>(value)
    }

    fn fluence(value: f64) -> aequitas::systems::si::quantities::EnergyPerArea<f64> {
        aequitas::systems::si::quantities::EnergyPerArea::from_base(value)
    }

    fn frame(gantry_angle_rad: f64, couch_mm: f64, w: f64) -> DeliveryFrame<f64> {
        DeliveryFrame {
            projection: 0,
            gantry_angle_rad: angle(gantry_angle_rad),
            couch: length(couch_mm),
            leaf_fluence: vec![fluence(w)],
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
        let dose = accumulate_delivered_dose(
            &[frame(0.0, 8.0, 2.0)],
            &mu,
            BeamGeometry::Parallel {
                standoff: length(500.0),
            },
            length(2.0),
            length(0.25),
        )
        .expect("valid attenuation volume");
        assert_relative_eq!(dose.sum(), expected_axial_energy(0.05, 2.0), epsilon = 1e-9);
    }

    #[test]
    fn zero_fluence_delivers_zero_dose() {
        let mu = uniform_cube(0.05);
        let dose = accumulate_delivered_dose(
            &[frame(0.0, 8.0, 0.0)],
            &mu,
            BeamGeometry::Parallel {
                standoff: length(500.0),
            },
            length(2.0),
            length(0.5),
        )
        .expect("valid attenuation volume");
        assert_relative_eq!(dose.sum(), 0.0, epsilon = 1e-15);
    }

    #[test]
    fn dose_is_linear_in_fluence() {
        // Doubling every leaf fluence doubles the dose voxelwise.
        let mu = uniform_cube(0.05);
        let d1 = accumulate_delivered_dose(
            &[frame(0.0, 8.0, 1.0)],
            &mu,
            BeamGeometry::Parallel {
                standoff: length(500.0),
            },
            length(2.0),
            length(0.25),
        )
        .expect("valid attenuation volume");
        let d2 = accumulate_delivered_dose(
            &[frame(0.0, 8.0, 2.0)],
            &mu,
            BeamGeometry::Parallel {
                standoff: length(500.0),
            },
            length(2.0),
            length(0.25),
        )
        .expect("valid attenuation volume");
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
        let together = accumulate_delivered_dose(
            &[a.clone(), b.clone()],
            &mu,
            BeamGeometry::Parallel {
                standoff: length(500.0),
            },
            length(2.0),
            length(0.25),
        )
        .expect("valid attenuation volume");
        let da = accumulate_delivered_dose(
            &[a],
            &mu,
            BeamGeometry::Parallel {
                standoff: length(500.0),
            },
            length(2.0),
            length(0.25),
        )
        .expect("valid attenuation volume");
        let db = accumulate_delivered_dose(
            &[b],
            &mu,
            BeamGeometry::Parallel {
                standoff: length(500.0),
            },
            length(2.0),
            length(0.25),
        )
        .expect("valid attenuation volume");
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
            gantry_angle_rad: angle(0.0),
            couch: length(8.0),
            leaf_fluence: vec![fluence(1.0), fluence(1.0), fluence(1.0)],
        };
        let dose = accumulate_delivered_dose(
            &[f],
            &mu,
            BeamGeometry::Parallel {
                standoff: length(500.0),
            },
            length(2.0),
            length(0.25),
        )
        .expect("valid attenuation volume");
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

    /// The absorbed fraction `1 - exp(-mu*L)` is a subtractive cancellation: at
    /// `mu*L = 0.08` the result is about 0.0769, so `exp`'s relative error is
    /// amplified roughly 13x. The comparison is also against a sum over the whole
    /// 9x9x9 dose grid, an independent second computation. Two hundred fifty-six
    /// ulps of `T` bounds both effects together.
    const ACCUMULATED_DOSE_ULPS: f64 = 256.0;

    /// Asserts accumulated delivered dose conserves absorbed energy in one width.
    fn accumulated_dose_matches_the_absorbed_fraction<T>()
    where
        T: ShippedScalar + helios_math::GeometryScalar,
    {
        // GeometryScalar and FloatElement both define `from_f64`; name the
        // one that performs the literal conversion so the call is unambiguous.
        let cast = <T as helios_math::FloatElement>::from_f64;
        let zero = cast(0.0);
        let spacing = cast(2.0);
        let grid =
            VoxelGrid::<T>::axis_aligned([9, 9, 9], [spacing; 3], Point3::new(zero, zero, zero))
                .expect("valid axis-aligned grid");

        let mu = cast(0.05);
        let attenuation = Volume::from_shape_fn(grid, |_| mu);
        let incident = cast(2.0);
        let frame = DeliveryFrame {
            projection: 0,
            gantry_angle_rad: Angle::from_unit::<Radian>(zero),
            couch: Length::from_unit::<Millimeter>(cast(8.0)),
            leaf_fluence: vec![aequitas::systems::si::quantities::EnergyPerArea::from_base(
                incident,
            )],
        };

        let dose = accumulate_delivered_dose(
            &[frame],
            &attenuation,
            BeamGeometry::Parallel {
                standoff: Length::from_unit::<Millimeter>(cast(500.0)),
            },
            Length::from_unit::<Millimeter>(spacing),
            Length::from_unit::<Millimeter>(cast(0.25)),
        )
        .expect("valid attenuation volume");

        // All energy removed from the beam over the 1.6 cm path is deposited.
        let expected = incident * (cast(1.0) - (-(mu * cast(1.6))).exp());
        assert_relative_eq!(
            dose.sum(),
            expected,
            max_relative = T::EPSILON * cast(ACCUMULATED_DOSE_ULPS)
        );
    }

    #[test]
    fn accumulated_dose_matches_the_absorbed_fraction_in_single_precision() {
        accumulated_dose_matches_the_absorbed_fraction::<f32>();
    }

    #[test]
    fn accumulated_dose_matches_the_absorbed_fraction_in_double_precision() {
        accumulated_dose_matches_the_absorbed_fraction::<f64>();
    }

    #[test]
    fn point_source_reduces_to_parallel_at_large_sad() {
        // As SAD → ∞ the divergent fan degenerates to parallel: the total dose of
        // a multi-leaf frame matches the parallel geometry.
        let mu = uniform_cube(0.05);
        let f = DeliveryFrame {
            projection: 0,
            gantry_angle_rad: angle(0.0),
            couch: length(8.0),
            leaf_fluence: vec![fluence(1.0), fluence(1.0), fluence(1.0)],
        };
        let parallel = accumulate_delivered_dose(
            std::slice::from_ref(&f),
            &mu,
            BeamGeometry::Parallel {
                standoff: length(500.0),
            },
            length(2.0),
            length(0.25),
        )
        .expect("valid attenuation volume");
        let far = accumulate_delivered_dose(
            &[f],
            &mu,
            BeamGeometry::PointSource {
                source_axis: length(1.0e6),
            },
            length(2.0),
            length(0.25),
        )
        .expect("valid attenuation volume");
        assert_relative_eq!(far.sum(), parallel.sum(), max_relative = 1e-4);
    }

    #[test]
    fn point_source_fan_diverges_across_rows() {
        // A far off-axis beamlet stays in a single detector row when parallel, but
        // sweeps several rows under a divergent point-source fan — the defining fan
        // property. Fine 1 mm grid so the divergence resolves past nearest-voxel
        // quantization.
        let grid =
            VoxelGrid::axis_aligned([31, 31, 1], [1.0, 1.0, 1.0], Point3::new(0.0, 0.0, 0.0))
                .unwrap();
        let mu = Volume::from_shape_fn(grid, |_| 0.05);
        // Single lit leaf at +6 mm offset (leaf 2 of 3 at 6 mm pitch).
        let f = DeliveryFrame {
            projection: 0,
            gantry_angle_rad: angle(0.0),
            couch: length(0.0),
            leaf_fluence: vec![fluence(0.0), fluence(0.0), fluence(1.0)],
        };
        let lit_rows = |dose: &Volume<f64>| -> usize {
            (0..31)
                .filter(|&j| (0..31).any(|i| dose.get(i, j, 0).unwrap() > 0.0))
                .count()
        };
        let par = accumulate_delivered_dose(
            std::slice::from_ref(&f),
            &mu,
            BeamGeometry::Parallel {
                standoff: length(500.0),
            },
            length(6.0),
            length(0.25),
        )
        .expect("valid attenuation volume");
        let pts = accumulate_delivered_dose(
            &[f],
            &mu,
            BeamGeometry::PointSource {
                source_axis: length(30.0),
            },
            length(6.0),
            length(0.25),
        )
        .expect("valid attenuation volume");
        assert_eq!(lit_rows(&par), 1, "parallel beamlet must stay in one row");
        assert!(
            lit_rows(&pts) >= 3,
            "divergent fan must sweep multiple rows, got {}",
            lit_rows(&pts)
        );
    }

    #[test]
    fn terma_then_scatter_produces_lateral_penumbra() {
        // End-to-end stage 1 → stage 2: a single central +x beamlet deposits terma
        // only on the y = z = centre line; the scatter kernel then spreads dose to
        // laterally-adjacent voxels that received *zero* primary terma (penumbra),
        // while the identity kernel leaves the terma unchanged.
        use helios_solver::{scatter_superposition, symmetric_deposition_kernel};
        let mu = uniform_cube(0.05);
        let terma = accumulate_delivered_dose(
            &[frame(0.0, 8.0, 1.0)],
            &mu,
            BeamGeometry::Parallel {
                standoff: length(500.0),
            },
            length(2.0),
            length(0.25),
        )
        .expect("valid attenuation volume");

        // Off-line voxel (mid-beam x=4, one voxel over in y) gets no primary terma.
        assert_relative_eq!(terma.get(4, 3, 4).unwrap(), 0.0, epsilon = 1e-15);

        // Identity kernel: dose == terma (differential vs the primary reference).
        let identity = scatter_superposition(&terma, &[1.0], &[1.0], &[1.0]);
        assert_relative_eq!(
            identity.get(4, 4, 4).unwrap(),
            terma.get(4, 4, 4).unwrap(),
            epsilon = 1e-15
        );

        // Spread kernel: the off-line voxel now receives scattered dose.
        let k = symmetric_deposition_kernel(length_cm(0.5), length_cm(0.2), 1);
        let dose = scatter_superposition(&terma, &k, &k, &k);
        assert!(
            dose.get(4, 3, 4).unwrap() > 0.0,
            "lateral neighbour must receive scattered penumbra dose"
        );
    }


mod tests_anisotropic;
