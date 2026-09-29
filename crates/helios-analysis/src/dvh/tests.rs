    #![expect(
        clippy::unwrap_used,
        reason = "ratchet HELIOS-UNWRAP-1: pre-existing debt"
    )]
    use super::*;
    use eunomia::assert_relative_eq;
    use helios_domain::VoxelGrid;
    use helios_math::Point3;
    use helios_math::ShippedScalar;

    fn dimensionless(value: f64) -> Dimensionless<f64> {
        Dimensionless::from_base(value)
    }

    fn grid(dims: [usize; 3]) -> VoxelGrid<f64> {
        VoxelGrid::axis_aligned(dims, [1.0, 1.0, 1.0], Point3::new(0.0, 0.0, 0.0)).expect("grid")
    }

    #[test]
    fn uniform_dose_is_a_step_histogram() {
        let dose = Volume::from_shape_fn(grid([4, 4, 4]), |_| 2.5);
        let dvh = Dvh::from_volume(&dose);
        assert_eq!(dvh.count(), 64);
        assert_relative_eq!(dvh.min().into_base(), 2.5, epsilon = 1e-15);
        assert_relative_eq!(dvh.max().into_base(), 2.5, epsilon = 1e-15);
        assert_relative_eq!(dvh.mean().into_base(), 2.5, epsilon = 1e-15);
        // At/below 2.5 → full volume; above → none.
        assert_relative_eq!(
            dvh.volume_fraction_at_dose(AbsorbedDose::from_base(2.5)),
            1.0,
            epsilon = 1e-15
        );
        assert_relative_eq!(
            dvh.volume_fraction_at_dose(AbsorbedDose::from_base(2.5001)),
            0.0,
            epsilon = 1e-15
        );
        // Every Dx equals the uniform dose.
        assert_relative_eq!(
            dvh.dose_at_volume_fraction(0.5).into_base(),
            2.5,
            epsilon = 1e-15
        );
        assert_relative_eq!(
            dvh.dose_at_volume_fraction(1.0).into_base(),
            2.5,
            epsilon = 1e-15
        );
    }

    #[test]
    fn linear_ramp_has_known_mean_and_quantiles() {
        // 100 voxels with dose = index 0..99.
        let g = VoxelGrid::axis_aligned([100, 1, 1], [1.0, 1.0, 1.0], Point3::new(0.0, 0.0, 0.0))
            .unwrap();
        let dose = Volume::from_shape_fn(g, |idx| idx[0] as f64);
        let dvh = Dvh::from_volume(&dose);
        assert_relative_eq!(dvh.min().into_base(), 0.0, epsilon = 1e-12);
        assert_relative_eq!(dvh.max().into_base(), 99.0, epsilon = 1e-12);
        assert_relative_eq!(dvh.mean().into_base(), 49.5, epsilon = 1e-12);
        // Half the volume (50 hottest voxels) receives ≥ dose 50 (values 50..99).
        assert_relative_eq!(
            dvh.dose_at_volume_fraction(0.5).into_base(),
            50.0,
            epsilon = 1e-12
        );
        // Exactly half of voxels have dose ≥ 50 (values 50..99 = 50 voxels).
        assert_relative_eq!(
            dvh.volume_fraction_at_dose(AbsorbedDose::from_base(50.0)),
            0.5,
            epsilon = 1e-12
        );
    }

    #[test]
    fn masked_dvh_restricts_to_the_structure() {
        // 4×4×4 dose: a "target" (i<2) at 2.0 Gy, surrounding "OAR" (i≥2) at 8.0.
        let dose = Volume::from_shape_fn(grid([4, 4, 4]), |idx| if idx[0] < 2 { 2.0 } else { 8.0 });

        let target = Dvh::from_volume_masked(&dose, |idx| idx[0] < 2);
        assert_eq!(target.count(), 32); // half of 64 voxels
        assert_relative_eq!(target.mean().into_base(), 2.0, epsilon = 1e-15);
        assert_relative_eq!(target.max().into_base(), 2.0, epsilon = 1e-15);

        let oar = Dvh::from_volume_masked(&dose, |idx| idx[0] >= 2);
        assert_eq!(oar.count(), 32);
        assert_relative_eq!(oar.mean().into_base(), 8.0, epsilon = 1e-15);

        // Whole-volume mean (5.0) differs from either structure — masking matters.
        assert_relative_eq!(
            Dvh::from_volume(&dose).mean().into_base(),
            5.0,
            epsilon = 1e-15
        );
    }

    #[test]
    fn single_voxel_mask_is_a_point_dvh() {
        let dose = Volume::from_shape_fn(grid([3, 3, 3]), |idx| (idx[0] + idx[1] + idx[2]) as f64);
        let point = Dvh::from_volume_masked(&dose, |idx| idx == [2, 2, 2]);
        assert_eq!(point.count(), 1);
        assert_relative_eq!(point.min().into_base(), 6.0, epsilon = 1e-15);
        assert_relative_eq!(
            point.dose_at_volume_fraction(1.0).into_base(),
            6.0,
            epsilon = 1e-15
        );
    }

    #[test]
    fn homogeneity_index_is_zero_for_uniform_and_known_for_a_ramp() {
        // Uniform dose → every Dx equal → HI = 0.
        let uniform = Dvh::from_volume(&Volume::from_shape_fn(grid([4, 4, 4]), |_| 3.0));
        assert_relative_eq!(uniform.homogeneity_index(), 0.0, epsilon = 1e-15);

        // Ramp 0..99 over 100 voxels: D2 = 98, D98 = 2, D50 = 50 → HI = 96/50 = 1.92.
        let g = VoxelGrid::axis_aligned([100, 1, 1], [1.0, 1.0, 1.0], Point3::new(0.0, 0.0, 0.0))
            .unwrap();
        let ramp = Dvh::from_volume(&Volume::from_shape_fn(g, |idx| idx[0] as f64));
        assert_relative_eq!(ramp.homogeneity_index(), 96.0 / 50.0, epsilon = 1e-12);
    }

    /// The mean of a uniform volume is a sum of eight identical values divided by
    /// eight, so at most nine roundings separate it from the constant. Sixteen
    /// ulps of `T` bounds that.
    const UNIFORM_MEAN_ULPS: f64 = 16.0;

    /// Asserts the DVH mean of a uniform dose in one scalar width.
    fn dvh_mean_of_a_uniform_dose_is_that_dose<T: ShippedScalar>() {
        let zero = T::from_f64(0.0);
        let unit = T::from_f64(1.0);
        let grid =
            VoxelGrid::<T>::axis_aligned([2, 2, 2], [unit; 3], Point3::new(zero, zero, zero))
                .expect("valid axis-aligned grid");

        let dose = Volume::from_shape_fn(grid, |_| unit);
        let dvh = Dvh::from_volume(&dose);

        assert_relative_eq!(
            dvh.mean().into_base(),
            unit,
            max_relative = T::EPSILON * T::from_f64(UNIFORM_MEAN_ULPS)
        );
    }

    #[test]
    fn dvh_mean_of_a_uniform_dose_is_that_dose_in_single_precision() {
        dvh_mean_of_a_uniform_dose_is_that_dose::<f32>();
    }

    #[test]
    fn dvh_mean_of_a_uniform_dose_is_that_dose_in_double_precision() {
        dvh_mean_of_a_uniform_dose_is_that_dose::<f64>();
    }

    #[test]
    fn dvh_geud_reuses_the_sample_and_matches_the_canonical_law() {
        // Heterogeneous dose; the DVH method must equal direct Asclepius
        // evaluation over the identical borrowed sample.
        let dose = Volume::from_shape_fn(grid([4, 4, 4]), |idx| {
            1.0 + idx[0] as f64 + 0.5 * idx[1] as f64
        });
        let dvh = Dvh::from_volume(&dose);
        assert_eq!(dvh.dose_sample().len(), 64);
        for a in [1.0, 2.0, -3.0, 8.0] {
            let model = GeneralizedEquivalentUniformDose::new(
                VolumeEffect::new(a).expect("finite non-zero exponent"),
            );
            let direct = model
                .evaluate(dvh.dose_sample())
                .expect("valid dose observation");
            assert_relative_eq!(
                dvh.generalized_eud(dimensionless(a))
                    .expect("valid response")
                    .into_base(),
                direct.into_base(),
                max_relative = 1e-13
            );
        }
    }

    #[test]
    fn threshold_query_preserves_boundary_and_nan_semantics() {
        let g = VoxelGrid::axis_aligned([4, 1, 1], [1.0, 1.0, 1.0], Point3::new(0.0, 0.0, 0.0))
            .expect("grid");
        let finite = Volume::from_shape_fn(g, |[i, _, _]| i as f64);
        let finite_dvh = Dvh::from_volume(&finite);
        assert_eq!(
            finite_dvh.volume_fraction_at_dose(AbsorbedDose::from_base(-1.0)),
            1.0
        );
        assert_eq!(
            finite_dvh.volume_fraction_at_dose(AbsorbedDose::from_base(0.0)),
            1.0
        );
        assert_eq!(
            finite_dvh.volume_fraction_at_dose(AbsorbedDose::from_base(3.0)),
            0.25
        );
        assert_eq!(
            finite_dvh.volume_fraction_at_dose(AbsorbedDose::from_base(4.0)),
            0.0
        );
        assert_eq!(
            finite_dvh.volume_fraction_at_dose(AbsorbedDose::from_base(f64::NAN)),
            0.0
        );

        let with_nan = Volume::from_shape_fn(
            VoxelGrid::axis_aligned([3, 1, 1], [1.0, 1.0, 1.0], Point3::new(0.0, 0.0, 0.0))
                .expect("grid"),
            |[i, _, _]| [1.0, f64::NAN, 3.0][i],
        );
        let nan_dvh = Dvh::from_volume(&with_nan);
        assert_eq!(
            nan_dvh.volume_fraction_at_dose(AbsorbedDose::from_base(1.0)),
            2.0 / 3.0
        );
        assert_eq!(
            nan_dvh.volume_fraction_at_dose(AbsorbedDose::from_base(2.0)),
            1.0 / 3.0
        );
        assert!(matches!(
            nan_dvh.generalized_eud(dimensionless(1.0)),
            Err(ResponseError::InvalidObservation { index: 2, .. })
        ));
    }

    #[test]
    fn outcome_parameters_preserve_asclepius_validation_errors() {
        let dvh = Dvh::from_volume(&Volume::from_shape_fn(grid([2, 1, 1]), |_| 1.0));
        assert!(matches!(
            dvh.generalized_eud(dimensionless(0.0)),
            Err(ResponseError::InvalidValue(source))
                if source.kind() == asclepius::ValueKind::VolumeEffect
                    && source.constraint() == asclepius::ValueConstraint::FiniteNonZero
        ));
        assert!(matches!(
            dvh.tcp_logistic(dimensionless(1.0), AbsorbedDose::from_base(1.0), -0.5),
            Err(ResponseError::InvalidValue(source))
                if source.kind() == asclepius::ValueKind::Gamma50
                    && source.constraint() == asclepius::ValueConstraint::FinitePositive
        ));
        assert!(matches!(
            dvh.ntcp_lkb(dimensionless(1.0), AbsorbedDose::from_base(1.0), 0.0),
            Err(ResponseError::InvalidValue(source))
                if source.kind() == asclepius::ValueKind::LymanSlope
                    && source.constraint() == asclepius::ValueConstraint::FinitePositive
        ));
    }

    #[test]
    fn uniform_structure_outcome_equals_the_pointwise_model() {
        // A uniform-dose structure has gEUD = that dose for any a, so its TCP/NTCP
        // reduce to the pointwise outcome models at that dose.
        let d = 62.0;
        let dvh = Dvh::from_volume(&Volume::from_shape_fn(grid([3, 3, 3]), |_| d));
        assert_relative_eq!(
            dvh.generalized_eud(dimensionless(-10.0))
                .expect("valid response")
                .into_base(),
            d,
            max_relative = 1e-12
        );
        // NTCP with TD50 = d ⇒ t = 0 ⇒ 0.5; TCP with TCD50 = d ⇒ 0.5.
        assert_relative_eq!(
            dvh.ntcp_lkb(dimensionless(1.0), AbsorbedDose::from_base(d), 0.2)
                .expect("valid response"),
            0.5,
            epsilon = 1e-12
        );
        assert_relative_eq!(
            dvh.tcp_logistic(dimensionless(1.0), AbsorbedDose::from_base(d), 2.0)
                .expect("valid response"),
            0.5,
            epsilon = 1e-12
        );
        let geud = dvh
            .generalized_eud(dimensionless(2.0))
            .expect("valid response");
        let ntcp_model = LymanComplicationProbability::new(
            AbsorbedDose::from_base(50.0),
            LymanSlope::new(0.2).expect("positive Lyman slope"),
        )
        .expect("positive midpoint");
        let tcp_model = LogisticControlProbability::new(
            AbsorbedDose::from_base(55.0),
            Gamma50::new(2.0).expect("positive gamma50"),
        )
        .expect("positive midpoint");
        assert_relative_eq!(
            dvh.ntcp_lkb(dimensionless(2.0), AbsorbedDose::from_base(50.0), 0.2)
                .expect("valid response"),
            ntcp_model.evaluate(geud).expect("valid response").get(),
            epsilon = 1e-14
        );
        assert_relative_eq!(
            dvh.tcp_logistic(dimensionless(2.0), AbsorbedDose::from_base(55.0), 2.0)
                .expect("valid response"),
            tcp_model.evaluate(geud).expect("valid response").get(),
            epsilon = 1e-14
        );
    }

    #[test]
    fn masked_structure_outcome_reflects_only_masked_voxels() {
        // Two half-slabs at 20 and 80 Gy; masking the hot half gives gEUD ≈ 80,
        // higher NTCP than masking the cold half — the per-structure evaluation.
        let dose =
            Volume::from_shape_fn(grid([4, 4, 4]), |idx| if idx[0] < 2 { 20.0 } else { 80.0 });
        let hot = Dvh::from_volume_masked(&dose, |idx| idx[0] >= 2);
        let cold = Dvh::from_volume_masked(&dose, |idx| idx[0] < 2);
        assert_relative_eq!(
            hot.generalized_eud(dimensionless(1.0))
                .expect("valid response")
                .into_base(),
            80.0,
            epsilon = 1e-12
        );
        assert_relative_eq!(
            cold.generalized_eud(dimensionless(1.0))
                .expect("valid response")
                .into_base(),
            20.0,
            epsilon = 1e-12
        );
        assert!(
            hot.ntcp_lkb(dimensionless(1.0), AbsorbedDose::from_base(50.0), 0.2)
                .expect("valid response")
                > cold
                    .ntcp_lkb(dimensionless(1.0), AbsorbedDose::from_base(50.0), 0.2)
                    .expect("valid response")
        );
    }
