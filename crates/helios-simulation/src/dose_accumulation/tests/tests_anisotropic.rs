use super::super::*;
use super::*;
use eunomia::assert_relative_eq;

    // Energy-weighted centroid along the beam axis (x, at θ=0).
    fn beam_axis_centroid_x(dose: &Volume<f64>) -> f64 {
        let [nx, ny, nz] = dose.grid().dims();
        let (mut num, mut den) = (0.0, 0.0);
        for i in 0..nx {
            for j in 0..ny {
                for k in 0..nz {
                    let d = dose.get(i, j, k).unwrap();
                    num += d * dose.grid().voxel_center(i, j, k).x;
                    den += d;
                }
            }
        }
        num / den
    }

    #[test]
    fn anisotropic_isotropic_cone_matches_the_separable_scatter_pipeline() {
        // Single frame at θ=0 (beam = +x), an ISOTROPIC cone (equal up/down
        // ranges) and a 2 mm sample step (= voxel pitch): every sample lands on a
        // node, so the per-frame oriented scatter reduces to scatter_superposition
        // of the delivered terma — the differential oracle tying the new path to
        // the verified isotropic pipeline.
        use helios_solver::scatter_superposition;
        let mu = uniform_cube(0.05);
        let geom = BeamGeometry::Parallel {
            standoff: length(500.0),
        };
        let frames = [frame(0.0, 8.0, 2.0)];
        let voxel_spacing = length_cm(0.2);

        let cone = CollapsedCone::forward_peaked(
            length_cm(0.4),
            length_cm(0.4),
            length_cm(0.3),
            voxel_spacing,
            1,
            1,
            1,
        );
        let (beam_k, _) =
            forward_peaked_kernel(length_cm(0.4), length_cm(0.4), voxel_spacing, 1, 1);
        let lat = symmetric_deposition_kernel(length_cm(0.3), voxel_spacing, 1);

        let terma = accumulate_delivered_dose(&frames, &mu, geom, length(2.0), length(0.25))
            .expect("valid attenuation volume");
        let expected = scatter_superposition(&terma, &beam_k, &lat, &lat);
        let got = accumulate_delivered_dose_anisotropic(
            &frames,
            &mu,
            geom,
            length(2.0),
            length(0.25),
            &cone,
        )
        .expect("valid attenuation volume");

        for i in 0..9 {
            for j in 0..9 {
                for k in 0..9 {
                    assert_relative_eq!(
                        got.get(i, j, k).unwrap(),
                        expected.get(i, j, k).unwrap(),
                        epsilon = 1e-10
                    );
                }
            }
        }
    }

    #[test]
    fn forward_peaked_cone_shifts_delivered_dose_downstream() {
        // Same θ=0 delivery under an isotropic vs a forward-peaked cone: the
        // forward-peaked kernel must move the beam-axis dose centroid downstream
        // (larger x) — proof the anisotropy reaches delivered dose, not just the
        // solver primitive.
        let mu = uniform_cube(0.05);
        let geom = BeamGeometry::Parallel {
            standoff: length(500.0),
        };
        let frames = [frame(0.0, 8.0, 2.0)];
        let voxel_spacing = length_cm(0.2);
        let iso = CollapsedCone::forward_peaked(
            length_cm(0.4),
            length_cm(0.4),
            length_cm(0.3),
            voxel_spacing,
            2,
            2,
            1,
        );
        let fwd = CollapsedCone::forward_peaked(
            length_cm(0.1),
            length_cm(1.0),
            length_cm(0.3),
            voxel_spacing,
            2,
            2,
            1,
        );

        let d_iso = accumulate_delivered_dose_anisotropic(
            &frames,
            &mu,
            geom,
            length(2.0),
            length(0.25),
            &iso,
        )
        .expect("valid attenuation volume");
        let d_fwd = accumulate_delivered_dose_anisotropic(
            &frames,
            &mu,
            geom,
            length(2.0),
            length(0.25),
            &fwd,
        )
        .expect("valid attenuation volume");
        let (c_iso, c_fwd) = (beam_axis_centroid_x(&d_iso), beam_axis_centroid_x(&d_fwd));
        assert!(
            c_fwd > c_iso,
            "forward-peaked centroid {c_fwd} must exceed isotropic {c_iso}"
        );
    }

    #[test]
    fn single_component_poly_cone_matches_the_monoenergetic_cone() {
        // A one-component poly-energetic cone delivers exactly the same dose as
        // the monoenergetic cone with the same ranges (the differential oracle
        // tying the poly path to the verified monoenergetic one).
        let mu = uniform_cube(0.05);
        let geom = BeamGeometry::Parallel {
            standoff: length(500.0),
        };
        let frames = [frame(0.0, 8.0, 2.0)];
        let voxel_spacing = length_cm(0.2);
        let mono = CollapsedCone::forward_peaked(
            length_cm(0.1),
            length_cm(1.0),
            length_cm(0.3),
            voxel_spacing,
            2,
            3,
            1,
        );
        let poly = CollapsedCone::poly_forward_peaked(
            &[SpectralComponent {
                range_up: length_cm(0.1),
                range_down: length_cm(1.0),
                weight: relative_weight(5.0),
            }],
            length_cm(0.3),
            voxel_spacing,
            2,
            3,
            1,
        );
        let d_mono = accumulate_delivered_dose_anisotropic(
            &frames,
            &mu,
            geom,
            length(2.0),
            length(0.25),
            &mono,
        )
        .expect("valid attenuation volume");
        let d_poly = accumulate_delivered_dose_anisotropic(
            &frames,
            &mu,
            geom,
            length(2.0),
            length(0.25),
            &poly,
        )
        .expect("valid attenuation volume");
        for i in 0..9 {
            for j in 0..9 {
                for k in 0..9 {
                    assert_relative_eq!(
                        d_poly.get(i, j, k).unwrap(),
                        d_mono.get(i, j, k).unwrap(),
                        epsilon = 1e-13
                    );
                }
            }
        }
    }

    #[test]
    fn harder_spectrum_shifts_delivered_dose_further_downstream() {
        // Two-component spectra weighted soft vs hard: the harder-weighted beam
        // pushes the delivered-dose beam-axis centroid further downstream.
        let mu = uniform_cube(0.05);
        let geom = BeamGeometry::Parallel {
            standoff: length(500.0),
        };
        let frames = [frame(0.0, 8.0, 2.0)];
        let voxel_spacing = length_cm(0.2);
        let soft = SpectralComponent {
            range_up: length_cm(0.2),
            range_down: length_cm(0.3),
            weight: relative_weight(1.0),
        };
        let hard = SpectralComponent {
            range_up: length_cm(0.05),
            range_down: length_cm(1.5),
            weight: relative_weight(1.0),
        };
        let mostly_soft = CollapsedCone::poly_forward_peaked(
            &[
                SpectralComponent {
                    weight: relative_weight(9.0),
                    ..soft
                },
                hard,
            ],
            length_cm(0.3),
            voxel_spacing,
            2,
            2,
            1,
        );
        let mostly_hard = CollapsedCone::poly_forward_peaked(
            &[
                soft,
                SpectralComponent {
                    weight: relative_weight(9.0),
                    ..hard
                },
            ],
            length_cm(0.3),
            voxel_spacing,
            2,
            2,
            1,
        );
        let d_soft = accumulate_delivered_dose_anisotropic(
            &frames,
            &mu,
            geom,
            length(2.0),
            length(0.25),
            &mostly_soft,
        )
        .expect("valid attenuation volume");
        let d_hard = accumulate_delivered_dose_anisotropic(
            &frames,
            &mu,
            geom,
            length(2.0),
            length(0.25),
            &mostly_hard,
        )
        .expect("valid attenuation volume");
        assert!(
            beam_axis_centroid_x(&d_hard) > beam_axis_centroid_x(&d_soft),
            "harder spectrum centroid {} must exceed softer {}",
            beam_axis_centroid_x(&d_hard),
            beam_axis_centroid_x(&d_soft)
        );
    }

    #[test]
    fn anisotropic_dose_is_linear_in_fluence() {
        let mu = uniform_cube(0.05);
        let geom = BeamGeometry::Parallel {
            standoff: length(500.0),
        };
        let cone = CollapsedCone::forward_peaked(
            length_cm(0.1),
            length_cm(1.0),
            length_cm(0.3),
            length_cm(0.2),
            2,
            2,
            1,
        );
        let d1 = accumulate_delivered_dose_anisotropic(
            &[frame(0.0, 8.0, 1.0)],
            &mu,
            geom,
            length(2.0),
            length(0.25),
            &cone,
        )
        .expect("valid attenuation volume");
        let d2 = accumulate_delivered_dose_anisotropic(
            &[frame(0.0, 8.0, 2.0)],
            &mu,
            geom,
            length(2.0),
            length(0.25),
            &cone,
        )
        .expect("valid attenuation volume");
        assert_relative_eq!(d2.sum(), 2.0 * d1.sum(), epsilon = 1e-12);
    }
