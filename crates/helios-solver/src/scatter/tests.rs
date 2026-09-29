    #![expect(
        clippy::unwrap_used,
        reason = "ratchet HELIOS-UNWRAP-1: pre-existing debt"
    )]
    use super::*;
    use eunomia::assert_relative_eq;
    use helios_math::Point3;
    use helios_math::ShippedScalar;

    fn reference_convolve_axis(
        vol: &Volume<f64>,
        kernel: &[f64],
        center: usize,
        axis: usize,
    ) -> Volume<f64> {
        Volume::from_shape_fn(*vol.grid(), |mut index| {
            let pos = index[axis] as isize;
            kernel
                .iter()
                .enumerate()
                .filter_map(|(tap, &weight)| {
                    let source = pos - (tap as isize - center as isize);
                    (source >= 0 && source < vol.grid().dims()[axis] as isize).then(|| {
                        index[axis] = source as usize;
                        vol.get(index[0], index[1], index[2])
                            .expect("reference index is inside the validated grid")
                            * weight
                    })
                })
                .sum()
        })
    }

    fn grid() -> VoxelGrid<f64> {
        VoxelGrid::axis_aligned([7, 7, 7], [2.0, 2.0, 2.0], Point3::new(0.0, 0.0, 0.0))
            .expect("grid")
    }

    // Terma concentrated in the single centre voxel (3,3,3).
    fn point_terma() -> Volume<f64> {
        Volume::from_shape_fn(grid(), |idx| if idx == [3, 3, 3] { 1.0 } else { 0.0 })
    }

    fn length_cm<T: ShippedScalar>(value: T) -> Length<T> {
        Length::from_base(value * T::from_f64(0.01))
    }

    fn relative_weight<T: ShippedScalar>(value: T) -> Dimensionless<T> {
        Dimensionless::from_base(value)
    }

    #[test]
    fn forward_peaked_reduces_to_symmetric_for_equal_ranges() {
        // Equal up/down ranges and radii → identical taps and centre to the
        // symmetric kernel (differential oracle for the asymmetric constructor).
        let sym = symmetric_deposition_kernel(length_cm(0.5), length_cm(0.2), 2);
        let (fp, centre) =
            forward_peaked_kernel(length_cm(0.5), length_cm(0.5), length_cm(0.2), 2, 2);
        assert_eq!(centre, 2);
        assert_eq!(fp.len(), sym.len());
        for (a, b) in fp.iter().zip(&sym) {
            assert_relative_eq!(a, b, epsilon = 1e-15);
        }
        // And the anisotropic superposition equals the symmetric one exactly.
        let iso = scatter_superposition(&point_terma(), &sym, &sym, &sym);
        let aniso = anisotropic_scatter_superposition(&point_terma(), 0, &fp, centre, &sym);
        for i in 0..7 {
            for j in 0..7 {
                for k in 0..7 {
                    assert_relative_eq!(
                        aniso.get(i, j, k).unwrap(),
                        iso.get(i, j, k).unwrap(),
                        epsilon = 1e-15
                    );
                }
            }
        }
    }

    #[test]
    fn forward_peaking_puts_more_energy_downstream_than_upstream() {
        // Long downstream range, short upstream: the voxel one step downstream
        // (+x) of the point source must receive strictly more dose than one step
        // upstream (−x) — the defining collapsed-cone anisotropy. Laterally the
        // spread stays symmetric.
        let (fp, centre) =
            forward_peaked_kernel(length_cm(0.1), length_cm(1.0), length_cm(0.2), 1, 3);
        let lat = symmetric_deposition_kernel(length_cm(0.3), length_cm(0.2), 1);
        let dose = anisotropic_scatter_superposition(&point_terma(), 0, &fp, centre, &lat);
        let down = dose.get(4, 3, 3).unwrap(); // +x of the source at (3,3,3)
        let up = dose.get(2, 3, 3).unwrap(); // −x
        assert!(down > up, "downstream {down} must exceed upstream {up}");
        // Lateral symmetry is preserved.
        assert_relative_eq!(
            dose.get(3, 2, 3).unwrap(),
            dose.get(3, 4, 3).unwrap(),
            epsilon = 1e-14
        );
    }

    #[test]
    fn anisotropic_kernel_conserves_interior_point_energy() {
        // Source far enough from every boundary that no tail truncates: the total
        // dose equals the total terma (Σ=1 normalization on every axis).
        let (fp, centre) =
            forward_peaked_kernel(length_cm(0.2), length_cm(0.6), length_cm(0.2), 1, 2);
        let lat = symmetric_deposition_kernel(length_cm(0.4), length_cm(0.2), 1);
        let dose = anisotropic_scatter_superposition(&point_terma(), 0, &fp, centre, &lat);
        assert_relative_eq!(dose.sum(), 1.0, epsilon = 1e-13);
    }

    fn spectral(range_up: f64, range_down: f64, weight: f64) -> SpectralComponent<f64> {
        SpectralComponent {
            range_up: length_cm(range_up),
            range_down: length_cm(range_down),
            weight: relative_weight(weight),
        }
    }

    #[test]
    fn single_component_poly_reduces_to_the_monoenergetic_kernel() {
        // A one-component spectrum (any positive weight) is the monoenergetic
        // kernel exactly — the weight cancels in the Σ=1 renormalization.
        let (mono, mc) =
            forward_peaked_kernel(length_cm(0.15), length_cm(0.6), length_cm(0.2), 1, 3);
        let (poly, pc) =
            poly_forward_peaked_kernel(&[spectral(0.15, 0.6, 3.7)], length_cm(0.2), 1, 3);
        assert_eq!(mc, pc);
        assert_eq!(mono.len(), poly.len());
        for (a, b) in poly.iter().zip(&mono) {
            assert_relative_eq!(a, b, epsilon = 1e-15);
        }
    }

    #[test]
    fn poly_kernel_is_invariant_to_weight_scaling_and_sums_to_one() {
        // Scaling every weight by a constant leaves the convex combination
        // unchanged, and the kernel is normalized.
        let comps = [spectral(0.1, 0.4, 1.0), spectral(0.2, 1.2, 2.0)];
        let scaled = [spectral(0.1, 0.4, 10.0), spectral(0.2, 1.2, 20.0)];
        let (ka, _) = poly_forward_peaked_kernel(&comps, length_cm(0.2), 1, 3);
        let (kb, _) = poly_forward_peaked_kernel(&scaled, length_cm(0.2), 1, 3);
        for (a, b) in ka.iter().zip(&kb) {
            assert_relative_eq!(a, b, epsilon = 1e-15);
        }
        assert_relative_eq!(ka.iter().sum::<f64>(), 1.0, epsilon = 1e-14);
    }

    #[test]
    fn harder_spectrum_is_more_forward_peaked() {
        // Downstream/upstream mass ratio grows as weight shifts to the harder
        // (longer downstream range) component — the beam-hardening signature.
        let soft = spectral(0.2, 0.3, 1.0); // near-isotropic
        let hard = spectral(0.05, 1.5, 1.0); // strongly forward
        let ratio = |k: &[f64], centre: usize| -> f64 {
            let down: f64 = k[centre + 1..].iter().sum();
            let up: f64 = k[..centre].iter().sum();
            down / up
        };
        let (mostly_soft, c) = poly_forward_peaked_kernel(
            &[spectral(0.2, 0.3, 9.0), spectral(0.05, 1.5, 1.0)],
            length_cm(0.2),
            2,
            2,
        );
        let (mostly_hard, _) = poly_forward_peaked_kernel(
            &[spectral(0.2, 0.3, 1.0), spectral(0.05, 1.5, 9.0)],
            length_cm(0.2),
            2,
            2,
        );
        let _ = (soft, hard);
        assert!(
            ratio(&mostly_hard, c) > ratio(&mostly_soft, c),
            "harder spectrum ratio {} must exceed softer {}",
            ratio(&mostly_hard, c),
            ratio(&mostly_soft, c)
        );
    }

    #[test]
    fn empty_spectrum_is_the_identity_kernel() {
        let (k, centre) = poly_forward_peaked_kernel::<f64>(&[], length_cm(0.2), 2, 1);
        assert_eq!(centre, 2);
        assert_eq!(k, vec![0.0, 0.0, 1.0, 0.0]); // centred delta at radius_up = 2
    }

    fn single_component_poly_matches_the_monoenergetic_kernel<T: ShippedScalar>() {
        let (mono, _) = forward_peaked_kernel(
            length_cm(T::from_f64(0.1)),
            length_cm(T::from_f64(0.5)),
            length_cm(T::from_f64(0.2)),
            1,
            2,
        );
        let (poly, _) = poly_forward_peaked_kernel(
            &[SpectralComponent {
                range_up: length_cm(T::from_f64(0.1)),
                range_down: length_cm(T::from_f64(0.5)),
                weight: relative_weight(T::from_f64(2.0)),
            }],
            length_cm(T::from_f64(0.2)),
            1,
            2,
        );

        // One spectral component reduces to the monoenergetic kernel: the weight
        // divides out in normalization, so only that division separates the two.
        for (from_poly, from_mono) in poly.iter().zip(&mono) {
            assert_relative_eq!(
                from_poly,
                from_mono,
                max_relative = T::EPSILON * T::from_f64(WEIGHT_NORMALIZATION_ULPS)
            );
        }
    }

    #[test]
    fn single_component_poly_matches_the_monoenergetic_kernel_in_single_precision() {
        single_component_poly_matches_the_monoenergetic_kernel::<f32>();
    }

    #[test]
    fn single_component_poly_matches_the_monoenergetic_kernel_in_double_precision() {
        single_component_poly_matches_the_monoenergetic_kernel::<f64>();
    }

    /// A single spectral component's weight cancels in normalization, so the poly
    /// kernel differs from the monoenergetic one only by that division. Eight ulps
    /// bounds a handful of roundings.
    const WEIGHT_NORMALIZATION_ULPS: f64 = 8.0;

    /// Mirror symmetry across the beam axis compares two sums built from the same
    /// kernel taps in mirrored order. Reordering a sum of `n` terms perturbs it by
    /// at most `n * T::EPSILON`; these kernels have single-digit tap counts, so 16
    /// ulps bounds the reordering with margin.
    const MIRROR_SYMMETRY_ULPS: f64 = 16.0;

    /// Asserts beam-axis selection and transverse symmetry in one scalar width.
    fn anisotropic_scatter_is_symmetric_across_the_beam_axis<T: ShippedScalar>() {
        let zero = T::from_f64(0.0);
        let one = T::from_f64(1.0);
        let spacing = T::from_f64(2.0);
        let grid =
            VoxelGrid::<T>::axis_aligned([7, 7, 7], [spacing; 3], Point3::new(zero, zero, zero))
                .expect("valid axis-aligned grid");

        let terma = Volume::from_shape_fn(grid, |idx| if idx == [3, 3, 3] { one } else { zero });
        let (forward, centre) = forward_peaked_kernel(
            length_cm(T::from_f64(0.1)),
            length_cm(one),
            length_cm(T::from_f64(0.2)),
            1,
            2,
        );
        let lateral = symmetric_deposition_kernel(
            length_cm(T::from_f64(0.3)),
            length_cm(T::from_f64(0.2)),
            1,
        );

        // Beam along y: anisotropy shows on the j axis, not i.
        let dose = anisotropic_scatter_superposition(&terma, 1, &forward, centre, &lateral);
        assert!(
            dose.get(3, 4, 3).unwrap() > dose.get(3, 2, 3).unwrap(),
            "forward peaking must favour +y over -y"
        );
        assert_relative_eq!(
            dose.get(2, 3, 3).unwrap(),
            dose.get(4, 3, 3).unwrap(),
            max_relative = T::EPSILON * T::from_f64(MIRROR_SYMMETRY_ULPS)
        );
    }

    #[test]
    fn anisotropic_scatter_is_symmetric_across_the_beam_axis_in_single_precision() {
        anisotropic_scatter_is_symmetric_across_the_beam_axis::<f32>();
    }

    #[test]
    fn anisotropic_scatter_is_symmetric_across_the_beam_axis_in_double_precision() {
        anisotropic_scatter_is_symmetric_across_the_beam_axis::<f64>();
    }

    #[test]
    fn const_axis_matches_bounds_checked_reference_bitwise() {
        let terma = Volume::from_shape_fn(grid(), |[i, j, k]| {
            (i * 49 + j * 7 + k) as f64 * 0.125 - 3.0
        });
        let kernel = [0.125, 0.25, 0.375, 0.25];
        let center = 1;

        for axis in 0..3 {
            let expected = reference_convolve_axis(&terma, &kernel, center, axis);
            let actual = match axis {
                0 => convolve_axis_at_const::<f64, 0>(&terma, &kernel, center),
                1 => convolve_axis_at_const::<f64, 1>(&terma, &kernel, center),
                2 => convolve_axis_at_const::<f64, 2>(&terma, &kernel, center),
                _ => unreachable!("invariant: loop range contains only the three axes"),
            };
            assert_eq!(actual.as_slice(), expected.as_slice());
        }
    }

    #[test]
    fn delta_kernel_is_identity() {
        // [1] on every axis deposits energy locally → dose == terma exactly.
        let terma = Volume::from_shape_fn(grid(), |idx| (idx[0] + 2 * idx[1] + 3 * idx[2]) as f64);
        let dose = scatter_superposition(&terma, &[1.0], &[1.0], &[1.0]);
        for i in 0..7 {
            for j in 0..7 {
                for k in 0..7 {
                    assert_relative_eq!(
                        dose.get(i, j, k).unwrap(),
                        terma.get(i, j, k).unwrap(),
                        epsilon = 1e-15
                    );
                }
            }
        }
    }

    #[test]
    fn symmetric_kernel_spreads_symmetrically() {
        // Spread the point terma along x only; the two x-neighbours of the centre
        // receive equal dose (kernel symmetry), and both are < the centre.
        let kx = symmetric_deposition_kernel(length_cm(0.4), length_cm(0.2), 2);
        let dose = scatter_superposition(&point_terma(), &kx, &[1.0], &[1.0]);
        let left = dose.get(2, 3, 3).unwrap();
        let right = dose.get(4, 3, 3).unwrap();
        assert_relative_eq!(left, right, epsilon = 1e-14);
        assert!(left > 0.0 && left < dose.get(3, 3, 3).unwrap());
    }

    #[test]
    fn normalized_kernel_conserves_point_energy_in_interior() {
        // The centre is >= radius from every boundary, so no tail is truncated:
        // the total dose must equal the total terma (energy conservation).
        let k = symmetric_deposition_kernel(length_cm(0.5), length_cm(0.2), 2);
        let dose = scatter_superposition(&point_terma(), &k, &k, &k);
        assert_relative_eq!(dose.sum(), 1.0, epsilon = 1e-13);
    }

    #[test]
    fn off_axis_neighbour_receives_lateral_penumbra() {
        // A diagonal neighbour of the point source gets non-zero dose only because
        // the kernel spreads on all three axes (penumbra).
        let k = symmetric_deposition_kernel(length_cm(0.5), length_cm(0.2), 1);
        let dose = scatter_superposition(&point_terma(), &k, &k, &k);
        assert!(
            dose.get(4, 4, 4).unwrap() > 0.0,
            "diagonal voxel must be lit"
        );
    }

    #[test]
    fn superposition_is_linear_in_terma() {
        let k = symmetric_deposition_kernel(length_cm(0.4), length_cm(0.2), 2);
        let terma = point_terma();
        let scaled = Volume::from_shape_fn(*terma.grid(), |idx| {
            3.0 * terma.get(idx[0], idx[1], idx[2]).unwrap()
        });
        let d1 = scatter_superposition(&terma, &k, &k, &k);
        let d3 = scatter_superposition(&scaled, &k, &k, &k);
        assert_relative_eq!(d3.sum(), 3.0 * d1.sum(), epsilon = 1e-13);
        assert_relative_eq!(
            d3.get(3, 3, 3).unwrap(),
            3.0 * d1.get(3, 3, 3).unwrap(),
            epsilon = 1e-13
        );
    }

    #[test]
    fn symmetric_kernel_is_normalized_and_peaked_at_centre() {
        let k = symmetric_deposition_kernel(length_cm(0.5), length_cm(0.2), 3);
        assert_eq!(k.len(), 7);
        let sum: f64 = k.iter().sum();
        assert_relative_eq!(sum, 1.0, epsilon = 1e-15);
        // Centre tap (index 3) is the maximum; symmetric about it.
        assert!(k[3] > k[2] && k[2] > k[1]);
        assert_relative_eq!(k[2], k[4], epsilon = 1e-15);
        assert_relative_eq!(k[0], k[6], epsilon = 1e-15);
    }

    /// Energy conservation sums the whole 5x5x5 dose grid, so 125 additions
    /// accumulate against a unit total; worst-case growth is `n * T::EPSILON` with
    /// `n = 125`. Two hundred fifty-six ulps bounds that, and is still an order of
    /// magnitude tighter than the 1e-5 absolute bound it replaces.
    const ENERGY_CONSERVATION_ULPS: f64 = 256.0;

    /// Asserts interior point-source energy conservation in one scalar width.
    fn scatter_superposition_conserves_interior_energy<T: ShippedScalar>() {
        let zero = T::from_f64(0.0);
        let one = T::from_f64(1.0);
        let spacing = T::from_f64(2.0);
        let grid =
            VoxelGrid::<T>::axis_aligned([5, 5, 5], [spacing; 3], Point3::new(zero, zero, zero))
                .expect("valid axis-aligned grid");

        let terma = Volume::from_shape_fn(grid, |idx| if idx == [2, 2, 2] { one } else { zero });
        let kernel = symmetric_deposition_kernel(
            length_cm(T::from_f64(0.5)),
            length_cm(T::from_f64(0.2)),
            2,
        );
        let dose = scatter_superposition(&terma, &kernel, &kernel, &kernel);

        // A normalized kernel redistributes a unit interior source without loss.
        assert_relative_eq!(
            dose.sum(),
            one,
            max_relative = T::EPSILON * T::from_f64(ENERGY_CONSERVATION_ULPS)
        );
    }

    #[test]
    fn scatter_superposition_conserves_interior_energy_in_single_precision() {
        scatter_superposition_conserves_interior_energy::<f32>();
    }

    #[test]
    fn scatter_superposition_conserves_interior_energy_in_double_precision() {
        scatter_superposition_conserves_interior_energy::<f64>();
    }
