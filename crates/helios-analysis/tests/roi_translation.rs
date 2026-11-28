//! Translation invariance of population variance at each shipped precision.
use aequitas::systems::si::{quantities::Length, units::Millimeter};
use helios_analysis::{dose_roi_statistics, roi_statistics};
use helios_domain::{Volume, VoxelGrid};
use helios_math::{Point3, ShippedScalar};

fn unit_grid<T: ShippedScalar>(dims: [usize; 3]) -> VoxelGrid<T> {
    let zero = T::from_count(0);
    let unit = Length::from_unit::<Millimeter>(T::from_count(1));
    VoxelGrid::axis_aligned(dims, [unit; 3], Point3::new(zero, zero, zero)).expect("valid grid")
}

fn translated_pair_preserves_noise<T: ShippedScalar>(offset: T, separation: T) {
    for dims in [[2, 1, 1], [2, 4, 8]] {
        let volume = Volume::from_shape_fn(unit_grid(dims), |[i, j, k]| {
            [offset, offset + separation][(i + j + k) % 2]
        });
        // A balanced sample of a and a+d has mean a+d/2 and population
        // standard deviation |d|/2. The chosen standard deviations are
        // binary-exact. A one-ULP midpoint rounds to the even endpoint.
        let half = T::from_f64(0.5);
        let expected_mean = offset + separation * half;
        let expected_std = separation.abs() * half;
        let actual = roi_statistics(&volume, [0; 3], dims);
        // Summing 64 observations can round the mean differently; this test
        // pins its exact midpoint only for the analytically exact pair sum.
        if dims == [2, 1, 1] {
            assert_eq!(actual.mean, expected_mean);
        }
        assert_eq!(actual.std, expected_std);
        let dose = dose_roi_statistics(&volume, [0; 3], dims);
        assert_eq!(dose.mean.into_base(), actual.mean);
        assert_eq!(dose.std.into_base(), expected_std);
    }
}

#[test]
fn translated_pair_keeps_noise() {
    translated_pair_preserves_noise::<f32>(1024.0, 0.25);
    translated_pair_preserves_noise::<f32>(1024.0, 1.0 / 8192.0);
    translated_pair_preserves_noise::<f32>(-1024.0, 0.25);
    translated_pair_preserves_noise::<f32>(0.0, 0.25);
    translated_pair_preserves_noise::<f32>(1024.25, -0.25);
    translated_pair_preserves_noise::<f32>(f32::from_bits(0x3f80_0001), f32::EPSILON);
    translated_pair_preserves_noise::<f64>(134_217_728.0, 2.0);
    translated_pair_preserves_noise::<f64>(134_217_728.0, 1.0 / 33_554_432.0);
    translated_pair_preserves_noise::<f64>(-134_217_728.0, 2.0);
    translated_pair_preserves_noise::<f64>(0.0, 2.0);
    translated_pair_preserves_noise::<f64>(134_217_730.0, -2.0);
    translated_pair_preserves_noise::<f64>(f64::from_bits(0x3ff0_0000_0000_0001), f64::EPSILON);
}

#[test]
fn roi_bounds_keep_empty_singleton_and_clipped_contracts() {
    let volume = Volume::from_shape_fn(unit_grid([2, 2, 2]), |[i, j, k]| {
        f64::from(u8::try_from(4 * i + 2 * j + k).expect("2x2x2 voxel index fits in u8"))
    });
    let empty = roi_statistics(&volume, [1, 0, 0], [1, 2, 2]);
    assert_eq!((empty.mean, empty.std), (0.0, 0.0));
    let reversed = roi_statistics(&volume, [2, 2, 2], [1, 1, 1]);
    assert_eq!((reversed.mean, reversed.std), (0.0, 0.0));
    let singleton = roi_statistics(&volume, [1; 3], [2; 3]);
    assert_eq!((singleton.mean, singleton.std), (7.0, 0.0));
    let pair = roi_statistics(&volume, [1, 1, 0], [3; 3]);
    assert_eq!((pair.mean, pair.std), (6.5, 0.5));
}

#[test]
fn adjacent_decimal_values_keep_noise_across_a_volume() {
    let grid = unit_grid([16; 3]);
    let step = f32::from_bits(0x3dcc_ccce) - f32::from_bits(0x3dcc_cccd);
    for offset in [0.0, f32::from_bits(0x3dcc_cccd)] {
        let volume = Volume::from_shape_fn(grid, |[_, _, k]| [offset, offset + step][k % 2]);
        let expected = f32::from_bits(0x3180_0000);
        let intensity = roi_statistics(&volume, [0; 3], [16; 3]);
        assert_eq!(intensity.std, expected);
        let dose = dose_roi_statistics(&volume, [0; 3], [16; 3]);
        assert_eq!(dose.std.into_base(), expected);
    }
}

fn symmetric_pair_keeps_finite_mean_and_noise<T: ShippedScalar>(magnitude: T) {
    let zero = T::from_count(0);
    let volume =
        Volume::from_shape_fn(unit_grid([2, 1, 1]), |[i, _, _]| [-magnitude, magnitude][i]);
    let stats = roi_statistics(&volume, [0; 3], [2, 1, 1]);
    assert_eq!(stats.mean, zero);
    assert_eq!(stats.std, magnitude);
    let dose = dose_roi_statistics(&volume, [0; 3], [2, 1, 1]);
    assert_eq!(dose.mean.into_base(), zero);
    assert_eq!(dose.std.into_base(), magnitude);
}

#[test]
fn symmetric_extremes_preserve_range() {
    // Bits of 2^127 (2^1023 for f64), MAX, MIN_POSITIVE, the smallest subnormal.
    for bits in [0x7f00_0000_u32, 0x7f7f_ffff, 0x0080_0000, 1] {
        symmetric_pair_keeps_finite_mean_and_noise(f32::from_bits(bits));
    }
    for bits in [
        0x7fe0_0000_0000_0000_u64,
        0x7fef_ffff_ffff_ffff,
        0x0010_0000_0000_0000,
        1,
    ] {
        symmetric_pair_keeps_finite_mean_and_noise(f64::from_bits(bits));
    }
}

/// Relative error bound for the population standard deviation of lattice data.
///
/// With `u = EPSILON / 2` and `gamma_2 = 2u / (1 - 2u)`, `fl(x - c)` and the
/// division by the endpoint scale `s` perturb each normalized `y` by at most
/// `gamma_2`, while `sigma(x) = s * sigma(y)` exactly. A perturbation `e` moves
/// a standard deviation by at most `rms(e)`: relative `gamma_2 * kappa`, with
/// `kappa = rms(y) / sigma(y)` taken from the anchor, within half a lattice
/// step of the range midpoint. The compensated center error enters `sigma`
/// squared. Residual, square, sum, division by `n`, root and scale add `5u`;
/// the reference (conversion, root, division) adds `2.5u`. Underflow and terms
/// of order `n u^2`, `u^2 kappa^2` are neglected.
fn derived_std_bound<T: ShippedScalar>(kappa: T) -> T {
    let (one, two, five) = (T::from_count(1), T::from_count(2), T::from_count(5));
    let u = T::EPSILON / two;
    let gamma_2 = two * u / (one - two * u);
    gamma_2 * kappa + five * u + five * u / two
}

/// Assert `std` of `lower + step * k` over an `n = prod(dims)` lattice against
/// the exact integer-moment reference within [`derived_std_bound`].
fn lattice_noise_stays_within_derived_bound<T: ShippedScalar>(
    lower: T,
    step: T,
    dims: [usize; 3],
    lattice: &[u64],
) {
    let count_of = |value: u64| T::from_count(usize::try_from(value).expect("moment fits usize"));
    let value = |k: u64| lower + step * count_of(k);
    let exact = |&k: &u64| value(k) - lower == step * count_of(k);
    assert!(lattice.iter().all(exact), "lattice values are exact");
    let volume = Volume::from_shape_fn(unit_grid(dims), |[i, j, k]| {
        value(lattice[(i * dims[1] + j) * dims[2] + k])
    });
    // Integer moments are exact: n * S2 - S1^2 = n^2 * variance in lattice units.
    let n = u64::try_from(lattice.len()).expect("count fits u64");
    let s1 = lattice.iter().sum::<u64>();
    let s2 = lattice.iter().map(|k| k * k).sum::<u64>();
    let extremes = lattice
        .iter()
        .fold((u64::MAX, 0), |(l, h), &k| (l.min(k), h.max(k)));
    let scaled_variance = count_of(n * s2 - s1 * s1);
    // shift = 2n * (|mean - midpoint| + 1/2), the bracket in lattice steps (1/2: anchor).
    let shift = count_of((2 * s1).abs_diff(n * (extremes.0 + extremes.1)) + n);
    let kappa = (T::from_count(1) + shift * shift / (T::from_count(4) * scaled_variance)).sqrt();
    let reference = step * scaled_variance.sqrt() / count_of(n);
    let actual = roi_statistics(&volume, [0; 3], dims).std;
    let bound = derived_std_bound(kappa);
    assert!(
        (actual - reference).abs() <= bound * reference,
        "std {actual:?} differs from {reference:?} by more than relative {:e}",
        bound.to_f64()
    );
}

#[test]
fn offset_lattice_noise_stays_within_derived_bound() {
    // 37 * 41 * 43 voxels: a count that is not a power of two.
    let dims = [37, 41, 43];
    let count = dims.iter().product::<usize>();
    let mut state = 0x9e37_79b9_7f4a_7c15_u64;
    let mut draw = move |modulus: u64| {
        state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        (state >> 33) % modulus
    };
    // A product of three uniform draws skews the values toward the range start.
    let skewed: Vec<u64> = (0..count)
        .map(|_| draw(1000) * draw(1000) * draw(1000) / 1_000_000)
        .collect();
    let mut outlier: Vec<u64> = (1..count).map(|_| draw(16)).collect();
    outlier.push(4000);
    for lattice in [&skewed, &outlier] {
        // Each offset sets the lattice step to the unit in the last place.
        lattice_noise_stays_within_derived_bound(2.0_f32.powi(30), 128.0, dims, lattice);
        lattice_noise_stays_within_derived_bound(2.0_f64.powi(52), 1.0, dims, lattice);
    }
}

#[test]
fn sparse_adjacent_outlier_keeps_bernoulli_noise() {
    let lattice: Vec<u64> = (0..4096).map(|k| u64::from(k == 4095)).collect();
    for lower in [0.1_f32, 1.0, 1024.0] {
        let step = f32::from_bits(lower.to_bits() + 1) - lower;
        lattice_noise_stays_within_derived_bound(lower, step, [16; 3], &lattice);
    }
    for lower in [0.1_f64, 1.0, 1024.0] {
        let step = f64::from_bits(lower.to_bits() + 1) - lower;
        lattice_noise_stays_within_derived_bound(lower, step, [16; 3], &lattice);
    }
}

fn sparse_subnormal_outlier_keeps_noise<T: ShippedScalar>(magnitude: T, expected: T) {
    let grid = unit_grid([16; 3]);
    for offset in [T::from_count(0), magnitude * T::from_count(2)] {
        for step in [magnitude, -magnitude] {
            for outlier in [[0; 3], [15; 3]] {
                let volume = Volume::from_shape_fn(grid, |index| {
                    [offset, offset + step][usize::from(index == outlier)]
                });
                // One of 4096 samples differs by 128 minimum subnormal units.
                // sigma = 128 * sqrt((1/4096)*(4095/4096)) units, strictly
                // between 1.5 and 2 units, hence rounds to exactly 2 units.
                assert_eq!(roi_statistics(&volume, [0; 3], [16; 3]).std, expected);
                assert_eq!(
                    dose_roi_statistics(&volume, [0; 3], [16; 3])
                        .std
                        .into_base(),
                    expected
                );
            }
        }
    }
}

#[test]
fn sparse_subnormal_outlier_preserves_noise() {
    sparse_subnormal_outlier_keeps_noise(f32::from_bits(0x80), f32::from_bits(0x02));
    sparse_subnormal_outlier_keeps_noise(f64::from_bits(0x80), f64::from_bits(0x02));
}

fn non_finite_sample_makes_noise_nan<T: ShippedScalar>() {
    let (nan, inf, c) = (T::nan(), T::infinity(), T::from_count);
    let cases: [&[T]; 10] = [
        &[c(1), nan, c(2)],
        &[c(1), inf, c(2)],
        &[c(1), -inf, c(2)],
        &[c(1), nan, c(1)],
        &[nan, c(5)],
        &[c(1_000_000), nan, c(1_000_000), c(1_000_000)],
        &[nan],
        &[inf],
        &[c(3), inf, c(3)],
        &[c(3), -inf, c(3)],
    ];
    for values in cases {
        let dims = [values.len(), 1, 1];
        let volume = Volume::from_shape_vec(unit_grid(dims), values.to_vec()).expect("shape");
        let stats = roi_statistics(&volume, [0; 3], dims);
        assert!(stats.std.is_nan(), "std of {values:?}");
        assert_eq!(stats.mean.is_nan(), values.iter().any(|v| v.is_nan()));
    }
}

#[test]
fn non_finite_samples_propagate_to_nan_noise() {
    non_finite_sample_makes_noise_nan::<f32>();
    non_finite_sample_makes_noise_nan::<f64>();
}

fn overflowing_span_with_nonzero_anchor_keeps_noise<T: ShippedScalar>(quantum: T) {
    // The span 8q is 2^128 (f32) or 2^1024 (f64) and overflows; the anchor is -2q.
    let values = vec![-T::from_count(6) * quantum, T::from_count(2) * quantum];
    let volume = Volume::from_shape_vec(unit_grid([2, 1, 1]), values).expect("shape");
    assert_eq!(
        roi_statistics(&volume, [0; 3], [2, 1, 1]).std,
        T::from_count(4) * quantum
    );
}

#[test]
fn overflowing_span_keeps_noise() {
    overflowing_span_with_nonzero_anchor_keeps_noise(2.0_f32.powi(125));
    overflowing_span_with_nonzero_anchor_keeps_noise(2.0_f64.powi(1021));
}
