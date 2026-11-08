//! Translation invariance of population variance at each shipped precision.
use aequitas::systems::si::{quantities::Length, units::Millimeter};
use helios_analysis::{dose_roi_statistics, roi_statistics};
use helios_domain::{Volume, VoxelGrid};
use helios_math::{Point3, ShippedScalar};

fn translated_pair_preserves_noise<T: ShippedScalar>(offset: T, separation: T) {
    let zero = T::from_f64(0.0);
    let one = T::from_f64(1.0);
    for dims in [[2, 1, 1], [2, 4, 8]] {
        let grid = VoxelGrid::axis_aligned(
            dims,
            [Length::from_unit::<Millimeter>(one); 3],
            Point3::new(zero, zero, zero),
        )
        .expect("valid grid");
        let volume = Volume::from_shape_fn(grid, |[i, j, k]| {
            if (i + j + k) % 2 == 0 {
                offset
            } else {
                offset + separation
            }
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
fn translated_pair_keeps_single_precision_noise() {
    translated_pair_preserves_noise::<f32>(1024.0, 0.25);
    translated_pair_preserves_noise::<f32>(1024.0, 1.0 / 8192.0);
    translated_pair_preserves_noise::<f32>(-1024.0, 0.25);
    translated_pair_preserves_noise::<f32>(0.0, 0.25);
    translated_pair_preserves_noise::<f32>(1024.25, -0.25);
    translated_pair_preserves_noise::<f32>(f32::from_bits(0x3f80_0001), f32::EPSILON);
}

#[test]
fn translated_pair_keeps_double_precision_noise() {
    translated_pair_preserves_noise::<f64>(134_217_728.0, 2.0);
    translated_pair_preserves_noise::<f64>(134_217_728.0, 1.0 / 33_554_432.0);
    translated_pair_preserves_noise::<f64>(-134_217_728.0, 2.0);
    translated_pair_preserves_noise::<f64>(0.0, 2.0);
    translated_pair_preserves_noise::<f64>(134_217_730.0, -2.0);
    translated_pair_preserves_noise::<f64>(f64::from_bits(0x3ff0_0000_0000_0001), f64::EPSILON);
}

#[test]
fn roi_bounds_keep_empty_singleton_and_clipped_contracts() {
    let grid = VoxelGrid::axis_aligned(
        [2, 2, 2],
        [Length::from_unit::<Millimeter>(1.0); 3],
        Point3::new(0.0, 0.0, 0.0),
    )
    .expect("valid grid");
    let volume = Volume::from_shape_fn(grid, |[i, j, k]| {
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
    let grid = VoxelGrid::axis_aligned(
        [16; 3],
        [Length::from_unit::<Millimeter>(1.0_f32); 3],
        Point3::new(0.0, 0.0, 0.0),
    )
    .expect("valid grid");
    let step = f32::from_bits(0x3dcc_ccce) - f32::from_bits(0x3dcc_cccd);
    for offset in [0.0, f32::from_bits(0x3dcc_cccd)] {
        let volume = Volume::from_shape_fn(
            grid,
            |[_, _, k]| {
                if k % 2 == 0 {
                    offset
                } else {
                    offset + step
                }
            },
        );
        let expected = f32::from_bits(0x3180_0000);
        let intensity = roi_statistics(&volume, [0; 3], [16; 3]);
        assert_eq!(intensity.std, expected);
        let dose = dose_roi_statistics(&volume, [0; 3], [16; 3]);
        assert_eq!(dose.std.into_base(), expected);
    }
}

fn symmetric_pair_keeps_finite_mean_and_noise<T: ShippedScalar>(magnitude: T) {
    let zero = T::from_f64(0.0);
    let grid = VoxelGrid::axis_aligned(
        [2, 1, 1],
        [Length::from_unit::<Millimeter>(T::from_f64(1.0)); 3],
        Point3::new(zero, zero, zero),
    )
    .expect("valid grid");
    let volume = Volume::from_shape_fn(
        grid,
        |[i, _, _]| {
            if i == 0 {
                -magnitude
            } else {
                magnitude
            }
        },
    );
    let stats = roi_statistics(&volume, [0; 3], [2, 1, 1]);
    assert_eq!(stats.mean, zero);
    assert_eq!(stats.std, magnitude);
    let dose = dose_roi_statistics(&volume, [0; 3], [2, 1, 1]);
    assert_eq!(dose.mean.into_base(), zero);
    assert_eq!(dose.std.into_base(), magnitude);
}

#[test]
fn symmetric_extremes_preserve_single_precision_range() {
    for magnitude in [
        f32::from_bits(0x7f00_0000),
        f32::MAX,
        f32::MIN_POSITIVE,
        f32::from_bits(1),
    ] {
        symmetric_pair_keeps_finite_mean_and_noise(magnitude);
    }
}

#[test]
fn symmetric_extremes_preserve_double_precision_range() {
    for magnitude in [
        f64::from_bits(0x7fe0_0000_0000_0000),
        f64::MAX,
        f64::MIN_POSITIVE,
        f64::from_bits(1),
    ] {
        symmetric_pair_keeps_finite_mean_and_noise(magnitude);
    }
}

#[test]
fn sparse_adjacent_outlier_keeps_bernoulli_noise() {
    let grid = VoxelGrid::axis_aligned(
        [16; 3],
        [Length::from_unit::<Millimeter>(1.0_f32); 3],
        Point3::new(0.0, 0.0, 0.0),
    )
    .expect("valid grid");
    for bits in [0x3dcc_cccd, 0x3f80_0000, 0x4480_0000] {
        let lower = f32::from_bits(bits);
        let upper = f32::from_bits(bits + 1);
        let volume =
            Volume::from_shape_fn(grid, |index| if index == [15; 3] { upper } else { lower });
        // One of N samples differs by d: variance = d^2 * (1/N) * (1-1/N).
        // Compensation contributes at most about 2u per sum; normalization,
        // centering and sqrt add fewer than 12 roundings for this binary input.
        let probability = 1.0_f32 / 4096.0;
        let expected = (upper - lower) * (probability * (1.0 - probability)).sqrt();
        let actual = roi_statistics(&volume, [0; 3], [16; 3]).std;
        assert!(
            (actual - expected).abs() <= 16.0 * f32::EPSILON * expected,
            "actual {actual}, expected {expected}"
        );
    }
}

fn sparse_subnormal_outlier_keeps_noise<T: ShippedScalar>(magnitude: T, expected: T) {
    let zero = T::from_f64(0.0);
    let grid = VoxelGrid::axis_aligned(
        [16; 3],
        [Length::from_unit::<Millimeter>(T::from_f64(1.0)); 3],
        Point3::new(zero, zero, zero),
    )
    .expect("valid grid");
    for offset in [zero, magnitude * T::from_f64(2.0)] {
        for step in [magnitude, -magnitude] {
            for outlier in [[0; 3], [15; 3]] {
                let volume = Volume::from_shape_fn(grid, |index| {
                    if index == outlier {
                        offset + step
                    } else {
                        offset
                    }
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
fn sparse_subnormal_outlier_preserves_single_precision() {
    sparse_subnormal_outlier_keeps_noise(f32::from_bits(0x80), f32::from_bits(0x02));
}

#[test]
fn sparse_subnormal_outlier_preserves_double_precision() {
    sparse_subnormal_outlier_keeps_noise(f64::from_bits(0x80), f64::from_bits(0x02));
}
