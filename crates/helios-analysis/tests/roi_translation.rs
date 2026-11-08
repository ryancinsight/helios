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
    let volume = Volume::from_shape_fn(grid, |[i, j, k]| (4 * i + 2 * j + k) as f64);
    let empty = roi_statistics(&volume, [1, 0, 0], [1, 2, 2]);
    assert_eq!((empty.mean, empty.std), (0.0, 0.0));
    let reversed = roi_statistics(&volume, [2, 2, 2], [1, 1, 1]);
    assert_eq!((reversed.mean, reversed.std), (0.0, 0.0));
    let singleton = roi_statistics(&volume, [1; 3], [2; 3]);
    assert_eq!((singleton.mean, singleton.std), (7.0, 0.0));
    let pair = roi_statistics(&volume, [1, 1, 0], [3; 3]);
    assert_eq!((pair.mean, pair.std), (6.5, 0.5));
}
