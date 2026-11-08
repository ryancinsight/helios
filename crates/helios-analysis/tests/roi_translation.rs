//! Translation invariance of population variance at each shipped precision.
use aequitas::systems::si::{quantities::Length, units::Millimeter};
use helios_analysis::{dose_roi_statistics, roi_statistics};
use helios_domain::{Volume, VoxelGrid};
use helios_math::{Point3, ShippedScalar};

fn translated_pair_preserves_noise<T: ShippedScalar>(offset: T, separation: T) {
    let zero = T::from_f64(0.0);
    let one = T::from_f64(1.0);
    let grid = VoxelGrid::axis_aligned(
        [2, 1, 1],
        [Length::from_unit::<Millimeter>(one); 3],
        Point3::new(zero, zero, zero),
    )
    .expect("valid grid");
    let volume = Volume::from_shape_fn(
        grid,
        |[i, _, _]| {
            if i == 0 {
                offset
            } else {
                offset + separation
            }
        },
    );
    // For the two-point sample [a, a + d], the exact mean is a + d/2,
    // and population standard deviation is |d|/2. All chosen values and
    // these results are binary-exact at the tested precision.
    let half = T::from_f64(0.5);
    let expected_mean = offset + separation * half;
    let expected_std = separation * half;
    let actual = roi_statistics(&volume, [0; 3], [2, 1, 1]);
    assert_eq!(actual.mean, expected_mean);
    assert_eq!(actual.std, expected_std);
    let dose = dose_roi_statistics(&volume, [0; 3], [2, 1, 1]);
    assert_eq!(dose.mean.into_base(), expected_mean);
    assert_eq!(dose.std.into_base(), expected_std);
}

#[test]
fn translated_pair_keeps_single_precision_noise() {
    translated_pair_preserves_noise::<f32>(1024.0, 0.25);
    translated_pair_preserves_noise::<f32>(-1024.0, 0.25);
    translated_pair_preserves_noise::<f32>(0.0, 0.25);
}

#[test]
fn translated_pair_keeps_double_precision_noise() {
    translated_pair_preserves_noise::<f64>(134_217_728.0, 2.0);
    translated_pair_preserves_noise::<f64>(-134_217_728.0, 2.0);
    translated_pair_preserves_noise::<f64>(0.0, 2.0);
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
