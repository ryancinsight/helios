//! MVCT image-quality metrics: reconstruction accuracy, noise, and contrast.
//!
//! Quantitative instruments for the imaging validation gate. Reconstruction
//! **accuracy** is measured against a ground-truth attenuation volume
//! ([`volume_rmse`], [`volume_relative_l2_error`]); **noise** is the intensity
//! standard deviation over a uniform region of interest ([`roi_statistics`]); and
//! **contrast** between two materials/ROIs is the Michelson contrast
//! ([`michelson_contrast`]) with its noise-normalized form, the contrast-to-noise
//! ratio ([`contrast_to_noise_ratio`]). All are generic over the [`Scalar`] seam.

use helios_core::HeliosError;
use helios_domain::Volume;
use helios_math::{NumericElement, Scalar};

/// Mean and (population) standard deviation of a region of interest.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RoiStats<T: Scalar> {
    /// Mean intensity over the ROI (the signal level).
    pub mean: T,
    /// Population standard deviation over the ROI (the noise level).
    pub std: T,
}

/// Absolute value via the ordered seam (`max(x, −x)`), avoiding a dependency on a
/// dedicated `abs` method.
#[inline]
fn abs<T: Scalar>(x: T) -> T {
    x.max_scalar(-x)
}

/// Mean and population standard deviation over the half-open index box
/// `[min, max)` (inclusive `min`, exclusive `max`) of `volume`.
///
/// The standard deviation over a *uniform* region is the MVCT noise metric (a
/// noiseless region yields `std = 0`). An empty box returns zero statistics.
#[must_use]
pub fn roi_statistics<T: Scalar>(
    volume: &Volume<T>,
    min: [usize; 3],
    max: [usize; 3],
) -> RoiStats<T> {
    let zero = <T as NumericElement>::ZERO;
    let mut sum = zero;
    let mut sum_sq = zero;
    let mut count = 0usize;
    for i in min[0]..max[0] {
        for j in min[1]..max[1] {
            for k in min[2]..max[2] {
                if let Some(v) = volume.get(i, j, k) {
                    sum += v;
                    sum_sq += v * v;
                    count += 1;
                }
            }
        }
    }
    if count == 0 {
        return RoiStats {
            mean: zero,
            std: zero,
        };
    }
    let inv_n = T::from_f64(count as f64).recip();
    let mean = sum * inv_n;
    // Var = E[x²] − E[x]²; clamp a tiny negative from rounding before sqrt.
    let var = (sum_sq * inv_n - mean * mean).max_scalar(zero);
    RoiStats {
        mean,
        std: var.sqrt(),
    }
}

/// Michelson contrast `(a − b) / (a + b)` between two intensity levels.
///
/// Defined for non-negative intensities not both zero (`a + b > 0`).
#[must_use]
pub fn michelson_contrast<T: Scalar>(a: T, b: T) -> T {
    (a - b) * (a + b).recip()
}

/// Contrast-to-noise ratio `|signal − background| / noise_std`.
///
/// The noise-normalized detectability of a feature against its background; the
/// standard MVCT low-contrast metric. `noise_std` must be positive.
#[must_use]
pub fn contrast_to_noise_ratio<T: Scalar>(signal_mean: T, background_mean: T, noise_std: T) -> T {
    abs(signal_mean - background_mean) * noise_std.recip()
}

/// Accumulate `(Σ(a−b)², Σb², n)` over two equally-shaped volumes.
fn error_accumulate<T: Scalar>(a: &Volume<T>, b: &Volume<T>) -> Result<(T, T, usize), HeliosError> {
    let da = a.grid().dims();
    let db = b.grid().dims();
    if da != db {
        return Err(HeliosError::InvalidDomainValue {
            field: "image_quality::volume error",
            value: (db[0] * db[1] * db[2]) as f64,
            reason: "the two volumes must have identical grid dimensions",
        });
    }
    let zero = <T as NumericElement>::ZERO;
    let (mut sq_diff, mut sq_ref) = (zero, zero);
    for i in 0..da[0] {
        for j in 0..da[1] {
            for k in 0..da[2] {
                let va = a.get(i, j, k).expect("index within grid");
                let vb = b.get(i, j, k).expect("index within grid");
                let d = va - vb;
                sq_diff += d * d;
                sq_ref += vb * vb;
            }
        }
    }
    Ok((sq_diff, sq_ref, da[0] * da[1] * da[2]))
}

/// Root-mean-square error between a reconstruction `recon` and ground truth
/// `truth` (identical grids). `0` for identical volumes.
///
/// # Errors
/// [`HeliosError::InvalidDomainValue`] if the grid dimensions differ.
pub fn volume_rmse<T: Scalar>(recon: &Volume<T>, truth: &Volume<T>) -> Result<T, HeliosError> {
    let (sq_diff, _sq_ref, n) = error_accumulate(recon, truth)?;
    Ok((sq_diff * T::from_f64(n as f64).recip()).sqrt())
}

/// Relative L2 error `‖recon − truth‖₂ / ‖truth‖₂` (identical grids). `0` for an
/// exact reconstruction; `1` when `recon` is uniformly zero against a non-zero
/// truth.
///
/// # Errors
/// [`HeliosError::InvalidDomainValue`] if the grid dimensions differ or `truth`
/// has zero norm (relative error undefined).
pub fn volume_relative_l2_error<T: Scalar>(
    recon: &Volume<T>,
    truth: &Volume<T>,
) -> Result<T, HeliosError> {
    let (sq_diff, sq_ref, _n) = error_accumulate(recon, truth)?;
    if sq_ref <= <T as NumericElement>::ZERO {
        return Err(HeliosError::InvalidDomainValue {
            field: "image_quality::volume_relative_l2_error",
            value: 0.0,
            reason: "truth has zero L2 norm; relative error is undefined",
        });
    }
    Ok((sq_diff * sq_ref.recip()).sqrt())
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use helios_domain::VoxelGrid;
    use helios_math::Point3;

    fn grid(dims: [usize; 3]) -> VoxelGrid<f64> {
        VoxelGrid::axis_aligned(dims, [1.0, 1.0, 1.0], Point3::new(0.0, 0.0, 0.0)).expect("grid")
    }

    #[test]
    fn roi_statistics_of_uniform_region_has_zero_noise() {
        let vol = Volume::from_shape_fn(grid([3, 3, 3]), |_| 7.0);
        let s = roi_statistics(&vol, [0, 0, 0], [3, 3, 3]);
        assert_relative_eq!(s.mean, 7.0, epsilon = 1e-14);
        assert_relative_eq!(s.std, 0.0, epsilon = 1e-14);
    }

    #[test]
    fn roi_statistics_matches_hand_computed_mean_and_std() {
        // Two voxels [2, 4]: mean 3, population var ((2−3)²+(4−3)²)/2 = 1, std 1.
        let vol = Volume::from_shape_vec(grid([1, 1, 2]), vec![2.0, 4.0]).unwrap();
        let s = roi_statistics(&vol, [0, 0, 0], [1, 1, 2]);
        assert_relative_eq!(s.mean, 3.0, epsilon = 1e-14);
        assert_relative_eq!(s.std, 1.0, epsilon = 1e-14);
    }

    #[test]
    fn michelson_contrast_is_the_normalized_difference() {
        assert_relative_eq!(michelson_contrast(3.0, 1.0), 0.5, epsilon = 1e-15);
        assert_relative_eq!(michelson_contrast(5.0, 5.0), 0.0, epsilon = 1e-15);
    }

    #[test]
    fn cnr_is_absolute_contrast_over_noise() {
        // |10 − 4| / 2 = 3; symmetric in the two levels.
        assert_relative_eq!(
            contrast_to_noise_ratio(10.0, 4.0, 2.0),
            3.0,
            epsilon = 1e-15
        );
        assert_relative_eq!(
            contrast_to_noise_ratio(4.0, 10.0, 2.0),
            3.0,
            epsilon = 1e-15
        );
    }

    #[test]
    fn rmse_is_zero_for_identical_and_equals_constant_offset() {
        let truth = Volume::from_shape_fn(grid([2, 2, 2]), |idx| idx[0] as f64 + idx[1] as f64);
        assert_relative_eq!(volume_rmse(&truth, &truth).unwrap(), 0.0, epsilon = 1e-15);
        // Uniform +0.5 offset → RMSE = 0.5.
        let shifted =
            Volume::from_shape_fn(grid([2, 2, 2]), |idx| idx[0] as f64 + idx[1] as f64 + 0.5);
        assert_relative_eq!(volume_rmse(&shifted, &truth).unwrap(), 0.5, epsilon = 1e-14);
    }

    #[test]
    fn relative_l2_error_matches_closed_form() {
        // truth ≡ 2, recon ≡ 2.5 over N voxels: ‖diff‖ = 0.5√N, ‖truth‖ = 2√N →
        // ratio = 0.25.
        let truth = Volume::from_shape_fn(grid([2, 2, 2]), |_| 2.0);
        let recon = Volume::from_shape_fn(grid([2, 2, 2]), |_| 2.5);
        assert_relative_eq!(
            volume_relative_l2_error(&recon, &truth).unwrap(),
            0.25,
            epsilon = 1e-14
        );
        assert_relative_eq!(
            volume_relative_l2_error(&truth, &truth).unwrap(),
            0.0,
            epsilon = 1e-15
        );
    }

    #[test]
    fn mismatched_grids_and_zero_truth_norm_error() {
        let a = Volume::from_shape_fn(grid([2, 2, 2]), |_| 1.0);
        let b = Volume::from_shape_fn(grid([2, 2, 3]), |_| 1.0);
        assert!(volume_rmse(&a, &b).is_err());
        let zero_truth = Volume::from_shape_fn(grid([2, 2, 2]), |_| 0.0);
        assert!(volume_relative_l2_error(&a, &zero_truth).is_err());
    }

    #[test]
    fn metrics_are_generic_over_scalar_f32() {
        let vol = Volume::from_shape_vec(
            VoxelGrid::<f32>::axis_aligned([1, 1, 2], [1.0, 1.0, 1.0], Point3::new(0.0, 0.0, 0.0))
                .unwrap(),
            vec![2.0_f32, 4.0],
        )
        .unwrap();
        let s = roi_statistics(&vol, [0, 0, 0], [1, 1, 2]);
        assert_relative_eq!(s.std, 1.0_f32, epsilon = 1e-6);
        assert_relative_eq!(michelson_contrast(3.0_f32, 1.0), 0.5, epsilon = 1e-7);
    }
}
