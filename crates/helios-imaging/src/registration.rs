//! IGRT rigid setup correction: integer-voxel translation registration.
//!
//! Aligns a daily image (e.g. an MVCT reconstruction) to a planning reference by
//! finding the whole-voxel displacement that minimizes the mean squared intensity
//! difference over their overlap — the setup-error / couch-shift estimate that
//! image-guided radiation therapy applies before delivery.
//!
//! The exhaustive search is `ritk-registration`'s
//! [`classical::translation`](ritk_registration::classical::translation) kernel:
//! two borrowed flat buffers, a zero-sized metric policy resolved statically, and
//! no allocation. Helios supplies the [`Volume`] views and the displacement sign
//! convention, and re-exports the kernel's typed error as
//! [`TranslationRegistrationError`] so callers need not name the provider.
//! Sub-voxel interpolation, rotation, and deformable registration (mutual
//! information, also `ritk`) extend it; the whole-voxel search here is the
//! deterministic, analytically verifiable base case (recovering a known applied
//! shift exactly).

use helios_domain::Volume;
use helios_math::GeometryScalar;
use ritk_registration::classical::translation::{
    self, MeanSquaredDifference, NormalizedCrossCorrelation,
};

/// Why a translation search produced no displacement.
///
/// Re-exported from `ritk-registration`, which owns the search and the checks
/// that precede it: a shape that does not match the buffers, a search radius that
/// is not representable as a signed offset, a non-finite voxel, or a search in
/// which every candidate had an undefined metric.
pub use ritk_registration::classical::translation::TranslationRegistrationError;

/// Estimate the integer-voxel displacement `s` of `moving` relative to `fixed`.
///
/// Returns the `s ∈ [−max_shift, max_shift]` (per axis) minimizing the mean
/// squared difference `mean_v (moving(v) − fixed(v − s))²` over the voxels where
/// both samples exist. For a `moving` image that is `fixed` translated by `s`,
/// the minimum (zero) is at exactly that `s` — so `s` is the setup displacement to
/// correct. `fixed` and `moving` are assumed to share a grid.
///
/// # Errors
///
/// [`TranslationRegistrationError::ShapeMismatch`] when the grid dimensions do not
/// account for both buffers (including a dimension product that overflows),
/// [`TranslationRegistrationError::SearchRadiusOverflow`] when a radius exceeds
/// `isize::MAX`, [`TranslationRegistrationError::NonFiniteInput`] when either
/// volume holds a NaN or infinite voxel — every candidate cost is then undefined,
/// so the volume is reported rather than silently scored — and
/// [`TranslationRegistrationError::UndefinedMetric`] when no candidate had a
/// non-empty overlap.
///
/// The mean-over-overlap SSD assumes **textured** images (real CT/MVCT), where any
/// misalignment leaves residual structure; on a near-flat image a large shift that
/// slides all features out of the overlap can tie the true minimum. Use
/// [`register_translation_ncc`] for that case.
///
/// Cost: `∏(2·max_shift + 1)` candidate shifts × overlap voxels (exhaustive); a
/// coarse-to-fine search and `ritk` mutual-information registration scale it up.
pub fn register_translation<T: GeometryScalar>(
    fixed: &Volume<T>,
    moving: &Volume<T>,
    max_shift: [usize; 3],
) -> Result<[isize; 3], TranslationRegistrationError> {
    translation::register_translation::<T, MeanSquaredDifference>(
        fixed.as_slice(),
        moving.as_slice(),
        fixed.grid().dims(),
        max_shift,
    )
}

/// Normalized-cross-correlation (NCC) variant of [`register_translation`],
/// **robust on low-texture images**.
///
/// Returns the integer-voxel displacement `s ∈ [−max_shift, max_shift]` maximizing
/// the NCC over the overlap:
/// `NCC(s) = Σ(m−m̄)(f−f̄) / (N·σ_m·σ_f)`, `m = moving(v)`, `f = fixed(v − s)`.
/// Because NCC measures *correlation* (invariant to intensity offset/scale), a
/// shift that slides all structure out of the overlap leaves a near-constant
/// (zero-variance) region — which the metric rejects rather than scoring as a
/// perfect match, curing the SSD false-minimum on near-flat images (the H-044
/// limitation). Candidates whose variance product is not positive are skipped; if
/// no shift yields a valid correlation the search reports
/// [`TranslationRegistrationError::UndefinedMetric`].
///
/// # Errors
///
/// As [`register_translation`], except that an all-zero-variance search reaches
/// [`TranslationRegistrationError::UndefinedMetric`] rather than returning a
/// displacement.
pub fn register_translation_ncc<T: GeometryScalar>(
    fixed: &Volume<T>,
    moving: &Volume<T>,
    max_shift: [usize; 3],
) -> Result<[isize; 3], TranslationRegistrationError> {
    translation::register_translation::<T, NormalizedCrossCorrelation>(
        fixed.as_slice(),
        moving.as_slice(),
        fixed.grid().dims(),
        max_shift,
    )
}

#[cfg(test)]
mod tests {
    #![expect(
        clippy::unwrap_used,
        reason = "ratchet HELIOS-UNWRAP-1: pre-existing debt"
    )]
    use super::*;
    use helios_domain::VoxelGrid;
    use helios_math::Point3;
    use helios_math::ShippedScalar;

    // A textured phantom: a paraboloid bowl centred at `(cx, cy)`. Every voxel
    // carries distinct signal (like a real image), so `moving = fixed` translated
    // by `s` has a unique SSD minimum at exactly `s` — no flat region to slide a
    // feature out of. Two independent quadratic terms → two independent linear
    // constraints → a unique two-axis minimum.
    fn bowl(cx: f64, cy: f64) -> Volume<f64> {
        let grid = VoxelGrid::axis_aligned([9, 9, 1], [2.0, 2.0, 2.0], Point3::new(0.0, 0.0, 0.0))
            .unwrap();
        Volume::from_shape_fn(grid, move |idx| {
            let (di, dj) = (idx[0] as f64 - cx, idx[1] as f64 - cy);
            di * di + dj * dj
        })
    }

    #[test]
    fn recovers_a_known_applied_shift() {
        // moving bowl centred at (4,4), fixed at (2,5) → displacement (2,−1,0).
        let fixed = bowl(2.0, 5.0);
        let moving = bowl(4.0, 4.0);
        assert_eq!(
            register_translation(&fixed, &moving, [3, 3, 0]),
            Ok([2, -1, 0])
        );
    }

    #[test]
    fn identical_images_register_to_zero() {
        let fixed = bowl(4.0, 4.0);
        assert_eq!(
            register_translation(&fixed, &fixed, [2, 2, 0]),
            Ok([0, 0, 0])
        );
    }

    #[test]
    fn recovers_a_negative_shift() {
        // moving centred at (3,4), fixed at (6,6) → displacement (−3,−2,0).
        let fixed = bowl(6.0, 6.0);
        let moving = bowl(3.0, 4.0);
        assert_eq!(
            register_translation(&fixed, &moving, [3, 3, 0]),
            Ok([-3, -2, 0])
        );
    }

    /// Asserts sum-of-squares translation recovery in one scalar width.
    ///
    /// The result is an integer voxel offset, so this assertion is exact and
    /// carries no tolerance: any bound would be meaningless on a discrete search
    /// result, and an off-by-one would be a defect rather than a rounding effect.
    fn translation_registration_recovers_a_known_shift<T>()
    where
        T: ShippedScalar + helios_math::GeometryScalar,
    {
        let cast = <T as helios_math::FloatElement>::from_f64;
        let zero = cast(0.0);
        let spacing = cast(2.0);
        let grid =
            VoxelGrid::<T>::axis_aligned([9, 9, 1], [spacing; 3], Point3::new(zero, zero, zero))
                .expect("valid axis-aligned grid");

        // Quadratic bowl centred at (cx, cy): its minimum locates the shift.
        let bowl = |cx: f64, cy: f64| {
            Volume::from_shape_fn(grid, move |idx| {
                let di = cast(idx[0] as f64 - cx);
                let dj = cast(idx[1] as f64 - cy);
                di * di + dj * dj
            })
        };

        assert_eq!(
            register_translation(&bowl(3.0, 3.0), &bowl(5.0, 2.0), [3, 3, 0]),
            Ok([2, -1, 0])
        );
    }

    #[test]
    fn translation_registration_recovers_a_known_shift_in_single_precision() {
        translation_registration_recovers_a_known_shift::<f32>();
    }

    #[test]
    fn translation_registration_recovers_a_known_shift_in_double_precision() {
        translation_registration_recovers_a_known_shift::<f64>();
    }

    // Low-texture phantom: flat background 1.0 with a single bright voxel at
    // `spike` — the SSD-over-overlap case where a large shift can slide the feature
    // out of overlap and tie the true minimum (H-044 note).
    fn spike(sx: usize, sy: usize) -> Volume<f64> {
        let grid = VoxelGrid::axis_aligned([9, 9, 1], [2.0, 2.0, 2.0], Point3::new(0.0, 0.0, 0.0))
            .unwrap();
        Volume::from_shape_fn(grid, move |idx| {
            if idx[0] == sx && idx[1] == sy {
                10.0
            } else {
                1.0
            }
        })
    }

    #[test]
    fn ncc_recovers_shift_on_a_low_texture_image() {
        // moving spike at (4,2), fixed at (2,3) → displacement (2,−1,0). NCC rejects
        // the flat-overlap shifts that make plain SSD ambiguous here.
        let fixed = spike(2, 3);
        let moving = spike(4, 2);
        assert_eq!(
            register_translation_ncc(&fixed, &moving, [3, 3, 0]),
            Ok([2, -1, 0])
        );
    }

    #[test]
    fn ncc_recovers_shift_on_a_textured_image_and_zero_for_identical() {
        assert_eq!(
            register_translation_ncc(&bowl(2.0, 5.0), &bowl(4.0, 4.0), [3, 3, 0]),
            Ok([2, -1, 0])
        );
        assert_eq!(
            register_translation_ncc(&bowl(4.0, 4.0), &bowl(4.0, 4.0), [2, 2, 0]),
            Ok([0, 0, 0])
        );
    }

    /// Asserts normalized cross-correlation translation recovery in one width.
    ///
    /// Exact for the same reason as the sum-of-squares case: the search returns an
    /// integer voxel offset.
    fn ncc_registration_recovers_a_known_shift<T>()
    where
        T: ShippedScalar + helios_math::GeometryScalar,
    {
        let cast = <T as helios_math::FloatElement>::from_f64;
        let zero = cast(0.0);
        let spacing = cast(2.0);
        let grid =
            VoxelGrid::<T>::axis_aligned([9, 9, 1], [spacing; 3], Point3::new(zero, zero, zero))
                .expect("valid axis-aligned grid");

        // Single bright voxel on a unit background; correlation peaks on the shift.
        let spike = |sx: usize, sy: usize| {
            Volume::from_shape_fn(grid, move |idx| {
                if idx[0] == sx && idx[1] == sy {
                    cast(10.0)
                } else {
                    cast(1.0)
                }
            })
        };

        assert_eq!(
            register_translation_ncc(&spike(3, 3), &spike(5, 2), [3, 3, 0]),
            Ok([2, -1, 0])
        );
    }

    #[test]
    fn ncc_registration_recovers_a_known_shift_in_single_precision() {
        ncc_registration_recovers_a_known_shift::<f32>();
    }

    #[test]
    fn ncc_registration_recovers_a_known_shift_in_double_precision() {
        ncc_registration_recovers_a_known_shift::<f64>();
    }

    // ── Delegated-kernel error surface ────────────────────────────────────────
    //
    // These pin the checks the `ritk-registration` kernel performs before the
    // search. They matter because the alternative — scoring a volume that cannot
    // be scored — returns a plausible-looking displacement, which on an IGRT
    // couch-shift estimate is a silently wrong correction.

    #[test]
    fn a_volume_that_does_not_match_the_reference_grid_is_reported() {
        let fixed = bowl(4.0, 4.0);
        let grid = VoxelGrid::axis_aligned([5, 5, 1], [2.0, 2.0, 2.0], Point3::new(0.0, 0.0, 0.0))
            .unwrap();
        let moving = Volume::from_shape_fn(grid, |_| 1.0);
        assert_eq!(
            register_translation(&fixed, &moving, [1, 1, 0]),
            Err(TranslationRegistrationError::ShapeMismatch {
                dimensions: [9, 9, 1],
                expected: Some(81),
                fixed_len: 81,
                moving_len: 25,
            })
        );
    }

    #[test]
    fn a_non_finite_voxel_is_reported_rather_than_scored() {
        let grid = VoxelGrid::axis_aligned([9, 9, 1], [2.0, 2.0, 2.0], Point3::new(0.0, 0.0, 0.0))
            .unwrap();
        // Flat index of voxel (3, 0, 0) in the (i·ny + j)·nz + k layout.
        let poisoned =
            Volume::from_shape_fn(grid, |idx| if idx == [3, 0, 0] { f64::NAN } else { 1.0 });
        assert_eq!(
            register_translation(&poisoned, &poisoned, [1, 1, 0]),
            Err(TranslationRegistrationError::NonFiniteInput {
                buffer: "fixed",
                index: 27,
            })
        );
        // The moving buffer is checked too, and the label says which one failed.
        let clean = bowl(4.0, 4.0);
        assert_eq!(
            register_translation(&clean, &poisoned, [1, 1, 0]),
            Err(TranslationRegistrationError::NonFiniteInput {
                buffer: "moving",
                index: 27,
            })
        );
    }

    #[test]
    fn ncc_over_a_constant_volume_has_no_defined_metric() {
        let grid = VoxelGrid::axis_aligned([4, 4, 1], [2.0, 2.0, 2.0], Point3::new(0.0, 0.0, 0.0))
            .unwrap();
        let flat = Volume::from_shape_fn(grid, |_| 1.0);
        assert_eq!(
            register_translation_ncc(&flat, &flat, [1, 1, 0]),
            Err(TranslationRegistrationError::UndefinedMetric)
        );
    }
}
