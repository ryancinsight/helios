//! IGRT rigid setup correction: integer-voxel translation registration.
//!
//! Aligns a daily image (e.g. an MVCT reconstruction) to a planning reference by
//! finding the whole-voxel displacement that minimizes the mean squared intensity
//! difference over their overlap — the setup-error / couch-shift estimate that
//! image-guided radiation therapy applies before delivery.
//!
//! This is the translation-only, whole-voxel core. Sub-voxel interpolation,
//! rotation, and deformable registration (mutual-information, via `ritk`) extend
//! it; the exhaustive whole-voxel search here is the deterministic, analytically
//! verifiable base case (recovering a known applied shift exactly).

use helios_domain::Volume;
use helios_math::GeometryScalar;

/// Estimate the integer-voxel displacement `s` of `moving` relative to `fixed`.
///
/// Returns the `s ∈ [−max_shift, max_shift]` (per axis) minimizing the mean
/// squared difference `mean_v (moving(v) − fixed(v − s))²` over the voxels where
/// both samples exist. For a `moving` image that is `fixed` translated by `s`,
/// the minimum (zero) is at exactly that `s` — so `s` is the setup displacement to
/// correct. `fixed` and `moving` are assumed to share a grid.
///
/// The mean-over-overlap SSD assumes **textured** images (real CT/MVCT), where any
/// misalignment leaves residual structure; on a near-flat image a large shift that
/// slides all features out of the overlap can tie the true minimum. A masked /
/// normalized-cross-correlation metric (H-044b) removes that assumption.
///
/// Cost: `∏(2·max_shift + 1)` candidate shifts × overlap voxels (exhaustive); a
/// coarse-to-fine search and `ritk` mutual-information registration scale it up.
#[must_use]
pub fn register_translation<T: GeometryScalar>(
    fixed: &Volume<T>,
    moving: &Volume<T>,
    max_shift: [usize; 3],
) -> [isize; 3] {
    let dims = fixed.grid().dims();
    let r = [
        max_shift[0] as isize,
        max_shift[1] as isize,
        max_shift[2] as isize,
    ];

    let mut best = [0isize; 3];
    let mut best_cost = f64::INFINITY;
    for s0 in -r[0]..=r[0] {
        for s1 in -r[1]..=r[1] {
            for s2 in -r[2]..=r[2] {
                let mut ssd = 0.0f64;
                let mut count = 0usize;
                for i in 0..dims[0] {
                    for j in 0..dims[1] {
                        for k in 0..dims[2] {
                            let (fi, fj, fk) = (i as isize - s0, j as isize - s1, k as isize - s2);
                            if fi < 0 || fj < 0 || fk < 0 {
                                continue;
                            }
                            let (Some(m), Some(f)) = (
                                moving.get(i, j, k),
                                fixed.get(fi as usize, fj as usize, fk as usize),
                            ) else {
                                continue;
                            };
                            let d = m.to_f64() - f.to_f64();
                            ssd += d * d;
                            count += 1;
                        }
                    }
                }
                if count == 0 {
                    continue;
                }
                let cost = ssd / count as f64;
                if cost < best_cost {
                    best_cost = cost;
                    best = [s0, s1, s2];
                }
            }
        }
    }
    best
}

#[cfg(test)]
mod tests {
    use super::*;
    use helios_domain::VoxelGrid;
    use helios_math::Point3;

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
        assert_eq!(register_translation(&fixed, &moving, [3, 3, 0]), [2, -1, 0]);
    }

    #[test]
    fn identical_images_register_to_zero() {
        let fixed = bowl(4.0, 4.0);
        assert_eq!(register_translation(&fixed, &fixed, [2, 2, 0]), [0, 0, 0]);
    }

    #[test]
    fn recovers_a_negative_shift() {
        // moving centred at (3,4), fixed at (6,6) → displacement (−3,−2,0).
        let fixed = bowl(6.0, 6.0);
        let moving = bowl(3.0, 4.0);
        assert_eq!(
            register_translation(&fixed, &moving, [3, 3, 0]),
            [-3, -2, 0]
        );
    }

    #[test]
    fn registration_is_generic_over_scalar_f32() {
        let grid =
            VoxelGrid::<f32>::axis_aligned([9, 9, 1], [2.0, 2.0, 2.0], Point3::new(0.0, 0.0, 0.0))
                .unwrap();
        let bowl_f32 = |cx: f32, cy: f32| {
            Volume::from_shape_fn(grid, move |idx| {
                let (di, dj) = (idx[0] as f32 - cx, idx[1] as f32 - cy);
                di * di + dj * dj
            })
        };
        assert_eq!(
            register_translation(&bowl_f32(3.0, 3.0), &bowl_f32(5.0, 2.0), [3, 3, 0]),
            [2, -1, 0]
        );
    }
}
