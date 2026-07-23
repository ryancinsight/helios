//! Geometric region-of-interest (ROI) masks for structure-masked analysis.
//!
//! Produce voxel-index predicates for common analytic ROI shapes (sphere, box)
//! over a [`VoxelGrid`], to feed [`Dvh::from_volume_masked`](crate::Dvh::from_volume_masked)
//! and other masked metrics. This gives per-structure DVH/statistics on
//! geometrically-defined structures (phantom QA inserts, simple planning targets)
//! without a hand-written closure at each call site; contour-defined ROIs (via a
//! `ritk` RT-struct rasterization) are a separate path (H-032d).

use helios_domain::VoxelGrid;
use helios_math::{Point3, Scalar};

/// A voxel-index predicate selecting voxels whose **centre** lies within `radius`
/// (world units) of `centre` on `grid` — a spherical ROI.
///
/// The returned closure is `Fn([usize; 3]) -> bool`, directly usable as the mask
/// for [`Dvh::from_volume_masked`](crate::Dvh::from_volume_masked).
pub fn spherical_mask<T: Scalar>(
    grid: VoxelGrid<T>,
    centre: Point3<T>,
    radius: T,
) -> impl Fn([usize; 3]) -> bool {
    let r_sq = radius * radius;
    move |idx| {
        let c = grid.voxel_center(idx[0], idx[1], idx[2]);
        (c - centre).norm_squared() <= r_sq
    }
}

/// A voxel-index predicate selecting voxels whose **centre** lies inside the closed
/// axis-aligned world box `[min, max]` on `grid` — a rectangular ROI.
pub fn box_mask<T: Scalar>(
    grid: VoxelGrid<T>,
    min: Point3<T>,
    max: Point3<T>,
) -> impl Fn([usize; 3]) -> bool {
    move |idx| {
        let c = grid.voxel_center(idx[0], idx[1], idx[2]);
        c.x >= min.x && c.x <= max.x && c.y >= min.y && c.y <= max.y && c.z >= min.z && c.z <= max.z
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Dvh;
    use eunomia::assert_relative_eq;
    use helios_domain::Volume;

    fn grid() -> VoxelGrid<f64> {
        // 5×5×1, 2 mm spacing → centre voxel (2,2,0) at world (4,4,0).
        VoxelGrid::axis_aligned([5, 5, 1], [2.0, 2.0, 2.0], Point3::new(0.0, 0.0, 0.0)).unwrap()
    }

    #[test]
    fn spherical_mask_selects_voxels_within_radius() {
        // Radius 2.5 mm about (4,4,0): centre + 4-connected neighbours (2 mm) in,
        // diagonals (2√2 ≈ 2.83 mm) and far voxels out.
        let mask = spherical_mask(grid(), Point3::new(4.0, 4.0, 0.0), 2.5);
        assert!(mask([2, 2, 0])); // distance 0
        assert!(mask([1, 2, 0])); // distance 2 mm
        assert!(mask([2, 3, 0])); // distance 2 mm
        assert!(!mask([1, 1, 0])); // 2.83 mm
        assert!(!mask([0, 0, 0])); // far
    }

    #[test]
    fn masked_dvh_over_spherical_roi_reads_the_roi_dose() {
        // Dose 5 inside the sphere, 1 outside; the sphere-masked DVH mean is 5.
        let g = grid();
        let mask = spherical_mask(g, Point3::new(4.0, 4.0, 0.0), 2.5);
        let dose = Volume::from_shape_fn(g, |idx| if mask(idx) { 5.0 } else { 1.0 });
        let roi =
            Dvh::from_volume_masked(&dose, spherical_mask(g, Point3::new(4.0, 4.0, 0.0), 2.5));
        assert_relative_eq!(roi.mean().into_base(), 5.0, epsilon = 1e-15);
        // 5 voxels (centre + 4 neighbours) are within the radius.
        assert_eq!(roi.count(), 5);
    }

    #[test]
    fn box_mask_selects_the_axis_aligned_region() {
        // World box [2,6]×[2,6]×[−1,1] → voxel indices 1..=3 in x and y (world 2..6).
        let g = grid();
        let mask = box_mask(g, Point3::new(2.0, 2.0, -1.0), Point3::new(6.0, 6.0, 1.0));
        assert!(mask([1, 1, 0]) && mask([3, 3, 0]) && mask([2, 2, 0]));
        assert!(!mask([0, 2, 0]) && !mask([4, 2, 0]));
        // 3×3 = 9 voxels selected.
        let dose = Volume::from_shape_fn(g, |_| 2.0);
        let roi = Dvh::from_volume_masked(&dose, mask);
        assert_eq!(roi.count(), 9);
    }

    #[test]
    fn roi_masks_are_generic_over_scalar_f32() {
        let g =
            VoxelGrid::<f32>::axis_aligned([5, 5, 1], [2.0, 2.0, 2.0], Point3::new(0.0, 0.0, 0.0))
                .unwrap();
        let mask = spherical_mask(g, Point3::new(4.0_f32, 4.0, 0.0), 2.5);
        assert!(mask([2, 2, 0]));
        assert!(!mask([0, 0, 0]));
    }
}
