//! Voxel-grid geometry: the discrete-index ↔ world/patient coordinate map.

use helios_core::HeliosError;
use helios_math::{Isometry3, NumericElement, Point3, Scalar, Translation3, UnitQuaternion};

/// A regular 3-D voxel grid with an anisotropic spacing and a rigid patient pose.
///
/// The grid defines a mapping between **continuous voxel-index coordinates**
/// `(i, j, k)` — where integer indices land on voxel-center sample nodes — and
/// **world/patient coordinates** in millimetres:
///
/// ```text
/// world = pose · (i·sx, j·sy, k·sz)
/// ```
///
/// `pose` is a rigid transform ([`Isometry3`]) carrying the grid's origin
/// (translation) and orientation (rotation / DICOM direction cosines). The
/// per-axis spacing `(sx, sy, sz)` is the anisotropic scale applied before the
/// rigid transform, so the full index→world map is affine.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct VoxelGrid<T: Scalar> {
    dims: [usize; 3],
    spacing: [T; 3],
    pose: Isometry3<T>,
}

const AXIS_FIELD: [&str; 3] = [
    "VoxelGrid::spacing.x",
    "VoxelGrid::spacing.y",
    "VoxelGrid::spacing.z",
];
const DIM_FIELD: [&str; 3] = [
    "VoxelGrid::dims.x",
    "VoxelGrid::dims.y",
    "VoxelGrid::dims.z",
];

impl<T: Scalar> VoxelGrid<T> {
    /// Construct a grid from dimensions, per-axis spacing (mm), and a rigid pose.
    ///
    /// # Errors
    /// Returns [`HeliosError::InvalidDomainValue`] if any dimension is zero, the
    /// total voxel count overflows `usize`, or any spacing is non-finite or not
    /// strictly positive.
    pub fn new(dims: [usize; 3], spacing: [T; 3], pose: Isometry3<T>) -> Result<Self, HeliosError> {
        for (axis, &d) in dims.iter().enumerate() {
            if d == 0 {
                return Err(HeliosError::InvalidDomainValue {
                    field: DIM_FIELD[axis],
                    value: 0.0,
                    reason: "voxel-grid dimension must be non-zero",
                });
            }
        }
        // Reject grids whose voxel count would overflow addressable memory.
        dims[0]
            .checked_mul(dims[1])
            .and_then(|p| p.checked_mul(dims[2]))
            .ok_or(HeliosError::InvalidDomainValue {
                field: "VoxelGrid::num_voxels",
                value: f64::INFINITY,
                reason: "total voxel count overflows usize",
            })?;
        for (axis, &s) in spacing.iter().enumerate() {
            if !s.is_finite() {
                return Err(HeliosError::InvalidDomainValue {
                    field: AXIS_FIELD[axis],
                    value: s.to_f64(),
                    reason: "spacing must be finite",
                });
            }
            if s <= <T as NumericElement>::ZERO {
                return Err(HeliosError::InvalidDomainValue {
                    field: AXIS_FIELD[axis],
                    value: s.to_f64(),
                    reason: "spacing must be strictly positive",
                });
            }
        }
        Ok(Self {
            dims,
            spacing,
            pose,
        })
    }

    /// Construct an axis-aligned grid whose axes coincide with the world axes and
    /// whose voxel `(0,0,0)` center sits at `origin`.
    ///
    /// # Errors
    /// As [`VoxelGrid::new`].
    pub fn axis_aligned(
        dims: [usize; 3],
        spacing: [T; 3],
        origin: Point3<T>,
    ) -> Result<Self, HeliosError> {
        let pose = Isometry3::from_parts(
            Translation3::new(origin.x, origin.y, origin.z),
            UnitQuaternion::identity(),
        );
        Self::new(dims, spacing, pose)
    }

    /// Grid dimensions `(nx, ny, nz)`.
    #[must_use]
    pub const fn dims(&self) -> [usize; 3] {
        self.dims
    }

    /// Total number of voxels `nx·ny·nz` (validated non-overflowing at construction).
    #[must_use]
    pub fn num_voxels(&self) -> usize {
        self.dims[0] * self.dims[1] * self.dims[2]
    }

    /// Per-axis spacing `(sx, sy, sz)` in millimetres.
    #[must_use]
    pub fn spacing(&self) -> [T; 3] {
        self.spacing
    }

    /// The rigid patient pose (origin translation + orientation).
    #[must_use]
    pub fn pose(&self) -> Isometry3<T> {
        self.pose
    }

    /// Map continuous voxel-index coordinates to world/patient coordinates (mm).
    #[must_use]
    pub fn index_to_world(&self, index: Point3<T>) -> Point3<T> {
        let metric = Point3::new(
            index.x * self.spacing[0],
            index.y * self.spacing[1],
            index.z * self.spacing[2],
        );
        self.pose.transform_point(metric)
    }

    /// Map world/patient coordinates (mm) to continuous voxel-index coordinates.
    #[must_use]
    pub fn world_to_index(&self, world: Point3<T>) -> Point3<T> {
        let local = self.pose.inverse().transform_point(world);
        Point3::new(
            local.x * self.spacing[0].recip(),
            local.y * self.spacing[1].recip(),
            local.z * self.spacing[2].recip(),
        )
    }

    /// World-space center of integer voxel `(i, j, k)`.
    #[must_use]
    pub fn voxel_center(&self, i: usize, j: usize, k: usize) -> Point3<T> {
        let index = Point3::new(
            T::from_f64(i as f64),
            T::from_f64(j as f64),
            T::from_f64(k as f64),
        );
        self.index_to_world(index)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    fn grid() -> VoxelGrid<f64> {
        VoxelGrid::axis_aligned([4, 5, 6], [2.0, 3.0, 4.0], Point3::new(10.0, 20.0, 30.0))
            .expect("valid grid")
    }

    #[test]
    fn rejects_zero_dimension_and_bad_spacing() {
        assert!(
            VoxelGrid::axis_aligned([0, 1, 1], [1.0, 1.0, 1.0], Point3::new(0.0, 0.0, 0.0))
                .is_err()
        );
        assert!(
            VoxelGrid::axis_aligned([1, 1, 1], [0.0, 1.0, 1.0], Point3::new(0.0, 0.0, 0.0))
                .is_err()
        );
        assert!(
            VoxelGrid::axis_aligned([1, 1, 1], [1.0, -1.0, 1.0], Point3::new(0.0, 0.0, 0.0))
                .is_err()
        );
        assert!(VoxelGrid::axis_aligned(
            [1, 1, 1],
            [1.0, 1.0, f64::NAN],
            Point3::new(0.0, 0.0, 0.0)
        )
        .is_err());
    }

    #[test]
    fn num_voxels_is_product_of_dims() {
        assert_eq!(grid().num_voxels(), 4 * 5 * 6);
    }

    #[test]
    fn axis_aligned_index_to_world_applies_spacing_and_origin() {
        let g = grid();
        // Voxel (1,1,1): world = origin + (1·2, 1·3, 1·4) = (12, 23, 34).
        let w = g.voxel_center(1, 1, 1);
        assert_eq!(w, Point3::new(12.0, 23.0, 34.0));
        // Origin voxel maps to the origin.
        assert_eq!(g.voxel_center(0, 0, 0), Point3::new(10.0, 20.0, 30.0));
    }

    #[test]
    fn world_index_round_trip_is_identity() {
        let g = grid();
        let idx = Point3::new(2.5, 1.25, 3.75);
        let back = g.world_to_index(g.index_to_world(idx));
        assert_relative_eq!(back.x, idx.x, epsilon = 1e-12);
        assert_relative_eq!(back.y, idx.y, epsilon = 1e-12);
        assert_relative_eq!(back.z, idx.z, epsilon = 1e-12);
    }

    #[test]
    fn rotated_pose_round_trips() {
        // 90° rotation about z: (x,y,z) -> (-y, x, z), plus a translation.
        let axis = helios_math::Vector3::new(0.0_f64, 0.0, 1.0);
        let rot = UnitQuaternion::from_axis_angle(
            helios_math::UnitVector3::new_normalize(axis),
            core::f64::consts::FRAC_PI_2,
        );
        let pose = Isometry3::from_parts(Translation3::new(1.0, 2.0, 3.0), rot);
        let g = VoxelGrid::new([3, 3, 3], [1.0, 1.0, 1.0], pose).expect("valid");
        let idx = Point3::new(2.0, 0.0, 1.0);
        let back = g.world_to_index(g.index_to_world(idx));
        assert_relative_eq!(back.x, idx.x, epsilon = 1e-12);
        assert_relative_eq!(back.y, idx.y, epsilon = 1e-12);
        assert_relative_eq!(back.z, idx.z, epsilon = 1e-12);
    }
}
