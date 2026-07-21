//! Voxel-grid geometry: the discrete-index ↔ world/patient coordinate map.

use helios_core::HeliosError;
use helios_math::{Isometry3, NumericElement, Point3, Scalar, Translation3, UnitQuaternion};

/// A regular, oriented 3-D voxel grid with anisotropic spacing.
///
/// The grid maps **continuous voxel-index coordinates** `(i, j, k)` — integer
/// indices land on voxel-center sample nodes — to **world/patient coordinates**
/// in millimetres:
///
/// ```text
/// world = pose · (i·sx, j·sy, k·sz)
/// ```
///
/// `pose` is an [`Isometry3`] from local scaled-index space to world/patient
/// space. [`Self::axis_aligned`] uses the identity rotation; [`Self::oriented`]
/// preserves a validated DICOM or other acquisition orientation without copying
/// volume data.
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
    /// Construct an axis-aligned grid whose voxel `(0,0,0)` center sits at
    /// `origin` (world/patient mm) with per-axis `spacing` (mm).
    ///
    /// # Errors
    /// Returns [`HeliosError::InvalidDomainValue`] if any dimension is zero, the
    /// total voxel count overflows `usize`, or any spacing is non-finite or not
    /// strictly positive.
    pub fn axis_aligned(
        dims: [usize; 3],
        spacing: [T; 3],
        origin: Point3<T>,
    ) -> Result<Self, HeliosError> {
        Self::oriented(dims, spacing, origin, UnitQuaternion::identity())
    }

    /// Construct an oriented grid whose voxel `(0,0,0)` center sits at
    /// `origin` (world/patient mm), with local-axis `spacing` (mm) rotated by
    /// `rotation` into world/patient space.
    ///
    /// The rotation is a Leto [`UnitQuaternion`], so affine scales, reflections,
    /// and non-orthogonal bases cannot enter a grid pose. External direction
    /// cosines are validated at their input boundary with
    /// [`UnitQuaternion::try_from_rotation_columns`].
    ///
    /// # Errors
    /// Returns [`HeliosError::InvalidDomainValue`] if any dimension is zero, the
    /// total voxel count overflows `usize`, or any spacing is non-finite or not
    /// strictly positive.
    pub fn oriented(
        dims: [usize; 3],
        spacing: [T; 3],
        origin: Point3<T>,
        rotation: UnitQuaternion<T>,
    ) -> Result<Self, HeliosError> {
        Self::validate_layout(dims, spacing)?;
        Ok(Self {
            dims,
            spacing,
            pose: Isometry3::from_parts(Translation3::from_vector(origin.coords), rotation),
        })
    }

    fn validate_layout(dims: [usize; 3], spacing: [T; 3]) -> Result<(), HeliosError> {
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
        Ok(())
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

    /// World/patient position of voxel `(0,0,0)`'s center.
    #[must_use]
    pub fn origin(&self) -> Point3<T> {
        Point3::from_coords(self.pose.translation)
    }

    /// Local-index-to-world/patient rigid pose.
    #[must_use]
    pub fn pose(&self) -> Isometry3<T> {
        self.pose
    }

    /// Map continuous voxel-index coordinates to world/patient coordinates (mm).
    #[must_use]
    pub fn index_to_world(&self, index: Point3<T>) -> Point3<T> {
        self.pose.transform_point(Point3::new(
            index.x * self.spacing[0],
            index.y * self.spacing[1],
            index.z * self.spacing[2],
        ))
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
        self.index_to_world(Point3::new(
            T::from_f64(i as f64),
            T::from_f64(j as f64),
            T::from_f64(k as f64),
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use eunomia::assert_relative_eq;
    use helios_math::Vector3;

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

    fn oriented_grid_preserves_index_world_contract<T: Scalar>(tolerance: T) {
        let zero = T::from_f64(0.0);
        let one = T::from_f64(1.0);
        let rotation = UnitQuaternion::try_from_rotation_columns(
            Vector3::new(zero, one, zero),
            Vector3::new(-one, zero, zero),
            Vector3::new(zero, zero, one),
            tolerance,
        )
        .expect("a right-handed quarter-turn basis is a valid grid pose");
        let grid = VoxelGrid::oriented(
            [4, 5, 6],
            [T::from_f64(2.0), T::from_f64(3.0), T::from_f64(4.0)],
            Point3::new(T::from_f64(10.0), T::from_f64(20.0), T::from_f64(30.0)),
            rotation,
        )
        .expect("valid oriented grid");

        let world = grid.voxel_center(1, 1, 1);
        let expected = [T::from_f64(7.0), T::from_f64(22.0), T::from_f64(34.0)];
        for (&actual, expected) in world.coords.data.iter().zip(expected) {
            assert!((actual - expected).abs() <= tolerance);
        }

        let index = Point3::new(T::from_f64(1.5), T::from_f64(0.75), T::from_f64(2.25));
        let round_trip = grid.world_to_index(grid.index_to_world(index));
        for (&actual, &expected) in round_trip.coords.data.iter().zip(index.coords.data.iter()) {
            assert!((actual - expected).abs() <= tolerance);
        }
    }

    #[test]
    fn oriented_grid_preserves_f32_index_world_contract() {
        oriented_grid_preserves_index_world_contract(1.0e-5_f32);
    }

    #[test]
    fn oriented_grid_preserves_f64_index_world_contract() {
        oriented_grid_preserves_index_world_contract(1.0e-12_f64);
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
}
