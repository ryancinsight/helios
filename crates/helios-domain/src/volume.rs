//! Dense scalar volume over a [`VoxelGrid`], with trilinear sampling.

use crate::grid::VoxelGrid;
use helios_core::HeliosError;
use helios_math::{NumericElement, Point3, Scalar};
use leto::Array3;

/// A dense scalar field sampled on a [`VoxelGrid`].
///
/// Backed by a C-contiguous leto [`Array3`] in `(i, j, k)` order (last axis
/// contiguous). Suitable for CT/MVCT densities, dose grids, and projection
/// stacks; the element type is the [`Scalar`] seam so the same volume works at
/// `f32` (GPU staging) or `f64` (reference) without a separate type.
#[derive(Debug, Clone)]
pub struct Volume<T: Scalar> {
    grid: VoxelGrid<T>,
    data: Array3<T>,
}

/// Linear interpolation `a + (b − a)·t`.
#[inline]
fn lerp<T: Scalar>(a: T, b: T, t: T) -> T {
    a + (b - a) * t
}

impl<T: Scalar> Volume<T> {
    /// Build a volume by evaluating `f` at each integer voxel index `[i, j, k]`.
    pub fn from_shape_fn<F>(grid: VoxelGrid<T>, f: F) -> Self
    where
        F: FnMut([usize; 3]) -> T,
    {
        let data = Array3::from_shape_fn(grid.dims(), f);
        Self { grid, data }
    }

    /// Build a zero-filled volume over `grid`.
    #[must_use]
    pub fn zeros(grid: VoxelGrid<T>) -> Self {
        let data = Array3::zeros(grid.dims());
        Self { grid, data }
    }

    /// Build a volume from a flat `(i, j, k)`-order (C-contiguous) vector.
    ///
    /// # Errors
    /// Returns [`HeliosError::InvalidDomainValue`] if `values.len()` does not
    /// equal the grid's voxel count.
    pub fn from_shape_vec(grid: VoxelGrid<T>, values: Vec<T>) -> Result<Self, HeliosError> {
        let expected = grid.num_voxels();
        if values.len() != expected {
            return Err(HeliosError::InvalidDomainValue {
                field: "Volume::from_shape_vec",
                value: values.len() as f64,
                reason: "value count does not match grid voxel count",
            });
        }
        let data = Array3::from_shape_vec(grid.dims(), values)
            .expect("length checked against grid dimensions");
        Ok(Self { grid, data })
    }

    /// The grid this volume is sampled on.
    #[must_use]
    pub fn grid(&self) -> &VoxelGrid<T> {
        &self.grid
    }

    /// C-contiguous flat offset of voxel `(i, j, k)` for dims `(nx, ny, nz)`.
    #[inline]
    fn flat_index(&self, i: usize, j: usize, k: usize) -> usize {
        let [_nx, ny, nz] = self.grid.dims();
        (i * ny + j) * nz + k
    }

    /// Zero-copy borrowed view of the underlying storage.
    ///
    /// Layout contract: C-contiguous `(i, j, k)` order with `k` fastest —
    /// `flat(i, j, k) = (i·ny + j)·nz + k` for dims `[nx, ny, nz]`. Hot kernels
    /// iterate this slice directly instead of per-voxel [`get`](Self::get)
    /// (which bounds-checks every axis); the layout is an API guarantee, locked
    /// by the `get`-vs-generator test.
    #[must_use]
    #[inline]
    pub fn as_slice(&self) -> &[T] {
        self.data
            .as_slice()
            .expect("invariant: volume storage is C-contiguous")
    }

    /// Value at integer voxel `(i, j, k)`, or `None` if out of bounds.
    #[must_use]
    pub fn get(&self, i: usize, j: usize, k: usize) -> Option<T> {
        let [nx, ny, nz] = self.grid.dims();
        if i >= nx || j >= ny || k >= nz {
            return None;
        }
        Some(self.as_slice()[self.flat_index(i, j, k)])
    }

    /// Add `value` to voxel `(i, j, k)` in place (dose / energy accumulation).
    ///
    /// Out-of-bounds indices are ignored (no-op), so a deposition scatter that
    /// rounds to a boundary voxel cannot panic or corrupt a neighbour.
    pub fn add_at(&mut self, i: usize, j: usize, k: usize, value: T) {
        let [nx, ny, nz] = self.grid.dims();
        if i >= nx || j >= ny || k >= nz {
            return;
        }
        let idx = (i * ny + j) * nz + k;
        let s = self
            .data
            .as_slice_mut()
            .expect("invariant: volume storage is C-contiguous");
        s[idx] += value;
    }

    /// Sum of all voxel values (e.g. total deposited energy / monitor units).
    #[must_use]
    pub fn sum(&self) -> T {
        self.as_slice()
            .iter()
            .copied()
            .fold(<T as NumericElement>::ZERO, |a, b| a + b)
    }

    /// Trilinear sample at continuous voxel-index coordinates.
    ///
    /// Returns `None` if the point is non-finite or lies outside the node range
    /// `[0, nx−1] × [0, ny−1] × [0, nz−1]`. Integer coordinates return the exact
    /// node value; a linear/affine field is reproduced exactly (the analytical
    /// oracle used in the tests).
    #[must_use]
    pub fn sample_trilinear(&self, index: Point3<T>) -> Option<T> {
        let [nx, ny, nz] = self.grid.dims();
        let coords = [index.x, index.y, index.z];
        let extents = [nx, ny, nz];

        let mut lo = [0usize; 3];
        let mut hi = [0usize; 3];
        let mut frac = [<T as NumericElement>::ZERO; 3];
        for axis in 0..3 {
            let c = coords[axis];
            if !c.is_finite() {
                return None;
            }
            let max_node = T::from_f64((extents[axis] - 1) as f64);
            if c < <T as NumericElement>::ZERO || c > max_node {
                return None;
            }
            let floor = c.floor();
            let i0 = floor.to_f64() as usize;
            lo[axis] = i0;
            hi[axis] = (i0 + 1).min(extents[axis] - 1);
            frac[axis] = c - floor;
        }

        let v = |i: usize, j: usize, k: usize| self.as_slice()[self.flat_index(i, j, k)];
        // Interpolate along x, then y, then z.
        let (tx, ty, tz) = (frac[0], frac[1], frac[2]);
        let c00 = lerp(v(lo[0], lo[1], lo[2]), v(hi[0], lo[1], lo[2]), tx);
        let c10 = lerp(v(lo[0], hi[1], lo[2]), v(hi[0], hi[1], lo[2]), tx);
        let c01 = lerp(v(lo[0], lo[1], hi[2]), v(hi[0], lo[1], hi[2]), tx);
        let c11 = lerp(v(lo[0], hi[1], hi[2]), v(hi[0], hi[1], hi[2]), tx);
        let c0 = lerp(c00, c10, ty);
        let c1 = lerp(c01, c11, ty);
        Some(lerp(c0, c1, tz))
    }

    /// Trilinear sample at a world/patient point (mm), or `None` if it maps
    /// outside the grid.
    #[must_use]
    pub fn sample_world(&self, world: Point3<T>) -> Option<T> {
        self.sample_trilinear(self.grid.world_to_index(world))
    }
}

#[cfg(test)]
mod tests {
    #![expect(
        clippy::unwrap_used,
        reason = "ratchet HELIOS-UNWRAP-1: pre-existing debt"
    )]
    use super::*;
    use eunomia::assert_relative_eq;
    use helios_math::ShippedScalar;

    fn affine_grid() -> VoxelGrid<f64> {
        VoxelGrid::axis_aligned([4, 5, 6], [2.0, 3.0, 4.0], Point3::new(10.0, 20.0, 30.0))
            .expect("valid grid")
    }

    // Affine field in index space: value = 2i + 3j + 5k + 7.
    fn affine_value(idx: [usize; 3]) -> f64 {
        2.0 * idx[0] as f64 + 3.0 * idx[1] as f64 + 5.0 * idx[2] as f64 + 7.0
    }

    #[test]
    fn get_matches_generator_locking_c_contiguous_layout() {
        // Distinct value per index locks flat_index against leto's C-order fill.
        let g = affine_grid();
        let vol = Volume::from_shape_fn(g, |idx| (idx[0] * 100 + idx[1] * 10 + idx[2]) as f64);
        for i in 0..4 {
            for j in 0..5 {
                for k in 0..6 {
                    assert_eq!(vol.get(i, j, k), Some((i * 100 + j * 10 + k) as f64));
                }
            }
        }
        assert_eq!(vol.get(4, 0, 0), None);
    }

    #[test]
    fn trilinear_reproduces_affine_field_exactly() {
        // Trilinear interpolation is exact for fields linear in each axis.
        let vol = Volume::from_shape_fn(affine_grid(), affine_value);
        let p = Point3::new(1.5, 2.25, 3.75);
        let expected = 2.0 * 1.5 + 3.0 * 2.25 + 5.0 * 3.75 + 7.0;
        assert_relative_eq!(vol.sample_trilinear(p).unwrap(), expected, epsilon = 1e-12);
    }

    #[test]
    fn trilinear_at_integer_node_returns_node_value() {
        let vol = Volume::from_shape_fn(affine_grid(), affine_value);
        assert_relative_eq!(
            vol.sample_trilinear(Point3::new(1.0, 2.0, 3.0)).unwrap(),
            affine_value([1, 2, 3]),
            epsilon = 1e-12
        );
    }

    #[test]
    fn sample_out_of_bounds_and_nonfinite_is_none() {
        let vol = Volume::from_shape_fn(affine_grid(), affine_value);
        assert_eq!(vol.sample_trilinear(Point3::new(-0.1, 0.0, 0.0)), None);
        assert_eq!(vol.sample_trilinear(Point3::new(3.5, 0.0, 0.0)), None); // nx-1 = 3
        assert_eq!(vol.sample_trilinear(Point3::new(f64::NAN, 0.0, 0.0)), None);
    }

    #[test]
    fn sample_world_hits_expected_voxel() {
        let g = affine_grid();
        let vol = Volume::from_shape_fn(g, affine_value);
        // World center of voxel (1,2,3) must sample that node's value.
        let world = vol.grid().voxel_center(1, 2, 3);
        assert_relative_eq!(
            vol.sample_world(world).unwrap(),
            affine_value([1, 2, 3]),
            epsilon = 1e-9
        );
    }

    /// Trilinear interpolation of an affine field is exact in exact arithmetic,
    /// so the error is only the eight weighted products and their summation.
    /// Sixteen ulps of `T` bounds that.
    const TRILINEAR_ULPS: f64 = 16.0;

    /// Asserts trilinear sampling reproduces an affine field in one width.
    fn trilinear_sampling_reproduces_an_affine_field<T>()
    where
        T: ShippedScalar + helios_math::GeometryScalar,
    {
        // GeometryScalar and FloatElement both define `from_f64`; name the
        // one that performs the literal conversion so the call is unambiguous.
        let cast = <T as helios_math::FloatElement>::from_f64;
        let zero = cast(0.0);
        let unit = cast(1.0);
        let grid =
            VoxelGrid::<T>::axis_aligned([2, 2, 2], [unit; 3], Point3::new(zero, zero, zero))
                .expect("valid axis-aligned grid");

        // f(i,j,k) = i + j + k is affine, so trilinear interpolation is exact.
        let volume = Volume::from_shape_fn(grid, |idx| {
            cast(idx[0] as f64) + cast(idx[1] as f64) + cast(idx[2] as f64)
        });

        let half = cast(0.5);
        assert_relative_eq!(
            volume
                .sample_trilinear(Point3::new(half, half, half))
                .expect("midpoint lies inside the grid"),
            cast(1.5),
            max_relative = T::EPSILON * cast(TRILINEAR_ULPS)
        );
    }

    #[test]
    fn trilinear_sampling_reproduces_an_affine_field_in_single_precision() {
        trilinear_sampling_reproduces_an_affine_field::<f32>();
    }

    #[test]
    fn trilinear_sampling_reproduces_an_affine_field_in_double_precision() {
        trilinear_sampling_reproduces_an_affine_field::<f64>();
    }

    #[test]
    fn add_at_accumulates_and_ignores_out_of_bounds() {
        let mut vol = Volume::zeros(affine_grid()); // dims [4,5,6]
        vol.add_at(1, 2, 3, 5.0);
        vol.add_at(1, 2, 3, 2.5);
        assert_relative_eq!(vol.get(1, 2, 3).unwrap(), 7.5, epsilon = 1e-15);
        // Out-of-bounds scatter is a silent no-op: total is unchanged.
        vol.add_at(4, 0, 0, 100.0);
        assert_relative_eq!(vol.sum(), 7.5, epsilon = 1e-15);
    }

    #[test]
    fn from_shape_vec_rejects_wrong_length() {
        let g = affine_grid();
        assert!(Volume::from_shape_vec(g, vec![0.0; 10]).is_err());
    }
}
