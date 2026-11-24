//! Map validated CT frames into oriented Helios volumes.

use crate::grid::VoxelGrid;
use crate::volume::Volume;
use aequitas::systems::si::{quantities::Length, units::Millimeter};
use eunomia::UnitScalar;
use helios_math::{Point3, Scalar};
use helios_quantities::HeliosError;

use super::frame::{read_slice, SliceRaw};
use super::geometry::{allclose3, dot3, rotation_from_axes, ORIENTATION_TOL};

fn oriented_grid_from_dicom<T: Scalar + UnitScalar>(
    dims: [usize; 3],
    spacing_mm: [f64; 3],
    origin: [f64; 3],
    row_dir: [f64; 3],
    col_dir: [f64; 3],
    normal_dir: [f64; 3],
) -> Result<VoxelGrid<T>, HeliosError> {
    // DICOM Pixel Spacing / Slice Thickness are millimetres by standard
    // (PS3.3 C.8); the typed boundary records that contract here.
    let spacing = spacing_mm.map(|s| Length::from_unit::<Millimeter>(T::from_f64(s)));
    let rotation = rotation_from_axes::<T>(row_dir, col_dir, normal_dir)?;
    VoxelGrid::oriented(
        dims,
        spacing,
        Point3::new(
            T::from_f64(origin[0]),
            T::from_f64(origin[1]),
            T::from_f64(origin[2]),
        ),
        rotation,
    )
}

/// Scatter one slice's row-major HU frame into a stacked C-contiguous
/// `(i = col, j = row, k)` buffer of shape `[cols, rows, nz]`:
/// `flat(i, j, k) = (i·rows + j)·nz + k`.
fn scatter_slice<T: Scalar>(dst: &mut [T], slice: &SliceRaw, k: usize, nz: usize) {
    let (rows, cols) = (slice.rows, slice.cols);
    for row in 0..rows {
        for col in 0..cols {
            dst[(col * rows + row) * nz + k] = T::from_f64(f64::from(slice.hu[row * cols + col]));
        }
    }
}

/// Load a single-slice DICOM CT/MVCT image into a [`Volume`] of Hounsfield units.
///
/// The pixel frame is decoded with the file's `RescaleSlope`/`RescaleIntercept`,
/// so the volume holds HU directly. Grid geometry: `dims = [Columns, Rows, 1]`
/// (voxel index `i = column`/x, `j = row`/y, `k = 0`); spacing
/// `[PixelSpacing_col, PixelSpacing_row, SliceThickness]` (mm); origin from
/// the required `ImagePositionPatient`; and orientation from the required
/// `ImageOrientationPatient` direction cosines.
///
/// # Errors
/// [`HeliosError::Dicom`] if the file cannot be parsed/decoded or a required
/// geometry attribute is missing or malformed; [`HeliosError::InvalidDomainValue`]
/// if the resulting grid dimensions/spacing are invalid.
pub fn load_ct_slice<T: Scalar + UnitScalar>(
    path: impl AsRef<std::path::Path>,
) -> Result<Volume<T>, HeliosError> {
    let slice = read_slice(path.as_ref())?;
    let grid = oriented_grid_from_dicom(
        [slice.cols, slice.rows, 1],
        [slice.col_spacing, slice.row_spacing, slice.thickness],
        slice.origin,
        slice.row_dir,
        slice.col_dir,
        slice.normal_dir,
    )?;
    let mut data = vec![T::from_f64(0.0); slice.rows * slice.cols];
    scatter_slice(&mut data, &slice, 0, 1);
    Volume::from_shape_vec(grid, data)
}

/// Consistency tolerance for in-plane geometry and slice spacing (mm).
///
/// DICOM stores positions/spacings as decimal strings; series slices share an
/// identical in-plane grid and a constant stack-axis step. 1 µm (`1e-3` mm) is
/// tight enough to catch a missing slice (gap = 2× spacing) or a mismatched grid
/// while tolerating decimal-string round-off.
const GEOMETRY_TOL_MM: f64 = 1.0e-3;

/// Load a multi-slice DICOM CT/MVCT **series** into a 3-D [`Volume`] of Hounsfield
/// units.
///
/// Every slice is parsed and decoded (HU); the slices are validated to share an
/// identical in-plane grid (`Rows`/`Columns`/`PixelSpacing`/orientation/in-plane
/// origin), sorted by stack position (`ImagePositionPatient` projected onto the
/// DICOM slice normal), and stacked along `k`. The z spacing is derived from the
/// (uniform) consecutive stack positions. Result geometry:
/// `dims = [Columns, Rows, nslices]`, spacing `[col, row, Δz]` (mm), origin at the
/// lowest-z slice.
///
/// # Errors
/// [`HeliosError::Dicom`] if `paths` is empty, any slice fails to load, the slices
/// disagree in in-plane geometry, or the z spacing is non-uniform (beyond
/// `GEOMETRY_TOL_MM`); [`HeliosError::InvalidDomainValue`] if the derived grid is
/// invalid (e.g. duplicate slice positions → zero spacing).
pub fn load_ct_series<T: Scalar + UnitScalar, P: AsRef<std::path::Path>>(
    paths: &[P],
) -> Result<Volume<T>, HeliosError> {
    if paths.is_empty() {
        return Err(HeliosError::Dicom {
            reason: "empty DICOM series (no slice paths)".to_owned(),
        });
    }
    let mut slices: Vec<SliceRaw> = paths
        .iter()
        .map(|p| read_slice(p.as_ref()))
        .collect::<Result<_, _>>()?;

    // In-plane geometry must be identical across the series.
    let (rows, cols) = (slices[0].rows, slices[0].cols);
    let (col_sp, row_sp) = (slices[0].col_spacing, slices[0].row_spacing);
    let row_dir = slices[0].row_dir;
    let col_dir = slices[0].col_dir;
    let normal_dir = slices[0].normal_dir;
    let in_plane_origin_row = dot3(slices[0].origin, row_dir);
    let in_plane_origin_col = dot3(slices[0].origin, col_dir);
    for s in &slices[1..] {
        let in_plane_row = dot3(s.origin, row_dir);
        let in_plane_col = dot3(s.origin, col_dir);
        let consistent = s.rows == rows
            && s.cols == cols
            && (s.col_spacing - col_sp).abs() <= GEOMETRY_TOL_MM
            && (s.row_spacing - row_sp).abs() <= GEOMETRY_TOL_MM
            && (in_plane_row - in_plane_origin_row).abs() <= GEOMETRY_TOL_MM
            && (in_plane_col - in_plane_origin_col).abs() <= GEOMETRY_TOL_MM
            && allclose3(s.row_dir, row_dir, ORIENTATION_TOL)
            && allclose3(s.col_dir, col_dir, ORIENTATION_TOL)
            && allclose3(s.normal_dir, normal_dir, ORIENTATION_TOL);
        if !consistent {
            return Err(HeliosError::Dicom {
                reason: "series slices have inconsistent in-plane geometry".to_owned(),
            });
        }
    }

    // Order along the stack axis and derive a uniform z spacing.
    slices.sort_by(|a, b| a.stack_position.total_cmp(&b.stack_position));
    let nz = slices.len();
    let z_spacing = if nz > 1 {
        slices[1].stack_position - slices[0].stack_position
    } else {
        slices[0].thickness
    };
    for w in slices.windows(2) {
        if ((w[1].stack_position - w[0].stack_position) - z_spacing).abs() > GEOMETRY_TOL_MM {
            return Err(HeliosError::Dicom {
                reason: "non-uniform slice spacing (missing or duplicate slice?)".to_owned(),
            });
        }
    }

    let grid = oriented_grid_from_dicom(
        [cols, rows, nz],
        [col_sp, row_sp, z_spacing],
        slices[0].origin,
        row_dir,
        col_dir,
        normal_dir,
    )?;

    let mut data = vec![T::from_f64(0.0); rows * cols * nz];
    for (k, slice) in slices.iter().enumerate() {
        scatter_slice(&mut data, slice, k, nz);
    }
    Volume::from_shape_vec(grid, data)
}
