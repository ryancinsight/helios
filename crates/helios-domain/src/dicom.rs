//! DICOM CT/MVCT load path (real-input boundary), built on `ritk-dicom`.
//!
//! Reads a single-slice DICOM image into a Helios [`Volume`] of Hounsfield units:
//! `ritk-dicom` parses the file and decodes the pixel frame (applying the
//! `RescaleSlope`/`RescaleIntercept` calibration), and this module maps the frame
//! plus the geometry attributes (`Rows`, `Columns`, `PixelSpacing`,
//! `ImagePositionPatient`, `ImageOrientationPatient`) into a typed [`Volume`] on
//! an oriented [`VoxelGrid`]. The pixel layout carries both `BitsAllocated` and
//! required `BitsStored`, so container padding cannot change decoded values. This
//! is the trust boundary: external file bytes become validated typed domain values
//! here.
//!
//! Both single-slice and multi-slice series paths are provided:
//! [`load_ct_slice`] loads one slice (`nz = 1`), and [`load_ct_series`] validates
//! and stacks a full slice set along the oriented stack axis.
//!
//! Feature-gated behind `dicom` so the RITK DICOM provider stays out of the core
//! build. The feature gates a complete implementation, not a stub.

use crate::grid::VoxelGrid;
use crate::volume::Volume;
use helios_core::HeliosError;
use helios_math::{Point3, Scalar, UnitQuaternion, Vector3};
use ritk_dicom::{
    decode_frame_with, parse_file_with, tags, DecodeFrameRequest, DicomAttributeRead,
    DicomRsBackend, DicomTag, PixelLayout, PixelSignedness, TransferSyntaxKind,
};

type Object = <DicomRsBackend as ritk_dicom::DicomParseBackend>::Object;
const IMAGE_ORIENTATION_PATIENT: DicomTag = DicomTag::new(0x0020, 0x0037);
const ORIENTATION_TOL: f64 = 1.0e-6;

fn dicom_err(step: &str, e: impl core::fmt::Display) -> HeliosError {
    HeliosError::Dicom {
        reason: format!("{step}: {e}"),
    }
}

/// Required unsigned-short attribute.
fn req_usize(obj: &Object, tag: DicomTag, name: &'static str) -> Result<usize, HeliosError> {
    obj.required_unsigned(tag, name)
        .map(usize::from)
        .map_err(|e| dicom_err(name, e))
}

/// Optional unsigned-short attribute with a default.
fn opt_u16(
    obj: &Object,
    tag: DicomTag,
    name: &'static str,
    default: u16,
) -> Result<u16, HeliosError> {
    obj.optional_unsigned(tag, name)
        .map(|value| value.unwrap_or(default))
        .map_err(|e| dicom_err(name, e))
}

/// Optional decimal-string scalar with a default.
fn opt_f64(
    obj: &Object,
    tag: DicomTag,
    name: &'static str,
    default: f64,
) -> Result<f64, HeliosError> {
    obj.optional_decimal(tag, name)
        .map(|value| value.unwrap_or(default))
        .map_err(|e| dicom_err(name, e))
}

/// Required multi-valued decimal string with an exact value count.
fn required_f64_array<const N: usize>(
    obj: &Object,
    tag: DicomTag,
    name: &'static str,
) -> Result<[f64; N], HeliosError> {
    let values = obj
        .optional_decimal_vec(tag, name)
        .map_err(|e| dicom_err(name, e))?
        .ok_or_else(|| HeliosError::Dicom {
            reason: format!("{name}: missing required attribute"),
        })?;
    values
        .try_into()
        .map_err(|values: Vec<f64>| HeliosError::Dicom {
            reason: format!("{name} expected {N} values, got {}", values.len()),
        })
}

/// One parsed+decoded DICOM slice in native (f64/mm/HU) form, before it is mapped
/// into a typed [`Volume`]. In-plane geometry is kept as `f64` for consistency
/// checks across a series; `hu` is row-major `[row·cols + col]`.
struct SliceRaw {
    rows: usize,
    cols: usize,
    col_spacing: f64,
    row_spacing: f64,
    thickness: f64,
    origin: [f64; 3],
    /// Direction cosine for the column index axis (`i`), from DICOM row direction.
    row_dir: [f64; 3],
    /// Direction cosine for the row index axis (`j`), from DICOM column direction.
    col_dir: [f64; 3],
    /// Slice-stack normal (`row_dir × col_dir`), right-handed.
    normal_dir: [f64; 3],
    /// Origin projected onto `normal_dir` (slice position along stack axis, mm).
    stack_position: f64,
    hu: Vec<f32>,
}

/// Parse and decode one DICOM slice into [`SliceRaw`] (HU + geometry).
fn read_slice(path: &std::path::Path) -> Result<SliceRaw, HeliosError> {
    let obj = parse_file_with::<DicomRsBackend, _>(path).map_err(|e| dicom_err("parse", e))?;

    let rows = req_usize(&obj, tags::ROWS, "Rows")?;
    let cols = req_usize(&obj, tags::COLUMNS, "Columns")?;
    let samples_per_pixel = usize::from(opt_u16(
        &obj,
        tags::SAMPLES_PER_PIXEL,
        "SamplesPerPixel",
        1,
    )?);
    let bits_allocated = opt_u16(&obj, tags::BITS_ALLOCATED, "BitsAllocated", 16)?;
    let bits_stored = obj
        .required_unsigned(tags::BITS_STORED, "BitsStored")
        .map_err(|e| dicom_err("BitsStored", e))?;
    let pixel_representation =
        if opt_u16(&obj, tags::PIXEL_REPRESENTATION, "PixelRepresentation", 0)? == 1 {
            PixelSignedness::Signed
        } else {
            PixelSignedness::Unsigned
        };
    let rescale_slope = opt_f64(&obj, tags::RESCALE_SLOPE, "RescaleSlope", 1.0)? as f32;
    let rescale_intercept = opt_f64(&obj, tags::RESCALE_INTERCEPT, "RescaleIntercept", 0.0)? as f32;

    // PixelSpacing is [row_spacing, col_spacing] (mm), and is required for a
    // patient-coordinate grid. Silent unit-spacing recovery would change the
    // physical meaning of every voxel without an error at the trust boundary.
    let [row_spacing, col_spacing] =
        required_f64_array::<2>(&obj, tags::PIXEL_SPACING, "PixelSpacing")?;
    let thickness = opt_f64(&obj, tags::SLICE_THICKNESS, "SliceThickness", 1.0)?;

    let origin =
        required_f64_array::<3>(&obj, tags::IMAGE_POSITION_PATIENT, "ImagePositionPatient")?;
    let orientation = required_f64_array::<6>(
        &obj,
        IMAGE_ORIENTATION_PATIENT,
        "ImageOrientationPatient (0020,0037)",
    )?;
    let row_dir = normalize3(
        [orientation[0], orientation[1], orientation[2]],
        "ImageOrientationPatient row direction",
    )?;
    let col_dir = normalize3(
        [orientation[3], orientation[4], orientation[5]],
        "ImageOrientationPatient column direction",
    )?;
    let normal_dir = normalize3(
        cross3(row_dir, col_dir),
        "ImageOrientationPatient normal direction",
    )?;
    let stack_position = dot3(origin, normal_dir);

    let transfer_syntax = TransferSyntaxKind::from_uid(obj.transfer_syntax_uid());
    let frame = decode_frame_with::<DicomRsBackend>(
        &obj,
        DecodeFrameRequest {
            frame_index: 0,
            transfer_syntax,
            layout: PixelLayout {
                rows,
                cols,
                samples_per_pixel,
                bits_allocated,
                bits_stored,
                pixel_representation,
                rescale_slope,
                rescale_intercept,
            },
        },
    )
    .map_err(|e| dicom_err("decode", e))?;

    if frame.pixels.len() != rows * cols {
        return Err(HeliosError::Dicom {
            reason: format!(
                "decoded pixel count {} != Rows·Columns {}",
                frame.pixels.len(),
                rows * cols
            ),
        });
    }

    Ok(SliceRaw {
        rows,
        cols,
        col_spacing,
        row_spacing,
        thickness,
        origin,
        row_dir,
        col_dir,
        normal_dir,
        stack_position,
        hu: frame.pixels,
    })
}

fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn normalize3(v: [f64; 3], component: &'static str) -> Result<[f64; 3], HeliosError> {
    let norm2 = dot3(v, v);
    if !norm2.is_finite() || norm2 <= ORIENTATION_TOL * ORIENTATION_TOL {
        return Err(HeliosError::Dicom {
            reason: format!("{component} is zero-length or non-finite"),
        });
    }
    let inv = norm2.sqrt().recip();
    Ok([v[0] * inv, v[1] * inv, v[2] * inv])
}

fn rotation_from_axes<T: Scalar>(
    row_dir: [f64; 3],
    col_dir: [f64; 3],
    normal_dir: [f64; 3],
) -> Result<UnitQuaternion<T>, HeliosError> {
    UnitQuaternion::try_from_rotation_columns(
        Vector3::new(
            T::from_f64(row_dir[0]),
            T::from_f64(row_dir[1]),
            T::from_f64(row_dir[2]),
        ),
        Vector3::new(
            T::from_f64(col_dir[0]),
            T::from_f64(col_dir[1]),
            T::from_f64(col_dir[2]),
        ),
        Vector3::new(
            T::from_f64(normal_dir[0]),
            T::from_f64(normal_dir[1]),
            T::from_f64(normal_dir[2]),
        ),
        T::from_f64(ORIENTATION_TOL),
    )
    .map_err(|e| HeliosError::Dicom {
        reason: format!("invalid ImageOrientationPatient basis: {e}"),
    })
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
pub fn load_ct_slice<T: Scalar>(
    path: impl AsRef<std::path::Path>,
) -> Result<Volume<T>, HeliosError> {
    let slice = read_slice(path.as_ref())?;
    let rotation = rotation_from_axes::<T>(slice.row_dir, slice.col_dir, slice.normal_dir)?;
    let grid = VoxelGrid::oriented(
        [slice.cols, slice.rows, 1],
        [
            T::from_f64(slice.col_spacing),
            T::from_f64(slice.row_spacing),
            T::from_f64(slice.thickness),
        ],
        Point3::new(
            T::from_f64(slice.origin[0]),
            T::from_f64(slice.origin[1]),
            T::from_f64(slice.origin[2]),
        ),
        rotation,
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
pub fn load_ct_series<T: Scalar, P: AsRef<std::path::Path>>(
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
        let same_row_dir =
            (0..3).all(|axis| (s.row_dir[axis] - row_dir[axis]).abs() <= ORIENTATION_TOL);
        let same_col_dir =
            (0..3).all(|axis| (s.col_dir[axis] - col_dir[axis]).abs() <= ORIENTATION_TOL);
        let same_normal_dir =
            (0..3).all(|axis| (s.normal_dir[axis] - normal_dir[axis]).abs() <= ORIENTATION_TOL);
        let in_plane_row = dot3(s.origin, row_dir);
        let in_plane_col = dot3(s.origin, col_dir);
        let consistent = s.rows == rows
            && s.cols == cols
            && (s.col_spacing - col_sp).abs() <= GEOMETRY_TOL_MM
            && (s.row_spacing - row_sp).abs() <= GEOMETRY_TOL_MM
            && (in_plane_row - in_plane_origin_row).abs() <= GEOMETRY_TOL_MM
            && (in_plane_col - in_plane_origin_col).abs() <= GEOMETRY_TOL_MM
            && same_row_dir
            && same_col_dir
            && same_normal_dir;
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

    let rotation = rotation_from_axes::<T>(row_dir, col_dir, normal_dir)?;
    let grid = VoxelGrid::oriented(
        [cols, rows, nz],
        [
            T::from_f64(col_sp),
            T::from_f64(row_sp),
            T::from_f64(z_spacing),
        ],
        Point3::new(
            T::from_f64(slices[0].origin[0]),
            T::from_f64(slices[0].origin[1]),
            T::from_f64(slices[0].origin[2]),
        ),
        rotation,
    )?;

    let mut data = vec![T::from_f64(0.0); rows * cols * nz];
    for (k, slice) in slices.iter().enumerate() {
        scatter_slice(&mut data, slice, k, nz);
    }
    Volume::from_shape_vec(grid, data)
}

#[cfg(test)]
mod tests;
