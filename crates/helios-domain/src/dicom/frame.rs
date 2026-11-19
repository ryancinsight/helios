//! Read one CT frame and the acquisition geometry required by Helios.

use helios_quantities::HeliosError;
use ritk_dicom::dataset::DicomFile;
use ritk_dicom::{
    DecodeFrameRequest, DicomAttributeRead, DicomTag, NativeBackend, PixelLayout, PixelSignedness,
    TransferSyntaxKind, decode_frame_with, tags,
};

use super::geometry::{cross3, dot3, normalize3};

fn dicom_err(step: &str, error: impl core::fmt::Display) -> HeliosError {
    HeliosError::Dicom {
        reason: format!("{step}: {error}"),
    }
}

/// Required unsigned-short attribute.
fn req_usize(obj: &DicomFile, tag: DicomTag, name: &'static str) -> Result<usize, HeliosError> {
    obj.required_unsigned(tag, name)
        .map(usize::from)
        .map_err(|e| dicom_err(name, e))
}

/// Optional unsigned-short attribute with a default.
fn opt_u16(
    obj: &DicomFile,
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
    obj: &DicomFile,
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
    obj: &DicomFile,
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
pub(super) struct SliceRaw {
    pub(super) rows: usize,
    pub(super) cols: usize,
    pub(super) col_spacing: f64,
    pub(super) row_spacing: f64,
    pub(super) thickness: f64,
    pub(super) origin: [f64; 3],
    /// Direction cosine for the column index axis (`i`), from DICOM row direction.
    pub(super) row_dir: [f64; 3],
    /// Direction cosine for the row index axis (`j`), from DICOM column direction.
    pub(super) col_dir: [f64; 3],
    /// Slice-stack normal (`row_dir × col_dir`), right-handed.
    pub(super) normal_dir: [f64; 3],
    /// Origin projected onto `normal_dir` (slice position along stack axis, mm).
    pub(super) stack_position: f64,
    pub(super) hu: Vec<f32>,
}

/// Parse and decode one DICOM slice into [`SliceRaw`] (HU + geometry).
pub(super) fn read_slice(path: &std::path::Path) -> Result<SliceRaw, HeliosError> {
    let obj = DicomFile::open(path).map_err(|e| dicom_err("parse", e))?;

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
        tags::IMAGE_ORIENTATION_PATIENT,
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
    let frame = decode_frame_with::<NativeBackend>(
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
