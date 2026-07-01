//! DICOM CT/MVCT load path (real-input boundary), built on `ritk-dicom`.
//!
//! Reads a single-slice DICOM image into a Helios [`Volume`] of Hounsfield units:
//! `ritk-dicom` parses the file and decodes the pixel frame (applying the
//! `RescaleSlope`/`RescaleIntercept` calibration), and this module maps the frame
//! plus the geometry attributes (`Rows`, `Columns`, `PixelSpacing`,
//! `ImagePositionPatient`) into a typed [`Volume`] on an axis-aligned
//! [`VoxelGrid`]. This is the trust boundary: external file bytes become validated
//! typed domain values here.
//!
//! Multi-slice series stacking (sort by `ImagePositionPatient`, derive the z
//! spacing) is a follow-up; this loads one slice (`nz = 1`).
//!
//! Feature-gated behind `dicom` so the heavy `dicom-rs` parser stays out of the
//! core build. The feature gates a complete implementation, not a stub.

use crate::grid::VoxelGrid;
use crate::volume::Volume;
use dicom::core::Tag;
use helios_core::HeliosError;
use helios_math::{Point3, Scalar};
use ritk_dicom::{
    decode_frame_with, parse_file_with, DecodeFrameRequest, DicomRsBackend, PixelLayout,
    PixelSignedness, TransferSyntaxKind,
};

// DICOM tags (group, element).
const ROWS: Tag = Tag(0x0028, 0x0010);
const COLUMNS: Tag = Tag(0x0028, 0x0011);
const SAMPLES_PER_PIXEL: Tag = Tag(0x0028, 0x0002);
const BITS_ALLOCATED: Tag = Tag(0x0028, 0x0100);
const PIXEL_REPRESENTATION: Tag = Tag(0x0028, 0x0103);
const RESCALE_INTERCEPT: Tag = Tag(0x0028, 0x1052);
const RESCALE_SLOPE: Tag = Tag(0x0028, 0x1053);
const PIXEL_SPACING: Tag = Tag(0x0028, 0x0030);
const SLICE_THICKNESS: Tag = Tag(0x0018, 0x0050);
const IMAGE_POSITION_PATIENT: Tag = Tag(0x0020, 0x0032);

type Object = <DicomRsBackend as ritk_dicom::DicomParseBackend>::Object;

fn dicom_err(step: &str, e: impl core::fmt::Display) -> HeliosError {
    HeliosError::Dicom {
        reason: format!("{step}: {e}"),
    }
}

/// Required unsigned-short attribute.
fn req_usize(obj: &Object, tag: Tag, name: &'static str) -> Result<usize, HeliosError> {
    let v: u16 = obj
        .element(tag)
        .map_err(|e| dicom_err(name, e))?
        .value()
        .to_int::<u16>()
        .map_err(|e| dicom_err(name, e))?;
    Ok(v as usize)
}

/// Optional unsigned-short attribute with a default.
fn opt_u16(obj: &Object, tag: Tag, default: u16) -> u16 {
    obj.element(tag)
        .ok()
        .and_then(|e| e.value().to_int::<u16>().ok())
        .unwrap_or(default)
}

/// Optional decimal-string scalar with a default.
fn opt_f64(obj: &Object, tag: Tag, default: f64) -> f64 {
    obj.element(tag)
        .ok()
        .and_then(|e| e.value().to_float64().ok())
        .unwrap_or(default)
}

/// Optional multi-valued decimal string.
fn multi_f64(obj: &Object, tag: Tag) -> Option<Vec<f64>> {
    obj.element(tag)
        .ok()
        .and_then(|e| e.value().to_multi_float64().ok())
}

/// Load a single-slice DICOM CT/MVCT image into a [`Volume`] of Hounsfield units.
///
/// The pixel frame is decoded with the file's `RescaleSlope`/`RescaleIntercept`,
/// so the volume holds HU directly. Grid geometry: `dims = [Columns, Rows, 1]`
/// (voxel index `i = column`/x, `j = row`/y, `k = 0`); spacing
/// `[PixelSpacing_col, PixelSpacing_row, SliceThickness]` (mm); origin from
/// `ImagePositionPatient` (defaulting to the origin when absent).
///
/// # Errors
/// [`HeliosError::Dicom`] if the file cannot be parsed/decoded or a required
/// geometry attribute is missing or malformed; [`HeliosError::InvalidDomainValue`]
/// if the resulting grid dimensions/spacing are invalid.
pub fn load_ct_slice<T: Scalar>(
    path: impl AsRef<std::path::Path>,
) -> Result<Volume<T>, HeliosError> {
    let obj =
        parse_file_with::<DicomRsBackend, _>(path.as_ref()).map_err(|e| dicom_err("parse", e))?;

    let rows = req_usize(&obj, ROWS, "Rows")?;
    let cols = req_usize(&obj, COLUMNS, "Columns")?;
    let samples_per_pixel = opt_u16(&obj, SAMPLES_PER_PIXEL, 1) as usize;
    let bits_allocated = opt_u16(&obj, BITS_ALLOCATED, 16);
    let pixel_representation = if opt_u16(&obj, PIXEL_REPRESENTATION, 0) == 1 {
        PixelSignedness::Signed
    } else {
        PixelSignedness::Unsigned
    };
    let rescale_slope = opt_f64(&obj, RESCALE_SLOPE, 1.0) as f32;
    let rescale_intercept = opt_f64(&obj, RESCALE_INTERCEPT, 0.0) as f32;

    // PixelSpacing is [row_spacing, col_spacing] (mm); default to 1 mm isotropic.
    let spacing = multi_f64(&obj, PIXEL_SPACING).unwrap_or_default();
    let row_spacing = spacing.first().copied().unwrap_or(1.0);
    let col_spacing = spacing.get(1).copied().unwrap_or(row_spacing);
    let thickness = opt_f64(&obj, SLICE_THICKNESS, 1.0);

    let ipp = multi_f64(&obj, IMAGE_POSITION_PATIENT).unwrap_or_default();
    let origin = Point3::new(
        T::from_f64(ipp.first().copied().unwrap_or(0.0)),
        T::from_f64(ipp.get(1).copied().unwrap_or(0.0)),
        T::from_f64(ipp.get(2).copied().unwrap_or(0.0)),
    );

    let transfer_syntax =
        TransferSyntaxKind::from_uid(obj.meta().transfer_syntax.trim_end_matches('\0'));

    let layout = PixelLayout {
        rows,
        cols,
        samples_per_pixel,
        bits_allocated,
        pixel_representation,
        rescale_slope,
        rescale_intercept,
    };
    let frame = decode_frame_with::<DicomRsBackend>(
        &obj,
        DecodeFrameRequest {
            frame_index: 0,
            transfer_syntax,
            layout,
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

    let grid = VoxelGrid::axis_aligned(
        [cols, rows, 1],
        [
            T::from_f64(col_spacing),
            T::from_f64(row_spacing),
            T::from_f64(thickness),
        ],
        origin,
    )?;

    // Re-order the DICOM row-major frame `[row·cols + col]` into the volume's
    // C-contiguous `(i, j, k)` layout with i = column, j = row, k = 0:
    // flat(i, j, 0) = i·rows + j.
    let mut data = vec![T::from_f64(0.0); rows * cols];
    for row in 0..rows {
        for col in 0..cols {
            data[col * rows + row] = T::from_f64(f64::from(frame.pixels[row * cols + col]));
        }
    }
    Volume::from_shape_vec(grid, data)
}

#[cfg(test)]
mod tests {
    use super::*;
    use dicom::core::smallvec::SmallVec;
    use dicom::core::{DataElement, PrimitiveValue, VR};
    use dicom::object::{FileMetaTableBuilder, InMemDicomObject};

    // Write a synthetic 2×2 unsigned-16 CT slice with a known HU pattern and
    // geometry, then load it back — a deterministic round-trip (no external
    // fixture). Raw pixels [10,20,30,40]; slope 2, intercept −10 → HU
    // [10,30,50,70]. PixelSpacing 0.8 (row) / 1.25 (col) mm; origin (5,7,9).
    fn write_slice(path: &std::path::Path) {
        let mut obj = InMemDicomObject::new_empty();
        let put = |obj: &mut InMemDicomObject, tag, vr, val: PrimitiveValue| {
            obj.put(DataElement::new(tag, vr, val));
        };
        put(
            &mut obj,
            Tag(0x0008, 0x0016),
            VR::UI,
            "1.2.840.10008.5.1.4.1.1.2".into(),
        );
        put(&mut obj, Tag(0x0008, 0x0018), VR::UI, "2.25.4242".into());
        put(&mut obj, ROWS, VR::US, PrimitiveValue::from(2_u16));
        put(&mut obj, COLUMNS, VR::US, PrimitiveValue::from(2_u16));
        put(
            &mut obj,
            SAMPLES_PER_PIXEL,
            VR::US,
            PrimitiveValue::from(1_u16),
        );
        put(
            &mut obj,
            BITS_ALLOCATED,
            VR::US,
            PrimitiveValue::from(16_u16),
        );
        put(
            &mut obj,
            PIXEL_REPRESENTATION,
            VR::US,
            PrimitiveValue::from(0_u16),
        );
        put(&mut obj, RESCALE_SLOPE, VR::DS, "2".into());
        put(&mut obj, RESCALE_INTERCEPT, VR::DS, "-10".into());
        put(&mut obj, PIXEL_SPACING, VR::DS, "0.8\\1.25".into());
        put(&mut obj, SLICE_THICKNESS, VR::DS, "3".into());
        put(&mut obj, IMAGE_POSITION_PATIENT, VR::DS, "5\\7\\9".into());
        put(
            &mut obj,
            Tag(0x7FE0, 0x0010),
            VR::OW,
            PrimitiveValue::U16(SmallVec::from_vec(vec![10_u16, 20, 30, 40])),
        );
        obj.with_meta(
            FileMetaTableBuilder::new()
                .media_storage_sop_class_uid("1.2.840.10008.5.1.4.1.1.2")
                .media_storage_sop_instance_uid("2.25.4242")
                .transfer_syntax("1.2.840.10008.1.2.1"),
        )
        .expect("valid meta")
        .write_to_file(path)
        .expect("write dicom");
    }

    #[test]
    fn round_trips_a_synthetic_ct_slice_to_hu_volume() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("slice.dcm");
        write_slice(&path);

        let vol: Volume<f64> = load_ct_slice(&path).expect("load");
        let grid = vol.grid();
        assert_eq!(grid.dims(), [2, 2, 1]);
        // spacing = [col, row, thickness] = [1.25, 0.8, 3.0].
        let sp = grid.spacing();
        assert!(
            (sp[0] - 1.25).abs() < 1e-12
                && (sp[1] - 0.8).abs() < 1e-12
                && (sp[2] - 3.0).abs() < 1e-12
        );

        // HU = raw·2 − 10; DICOM row-major [10,20,30,40] → HU [10,30,50,70].
        // Volume index (i=col, j=row): (row0,col0)=10, (row0,col1)=30,
        // (row1,col0)=50, (row1,col1)=70.
        assert_eq!(vol.get(0, 0, 0), Some(10.0)); // col0,row0
        assert_eq!(vol.get(1, 0, 0), Some(30.0)); // col1,row0
        assert_eq!(vol.get(0, 1, 0), Some(50.0)); // col0,row1
        assert_eq!(vol.get(1, 1, 0), Some(70.0)); // col1,row1
    }

    #[test]
    fn missing_file_is_a_dicom_error_not_a_panic() {
        let err = load_ct_slice::<f64>("does_not_exist.dcm").unwrap_err();
        assert!(matches!(err, HeliosError::Dicom { .. }));
    }
}
