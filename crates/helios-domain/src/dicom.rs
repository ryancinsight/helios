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

/// One parsed+decoded DICOM slice in native (f64/mm/HU) form, before it is mapped
/// into a typed [`Volume`]. In-plane geometry is kept as `f64` for consistency
/// checks across a series; `hu` is row-major `[row·cols + col]`.
struct SliceRaw {
    rows: usize,
    cols: usize,
    col_spacing: f64,
    row_spacing: f64,
    thickness: f64,
    origin_x: f64,
    origin_y: f64,
    /// `ImagePositionPatient` z (slice position along the stack axis, mm).
    z: f64,
    hu: Vec<f32>,
}

/// Parse and decode one DICOM slice into [`SliceRaw`] (HU + geometry).
fn read_slice(path: &std::path::Path) -> Result<SliceRaw, HeliosError> {
    let obj = parse_file_with::<DicomRsBackend, _>(path).map_err(|e| dicom_err("parse", e))?;

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
    let origin_x = ipp.first().copied().unwrap_or(0.0);
    let origin_y = ipp.get(1).copied().unwrap_or(0.0);
    let z = ipp.get(2).copied().unwrap_or(0.0);

    let transfer_syntax =
        TransferSyntaxKind::from_uid(obj.meta().transfer_syntax.trim_end_matches('\0'));
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
        origin_x,
        origin_y,
        z,
        hu: frame.pixels,
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
/// `ImagePositionPatient` (defaulting to the origin when absent).
///
/// # Errors
/// [`HeliosError::Dicom`] if the file cannot be parsed/decoded or a required
/// geometry attribute is missing or malformed; [`HeliosError::InvalidDomainValue`]
/// if the resulting grid dimensions/spacing are invalid.
pub fn load_ct_slice<T: Scalar>(
    path: impl AsRef<std::path::Path>,
) -> Result<Volume<T>, HeliosError> {
    let slice = read_slice(path.as_ref())?;
    let grid = VoxelGrid::axis_aligned(
        [slice.cols, slice.rows, 1],
        [
            T::from_f64(slice.col_spacing),
            T::from_f64(slice.row_spacing),
            T::from_f64(slice.thickness),
        ],
        Point3::new(
            T::from_f64(slice.origin_x),
            T::from_f64(slice.origin_y),
            T::from_f64(slice.z),
        ),
    )?;
    let mut data = vec![T::from_f64(0.0); slice.rows * slice.cols];
    scatter_slice(&mut data, &slice, 0, 1);
    Volume::from_shape_vec(grid, data)
}

/// Consistency tolerance for in-plane geometry and slice spacing (mm).
///
/// DICOM stores positions/spacings as decimal strings; axial-series slices share
/// an identical in-plane grid and a constant `ImagePositionPatient` z step. 1 µm
/// (`1e-3` mm) is tight enough to catch a missing slice (gap = 2× spacing) or a
/// mismatched grid while tolerating decimal-string round-off.
const GEOMETRY_TOL_MM: f64 = 1.0e-3;

/// Load a multi-slice DICOM CT/MVCT **series** into a 3-D [`Volume`] of Hounsfield
/// units.
///
/// Every slice is parsed and decoded (HU); the slices are validated to share an
/// identical in-plane grid (`Rows`/`Columns`/`PixelSpacing`/in-plane origin),
/// sorted by `ImagePositionPatient` z, and stacked along `k`. The z spacing is
/// derived from the (uniform) consecutive slice positions. Result geometry:
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
    let (ox, oy) = (slices[0].origin_x, slices[0].origin_y);
    for s in &slices[1..] {
        let consistent = s.rows == rows
            && s.cols == cols
            && (s.col_spacing - col_sp).abs() <= GEOMETRY_TOL_MM
            && (s.row_spacing - row_sp).abs() <= GEOMETRY_TOL_MM
            && (s.origin_x - ox).abs() <= GEOMETRY_TOL_MM
            && (s.origin_y - oy).abs() <= GEOMETRY_TOL_MM;
        if !consistent {
            return Err(HeliosError::Dicom {
                reason: "series slices have inconsistent in-plane geometry".to_owned(),
            });
        }
    }

    // Order along the stack axis and derive a uniform z spacing.
    slices.sort_by(|a, b| a.z.total_cmp(&b.z));
    let nz = slices.len();
    let z_spacing = if nz > 1 {
        slices[1].z - slices[0].z
    } else {
        slices[0].thickness
    };
    for w in slices.windows(2) {
        if ((w[1].z - w[0].z) - z_spacing).abs() > GEOMETRY_TOL_MM {
            return Err(HeliosError::Dicom {
                reason: "non-uniform slice spacing (missing or duplicate slice?)".to_owned(),
            });
        }
    }

    let grid = VoxelGrid::axis_aligned(
        [cols, rows, nz],
        [
            T::from_f64(col_sp),
            T::from_f64(row_sp),
            T::from_f64(z_spacing),
        ],
        Point3::new(T::from_f64(ox), T::from_f64(oy), T::from_f64(slices[0].z)),
    )?;

    let mut data = vec![T::from_f64(0.0); rows * cols * nz];
    for (k, slice) in slices.iter().enumerate() {
        scatter_slice(&mut data, slice, k, nz);
    }
    Volume::from_shape_vec(grid, data)
}

#[cfg(test)]
mod tests {
    use super::*;
    use dicom::core::smallvec::SmallVec;
    use dicom::core::{DataElement, PrimitiveValue, VR};
    use dicom::object::{FileMetaTableBuilder, InMemDicomObject};

    // Write a synthetic 2×2 unsigned-16 CT slice at position `z_mm` with a known
    // HU pattern and geometry (no external fixture). Slope 2, intercept −10;
    // PixelSpacing 0.8 (row) / 1.25 (col) mm; in-plane origin (5,7); a unique SOP
    // instance UID per file.
    fn write_slice_at(path: &std::path::Path, pixels: [u16; 4], z_mm: f64, uid: &str) {
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
        put(&mut obj, Tag(0x0008, 0x0018), VR::UI, uid.into());
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
        put(
            &mut obj,
            IMAGE_POSITION_PATIENT,
            VR::DS,
            format!("5\\7\\{z_mm}").into(),
        );
        put(
            &mut obj,
            Tag(0x7FE0, 0x0010),
            VR::OW,
            PrimitiveValue::U16(SmallVec::from_vec(pixels.to_vec())),
        );
        obj.with_meta(
            FileMetaTableBuilder::new()
                .media_storage_sop_class_uid("1.2.840.10008.5.1.4.1.1.2")
                .media_storage_sop_instance_uid(uid)
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
        write_slice_at(&path, [10, 20, 30, 40], 9.0, "2.25.4242");

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

    // Write a 3-slice series (deliberately out of z-order on disk) at z = 0, 4, 8
    // mm with a distinct HU tag per slice. Returns the temp dir (kept alive) + the
    // shuffled paths.
    fn write_series(dir: &std::path::Path) -> Vec<std::path::PathBuf> {
        // (z, base-pixels): HU = base·2 − 10.
        let specs = [
            (8.0_f64, [100_u16, 100, 100, 100], "2.25.3"),
            (0.0_f64, [10_u16, 10, 10, 10], "2.25.1"),
            (4.0_f64, [55_u16, 55, 55, 55], "2.25.2"),
        ];
        let mut paths = Vec::new();
        for (i, (z, px, uid)) in specs.iter().enumerate() {
            let p = dir.join(format!("s{i}.dcm"));
            write_slice_at(&p, *px, *z, uid);
            paths.push(p);
        }
        paths
    }

    #[test]
    fn series_stacks_sorted_by_position_with_derived_spacing() {
        let dir = tempfile::tempdir().unwrap();
        let paths = write_series(dir.path());
        let vol: Volume<f64> = load_ct_series(&paths).expect("series load");

        // dims = [cols, rows, nz] = [2, 2, 3]; z spacing derived from 0,4,8 → 4.
        assert_eq!(vol.grid().dims(), [2, 2, 3]);
        assert!((vol.grid().spacing()[2] - 4.0).abs() < 1e-12);
        // Origin z is the lowest slice position (0), regardless of input order.
        assert!((vol.grid().origin().z - 0.0).abs() < 1e-12);

        // k is sorted by z: k0 (z=0) HU=10·2−10=10, k1 (z=4) HU=55·2−10=100,
        // k2 (z=8) HU=100·2−10=190. Uniform in-plane, so any (i,j) matches.
        assert_eq!(vol.get(0, 0, 0), Some(10.0));
        assert_eq!(vol.get(1, 1, 1), Some(100.0));
        assert_eq!(vol.get(0, 1, 2), Some(190.0));
    }

    #[test]
    fn single_path_series_equals_single_slice_load() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("one.dcm");
        write_slice_at(&path, [10, 20, 30, 40], 9.0, "2.25.9");
        let series: Volume<f64> = load_ct_series(std::slice::from_ref(&path)).unwrap();
        let single: Volume<f64> = load_ct_slice(&path).unwrap();
        assert_eq!(series.grid().dims(), single.grid().dims());
        for j in 0..2 {
            for i in 0..2 {
                assert_eq!(series.get(i, j, 0), single.get(i, j, 0));
            }
        }
    }

    #[test]
    fn empty_and_non_uniform_series_error() {
        let empty: &[std::path::PathBuf] = &[];
        assert!(matches!(
            load_ct_series::<f64, _>(empty),
            Err(HeliosError::Dicom { .. })
        ));

        // Slices at z = 0, 4, 10 → gap (missing slice) → non-uniform spacing.
        let dir = tempfile::tempdir().unwrap();
        let mut paths = Vec::new();
        for (i, z) in [0.0, 4.0, 10.0].iter().enumerate() {
            let p = dir.path().join(format!("g{i}.dcm"));
            write_slice_at(&p, [10, 10, 10, 10], *z, &format!("2.25.{}", 20 + i));
            paths.push(p);
        }
        assert!(matches!(
            load_ct_series::<f64, _>(&paths),
            Err(HeliosError::Dicom { .. })
        ));
    }
}
