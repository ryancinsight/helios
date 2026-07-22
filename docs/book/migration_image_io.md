# Chapter 36 — Ritk: Image I/O — DICOM, NIfTI, PNG

Helios ingests CT, MVCT, and dose-evaluation images in three major
formats: **DICOM** (clinical standard), **NIfTI** (research), and PNG
(visualization).  These migrate to the shared **Ritk** image I/O crate,
which itself depends on `eunomia::RealField` and `leto::NdArray` for its
storage representation.

## The Ritk Surface

```rust
pub enum ImageFormat {
    Dicom,
    Nifti,
    Png,
}

pub struct ImageHandle<F: FloatElement> {
    pub voxels:    NdArray<F, Ix3>,
    pub spacing_mm: [F; 3],
    pub origin_mm:  Point3<F>,
    pub format:    ImageFormat,
    pub metadata:  ImageMetadata,
}

pub trait ImageReader {
    fn read<F: FloatElement>(path: &Path) -> Result<ImageHandle<F>, RitkError>;
}

pub trait ImageWriter {
    fn write<F: FloatElement>(handle: &ImageHandle<F>, path: &Path) -> Result<(), RitkError>;
}
```

A typical helios CT ingestion becomes:

```rust
use ritk_dicom::DicomReader;
use ritk_nifti::NiftiReader;

let handle = DicomReader::read::<f64>("/path/to/ct.dcm")?;
// or:
let handle = NiftiReader::read::<f64>("/path/to/ct.nii.gz")?;

let voxel_grid = VoxelGrid::from_handle(&handle);
```

The reader returns an `ImageHandle<F>` whose `voxels` slot is a
`NdArray<F, Ix3>` — the same shape that the dose and gamma-index code
already consumes.

## Migration Procedure

| Legacy | Atlas |
|---|---|
| `dicom` crate's `InMemoryDicomObject` | `ritk_dicom::DicomElement` |
| `nifti` crate's `NiftiObject` | `ritk_nifti::NiftiElement` (lazy archive) |
| `image` crate's `DynamicImage` | `ritk_png::PngReader` |
| hand-parsed `.dcm` header | `DicomReader::read` (one call) |
| runtime `as f64` conversion of HU | typed `HounsfieldUnit::try_from((handle, i, j, k))` |

## DICOM Tag Module

```rust
pub mod tags {
    pub const PATIENT_NAME:           &str = "0010:0010";
    pub const PIXEL_SPACING:          &str = "0028:0030";
    pub const SLICE_THICKNESS:        &str = "0018:0050";
    pub const MODALITY:               &str = "0008:0060";
    pub const RT_DOSE:                &str = "3004:0002";
}
```

DICOM tag constants are Atlas-wide — the same `PATIENT_NAME` works for
both CT ingestion and RT dose output.

## Helios-Specific Use Cases

| Use case | Where |
|---|---|
| CT ingestion | `helios-domain::dicom.rs` (DICOM → VoxelGrid) |
| MVCT registration | `helios-imaging::mvct.rs` (DICOM → ImageHandle → Coeus) |
| RT dose output | `helios-analysis::gamma.rs` (DoseGrid → ImageHandle → DICOM) |
| Visualization | `helios-python` (PNG) |

## Validation Examples

- [`tomotherapy_workflow`](examples/tomotherapy_workflow.md) — DICOM
  ingestion via Ritk.
- [`adaptive_rt_workflow`](examples/adaptive_rt_workflow.md) —
  MVCT acquisistion via Ritk.
- [`linac_dose_accumulation`](examples/linac_dose_accumulation.md) —
  RT dose output via Ritk.

## Further Reading

- [`ritk-image` source](../../../ritk/crates/ritk-image/)
- [`ritk-dicom` source](../../../ritk/crates/ritk-dicom/)
- [`ritk` source](../../../ritk/)
