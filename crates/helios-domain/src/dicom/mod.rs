//! DICOM CT/MVCT load path (real-input boundary), built on `ritk-dicom`.
//!
//! Reads a single-slice DICOM image into a Helios [`Volume`](crate::volume::Volume) of Hounsfield units:
//! `ritk-dicom` parses the file and decodes the pixel frame (applying the
//! `RescaleSlope`/`RescaleIntercept` calibration), and this module maps the frame
//! plus the geometry attributes (`Rows`, `Columns`, `PixelSpacing`,
//! `ImagePositionPatient`, `ImageOrientationPatient`) into a typed [`Volume`](crate::volume::Volume) on
//! an oriented [`VoxelGrid`](crate::grid::VoxelGrid). The pixel layout carries both `BitsAllocated` and
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

mod frame;
mod geometry;
mod volume;

pub use volume::{load_ct_series, load_ct_slice};

#[cfg(test)]
mod tests;
