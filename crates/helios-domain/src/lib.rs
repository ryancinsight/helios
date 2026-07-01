//! Helios domain layer.
//!
//! Patient and imaging geometry: the [`VoxelGrid`] that maps discrete voxel
//! indices to world/patient coordinates, and the dense [`Volume`] scalar field
//! (CT/MVCT densities, dose grids, projection stacks) sampled over that grid.
//!
//! The pure domain types depend only on `helios-core` (errors, validated units)
//! and `helios-math` (the `Scalar` seam and leto substrate). DICOM I/O
//! (`ritk-dicom`) is wired in at the boundary behind the `dicom` feature
//! ([`load_ct_slice`]): external file bytes become validated typed [`Volume`]s
//! there, keeping the core types infrastructure-agnostic.
#![forbid(unsafe_code)]
#![deny(missing_docs)]

#[cfg(feature = "dicom")]
mod dicom;
mod grid;
mod helical;
mod mlc;
mod volume;

#[cfg(feature = "dicom")]
pub use dicom::{load_ct_series, load_ct_slice};
pub use grid::VoxelGrid;
pub use helical::HelicalDelivery;
pub use mlc::{LeafOpenTimeSinogram, MlcModel};
pub use volume::Volume;
