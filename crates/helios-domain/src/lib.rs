//! Helios domain layer.
//!
//! Patient and imaging geometry: the [`VoxelGrid`] that maps discrete voxel
//! indices to world/patient coordinates, and the dense [`Volume`] scalar field
//! (CT/MVCT densities, dose grids, projection stacks) sampled over that grid.
//!
//! This layer is pure and infrastructure-agnostic: it depends only on
//! `helios-core` (errors, validated units) and `helios-math` (the `Scalar` seam
//! and leto substrate). DICOM I/O (`ritk-io`) and geometry kernels (`gaia`) are
//! wired in at the boundary in later increments; the types here are the domain
//! representation those loaders populate.
#![forbid(unsafe_code)]
#![deny(missing_docs)]

mod grid;
mod helical;
mod volume;

pub use grid::VoxelGrid;
pub use helical::HelicalDelivery;
pub use volume::Volume;
