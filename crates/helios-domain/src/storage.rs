//! Volumetric storage boundary: HDF5 persistence of [`Volume`]s via consus.
//!
//! Archives a dose/CT/MVCT [`Volume`] — data **and** grid geometry — to a
//! standard HDF5 file using the pure-Rust consus stack (the mandated Atlas
//! volumetric-storage component), and loads it back into a typed [`Volume`].
//! The file holds two root datasets:
//!
//! - `volume`: the 3-D scalar field, shape `[nx, ny, nz]`, IEEE-754 f64 LE, in
//!   the same C-contiguous `(i, j, k)` order the in-memory [`Volume`] uses.
//! - `geometry`: `[spacing_x, spacing_y, spacing_z, origin_x, origin_y,
//!   origin_z, x_axis.xyz, y_axis.xyz, z_axis.xyz]` (mm and dimensionless
//!   direction cosines), f64 LE. The final nine values are the world-space
//!   columns of the grid rotation, so the archive reconstructs the full rigid
//!   [`VoxelGrid`] pose without losing DICOM orientation.
//!
//! Values are serialized as `f64` at this boundary (the archive's fixed on-disk
//! precision, like the DICOM boundary's HU semantics); the in-memory `T: Scalar`
//! converts through `to_f64`/`from_f64`, so an `f64` volume round-trips bitwise.
//!
//! Feature-gated behind `storage` so the HDF5 machinery stays out of the core
//! build. The feature gates a complete implementation, not a stub.

use crate::grid::VoxelGrid;
use crate::volume::Volume;
use aequitas::systems::si::{quantities::Length, units::Millimeter};
use core::num::NonZeroUsize;
use eunomia::UnitScalar;
use helios_math::{Point3, Scalar, Vector3};
use helios_quantities::HeliosError;
use std::fs::File;

use consus_core::{ByteOrder, Datatype, LinkType, Shape};
use consus_hdf5::dataset::{Hdf5Dataset, StorageLayout};
use consus_hdf5::datatype::IeeeFloatEncoding;
use consus_hdf5::file::writer::{DatasetCreationProps, FileCreationProps, Hdf5FileBuilder};
use consus_hdf5::file::Hdf5File;

/// Root dataset name for the scalar field.
const VOLUME_DATASET: &str = "volume";
/// Root dataset name for the grid geometry (spacing, origin, and rotation).
const GEOMETRY_DATASET: &str = "geometry";
/// Geometry values: spacing, origin, then three rotation-matrix columns.
const GEOMETRY_VALUE_COUNT: usize = 15;

fn storage_err(step: &str, e: impl core::fmt::Display) -> HeliosError {
    HeliosError::Storage {
        reason: format!("{step}: {e}"),
    }
}

/// IEEE-754 binary64, little-endian — the archive datatype.
fn f64_le() -> Datatype {
    Datatype::Float {
        bits: NonZeroUsize::new(64).expect("invariant: 64 != 0"),
        byte_order: ByteOrder::LittleEndian,
    }
}

/// Serialize `values` as little-endian f64 bytes.
fn to_le_bytes(values: impl Iterator<Item = f64>) -> Vec<u8> {
    let mut bytes = Vec::new();
    for v in values {
        bytes.extend_from_slice(&v.to_le_bytes());
    }
    bytes
}

/// Convert already admitted binary64 archive bytes directly to the destination scalar.
fn decode_archive<T: Scalar>(bytes: &[u8]) -> impl ExactSizeIterator<Item = T> + '_ {
    let (words, remainder) = bytes.as_chunks::<{ size_of::<f64>() }>();
    debug_assert!(remainder.is_empty(), "validated archive byte count");
    words
        .iter()
        .map(|word| T::from_f64(f64::from_le_bytes(*word)))
}

/// Archive `volume` (data + rigid grid geometry) as an HDF5 file at `path`.
///
/// # Errors
/// [`HeliosError::Storage`] if HDF5 encoding or the filesystem write fails.
pub fn save_volume_hdf5<T: Scalar>(
    volume: &Volume<T>,
    path: impl AsRef<std::path::Path>,
) -> Result<(), HeliosError> {
    let grid = volume.grid();
    let [nx, ny, nz] = grid.dims();

    // Field bytes straight from the zero-copy slice view — already the archive's
    // C-contiguous (i, j, k) order (the Volume layout contract).
    let field: Vec<f64> = volume.as_slice().iter().map(|v| v.to_f64()).collect();
    let spacing = grid.spacing();
    let origin = grid.origin();
    let rotation = grid.pose().rotation;
    let zero = T::ZERO;
    let one = T::ONE;
    let x_axis = rotation.transform_vector(Vector3::new(one, zero, zero));
    let y_axis = rotation.transform_vector(Vector3::new(zero, one, zero));
    let z_axis = rotation.transform_vector(Vector3::new(zero, zero, one));
    let geometry = [
        spacing[0].to_f64(),
        spacing[1].to_f64(),
        spacing[2].to_f64(),
        origin.x.to_f64(),
        origin.y.to_f64(),
        origin.z.to_f64(),
        x_axis.x.to_f64(),
        x_axis.y.to_f64(),
        x_axis.z.to_f64(),
        y_axis.x.to_f64(),
        y_axis.y.to_f64(),
        y_axis.z.to_f64(),
        z_axis.x.to_f64(),
        z_axis.y.to_f64(),
        z_axis.z.to_f64(),
    ];

    let mut builder = Hdf5FileBuilder::new(FileCreationProps::default());
    let dcpl = DatasetCreationProps::default();
    builder
        .add_dataset(
            VOLUME_DATASET,
            &f64_le(),
            &Shape::fixed(&[nx, ny, nz]),
            &to_le_bytes(field.into_iter()),
            &dcpl,
        )
        .map_err(|e| storage_err("encode volume dataset", e))?;
    builder
        .add_dataset(
            GEOMETRY_DATASET,
            &f64_le(),
            &Shape::fixed(&[GEOMETRY_VALUE_COUNT]),
            &to_le_bytes(geometry.into_iter()),
            &dcpl,
        )
        .map_err(|e| storage_err("encode geometry dataset", e))?;
    let bytes = builder
        .finish()
        .map_err(|e| storage_err("finalize HDF5 file", e))?;
    std::fs::write(path.as_ref(), bytes).map_err(|e| storage_err("write file", e))
}

/// Admit an archive dataset using its original HDF5 datatype message.
fn admit_root_dataset(
    file: &Hdf5File<File>,
    children: &[(String, u64, LinkType)],
    name: &str,
) -> Result<Hdf5Dataset, HeliosError> {
    let (_, header_addr, _) = children
        .iter()
        .find(|(entry, _, link)| entry == name && *link == LinkType::Hard)
        .ok_or_else(|| storage_err(name, "expected a named root dataset with a hard link"))?;
    let dataset = file
        .dataset_at_with_ieee_float_encoding(*header_addr, IeeeFloatEncoding::Binary64LittleEndian)
        .map_err(|error| {
            storage_err(&format!("dataset '{name}' IEEE datatype admission"), error)
        })?;
    if dataset.shape.has_unlimited()
        || dataset.layout != StorageLayout::Contiguous
        || !dataset.filters.is_empty()
    {
        return Err(storage_err(
            name,
            "expected fixed, unfiltered contiguous storage",
        ));
    }
    Ok(dataset)
}

/// Bound the contiguous transfer by both the shape and the underlying file length.
fn contiguous_region(
    dataset: &Hdf5Dataset,
    name: &str,
    file_len: u64,
) -> Result<(u64, usize), HeliosError> {
    let byte_count = dataset
        .shape
        .current_dims()
        .iter()
        .try_fold(size_of::<f64>(), |count, &dim| count.checked_mul(dim))
        .ok_or_else(|| storage_err(name, "dataset byte count overflows usize"))?;
    let address = dataset
        .data_address
        .ok_or_else(|| storage_err(name, "missing contiguous data address"))?;
    let end = u64::try_from(byte_count)
        .ok()
        .and_then(|count| address.checked_add(count));
    if end.is_none_or(|end| end > file_len) {
        return Err(storage_err(
            name,
            "contiguous data exceeds the archive length",
        ));
    }
    Ok((address, byte_count))
}

/// Load a [`Volume`] previously archived by [`save_volume_hdf5`] from `path`.
///
/// The source file must remain stable during the synchronous read.
///
/// # Errors
/// [`HeliosError::Storage`] if the archive cannot be read or its datasets are
/// missing, malformed, or incompatible with the binary64 little-endian schema.
/// [`HeliosError::InvalidDomainValue`] if the geometry does not form a valid grid.
pub fn load_volume_hdf5<T: Scalar + UnitScalar>(
    path: impl AsRef<std::path::Path>,
) -> Result<Volume<T>, HeliosError> {
    let input = File::open(path.as_ref()).map_err(|error| storage_err("open file", error))?;
    let file_len =
        consus_io::Length::len(&input).map_err(|error| storage_err("file length", error))?;
    let file = Hdf5File::open(input).map_err(|error| storage_err("open HDF5", error))?;
    let children = file
        .list_root_group()
        .map_err(|error| storage_err("list root group", error))?;

    // Both datatypes must be admitted before either payload is read.
    let geometry = admit_root_dataset(&file, &children, GEOMETRY_DATASET)?;
    let volume = admit_root_dataset(&file, &children, VOLUME_DATASET)?;
    let geometry_dims = geometry.shape.current_dims();
    if geometry_dims.as_slice() != [GEOMETRY_VALUE_COUNT] {
        return Err(storage_err(
            GEOMETRY_DATASET,
            format_args!("expected shape [{GEOMETRY_VALUE_COUNT}], got {geometry_dims:?}"),
        ));
    }
    let volume_dims = volume.shape.current_dims();
    let dims = <[usize; 3]>::try_from(volume_dims.as_slice()).map_err(|_| {
        storage_err(
            VOLUME_DATASET,
            format_args!("expected rank 3, got {}", volume_dims.len()),
        )
    })?;
    let (geometry_addr, _) = contiguous_region(&geometry, GEOMETRY_DATASET, file_len)?;
    let (volume_addr, volume_bytes) = contiguous_region(&volume, VOLUME_DATASET, file_len)?;

    let mut geometry_bytes = [0; GEOMETRY_VALUE_COUNT * size_of::<f64>()];
    file.read_contiguous_dataset_bytes(geometry_addr, 0, &mut geometry_bytes)
        .map_err(|error| storage_err("read geometry dataset bytes", error))?;
    let mut geometry_values = decode_archive::<T>(&geometry_bytes);
    let geom: [T; GEOMETRY_VALUE_COUNT] =
        core::array::from_fn(|_| geometry_values.next().expect("one word per geometry value"));

    let rotation = helios_math::UnitQuaternion::try_from_rotation_columns(
        Vector3::new(geom[6], geom[7], geom[8]),
        Vector3::new(geom[9], geom[10], geom[11]),
        Vector3::new(geom[12], geom[13], geom[14]),
        T::EPSILON.sqrt(),
    )
    .map_err(|error| HeliosError::Storage {
        reason: format!("geometry rotation is not a rigid right-handed basis: {error}"),
    })?;
    // Archived geometry is millimetres, matching the grid's storage units.
    let grid = VoxelGrid::oriented(
        dims,
        [
            Length::from_unit::<Millimeter>(geom[0]),
            Length::from_unit::<Millimeter>(geom[1]),
            Length::from_unit::<Millimeter>(geom[2]),
        ],
        Point3::new(geom[3], geom[4], geom[5]),
        rotation,
    )?;
    // Read the admitted contiguous region through the caller-buffer API. The
    // provider's allocating raw reader limits a single allocation to 64 MiB.
    let mut bytes = Vec::new();
    bytes
        .try_reserve_exact(volume_bytes)
        .map_err(|error| storage_err("allocate volume dataset bytes", error))?;
    bytes.resize(volume_bytes, 0);
    file.read_contiguous_dataset_bytes(volume_addr, 0, &mut bytes)
        .map_err(|error| storage_err("read volume dataset bytes", error))?;
    let values = decode_archive::<T>(&bytes);
    let mut field = Vec::new();
    field
        .try_reserve_exact(values.len())
        .map_err(|error| storage_err("allocate volume values", error))?;
    field.extend(values);
    Volume::from_shape_vec(grid, field)
}

#[cfg(test)]
mod tests {
    use super::*;
    use aequitas::systems::si::{quantities::Length, units::Millimeter};

    fn test_volume() -> Volume<f64> {
        // Distinct value per voxel locks the (i, j, k) storage order; non-trivial
        // spacing/origin lock the geometry round-trip.
        let grid = VoxelGrid::axis_aligned(
            [4, 3, 2],
            [
                Length::from_unit::<Millimeter>(1.25),
                Length::from_unit::<Millimeter>(2.0),
                Length::from_unit::<Millimeter>(3.5),
            ],
            Point3::new(-5.0, 7.5, 11.0),
        )
        .expect("grid");
        Volume::from_shape_fn(grid, |idx| {
            100.0 * idx[0] as f64 + 10.0 * idx[1] as f64 + idx[2] as f64 + 0.125
        })
    }

    #[test]
    fn volume_round_trips_bitwise_through_hdf5() {
        let dir = tempfile::tempdir().expect("synthetic storage fixture round-trips");
        let path = dir.path().join("dose.h5");
        let original = test_volume();
        save_volume_hdf5(&original, &path).expect("save");

        let loaded: Volume<f64> = load_volume_hdf5(&path).expect("load");
        // Geometry is reconstructed exactly (f64 → f64 is bitwise).
        assert_eq!(loaded.grid().dims(), [4, 3, 2]);
        assert_eq!(loaded.grid().spacing(), [1.25, 2.0, 3.5]);
        assert_eq!(loaded.grid().origin().x, -5.0);
        assert_eq!(loaded.grid().origin().z, 11.0);
        // Every voxel is bitwise-identical (no evaluation-order change: pure
        // serialize/deserialize).
        for i in 0..4 {
            for j in 0..3 {
                for k in 0..2 {
                    assert_eq!(loaded.get(i, j, k), original.get(i, j, k), "({i},{j},{k})");
                }
            }
        }
    }

    #[test]
    fn hdf5_file_is_readable_as_standard_hdf5() {
        // The archive starts with the HDF5 superblock signature, so external
        // tools (h5py etc.) can open it — the interoperability contract.
        let dir = tempfile::tempdir().expect("synthetic storage fixture round-trips");
        let path = dir.path().join("sig.h5");
        save_volume_hdf5(&test_volume(), &path).expect("save");
        let bytes = std::fs::read(&path).expect("synthetic storage fixture round-trips");
        assert_eq!(&bytes[0..8], b"\x89HDF\r\n\x1a\n", "HDF5 signature");
    }

    #[test]
    fn f32_volume_round_trips_through_the_f64_archive() {
        let dir = tempfile::tempdir().expect("synthetic storage fixture round-trips");
        let path = dir.path().join("f32.h5");
        let grid = VoxelGrid::<f32>::axis_aligned(
            [2, 2, 2],
            [Length::from_unit::<Millimeter>(1.0); 3],
            Point3::new(0.0, 0.0, 0.0),
        )
        .expect("synthetic storage fixture round-trips");
        let original = Volume::from_shape_fn(grid, |idx| (idx[0] + 2 * idx[1] + 4 * idx[2]) as f32);
        save_volume_hdf5(&original, &path).expect("save");
        let loaded: Volume<f32> = load_volume_hdf5(&path).expect("load");
        for i in 0..2 {
            for j in 0..2 {
                for k in 0..2 {
                    // f32 → f64 → f32 is exact (f64 represents every f32).
                    assert_eq!(loaded.get(i, j, k), original.get(i, j, k));
                }
            }
        }
    }

    #[test]
    fn oriented_grid_pose_round_trips_through_hdf5() {
        let dir = tempfile::tempdir().expect("synthetic storage fixture round-trips");
        let path = dir.path().join("oriented.h5");
        let rotation = helios_math::UnitQuaternion::try_from_rotation_columns(
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(-1.0, 0.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
            1.0e-12,
        )
        .expect("right-handed quarter-turn basis");
        let grid = VoxelGrid::oriented(
            [4, 3, 2],
            [
                Length::from_unit::<Millimeter>(1.25),
                Length::from_unit::<Millimeter>(2.0),
                Length::from_unit::<Millimeter>(3.5),
            ],
            Point3::new(-5.0, 7.5, 11.0),
            rotation,
        )
        .expect("grid");
        let original = Volume::from_shape_fn(grid, |idx| {
            100.0 * idx[0] as f64 + 10.0 * idx[1] as f64 + idx[2] as f64
        });
        save_volume_hdf5(&original, &path).expect("save");

        let loaded: Volume<f64> = load_volume_hdf5(&path).expect("load");
        let expected = original.grid().index_to_world(Point3::new(2.0, 1.0, 1.0));
        let actual = loaded.grid().index_to_world(Point3::new(2.0, 1.0, 1.0));
        assert!((actual.x - expected.x).abs() <= 1.0e-12);
        assert!((actual.y - expected.y).abs() <= 1.0e-12);
        assert!((actual.z - expected.z).abs() <= 1.0e-12);
        assert_eq!(loaded.get(2, 1, 1), original.get(2, 1, 1));
    }

    #[test]
    fn missing_file_and_garbage_bytes_are_storage_errors() {
        assert!(matches!(
            load_volume_hdf5::<f64>("does_not_exist.h5"),
            Err(HeliosError::Storage { .. })
        ));
        let dir = tempfile::tempdir().expect("synthetic storage fixture round-trips");
        let path = dir.path().join("garbage.h5");
        std::fs::write(&path, b"not an hdf5 file").expect("synthetic storage fixture round-trips");
        assert!(matches!(
            load_volume_hdf5::<f64>(&path),
            Err(HeliosError::Storage { .. })
        ));
    }
}
