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
use core::num::NonZeroUsize;
use helios_core::HeliosError;
use helios_math::{Point3, Scalar, Vector3};

use consus_core::{ByteOrder, Datatype, Shape};
use consus_hdf5::file::writer::{DatasetCreationProps, FileCreationProps, Hdf5FileBuilder};
use consus_hdf5::file::Hdf5File;
use consus_io::MemCursor;

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

/// Parse little-endian f64 bytes; errors if the length is not a multiple of 8.
fn from_le_bytes(bytes: &[u8]) -> Result<Vec<f64>, HeliosError> {
    if !bytes.len().is_multiple_of(8) {
        return Err(HeliosError::Storage {
            reason: format!("dataset byte length {} is not a multiple of 8", bytes.len()),
        });
    }
    Ok(bytes
        .chunks_exact(8)
        .map(|c| f64::from_le_bytes(c.try_into().expect("chunks_exact(8) yields 8 bytes")))
        .collect())
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

/// Read the full contiguous payload of the named root dataset.
fn read_root_dataset(
    file: &Hdf5File<MemCursor>,
    name: &str,
) -> Result<(Vec<u8>, Vec<usize>), HeliosError> {
    let children = file
        .list_root_group()
        .map_err(|e| storage_err("list root group", e))?;
    let (_, header_addr, _) = children
        .into_iter()
        .find(|(n, _, _)| n == name)
        .ok_or_else(|| HeliosError::Storage {
            reason: format!("dataset '{name}' not found in file"),
        })?;
    let ds = file
        .dataset_at(header_addr)
        .map_err(|e| storage_err("read dataset metadata", e))?;
    if ds.shape.has_unlimited() {
        return Err(HeliosError::Storage {
            reason: format!("dataset '{name}' has a non-fixed shape"),
        });
    }
    let dims: Vec<usize> = ds.shape.current_dims().to_vec();
    let n_elems: usize = dims.iter().product();
    let data_addr = ds.data_address.ok_or_else(|| HeliosError::Storage {
        reason: format!("dataset '{name}' has no contiguous data address"),
    })?;
    let mut buf = vec![0u8; n_elems * 8];
    file.read_contiguous_dataset_bytes(data_addr, 0, &mut buf)
        .map_err(|e| storage_err("read dataset bytes", e))?;
    Ok((buf, dims))
}

/// Load a [`Volume`] previously archived by [`save_volume_hdf5`] from `path`.
///
/// # Errors
/// [`HeliosError::Storage`] if the file cannot be read/parsed or its datasets are
/// missing/malformed; [`HeliosError::InvalidDomainValue`] if the stored geometry
/// does not form a valid grid.
pub fn load_volume_hdf5<T: Scalar>(
    path: impl AsRef<std::path::Path>,
) -> Result<Volume<T>, HeliosError> {
    let bytes = std::fs::read(path.as_ref()).map_err(|e| storage_err("read file", e))?;
    let file =
        Hdf5File::open(MemCursor::from_bytes(bytes)).map_err(|e| storage_err("open HDF5", e))?;

    let (geom_bytes, geom_dims) = read_root_dataset(&file, GEOMETRY_DATASET)?;
    if geom_dims != [GEOMETRY_VALUE_COUNT] {
        return Err(HeliosError::Storage {
            reason: format!(
                "geometry dataset has shape {geom_dims:?}, expected [{GEOMETRY_VALUE_COUNT}]"
            ),
        });
    }
    let geom = from_le_bytes(&geom_bytes)?;
    if geom.len() != GEOMETRY_VALUE_COUNT {
        return Err(HeliosError::Storage {
            reason: format!(
                "geometry dataset has {} values, expected {GEOMETRY_VALUE_COUNT}",
                geom.len()
            ),
        });
    }

    let (field_bytes, field_dims) = read_root_dataset(&file, VOLUME_DATASET)?;
    if field_dims.len() != 3 {
        return Err(HeliosError::Storage {
            reason: format!("volume dataset has rank {}, expected 3", field_dims.len()),
        });
    }
    let field = from_le_bytes(&field_bytes)?;

    let rotation = helios_math::UnitQuaternion::try_from_rotation_columns(
        Vector3::new(
            T::from_f64(geom[6]),
            T::from_f64(geom[7]),
            T::from_f64(geom[8]),
        ),
        Vector3::new(
            T::from_f64(geom[9]),
            T::from_f64(geom[10]),
            T::from_f64(geom[11]),
        ),
        Vector3::new(
            T::from_f64(geom[12]),
            T::from_f64(geom[13]),
            T::from_f64(geom[14]),
        ),
        T::EPSILON.sqrt(),
    )
    .map_err(|error| HeliosError::Storage {
        reason: format!("geometry rotation is not a rigid right-handed basis: {error}"),
    })?;
    let grid = VoxelGrid::oriented(
        [field_dims[0], field_dims[1], field_dims[2]],
        [
            T::from_f64(geom[0]),
            T::from_f64(geom[1]),
            T::from_f64(geom[2]),
        ],
        Point3::new(
            T::from_f64(geom[3]),
            T::from_f64(geom[4]),
            T::from_f64(geom[5]),
        ),
        rotation,
    )?;
    Volume::from_shape_vec(grid, field.into_iter().map(T::from_f64).collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn test_volume() -> Volume<f64> {
        // Distinct value per voxel locks the (i, j, k) storage order; non-trivial
        // spacing/origin lock the geometry round-trip.
        let grid =
            VoxelGrid::axis_aligned([4, 3, 2], [1.25, 2.0, 3.5], Point3::new(-5.0, 7.5, 11.0))
                .expect("grid");
        Volume::from_shape_fn(grid, |idx| {
            100.0 * idx[0] as f64 + 10.0 * idx[1] as f64 + idx[2] as f64 + 0.125
        })
    }

    #[test]
    fn volume_round_trips_bitwise_through_hdf5() {
        let dir = tempfile::tempdir().unwrap();
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
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("sig.h5");
        save_volume_hdf5(&test_volume(), &path).expect("save");
        let bytes = std::fs::read(&path).unwrap();
        assert_eq!(&bytes[0..8], b"\x89HDF\r\n\x1a\n", "HDF5 signature");
    }

    #[test]
    fn f32_volume_round_trips_through_the_f64_archive() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("f32.h5");
        let grid =
            VoxelGrid::<f32>::axis_aligned([2, 2, 2], [1.0, 1.0, 1.0], Point3::new(0.0, 0.0, 0.0))
                .unwrap();
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
        let dir = tempfile::tempdir().unwrap();
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
            [1.25, 2.0, 3.5],
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
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("garbage.h5");
        std::fs::write(&path, b"not an hdf5 file").unwrap();
        assert!(matches!(
            load_volume_hdf5::<f64>(&path),
            Err(HeliosError::Storage { .. })
        ));
    }
}
