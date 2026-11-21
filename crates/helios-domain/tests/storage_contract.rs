//! Compatibility tests for the Helios HDF5 archive schema.

#![cfg(feature = "storage")]

use std::fs;
use std::num::NonZeroUsize;
use std::path::Path;

use consus_core::{ByteOrder, Datatype, Shape};
use consus_hdf5::file::{reader, writer, Hdf5File};
use consus_hdf5::object_header::message_types;
use consus_io::MemCursor;
use helios_domain::{load_volume_hdf5, save_volume_hdf5};
use helios_quantities::HeliosError;
use tracing::{info_span, subscriber::DefaultGuard};
use tracing_subscriber::fmt::format::FmtSpan;

const DIMS: [usize; 3] = [2, 2, 2];
const GEOMETRY: [f64; 15] = [
    1.25, 2.0, 3.5, -5.0, 7.5, 11.0, 0.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0,
];
const WORDS: [u64; 8] = [
    0x3ff0_0000_0000_0000, // 1
    0x8000_0000_0000_0000, // negative zero
    0x0000_0000_0000_0001, // smallest positive binary64
    0x3ff0_0000_1000_0000, // 1 + 2^-24, halfway between adjacent binary32 values
    0x3ff0_0000_1000_0001, // one binary64 ULP above that tie
    0xc014_0000_0000_0000, // -5
    0x400a_0000_0000_0000, // 3.25
    0x4026_0000_0000_0000, // 11
];
const EXPECTED_SINGLE: [u32; 8] = [
    0x3f80_0000,
    0x8000_0000,
    0,
    0x3f80_0000,
    0x3f80_0001,
    0xc0a0_0000,
    0x4050_0000,
    0x4130_0000,
];

fn trace() -> DefaultGuard {
    tracing::subscriber::set_default(
        tracing_subscriber::fmt()
            .with_test_writer()
            .with_span_events(FmtSpan::FULL)
            .finish(),
    )
}

fn datatype(order: ByteOrder) -> Datatype {
    Datatype::Float {
        bits: NonZeroUsize::new(64).expect("binary64 has nonzero width"),
        byte_order: order,
    }
}

fn payload(words: impl Iterator<Item = u64>, order: ByteOrder) -> Vec<u8> {
    words
        .flat_map(|bits| match order {
            ByteOrder::LittleEndian => bits.to_le_bytes(),
            ByteOrder::BigEndian => bits.to_be_bytes(),
        })
        .collect()
}

fn archive(
    dims: &[usize],
    words: impl Iterator<Item = u64>,
    volume_order: ByteOrder,
    geometry: &[f64],
    geometry_order: ByteOrder,
) -> Vec<u8> {
    let mut builder = writer::Hdf5FileBuilder::new(writer::FileCreationProps::default());
    let props = writer::DatasetCreationProps::default();
    builder
        .add_dataset(
            "volume",
            &datatype(volume_order),
            &Shape::fixed(dims),
            &payload(words, volume_order),
            &props,
        )
        .expect("write the volume dataset");
    builder
        .add_dataset(
            "geometry",
            &datatype(geometry_order),
            &Shape::fixed(&[geometry.len()]),
            &payload(geometry.iter().map(|value| value.to_bits()), geometry_order),
            &props,
        )
        .expect("write the geometry dataset");
    builder.finish().expect("finish the HDF5 fixture")
}

fn canonical() -> Vec<u8> {
    archive(
        &DIMS,
        WORDS.into_iter(),
        ByteOrder::LittleEndian,
        &GEOMETRY,
        ByteOrder::LittleEndian,
    )
}

fn edit_header(
    bytes: Vec<u8>,
    name: &str,
    message_type: u16,
    edit: impl FnOnce(&mut [u8]),
) -> Vec<u8> {
    let file = Hdf5File::open(MemCursor::from_bytes(bytes))
        .expect("open the canonical fixture for header editing");
    let (_, address, _) = file
        .list_root_group()
        .expect("list the fixture's root datasets")
        .into_iter()
        .find(|(entry, _, _)| entry == name)
        .expect("find the dataset whose header is being edited");
    let mut header = reader::read_object_header(file.source(), address, file.context())
        .expect("read the original object header");
    let message = header
        .messages
        .iter_mut()
        .find(|message| message.message_type == message_type)
        .expect("find the requested message in the object header");
    edit(&mut message.data);
    let messages: Vec<_> = header
        .messages
        .iter()
        .map(|message| (message.message_type, message.data.as_slice()))
        .collect();
    let mut state = writer::WriteState::new(writer::FileCreationProps::default());
    state.eof = address;
    let mut sink = file.source().clone();
    assert_eq!(
        writer::write_object_header_v2(&mut sink, &mut state, &messages)
            .expect("rewrite the fixture header and checksum"),
        address
    );
    let file = Hdf5File::open(sink).expect("reopen the fixture with the edited header");
    // A well-formed header must reach the ordinary metadata lookup.
    let _ = file
        .dataset_at(address)
        .expect("the edited object header remains structurally readable");
    file.source().clone().into_bytes()
}

fn check_values(path: &Path) {
    let (double, single) = {
        let _act = info_span!("act").entered();
        (
            load_volume_hdf5::<f64>(path).expect("load the canonical archive into binary64"),
            load_volume_hdf5::<f32>(path).expect("load the canonical archive into binary32"),
        )
    };
    let _assert = info_span!("assert").entered();
    assert_eq!(double.grid().dims(), DIMS);
    assert_eq!(single.grid().dims(), DIMS);
    assert!(double.as_slice().iter().map(|x| x.to_bits()).eq(WORDS));
    assert!(single
        .as_slice()
        .iter()
        .map(|x| x.to_bits())
        .eq(EXPECTED_SINGLE));
    assert_eq!(double.get(1, 0, 1), Some(-5.0));
    assert_eq!(single.get(1, 0, 1), Some(-5.0));
    assert_eq!(double.grid().spacing(), [1.25, 2.0, 3.5]);
    assert_eq!(single.grid().spacing(), [1.25, 2.0, 3.5]);
    let double_origin = double.grid().origin();
    let single_origin = single.grid().origin();
    assert_eq!(
        [double_origin.x, double_origin.y, double_origin.z],
        [-5.0, 7.5, 11.0]
    );
    assert_eq!(
        [single_origin.x, single_origin.y, single_origin.z],
        [-5.0, 7.5, 11.0]
    );
}

fn storage_error(bytes: &[u8]) -> String {
    let dir = tempfile::tempdir().expect("create a directory for the rejection fixture");
    let path = dir.path().join("archive.h5");
    fs::write(&path, bytes).expect("write the rejection fixture");
    let result = {
        let _act = info_span!("act").entered();
        load_volume_hdf5::<f64>(&path)
    };
    let _assert = info_span!("assert").entered();
    match result {
        Err(HeliosError::Storage { reason }) => reason,
        other => panic!("expected storage error, got {other:?}"),
    }
}

#[test]
fn canonical_archive_preserves_bits_pose_and_source_bytes() {
    let _trace = trace();
    let (dir, bytes) = {
        let _arrange = info_span!("arrange").entered();
        (
            tempfile::tempdir().expect("create a directory for the canonical archive"),
            canonical(),
        )
    };
    let path = dir.path().join("canonical.h5");
    fs::write(&path, &bytes).expect("write the canonical fixture");
    check_values(&path);

    let volume =
        load_volume_hdf5::<f64>(&path).expect("load the binary64 volume for pose verification");
    let single =
        load_volume_hdf5::<f32>(&path).expect("load the binary32 volume for pose verification");
    let world = single.grid().index_to_world([1.0, 1.0, 1.0].into());
    // Origin plus the three quarter-turned spacing vectors.
    for (value, target) in [world.x, world.y, world.z]
        .into_iter()
        .zip([-7.0, 8.75, 14.5])
    {
        assert!((value - target).abs() <= 8.0 * f32::EPSILON * target.abs());
    }
    let rotation = volume.grid().pose().rotation;
    let x_axis = rotation.transform_vector([1.0, 0.0, 0.0].into());
    let y_axis = rotation.transform_vector([0.0, 1.0, 0.0].into());
    let z_axis = rotation.transform_vector([0.0, 0.0, 1.0].into());
    // sqrt(2)/2 reconstruction and quaternion rotation involve rounded
    // products and sums. The fixed quarter-turn permits a small ULP envelope.
    for (actual, expected) in [
        (x_axis, [0.0, 1.0, 0.0]),
        (y_axis, [-1.0, 0.0, 0.0]),
        (z_axis, [0.0, 0.0, 1.0]),
    ] {
        for (value, target) in [actual.x, actual.y, actual.z].into_iter().zip(expected) {
            assert!((value - target).abs() <= 8.0 * f64::EPSILON);
        }
    }
    assert_eq!(
        fs::read(&path).expect("read the original canonical archive after loading"),
        bytes
    );
    let saved = dir.path().join("saved.h5");
    save_volume_hdf5(&volume, &saved).expect("write the canonical Helios archive");
    check_values(&saved);
    // Optional fixed output for the cross-revision SHA-256 comparison.
    if let Some(output) = std::env::var_os("HELIOS_CANONICAL_ARCHIVE_OUTPUT") {
        fs::copy(saved, output).expect("copy the Helios archive for cross-revision comparison");
    }
}

fn expect_bad_bias(name: &str, bias: u32) {
    let _trace = trace();
    let bytes = {
        let _arrange = info_span!("arrange", name, bias).entered();
        edit_header(canonical(), name, message_types::DATATYPE, |data| {
            assert_eq!(&data[16..20], &1023u32.to_le_bytes());
            data[16..20].copy_from_slice(&bias.to_le_bytes());
        })
    };
    let reason = storage_error(&bytes);
    assert!(
        reason.contains(name) && reason.contains("Binary64LittleEndian"),
        "{reason}"
    );
}

#[test]
fn rejects_custom_bias_one_below_in_volume() {
    expect_bad_bias("volume", 1022);
}

#[test]
fn rejects_custom_bias_one_below_in_geometry() {
    expect_bad_bias("geometry", 1022);
}

#[test]
fn rejects_custom_bias_one_above_in_volume() {
    expect_bad_bias("volume", 1024);
}

#[test]
fn rejects_custom_bias_one_above_in_geometry() {
    expect_bad_bias("geometry", 1024);
}

fn expect_big_endian_rejection(name: &str, volume_order: ByteOrder, geometry_order: ByteOrder) {
    let _trace = trace();
    let bytes = {
        let _arrange = info_span!("arrange", name).entered();
        archive(
            &DIMS,
            WORDS.into_iter(),
            volume_order,
            &GEOMETRY,
            geometry_order,
        )
    };
    let reason = storage_error(&bytes);
    assert!(
        reason.contains(name) && reason.contains("Binary64LittleEndian"),
        "{reason}"
    );
}

#[test]
fn rejects_big_endian_volume() {
    expect_big_endian_rejection("volume", ByteOrder::BigEndian, ByteOrder::LittleEndian);
}

#[test]
fn rejects_big_endian_geometry() {
    expect_big_endian_rejection("geometry", ByteOrder::LittleEndian, ByteOrder::BigEndian);
}

#[test]
fn full_width_padding_preserves_decoded_values() {
    let _trace = trace();
    let dir = tempfile::tempdir().expect("create a directory for the padding fixtures");
    let path = dir.path().join("padding.h5");
    for padding in 0u8..8 {
        let bytes = {
            let _arrange = info_span!("arrange", padding).entered();
            ["volume", "geometry"]
                .into_iter()
                .fold(canonical(), |bytes, name| {
                    edit_header(bytes, name, message_types::DATATYPE, |data| {
                        assert_eq!(data[1], 0x20);
                        data[1] |= padding << 1;
                    })
                })
        };
        fs::write(&path, &bytes).expect("write the selected padding fixture");
        check_values(&path);
        assert_eq!(
            fs::read(&path).expect("read the padding fixture after loading"),
            bytes
        );
    }
}

#[test]
fn both_datatypes_are_admitted_before_geometry_payload_io() {
    let _trace = trace();
    let bytes = {
        let _arrange = info_span!("arrange").entered();
        let bytes = edit_header(
            canonical(),
            "geometry",
            message_types::DATA_LAYOUT,
            |data| {
                assert_eq!(&data[..2], &[3, 1]); // writer's v3 contiguous layout
                data[2..10].copy_from_slice(&(u64::MAX - 4096).to_le_bytes());
            },
        );
        edit_header(bytes, "volume", message_types::DATATYPE, |data| {
            data[16..20].copy_from_slice(&1022u32.to_le_bytes());
        })
    };
    let reason = storage_error(&bytes);
    assert!(
        reason.contains("volume") && reason.contains("Binary64LittleEndian"),
        "{reason}"
    );
}

#[test]
fn rejects_wrong_geometry_shape_and_nonrigid_basis() {
    let _trace = trace();
    let le = ByteOrder::LittleEndian;
    let (wrong_shape, nonrigid) = {
        let _arrange = info_span!("arrange").entered();
        let mut reflected = GEOMETRY;
        reflected[14] = -1.0;
        (
            archive(&DIMS, WORDS.into_iter(), le, &GEOMETRY[..14], le),
            archive(&DIMS, WORDS.into_iter(), le, &reflected, le),
        )
    };
    assert!(storage_error(&wrong_shape).contains("geometry"));
    assert!(storage_error(&nonrigid).contains("right-handed"));
}

#[test]
fn reads_contiguous_dataset_above_provider_allocation_ceiling() {
    const LARGE_DIMS: [usize; 3] = [257, 256, 128];
    const COUNT: usize = LARGE_DIMS[0] * LARGE_DIMS[1] * LARGE_DIMS[2];
    fn voxel(index: usize) -> u16 {
        u16::try_from(index % 65_521).expect("the fixture modulus fits in a 16-bit word")
    }
    let _trace = trace();
    let (dir, path) = {
        let _arrange = info_span!("arrange").entered();
        let dir = tempfile::tempdir().expect("create a directory for the large archive");
        let path = dir.path().join("large.h5");
        let le = ByteOrder::LittleEndian;
        fs::write(
            &path,
            archive(
                &LARGE_DIMS,
                (0..COUNT).map(|i| f64::from(voxel(i)).to_bits()),
                le,
                &GEOMETRY,
                le,
            ),
        )
        .expect("write the large contiguous fixture");
        (dir, path)
    };
    {
        let double = {
            let _act = info_span!("act", destination = "binary64").entered();
            load_volume_hdf5::<f64>(&path).expect("load the large archive into binary64")
        };
        let _assert = info_span!("assert").entered();
        assert_eq!(double.grid().dims(), LARGE_DIMS);
        assert_eq!(double.as_slice().len(), COUNT);
        assert!(double
            .as_slice()
            .iter()
            .enumerate()
            .all(|(index, &value)| value.to_bits() == f64::from(voxel(index)).to_bits()));
    }
    {
        let single = {
            let _act = info_span!("act", destination = "binary32").entered();
            load_volume_hdf5::<f32>(&path).expect("load the large archive into binary32")
        };
        let _assert = info_span!("assert").entered();
        assert_eq!(single.grid().dims(), LARGE_DIMS);
        assert_eq!(single.as_slice().len(), COUNT);
        assert!(single
            .as_slice()
            .iter()
            .enumerate()
            .all(|(index, &value)| value.to_bits() == f32::from(voxel(index)).to_bits()));
    }
    drop(dir);
}
