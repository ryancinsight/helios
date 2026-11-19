use super::{load_ct_series, load_ct_slice};
use crate::volume::Volume;
use helios_quantities::HeliosError;

mod fixtures;
mod native_writer;

use fixtures::{
    BitsStoredFixture, OmittedGeometry, write_slice_at, write_slice_at_with_geometry,
    write_slice_at_with_orientation,
};

#[test]
fn round_trips_a_synthetic_ct_slice_to_hu_volume() {
    let dir = tempfile::tempdir().expect("test scratch dir is creatable");
    let path = dir.path().join("slice.dcm");
    write_slice_at(&path, [10, 20, 30, 40], 9.0, "2.25.4242");

    let vol: Volume<f64> = load_ct_slice(&path).expect("load");
    let grid = vol.grid();
    assert_eq!(grid.dims(), [2, 2, 1]);
    // spacing = [col, row, thickness] = [1.25, 0.8, 3.0].
    let sp = grid.spacing();
    assert!(
        (sp[0] - 1.25).abs() < 1e-12 && (sp[1] - 0.8).abs() < 1e-12 && (sp[2] - 3.0).abs() < 1e-12
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
fn twelve_bit_stored_samples_ignore_container_padding() {
    let dir = tempfile::tempdir().expect("test scratch dir is creatable");
    let path = dir.path().join("twelve-bit.dcm");
    write_slice_at_with_geometry(
        &path,
        [0xF000, 0xA001, 0xB800, 0xFFFF],
        9.0,
        "2.25.4244",
        "1\\0\\0\\0\\1\\0",
        OmittedGeometry::None,
        BitsStoredFixture::Value(12),
    );

    let volume: Volume<f64> = load_ct_slice(&path).expect("load twelve-bit slice");
    assert_eq!(volume.grid().dims(), [2, 2, 1]);
    assert_eq!(volume.get(0, 0, 0), Some(-10.0));
    assert_eq!(volume.get(1, 0, 0), Some(-8.0));
    assert_eq!(volume.get(0, 1, 0), Some(4086.0));
    assert_eq!(volume.get(1, 1, 0), Some(8180.0));
}

#[test]
fn missing_file_is_a_dicom_error_not_a_panic() {
    let err = load_ct_slice::<f64>("does_not_exist.dcm")
        .expect_err("a missing file is a Dicom error, not a panic");
    assert!(matches!(err, HeliosError::Dicom { .. }));
}

#[test]
fn missing_required_geometry_is_rejected_at_dicom_boundary() {
    let cases = [
        (OmittedGeometry::PixelSpacing, "PixelSpacing"),
        (OmittedGeometry::ImagePosition, "ImagePositionPatient"),
        (
            OmittedGeometry::ImageOrientation,
            "ImageOrientationPatient (0020,0037)",
        ),
    ];

    for (index, (omitted, name)) in cases.into_iter().enumerate() {
        let dir = tempfile::tempdir().expect("test scratch dir is creatable");
        let path = dir.path().join("missing-geometry.dcm");
        write_slice_at_with_geometry(
            &path,
            [10, 20, 30, 40],
            9.0,
            &format!("2.25.4242.{index}"),
            "1\\0\\0\\0\\1\\0",
            omitted,
            BitsStoredFixture::Value(16),
        );

        match load_ct_slice::<f64>(&path) {
            Err(HeliosError::Dicom { reason }) => {
                assert_eq!(reason, format!("{name}: missing required attribute"));
            }
            Err(other) => panic!("unexpected DICOM error: {other}"),
            Ok(_) => panic!("missing {name} was accepted"),
        }
    }
}

#[test]
fn missing_and_malformed_bits_stored_are_rejected_at_dicom_boundary() {
    let cases = [
        (BitsStoredFixture::Missing, "BitsStored: missing BitsStored"),
        (
            BitsStoredFixture::Value(17),
            "decode: bits_stored=17 is outside 1..=bits_allocated (16)",
        ),
        (
            BitsStoredFixture::Values([12, 16]),
            "BitsStored: BitsStored (0028,0101) has cardinality=2; expected exactly 1",
        ),
    ];

    for (index, (case, expected)) in cases.into_iter().enumerate() {
        let dir = tempfile::tempdir().expect("test scratch dir is creatable");
        let path = dir.path().join("invalid-bits-stored.dcm");
        write_slice_at_with_geometry(
            &path,
            [10, 20, 30, 40],
            9.0,
            &format!("2.25.4243.{index}"),
            "1\\0\\0\\0\\1\\0",
            OmittedGeometry::None,
            case,
        );

        match load_ct_slice::<f64>(&path) {
            Err(HeliosError::Dicom { reason }) => assert_eq!(reason, expected),
            Err(other) => panic!("unexpected DICOM error: {other}"),
            Ok(_) => panic!("missing or malformed BitsStored metadata was accepted"),
        }
    }
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
    let dir = tempfile::tempdir().expect("test scratch dir is creatable");
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
    let dir = tempfile::tempdir().expect("test scratch dir is creatable");
    let path = dir.path().join("one.dcm");
    write_slice_at(&path, [10, 20, 30, 40], 9.0, "2.25.9");
    let series: Volume<f64> =
        load_ct_series(std::slice::from_ref(&path)).expect("synthetic series round-trips");
    let single: Volume<f64> = load_ct_slice(&path).expect("synthetic DICOM fixture round-trips");
    assert_eq!(series.grid().dims(), single.grid().dims());
    for j in 0..2 {
        for i in 0..2 {
            assert_eq!(series.get(i, j, 0), single.get(i, j, 0));
        }
    }
}

#[test]
fn load_slice_preserves_oriented_iop_pose() {
    let dir = tempfile::tempdir().expect("test scratch dir is creatable");
    let path = dir.path().join("oriented.dcm");
    // row_dir = +Y, col_dir = -X, normal = +Z (right-handed).
    write_slice_at_with_orientation(
        &path,
        [10, 20, 30, 40],
        9.0,
        "2.25.4242.42",
        "0\\1\\0\\-1\\0\\0",
    );

    let vol: Volume<f64> = load_ct_slice(&path).expect("load");
    let g = vol.grid();
    // Voxel (1,0,0): origin + 1*col_spacing*row_dir = (5, 7+1.25, 9).
    let p_i = g.voxel_center(1, 0, 0);
    assert!((p_i.x - 5.0).abs() < 1e-12);
    assert!((p_i.y - 8.25).abs() < 1e-12);
    assert!((p_i.z - 9.0).abs() < 1e-12);
    // Voxel (0,1,0): origin + 1*row_spacing*col_dir = (5-0.8, 7, 9).
    let p_j = g.voxel_center(0, 1, 0);
    assert!((p_j.x - 4.2).abs() < 1e-12);
    assert!((p_j.y - 7.0).abs() < 1e-12);
    assert!((p_j.z - 9.0).abs() < 1e-12);
}

#[test]
fn empty_and_non_uniform_series_error() {
    let empty: &[std::path::PathBuf] = &[];
    assert!(matches!(
        load_ct_series::<f64, _>(empty),
        Err(HeliosError::Dicom { .. })
    ));

    // Slices at z = 0, 4, 10 → gap (missing slice) → non-uniform spacing.
    let dir = tempfile::tempdir().expect("test scratch dir is creatable");
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
