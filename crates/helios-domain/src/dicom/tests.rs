    use super::*;
    use ritk_dicom::{tags, DicomTag};

    fn append_element(output: &mut Vec<u8>, tag: DicomTag, vr: [u8; 2], value: &[u8]) {
        output.extend_from_slice(&tag.group.to_le_bytes());
        output.extend_from_slice(&tag.element.to_le_bytes());
        output.extend_from_slice(&vr);
        if matches!(
            vr,
            [b'O', b'B' | b'W' | b'F'] | [b'S', b'Q'] | [b'U', b'T' | b'N']
        ) {
            output.extend_from_slice(&[0, 0]);
            output.extend_from_slice(
                &u32::try_from(value.len())
                    .expect("synthetic DICOM element length fits u32")
                    .to_le_bytes(),
            );
        } else {
            output.extend_from_slice(
                &u16::try_from(value.len())
                    .expect("synthetic DICOM element length fits u16")
                    .to_le_bytes(),
            );
        }
        output.extend_from_slice(value);
    }

    fn text_value(vr: [u8; 2], value: &str) -> Vec<u8> {
        let mut bytes = value.as_bytes().to_vec();
        if !bytes.len().is_multiple_of(2) {
            bytes.push(if vr == *b"UI" { 0 } else { b' ' });
        }
        bytes
    }

    fn unsigned_value(value: u16) -> [u8; 2] {
        value.to_le_bytes()
    }

    #[derive(Clone, Copy, PartialEq, Eq)]
    enum OmittedGeometry {
        None,
        PixelSpacing,
        ImagePosition,
        ImageOrientation,
    }

    #[derive(Clone, Copy)]
    enum BitsStoredFixture {
        Missing,
        Value(u16),
        Values([u16; 2]),
    }

    // Write a synthetic 2×2 unsigned-16 CT slice at position `z_mm` with a known
    // HU pattern and geometry (no external fixture). Slope 2, intercept −10;
    // PixelSpacing 0.8 (row) / 1.25 (col) mm; in-plane origin (5,7); configurable
    // `ImageOrientationPatient`; a unique SOP instance UID per file.
    fn write_slice_at_with_geometry(
        path: &std::path::Path,
        pixels: [u16; 4],
        z_mm: f64,
        uid: &str,
        image_orientation_patient: &str,
        omitted_geometry: OmittedGeometry,
        bits_stored: BitsStoredFixture,
    ) {
        let sop_class = "1.2.840.10008.5.1.4.1.1.2";
        let transfer_syntax = "1.2.840.10008.1.2.1";

        // This test-only fixture uses the explicit-VR little-endian Part 10
        // wire format. Parsing and pixel decoding remain exclusively owned by
        // ritk-dicom; no consumer-side DICOM object model is constructed.
        let mut meta_body = Vec::new();
        append_element(
            &mut meta_body,
            DicomTag::new(0x0002, 0x0001),
            *b"OB",
            &[0, 1],
        );
        append_element(
            &mut meta_body,
            DicomTag::new(0x0002, 0x0002),
            *b"UI",
            &text_value(*b"UI", sop_class),
        );
        append_element(
            &mut meta_body,
            DicomTag::new(0x0002, 0x0003),
            *b"UI",
            &text_value(*b"UI", uid),
        );
        append_element(
            &mut meta_body,
            DicomTag::new(0x0002, 0x0010),
            *b"UI",
            &text_value(*b"UI", transfer_syntax),
        );
        append_element(
            &mut meta_body,
            DicomTag::new(0x0002, 0x0012),
            *b"UI",
            &text_value(*b"UI", "2.25.4242.1"),
        );

        let mut bytes = vec![0; 128];
        bytes.extend_from_slice(b"DICM");
        append_element(
            &mut bytes,
            DicomTag::new(0x0002, 0x0000),
            *b"UL",
            &u32::try_from(meta_body.len())
                .expect("synthetic DICOM meta length fits u32")
                .to_le_bytes(),
        );
        bytes.extend_from_slice(&meta_body);

        append_element(
            &mut bytes,
            DicomTag::new(0x0008, 0x0016),
            *b"UI",
            &text_value(*b"UI", sop_class),
        );
        append_element(
            &mut bytes,
            DicomTag::new(0x0008, 0x0018),
            *b"UI",
            &text_value(*b"UI", uid),
        );
        append_element(
            &mut bytes,
            DicomTag::new(0x0008, 0x0008),
            *b"CS",
            &text_value(*b"CS", "ORIGINAL\\PRIMARY\\AXIAL"),
        );
        append_element(&mut bytes, tags::ROWS, *b"US", &unsigned_value(2));
        append_element(&mut bytes, tags::COLUMNS, *b"US", &unsigned_value(2));
        append_element(
            &mut bytes,
            tags::SAMPLES_PER_PIXEL,
            *b"US",
            &unsigned_value(1),
        );
        append_element(
            &mut bytes,
            tags::BITS_ALLOCATED,
            *b"US",
            &unsigned_value(16),
        );
        let bits_stored_value = match bits_stored {
            BitsStoredFixture::Missing => 16,
            BitsStoredFixture::Value(value) => {
                append_element(
                    &mut bytes,
                    tags::BITS_STORED,
                    *b"US",
                    &unsigned_value(value),
                );
                value
            }
            BitsStoredFixture::Values([first, second]) => {
                let [first_low, first_high] = first.to_le_bytes();
                let [second_low, second_high] = second.to_le_bytes();
                append_element(
                    &mut bytes,
                    tags::BITS_STORED,
                    *b"US",
                    &[first_low, first_high, second_low, second_high],
                );
                first
            }
        };
        append_element(
            &mut bytes,
            DicomTag::new(0x0028, 0x0102),
            *b"US",
            &unsigned_value(bits_stored_value.saturating_sub(1)),
        );
        append_element(
            &mut bytes,
            tags::PIXEL_REPRESENTATION,
            *b"US",
            &unsigned_value(0),
        );
        append_element(
            &mut bytes,
            DicomTag::new(0x0028, 0x0004),
            *b"CS",
            &text_value(*b"CS", "MONOCHROME2"),
        );
        append_element(
            &mut bytes,
            tags::RESCALE_SLOPE,
            *b"DS",
            &text_value(*b"DS", "2"),
        );
        append_element(
            &mut bytes,
            tags::RESCALE_INTERCEPT,
            *b"DS",
            &text_value(*b"DS", "-10"),
        );
        if omitted_geometry != OmittedGeometry::PixelSpacing {
            append_element(
                &mut bytes,
                tags::PIXEL_SPACING,
                *b"DS",
                &text_value(*b"DS", "0.8\\1.25"),
            );
        }
        append_element(
            &mut bytes,
            tags::SLICE_THICKNESS,
            *b"DS",
            &text_value(*b"DS", "3"),
        );
        if omitted_geometry != OmittedGeometry::ImagePosition {
            append_element(
                &mut bytes,
                tags::IMAGE_POSITION_PATIENT,
                *b"DS",
                &text_value(*b"DS", &format!("5\\7\\{z_mm}")),
            );
        }
        if omitted_geometry != OmittedGeometry::ImageOrientation {
            append_element(
                &mut bytes,
                IMAGE_ORIENTATION_PATIENT,
                *b"DS",
                &text_value(*b"DS", image_orientation_patient),
            );
        }

        let mut pixel_bytes = Vec::with_capacity(pixels.len() * 2);
        for pixel in pixels {
            pixel_bytes.extend_from_slice(&pixel.to_le_bytes());
        }
        append_element(
            &mut bytes,
            DicomTag::new(0x7FE0, 0x0010),
            *b"OW",
            &pixel_bytes,
        );
        std::fs::write(path, bytes).expect("write synthetic DICOM");
    }

    fn write_slice_at_with_orientation(
        path: &std::path::Path,
        pixels: [u16; 4],
        z_mm: f64,
        uid: &str,
        image_orientation_patient: &str,
    ) {
        write_slice_at_with_geometry(
            path,
            pixels,
            z_mm,
            uid,
            image_orientation_patient,
            OmittedGeometry::None,
            BitsStoredFixture::Value(16),
        );
    }

    fn write_slice_at(path: &std::path::Path, pixels: [u16; 4], z_mm: f64, uid: &str) {
        write_slice_at_with_orientation(path, pixels, z_mm, uid, "1\\0\\0\\0\\1\\0");
    }

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
        let single: Volume<f64> =
            load_ct_slice(&path).expect("synthetic DICOM fixture round-trips");
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
