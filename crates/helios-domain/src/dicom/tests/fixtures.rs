use ritk_dicom::{DicomTag, tags};

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
pub(super) enum OmittedGeometry {
    None,
    PixelSpacing,
    ImagePosition,
    ImageOrientation,
}

#[derive(Clone, Copy)]
pub(super) enum BitsStoredFixture {
    Missing,
    Value(u16),
    Values([u16; 2]),
}

// Write a synthetic 2×2 unsigned-16 CT slice at position `z_mm` with a known
// HU pattern and geometry (no external fixture). Slope 2, intercept −10;
// PixelSpacing 0.8 (row) / 1.25 (col) mm; in-plane origin (5,7); configurable
// `ImageOrientationPatient`; a unique SOP instance UID per file.
pub(super) fn write_slice_at_with_geometry(
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
            tags::IMAGE_ORIENTATION_PATIENT,
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

pub(super) fn write_slice_at_with_orientation(
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

pub(super) fn write_slice_at(path: &std::path::Path, pixels: [u16; 4], z_mm: f64, uid: &str) {
    write_slice_at_with_orientation(path, pixels, z_mm, uid, "1\\0\\0\\0\\1\\0");
}
