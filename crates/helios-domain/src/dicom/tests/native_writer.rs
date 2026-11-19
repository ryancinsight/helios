use crate::volume::Volume;
use ritk_dicom::TransferSyntaxKind;
use ritk_dicom::dataset::DicomFile;
use tracing_subscriber::fmt::format::FmtSpan;

use super::{load_ct_slice, write_slice_at};

#[test]
fn native_writer_deflate_preserves_helios_geometry_and_calibration() {
    let _trace = tracing::subscriber::set_default(
        tracing_subscriber::fmt()
            .with_span_events(FmtSpan::FULL)
            .with_test_writer()
            .finish(),
    );
    let (_scratch, plain, deflated) = tracing::info_span!("arrange").in_scope(|| {
        let scratch = tempfile::tempdir().expect("test scratch");
        let plain = scratch.path().join("plain.dcm");
        let deflated = scratch.path().join("deflated.dcm");
        write_slice_at(&plain, [0, 10, 100, 1023], 9.0, "2.25.4242");
        (scratch, plain, deflated)
    });
    let (expected, actual): (Volume<f64>, Volume<f64>) =
        tracing::info_span!("act").in_scope(|| {
            let mut file = DicomFile::open(&plain).expect("independently encoded fixture loads");
            file.update_meta(|meta| {
                meta.transfer_syntax = TransferSyntaxKind::DeflatedExplicitVrLittleEndian
                    .uid()
                    .to_owned();
            });
            file.write_to_file(&deflated)
                .expect("native deflate writer");
            (
                load_ct_slice(&plain).expect("load plain CT"),
                load_ct_slice(&deflated).expect("load deflated CT"),
            )
        });
    tracing::info_span!("assert").in_scope(|| {
        assert_eq!(actual.grid().dims(), [2, 2, 1]);
        assert_eq!(actual.grid().dims(), expected.grid().dims());
        assert_eq!(actual.grid().spacing(), expected.grid().spacing());
        assert_eq!(
            actual.grid().origin(),
            helios_math::Point3::new(5.0, 7.0, 9.0)
        );
        assert_eq!(actual.grid().pose(), expected.grid().pose());
        for (col, row, hu) in [(0, 0, -10.0), (1, 0, 10.0), (0, 1, 190.0), (1, 1, 2036.0)] {
            assert_eq!(actual.get(col, row, 0), Some(hu));
            assert_eq!(expected.get(col, row, 0), Some(hu));
        }
    });
}
