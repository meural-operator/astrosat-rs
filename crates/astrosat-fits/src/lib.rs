//! Pure-Rust zero-copy FITS 4.0 and OGIP science product engine for AstroSat.

pub mod bintable;
pub mod card;
pub mod checksum;
pub mod error;
pub mod image;
pub mod io;
pub mod ogip;

pub use bintable::{BinTableBuilder, BinTableRef, ColumnData, ColumnSpec, ColumnType, VarDescriptor};
pub use card::{Card, CardValue, Header};
pub use checksum::{
    accumulate_checksum, compute_hdu_checksums, datasum, decode_checksum, encode_checksum,
    verify_hdu_checksum, FITS_BLOCK_SIZE,
};
pub use error::{FitsError, Result};
pub use image::{Bitpix, ImageBuilder, ImageRef};
pub use io::{write_fits, FitsFile, FitsWriter, Hdu, HduKind, HduRef};
pub use ogip::{
    apply_ogip_meta, apply_wcs, build_gti_hdu, build_lightcurve_hdu, build_orbit_hdu,
    build_response_hdus, build_spectrum_hdu, validate_ogip, IssueSeverity, OgipMeta,
    ResponseMatrixData, SpectrumData, SpectrumMeta, ValidationIssue, WcsParams,
};

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn checksum_encode_decode_roundtrip() {
        for &v in &[0u32, 1, 42, 0x1234_5678, 0xDEAD_BEEF, 0xFFFF_FFFE, 0xFFFF_FFFF] {
            let enc = encode_checksum(v);
            assert_eq!(enc.len(), 16);
            assert!(enc.bytes().all(|b| b.is_ascii_alphanumeric()));
            let dec = decode_checksum(&enc).unwrap();
            let expected = if v == 0xFFFF_FFFF { 0 } else { v };
            assert_eq!(dec, expected, "failed roundtrip for {v:#X}, enc={enc}");
        }
    }

    #[test]
    fn continue_long_string_roundtrip() {
        let long_val = "/data/caldb/czti/data/astrosat/czti/bcf/AS1cztbadpix20160526v01.fits_extra_long_comment_path_testing_continue_cards";
        assert!(long_val.len() > 80);

        let mut primary = Hdu::empty_primary();
        primary
            .header
            .set_str("LONGPATH", long_val, "very long CALDB path");
        primary
            .header
            .set_str("QUOTED", "O'Brien's Test", "escaped quotes");

        let mut buf = Vec::new();
        let mut writer = FitsWriter::new(&mut buf);
        writer.write_hdu_with_checksum(&mut primary).unwrap();
        writer.finish().unwrap();

        assert_eq!(buf.len() % FITS_BLOCK_SIZE, 0);

        let fits = FitsFile::from_bytes(buf).unwrap();
        let h0 = fits.hdu(0).unwrap();
        assert_eq!(h0.header().get_str("LONGPATH"), Some(long_val));
        assert_eq!(h0.header().get_str("QUOTED"), Some("O'Brien's Test"));
        assert_eq!(h0.verify_checksum(), (true, true));
    }

    #[test]
    fn bintable_all_types_and_vectors_roundtrip() {
        let mut primary = Hdu::empty_primary();
        let vec_u8: Vec<u8> = (0..(3 * 2242)).map(|i| (i % 251) as u8).collect();
        let vec_i16: Vec<i16> = (0..(3 * 16)).map(|i| (i as i16) - 20).collect();
        let vec_f32: Vec<f32> = (0..(3 * 8)).map(|i| (i as f32) * 0.25).collect();

        let mut bt = BinTableBuilder::new("TEST_EVENTS")
            .add_f64("TIME", Some("s"), vec![100.5, 200.25, 300.75])
            .add_f32("ENERGY", Some("keV"), vec![12.5, 35.2, 102.8])
            .add_u16("PHA", Some("chan"), vec![0, 512, 65535])
            .add_i16("PI", Some("chan"), vec![-1, 100, 511])
            .add_u32("PKTNO", None, vec![0, 100_000, 4_000_000_000])
            .add_i64("BIGCNT", None, vec![-99, 0, 9_000_000_000])
            .add_u8("QUADID", None, vec![0, 1, 3])
            .add_bool("FLAG", vec![true, false, true])
            .add_strings("DETNAM", 8, vec!["Q0".into(), "LAXPC20".into(), "CZTI".into()])
            .add_vec_u8("FRAME", None, 2242, vec_u8.clone())
            .add_vec_i16("SPEC16", None, 16, vec_i16.clone())
            .add_vec_f32("WEIGHTS", None, 8, vec_f32.clone())
            .build()
            .unwrap();

        let mut buf = Vec::new();
        let mut w = FitsWriter::new(&mut buf);
        w.write_hdu_with_checksum(&mut primary).unwrap();
        w.write_hdu_with_checksum(&mut bt).unwrap();
        w.finish().unwrap();

        let fits = FitsFile::from_bytes(buf).unwrap();
        assert_eq!(fits.len(), 2);
        let h1 = fits.hdu_by_name("TEST_EVENTS").unwrap();
        assert_eq!(h1.verify_checksum(), (true, true));
        let tbl = h1.as_bintable().unwrap();
        assert_eq!(tbl.nrows(), 3);
        assert_eq!(tbl.read_col_f64("TIME").unwrap(), vec![100.5, 200.25, 300.75]);
        assert_eq!(tbl.read_col_u16("PHA").unwrap(), vec![0, 512, 65535]);
        assert_eq!(tbl.read_col_i16("PI").unwrap(), vec![-1, 100, 511]);
        assert_eq!(
            tbl.read_col_u32("PKTNO").unwrap(),
            vec![0, 100_000, 4_000_000_000]
        );
        assert_eq!(
            tbl.read_col_i64("BIGCNT").unwrap(),
            vec![-99, 0, 9_000_000_000]
        );
        assert_eq!(tbl.read_col_u8("QUADID").unwrap(), vec![0, 1, 3]);
        assert_eq!(tbl.read_col_bool("FLAG").unwrap(), vec![true, false, true]);
        assert_eq!(
            tbl.read_col_strings("DETNAM").unwrap(),
            vec!["Q0", "LAXPC20", "CZTI"]
        );

        let (rep_b, got_u8) = tbl.read_vec_u8("FRAME").unwrap();
        assert_eq!(rep_b, 2242);
        assert_eq!(got_u8, vec_u8);
        assert_eq!(tbl.cell_bytes_by_name(1, "FRAME").unwrap(), &vec_u8[2242..4484]);

        let (rep_i, got_i16) = tbl.read_vec_i16("SPEC16").unwrap();
        assert_eq!(rep_i, 16);
        assert_eq!(got_i16, vec_i16);

        let (rep_f, got_f32) = tbl.read_vec_f32("WEIGHTS").unwrap();
        assert_eq!(rep_f, 8);
        assert_eq!(got_f32, vec_f32);
    }

    #[test]
    fn image_hdu_and_wcs_roundtrip() {
        let pixels: Vec<f32> = (0..(16 * 12)).map(|i| (i as f32) * 1.5 - 10.0).collect();
        let wcs = WcsParams {
            crpix1: 8.5,
            crpix2: 6.5,
            crval1: 83.6331,
            crval2: 22.0145,
            ..Default::default()
        };
        let mut h0 = ImageBuilder::from_f32(None, vec![16, 12], &pixels)
            .unwrap()
            .with_wcs(&wcs)
            .build_primary();
        let mut h1 = ImageBuilder::from_i16(Some("DPH"), vec![4, 4], &(0..16).collect::<Vec<_>>())
            .unwrap()
            .with_bzero_bscale(32768.0, 1.0)
            .build_extension();

        let mut buf = Vec::new();
        let mut w = FitsWriter::new(&mut buf);
        w.write_hdu_with_checksum(&mut h0).unwrap();
        w.write_hdu_with_checksum(&mut h1).unwrap();
        w.finish().unwrap();

        let fits = FitsFile::from_bytes(buf).unwrap();
        assert_eq!(fits.hdu(0).unwrap().verify_checksum(), (true, true));
        assert_eq!(fits.hdu(1).unwrap().verify_checksum(), (true, true));

        let img0 = fits.hdu(0).unwrap().as_image().unwrap();
        assert_eq!(img0.shape(), &[16, 12]);
        assert_eq!(img0.read_f32().unwrap(), pixels);
        assert!((fits.hdu(0).unwrap().header().get_f64("CRVAL1").unwrap() - 83.6331).abs() < 1e-6);

        let img1 = fits.hdu_by_name("DPH").unwrap().as_image().unwrap();
        let vals = img1.read_f64().unwrap();
        assert_eq!(vals[0], 32768.0);
        assert_eq!(vals[15], 32783.0);
    }

    #[test]
    fn ogip_products_build_and_validate() {
        let meta = OgipMeta {
            instrume: "CZTI".into(),
            detnam: Some("Q0".into()),
            object: Some("CRAB".into()),
            ra_obj: Some(83.6331),
            dec_obj: Some(22.0145),
            tstart: 190_000_000.0,
            tstop: 190_001_000.0,
            exposure: 950.0,
            ontime: Some(1000.0),
            ..Default::default()
        };

        let mut primary = Hdu::empty_primary();
        let mut gti = build_gti_hdu("GTI", &[190_000_000.0], &[190_001_000.0], &meta).unwrap();
        let spec_meta = SpectrumMeta {
            detchans: 4,
            tlmin1: 0,
            tlmax1: 3,
            ..Default::default()
        };
        let mut spec = build_spectrum_hdu(
            "SPECTRUM",
            &[0, 1, 2, 3],
            SpectrumData::Counts(vec![10, 25, 18, 4]),
            Some(&[3.16, 5.0, 4.24, 2.0]),
            None,
            None,
            &meta,
            &spec_meta,
        )
        .unwrap();
        let mut lc = build_lightcurve_hdu(
            "RATE",
            &[190_000_005.0, 190_000_015.0],
            &[12.3, 11.8],
            &[1.1, 1.05],
            Some(&[1.0, 0.98]),
            10.0,
            &meta,
        )
        .unwrap();
        let rmf_data = ResponseMatrixData {
            chantype: "PI".into(),
            detchans: 4,
            tlmin4: 0,
            energ_lo: vec![20.0, 40.0],
            energ_hi: vec![40.0, 60.0],
            n_grp: vec![1, 1],
            f_chan: vec![vec![0], vec![2]],
            n_chan: vec![vec![2], vec![2]],
            matrix: vec![vec![0.6, 0.4], vec![0.3, 0.7]],
            channels: vec![0, 1, 2, 3],
            e_min: vec![20.0, 30.0, 40.0, 50.0],
            e_max: vec![30.0, 40.0, 50.0, 60.0],
            lo_thresh: 1e-6,
        };
        let (mut mat, mut eb) = build_response_hdus(&rmf_data, &meta).unwrap();

        let mut buf = Vec::new();
        let mut w = FitsWriter::new(&mut buf);
        for h in [&mut primary, &mut gti, &mut spec, &mut lc, &mut mat, &mut eb] {
            w.write_hdu_with_checksum(h).unwrap();
        }
        w.finish().unwrap();

        let fits = FitsFile::from_bytes(buf).unwrap();
        let issues = validate_ogip(&fits);
        assert!(issues.is_empty(), "unexpected OGIP issues: {issues:?}");

        // Verify variable-length RMF matrix readback.
        let mat_tbl = fits.hdu_by_name("MATRIX").unwrap().as_bintable().unwrap();
        assert_eq!(mat_tbl.read_var_i16("F_CHAN", 0).unwrap(), vec![0]);
        assert_eq!(mat_tbl.read_var_i16("F_CHAN", 1).unwrap(), vec![2]);
        assert_eq!(mat_tbl.read_var_f32("MATRIX", 0).unwrap(), vec![0.6, 0.4]);
        assert_eq!(mat_tbl.read_var_f32("MATRIX", 1).unwrap(), vec![0.3, 0.7]);
    }

    #[test]
    fn nasa_fitsverify_zero_errors_zero_warnings() {
        let fitsverify = std::path::Path::new("/mnt/d/ISRO/.toolchain/bin/fitsverify");
        if !fitsverify.exists() {
            return;
        }
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("verify_all.fits");

        let meta = OgipMeta {
            instrume: "CZTI".into(),
            detnam: Some("Q0".into()),
            object: Some("CRAB".into()),
            ra_obj: Some(83.6331),
            dec_obj: Some(22.0145),
            tstart: 190_000_000.0,
            tstop: 190_001_000.0,
            exposure: 950.0,
            ontime: Some(1000.0),
            ..Default::default()
        };

        let mut p0 = Hdu::empty_primary();
        p0.header.set_f64("TSTART", 197_325_364.264_079_2, "exact sub-us start time");
        p0.header.set_f64("BIGSCI", 1.234_567_890_123_456e15, "16-sig-fig scientific float");
        p0.header.set_str(
            "LONGSTR",
            "a_very_long_string_value_exceeding_eighty_bytes_to_verify_fits_4_0_continue_card_support_cleanly",
            "long string",
        );
        let pixels: Vec<f32> = (0..64).map(|i| i as f32).collect();
        let img = ImageBuilder::from_f32(Some("SKYIMG"), vec![8, 8], &pixels)
            .unwrap()
            .with_wcs(&WcsParams::default())
            .build_extension();
        let gti = build_gti_hdu("GTI", &[190_000_000.0], &[190_001_000.0], &meta).unwrap();
        let spec_meta = SpectrumMeta {
            detchans: 4,
            tlmin1: 0,
            tlmax1: 3,
            ..Default::default()
        };
        let spec = build_spectrum_hdu(
            "SPECTRUM",
            &[0, 1, 2, 3],
            SpectrumData::Counts(vec![10, 25, 18, 4]),
            Some(&[3.16, 5.0, 4.24, 2.0]),
            None,
            None,
            &meta,
            &spec_meta,
        )
        .unwrap();
        let lc = build_lightcurve_hdu(
            "RATE",
            &[190_000_005.0, 190_000_015.0],
            &[12.3, 11.8],
            &[1.1, 1.05],
            Some(&[1.0, 0.98]),
            10.0,
            &meta,
        )
        .unwrap();
        let orb = build_orbit_hdu(
            &[190_000_000.0, 190_000_010.0],
            &[[7000.0, 0.0, 0.0], [6999.0, 75.0, 0.0]],
            &[[0.0, 7.5, 0.0], [-0.08, 7.5, 0.0]],
            &meta,
        )
        .unwrap();

        let rmf_data = ResponseMatrixData {
            chantype: "PI".into(),
            detchans: 4,
            tlmin4: 0,
            energ_lo: vec![20.0, 40.0],
            energ_hi: vec![40.0, 60.0],
            n_grp: vec![1, 1],
            f_chan: vec![vec![0], vec![2]],
            n_chan: vec![vec![2], vec![2]],
            matrix: vec![vec![0.6, 0.4], vec![0.3, 0.7]],
            channels: vec![0, 1, 2, 3],
            e_min: vec![20.0, 30.0, 40.0, 50.0],
            e_max: vec![30.0, 40.0, 50.0, 60.0],
            lo_thresh: 1e-6,
        };
        let (mat, eb) = build_response_hdus(&rmf_data, &meta).unwrap();

        let mut hdus = vec![p0, img, gti, spec, lc, orb, mat, eb];
        write_fits(&path, &mut hdus, true).unwrap();

        // Verify 16-significant-digit f64 header precision roundtrip (< 0.1 us on ~2e8 s timestamps).
        let reopened = FitsFile::open(&path).unwrap();
        let got_ts = reopened.hdu(0).unwrap().header().get_f64("TSTART").unwrap();
        assert!(
            (got_ts - 197_325_364.264_079_2).abs() < 1e-7,
            "TSTART header precision lost: got {got_ts}"
        );
        let got_sci = reopened.hdu(0).unwrap().header().get_f64("BIGSCI").unwrap();
        assert!(
            ((got_sci - 1.234_567_890_123_456e15) / 1.234_567_890_123_456e15).abs() < 1e-14,
            "BIGSCI header precision lost: got {got_sci}"
        );

        let out = std::process::Command::new(fitsverify)
            .arg(&path)
            .output()
            .expect("failed to execute fitsverify");
        let stdout = String::from_utf8_lossy(&out.stdout);
        let stderr = String::from_utf8_lossy(&out.stderr);
        assert!(
            out.status.success() && stdout.contains("0 warning(s) and 0 error(s)"),
            "fitsverify reported issues:\nSTDOUT:\n{stdout}\nSTDERR:\n{stderr}"
        );
    }

    #[test]
    fn malformed_fits_headers_and_heap_rejected_without_panic() {
        // 1. Invalid BITPIX = 9 must fail cleanly.
        let mut h_bad_bp = Header::new();
        h_bad_bp.set_bool("SIMPLE", true, "");
        h_bad_bp.set_i64("BITPIX", 9, "");
        h_bad_bp.set_i64("NAXIS", 0, "");
        assert!(FitsFile::from_bytes(h_bad_bp.to_fits_bytes()).is_err());

        // 2. Negative NAXIS must fail cleanly (no Vec::with_capacity panic).
        let mut h_neg_naxis = Header::new();
        h_neg_naxis.set_bool("SIMPLE", true, "");
        h_neg_naxis.set_i64("BITPIX", 8, "");
        h_neg_naxis.set_i64("NAXIS", -5, "");
        assert!(FitsFile::from_bytes(h_neg_naxis.to_fits_bytes()).is_err());

        // 3. Overflowing NAXIS1 * NAXIS2 must fail cleanly without wrapping.
        let mut h_ovf = Header::new();
        h_ovf.set_bool("SIMPLE", true, "");
        h_ovf.set_i64("BITPIX", 8, "");
        h_ovf.set_i64("NAXIS", 2, "");
        h_ovf.set_i64("NAXIS1", i64::MAX / 2, "");
        h_ovf.set_i64("NAXIS2", 10, "");
        assert!(FitsFile::from_bytes(h_ovf.to_fits_bytes()).is_err());

        // 4. Negative TFIELDS / NAXIS1 / PCOUNT in BINTABLE must return Err, never abort.
        let p0 = Hdu::empty_primary();
        let mut h_bad_bt = Header::new();
        h_bad_bt.set_str("XTENSION", "BINTABLE", "");
        h_bad_bt.set_i64("BITPIX", 8, "");
        h_bad_bt.set_i64("NAXIS", 2, "");
        h_bad_bt.set_i64("NAXIS1", 8, "");
        h_bad_bt.set_i64("NAXIS2", 0, "");
        h_bad_bt.set_i64("PCOUNT", 0, "");
        h_bad_bt.set_i64("GCOUNT", 1, "");
        h_bad_bt.set_i64("TFIELDS", -1, "");
        let mut buf = p0.header.to_fits_bytes();
        buf.extend_from_slice(&h_bad_bt.to_fits_bytes());
        let fits = FitsFile::from_bytes(buf).unwrap();
        assert!(fits.hdu(1).unwrap().as_bintable().is_err());

        // 5. Out-of-range float in get_i64 must return None (not saturate to i64::MAX).
        let mut hdr = Header::new();
        hdr.set_f64("HUGE", 1.0e30, "");
        assert_eq!(hdr.get_i64("HUGE"), None);

        // 6. Malformed P32 heap descriptor with u32::MAX count/offset must return Err.
        let mut bt = BinTableBuilder::new("VAR")
            .add_var_f32("ARR", None, vec![vec![1.0, 2.0]])
            .build()
            .unwrap();
        // Corrupt the descriptor in main table row 0 to (u32::MAX, u32::MAX).
        bt.data[0..8].copy_from_slice(&[0xFF; 8]);
        let mut buf2 = Vec::new();
        let mut w = FitsWriter::new(&mut buf2);
        w.write_hdu(&p0).unwrap();
        w.write_hdu(&bt).unwrap();
        w.finish().unwrap();
        let fits2 = FitsFile::from_bytes(buf2).unwrap();
        let tbl2 = fits2.hdu(1).unwrap().as_bintable().unwrap();
        assert!(tbl2.read_var_f32("ARR", 0).is_err());

        // 7. Extreme TZERO values (e.g. 1e30 or 2^63) must never panic with integer overflow in debug mode.
        let mut bt_tz = BinTableBuilder::new("TZ")
            .add_i32("C32", None, vec![10])
            .add_i16("C16", None, vec![10])
            .add_u8("C8", None, vec![10])
            .add_vec_i16("V16", None, 2, vec![10, 20])
            .build()
            .unwrap();
        bt_tz.header.set_f64("TZERO1", 1.0e30, "");
        bt_tz.header.set_f64("TZERO2", 9_223_372_036_854_775_808.0, "");
        bt_tz.header.set_f64("TZERO3", 9_223_372_036_854_775_808.0, "");
        bt_tz.header.set_f64("TZERO4", 1.0e30, "");
        let mut buf3 = Vec::new();
        let mut w3 = FitsWriter::new(&mut buf3);
        w3.write_hdu(&p0).unwrap();
        w3.write_hdu(&bt_tz).unwrap();
        w3.finish().unwrap();
        let fits3 = FitsFile::from_bytes(buf3).unwrap();
        let tbl3 = fits3.hdu(1).unwrap().as_bintable().unwrap();
        let _ = tbl3.read_col_i64("C32").unwrap();
        let _ = tbl3.read_col_i64("C16").unwrap();
        let _ = tbl3.read_col_i64("C8").unwrap();
        let _ = tbl3.read_vec_u16("V16").unwrap();
    }
}
