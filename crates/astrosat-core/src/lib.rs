//! Shared AstroSat timing, TCT, GTI, MKF, spacecraft orbit, Solar System Barycentric correction,
//! and CALDB CIF query engine.

pub mod bary;
pub mod caldb;
pub mod error;
pub mod gti;
pub mod mkf;
pub mod orbit;
pub mod tct;
pub mod time;
pub mod toml_lite;

pub use bary::{
    barycenter_fits, fk5_to_icrs, source_unit_vector_icrs, BaryCalculator, BaryCorrection,
    BaryEphemeris, BaryFitsOptions, BaryFitsSummary, CoordFrame, EphemerisKind, JplDeEphemeris,
    SolarSystemState, AU_M, C_LIGHT, GM_SUN_OVER_C3,
};
pub use caldb::{CaldbEntry, CaldbIndex};
pub use error::{CoreError, Result};
pub use gti::{GtiInterval, GtiSet};
pub use mkf::MkfTable;
pub use orbit::{OrbitState, OrbitTable};
pub use tct::{
    CztiTct, LaxpcTct, TctTable, CZTI_CLOCK_WRAP_SEC, CZTI_DEFAULT_FEB_LAG, LAXPC_TCT_RATIO_MAX,
    LAXPC_TCT_RATIO_MIN,
};
pub use time::{
    astrosat_utc_to_iso8601, astrosat_utc_to_jd_tt, astrosat_utc_to_mjd_tdb,
    astrosat_utc_to_mjd_tt, astrosat_utc_to_mjd_utc, astrosat_utc_to_tt_sec,
    leap_seconds_since_epoch, mjd_utc_to_astrosat_utc, parse_date_to_mjd, tai_minus_utc_astrosat,
    tai_minus_utc_mjd, tdb_minus_tt, tt_sec_to_astrosat_utc, ASTROSAT_MJDREFF, ASTROSAT_MJDREFI,
    JD_J2000, JD_MINUS_MJD, LEAP_SECONDS_TABLE, SECS_PER_DAY, TT_MINUS_TAI,
};
pub use toml_lite::{
    parse_toml_to_json, to_canonical_toml, toml_quote, AppConfig, BaryFileConfig, CztiFileConfig,
    LaxpcFileConfig, TimingFileConfig, UiFileConfig,
};

#[cfg(test)]
mod tests {
    use super::*;
    use astrosat_fits::{build_orbit_hdu, write_fits, BinTableBuilder, FitsFile, Hdu, OgipMeta};
    use tempfile::tempdir;

    #[test]
    fn leap_second_boundaries_and_roundtrip() {
        // AstroSat epoch: 2010-01-01 (MJD 55197) -> TAI - UTC = 34 s.
        assert_eq!(tai_minus_utc_mjd(55197.0), 34.0);
        // Just before & at 2012-07-01 (MJD 56109): 34 -> 35 s.
        assert_eq!(tai_minus_utc_mjd(56108.9999), 34.0);
        assert_eq!(tai_minus_utc_mjd(56109.0), 35.0);
        // Just before & at 2015-07-01 (MJD 57204): 35 -> 36 s.
        assert_eq!(tai_minus_utc_mjd(57203.9999), 35.0);
        assert_eq!(tai_minus_utc_mjd(57204.0), 36.0);
        // Just before & at 2017-01-01 (MJD 57754): 36 -> 37 s.
        assert_eq!(tai_minus_utc_mjd(57753.9999), 36.0);
        assert_eq!(tai_minus_utc_mjd(57754.0), 37.0);

        let t_2016 = mjd_utc_to_astrosat_utc(57394.0); // 2016-01-07 (Crab ObsID 9000000252)
        assert_eq!(tai_minus_utc_astrosat(t_2016), 36.0);
        assert_eq!(leap_seconds_since_epoch(t_2016), 2.0);
        assert_eq!(astrosat_utc_to_iso8601(t_2016), "2016-01-07T00:00:00.000");

        let t_tt = astrosat_utc_to_tt_sec(t_2016);
        assert!((t_tt - (t_2016 + 36.0 + 32.184)).abs() < 1e-9);
        assert!((tt_sec_to_astrosat_utc(t_tt) - t_2016).abs() < 1e-9);

        // Fairhead & Bretagnon TDB-TT amplitude is ~1.657 ms with annual period.
        let (jd_hi, jd_lo) = astrosat_utc_to_jd_tt(t_2016);
        let einstein = tdb_minus_tt(jd_hi, jd_lo);
        assert!(
            einstein.abs() < 2.0e-3,
            "unexpected TDB-TT amplitude: {einstein}"
        );
        // Check peak amplitude across a year (~1.66 ms).
        let max_einstein = (0..365)
            .map(|d| tdb_minus_tt(2_457_394.5 + d as f64, 0.0).abs())
            .fold(0.0f64, f64::max);
        assert!(
            (max_einstein - 1.66e-3).abs() < 1.0e-4,
            "annual max TDB-TT should be ~1.66 ms, got {max_einstein}"
        );
    }

    #[test]
    fn czti_and_laxpc_tct_interpolation() {
        // CZTI rollover across 2^20 = 1048576 s with -2.0 s FEB lag.
        let sps = vec![200_000_000.0, 200_000_100.0, 200_000_200.0, 200_000_300.0];
        let inst = vec![1_048_500.0, 24.0, 124.0, 224.0]; // wraps past 1048576 between idx 0 and 1!
        let czti = CztiTct::new(TctTable::from_pairs(sps, inst).unwrap()).unwrap();
        assert_eq!(czti.num_segments(), 1);
        // cztseccnt = 1048552.0 -> minus 2.0s lag = 1048550.0 -> 50s after idx 0 -> 200_000_050.0
        let t1 = czti.interpolate(1_048_552.0, 2.0, None).unwrap();
        assert!((t1 - 200_000_050.0).abs() < 1e-6, "got {t1}");
        // cztseccnt = 76.0 -> minus 2.0s lag = 74.0 -> 50s after idx 1 -> 200_000_150.0
        let t2 = czti.interpolate(76.0, 2.0, None).unwrap();
        assert!((t2 - 200_000_150.0).abs() < 1e-6, "got {t2}");

        // LAXPC glitch rejection: true clock ratio = 1.0005, inject a bad glitch row.
        let l_inst = vec![1000.0, 1100.0, 1150.0, 1200.0, 1300.0];
        let l_sps = vec![
            190_000_000.0,
            190_000_100.05,
            190_000_500.0, // GLITCH: dt_sps/dt_ins = 8.0!
            190_000_200.10,
            190_000_300.15,
        ];
        let laxpc = LaxpcTct::new(TctTable::from_pairs(l_sps, l_inst).unwrap()).unwrap();
        assert_eq!(laxpc.rejected_glitch_count(), 1);
        assert_eq!(laxpc.len(), 4);
        let t_mid = laxpc.interpolate(1150.0);
        assert!((t_mid - 190_000_150.075).abs() < 1e-6, "got {t_mid}");
        let t_nev = laxpc.interpolate_neville(1150.0, 4);
        assert!((t_nev - 190_000_150.075).abs() < 1e-6, "got {t_nev}");
    }

    #[test]
    fn gti_algebra_and_filtering() {
        let a = GtiSet::from_pairs(&[(20.0, 30.0), (0.0, 10.0), (9.5, 15.0), (50.0, 40.0)]);
        assert_eq!(a.len(), 2);
        assert_eq!(a.intervals()[0], GtiInterval::new(0.0, 15.0));
        assert_eq!(a.intervals()[1], GtiInterval::new(20.0, 30.0));
        assert!((a.exposure() - 25.0).abs() < 1e-12);

        // Merge gaps < 6.0 s merges [0,15] and [20,30] (gap=5.0).
        let merged = a.merge_gaps(6.0);
        assert_eq!(merged.len(), 1);
        assert!((merged.exposure() - 30.0).abs() < 1e-12);

        // Intersect & subtract.
        let b = GtiSet::from_pairs(&[(5.0, 25.0)]);
        let inter = a.intersect(&b);
        assert_eq!(
            inter.intervals(),
            &[GtiInterval::new(5.0, 15.0), GtiInterval::new(20.0, 25.0)]
        );
        let sub = a.subtract(&b);
        assert_eq!(
            sub.intervals(),
            &[GtiInterval::new(0.0, 5.0), GtiInterval::new(25.0, 30.0)]
        );
        assert_eq!(a.drop_short(12.0).intervals(), &[GtiInterval::new(0.0, 15.0)]);

        // Overlap & filtering.
        assert!((a.overlap_duration(10.0, 25.0) - 10.0).abs() < 1e-12);
        let mask = a.filter_mask(&[-1.0, 0.0, 12.0, 17.5, 25.0, 35.0]);
        assert_eq!(mask, vec![false, true, true, false, true, false]);
    }

    #[test]
    fn orbit_taylor_and_analytical_barycentric_physics() {
        // Circular 650 km equatorial LEO orbit (R = 7028 km, v = 7.53 km/s, omega = v/R).
        let r_km = 7028.0f64;
        let v_kms = 7.53f64;
        let omega = v_kms / r_km;
        let t0 = mjd_utc_to_astrosat_utc(57394.0); // 2016-01-07

        let mut times = Vec::new();
        let mut pos_km = Vec::new();
        let mut vel_kms = Vec::new();
        for i in 0..100 {
            let dt = (i as f64) * 10.0;
            let th = omega * dt;
            times.push(t0 + dt);
            pos_km.push([r_km * th.cos(), r_km * th.sin(), 0.0]);
            vel_kms.push([-v_kms * th.sin(), v_kms * th.cos(), 0.0]);
        }

        let orb = OrbitTable::from_arrays(times.clone(), pos_km.clone(), vel_kms.clone()).unwrap();
        let st = orb.interpolate(t0 + 25.0).unwrap();
        let th_true = omega * 25.0;
        let err_m = ((st.pos_m[0] - r_km * 1000.0 * th_true.cos()).powi(2)
            + (st.pos_m[1] - r_km * 1000.0 * th_true.sin()).powi(2))
        .sqrt();
        assert!(err_m < 2.0, "orbit Taylor interp error too large: {err_m} m");

        // Crab pulsar (RA=83.6331, DEC=22.0145) on 2016-01-07 (near opposition -> large positive Roemer!).
        let calc = BaryCalculator::new(
            83.6331,
            22.0145,
            CoordFrame::Icrs,
            orb,
            BaryEphemeris::analytical(),
        );
        let bc = calc.correct_time(t0 + 25.0).unwrap();
        // On Jan 7, Earth is near RA ~ 105 deg (close to Crab RA 83.6 deg), so Earth Roemer delay is ~+460 s!
        assert!(
            bc.roemer_earth_s > 400.0 && bc.roemer_earth_s < 510.0,
            "unexpected Crab Jan-7 Roemer delay: {} s",
            bc.roemer_earth_s
        );
        // Spacecraft Roemer delay <= R_orbit / c ~ 23.4 ms.
        assert!(bc.roemer_sc_s.abs() <= 0.025);
        // Spacecraft Lorentz term <= R * v_earth / c^2 ~ 2.4 us.
        assert!(bc.lorentz_sc_s.abs() <= 3.0e-6);
        // Shapiro term is O(10 us) near opposition.
        assert!(bc.shapiro_s.abs() < 20.0e-6);

        // Verify full FITS barycenter round-trip (including LAXPC Format A without MJDREF in HDU 0).
        let dir = tempdir().unwrap();
        let orb_path = dir.path().join("test.orb");
        let in_evt = dir.path().join("in_evt.fits");
        let out_evt = dir.path().join("out_bary.fits");

        let meta = OgipMeta {
            instrume: "LAXPC".into(),
            tstart: t0,
            tstop: t0 + 900.0,
            exposure: 500.0,
            ..Default::default()
        };
        let p0 = Hdu::empty_primary();
        let orb_hdu = build_orbit_hdu(&times, &pos_km, &vel_kms, &meta).unwrap();
        write_fits(&orb_path, &mut [p0.clone(), orb_hdu], true).unwrap();

        let evt_times = vec![t0 + 10.0, t0 + 25.0, t0 + 50.0];
        let evt_hdu = BinTableBuilder::new("EVENTS")
            .add_f64("TIME", Some("s"), evt_times.clone())
            .add_i16("CHANNEL", None, vec![10, 20, 30])
            .header_i64("MJDREFI", 55197, "MJD ref")
            .header_f64("MJDREFF", 0.0, "MJD ref frac")
            .header_f64("TSTART", t0 + 10.0, "start")
            .header_f64("TSTOP", t0 + 50.0, "stop")
            .build()
            .unwrap();
        let gti = GtiSet::from_pairs(&[(t0 + 10.0, t0 + 50.0)]);
        let gti_hdu = gti.to_fits_hdu("GTI", &meta).unwrap();
        write_fits(&in_evt, &mut [p0, evt_hdu, gti_hdu], true).unwrap();

        let summary = barycenter_fits(
            &in_evt,
            &out_evt,
            &orb_path,
            83.6331,
            22.0145,
            &BaryFitsOptions::default(),
        )
        .unwrap();
        assert_eq!(summary.events_corrected, 3);
        assert_eq!(summary.gti_rows_corrected, 1);

        let out_fits = FitsFile::open(&out_evt).unwrap();
        assert_eq!(out_fits.hdu(0).unwrap().header().get_str("TIMESYS"), Some("TDB"));
        assert_eq!(
            out_fits.hdu(1).unwrap().header().get_str("TIMEREF"),
            Some("SOLARSYSTEM")
        );
        assert_eq!(out_fits.hdu(1).unwrap().verify_checksum(), (true, true));
        let out_times = out_fits
            .hdu_by_name("EVENTS")
            .unwrap()
            .as_bintable()
            .unwrap()
            .read_col_f64("TIME")
            .unwrap();
        assert!((out_times[1] - bc.t_tdb_s).abs() < 1e-9);
    }

    #[test]
    fn caldb_index_query_works() {
        let dir = tempdir().unwrap();
        let indx_path = dir.path().join("caldb.indx");
        let p0 = Hdu::empty_primary();
        let cif = BinTableBuilder::new("CIF")
            .add_strings("TELESCOP", 10, vec!["ASTROSAT".into(), "ASTROSAT".into()])
            .add_strings("INSTRUME", 10, vec!["CZTI".into(), "CZTI".into()])
            .add_strings("DETNAM", 10, vec!["NONE".into(), "NONE".into()])
            .add_strings("FILTER", 10, vec!["NONE".into(), "NONE".into()])
            .add_strings("CAL_CNAM", 12, vec!["GAIN".into(), "GAIN".into()])
            .add_strings("CAL_DIR", 20, vec!["bcf".into(), "bcf".into()])
            .add_strings(
                "CAL_FILE",
                32,
                vec!["gain_v01.fits".into(), "gain_v02.fits".into()],
            )
            .add_i16("CAL_XNO", None, vec![1, 1])
            .add_strings("CAL_VSD", 12, vec!["2015-09-28".into(), "2016-05-01".into()])
            .add_strings("CAL_VST", 10, vec!["00:00:00".into(), "00:00:00".into()])
            .add_f64("REF_TIME", None, vec![57293.0, 57509.0])
            .add_i16("CAL_QUAL", None, vec![0, 0])
            .build()
            .unwrap();
        write_fits(&indx_path, &mut [p0, cif], true).unwrap();

        let caldb = CaldbIndex::open(dir.path(), None).unwrap();
        // Observation on MJD 57400 (Jan 2016) should pick v01; MJD 57600 (Jul 2016) should pick v02!
        let best_early = caldb
            .query_best("ASTROSAT", "CZTI", Some("Q0"), "GAIN", Some(57400.0), None)
            .unwrap();
        assert_eq!(best_early.cal_file, "gain_v01.fits");

        let best_late = caldb
            .query_best("ASTROSAT", "CZTI", Some("Q0"), "GAIN", Some(57600.0), None)
            .unwrap();
        assert_eq!(best_late.cal_file, "gain_v02.fits");
    }

    fn test_data_root() -> std::path::PathBuf {
        std::env::var("ASTROSAT_TEST_DATA")
            .map(std::path::PathBuf::from)
            .unwrap_or_else(|_| std::path::PathBuf::from("/mnt/d/ISRO/data"))
    }

    #[test]
    fn real_astrosat_datasets_integration() {
        let data_root = test_data_root();
        let caldb_root = data_root.join("caldb/czti");
        if !caldb_root.exists() {
            eprintln!(
                "Skipping real AstroSat datasets test (data dir not present: {})",
                caldb_root.display()
            );
            return;
        }

        // 1. Real CZTI CALDB CIF query.
        let caldb = CaldbIndex::open(&caldb_root, None).expect("open real CZTI caldb.indx");
        assert!(!caldb.entries().is_empty());
        let gain = caldb
            .query_best("ASTROSAT", "CZTI", None, "GAIN", Some(197_000_000.0), None)
            .expect("find CZTI GAIN in real CALDB");
        assert!(gain.full_path.exists(), "GAIN path missing: {:?}", gain.full_path);
        let badpix = caldb
            .query_best("ASTROSAT", "CZTI", None, "BADPIX", Some(197_000_000.0), None)
            .expect("find CZTI BADPIX in real CALDB");
        assert!(badpix.full_path.exists(), "BADPIX path missing: {:?}", badpix.full_path);

        // 2. Real CZTI ObsID 9000000406 TCT, MKF, ORB.
        let czti_l1 = data_root.join("sample/czti/czti_sample_data/20160331_T01_112T01_9000000406_level1");
        let czti_tct = CztiTct::from_fits(
            czti_l1.join("aux/AS1T01_112T01_9000000406czt_level1.tct"),
        )
        .expect("read real CZTI TCT");
        assert!(czti_tct.num_segments() >= 1);

        let czti_mkf = MkfTable::from_fits(
            czti_l1.join("AS1T01_112T01_9000000406czt_level2.mkf"),
        )
        .expect("read real CZTI MKF");
        assert!(czti_mkf.len() > 100);
        assert!(czti_mkf.column("ELV").is_some());

        let czti_orb = OrbitTable::from_fits(
            data_root.join("sample/czti_orb/AS1T01_112T01_9000000406czt_level1.orb"),
        )
        .expect("read real CZTI ORB");
        assert!(czti_orb.len() > 100);
        let mid_t = 0.5 * (czti_orb.tstart() + czti_orb.tstop());
        let st = czti_orb.interpolate(mid_t).expect("interp real CZTI ORB");
        let r_km = (st.pos_m[0].powi(2) + st.pos_m[1].powi(2) + st.pos_m[2].powi(2)).sqrt() * 1e-3;
        assert!(
            (6900.0..7200.0).contains(&r_km),
            "AstroSat LEO orbital radius should be ~7028 km, got {r_km:.1} km"
        );

        // 3. Real LAXPC ObsID 9000000252 (Crab) TCT, MKF, ORB + real barycentric FITS correction.
        let laxpc_aux = data_root.join("sample/laxpc_l1/20160107_G02_010T01_9000000252_level1/laxpc/01509/aux");
        let laxpc_tct_path = laxpc_aux.join("AS1G02_010T01_9000000252lxp_level1.tct");
        let laxpc_orb_path = laxpc_aux.join("aux1/AS1G02_010T01_9000000252lxp_level1.orb");

        let laxpc_tct = LaxpcTct::from_fits(&laxpc_tct_path).expect("read real LAXPC TCT");
        assert!(laxpc_tct.len() > 10);

        let laxpc_orb = OrbitTable::from_fits(&laxpc_orb_path).expect("read real LAXPC ORB");
        assert!(laxpc_orb.len() > 10);
        let orb_mid = 0.5 * (laxpc_orb.tstart() + laxpc_orb.tstop());
        let calc = BaryCalculator::new(
            83.6331,
            22.0145,
            CoordFrame::Icrs,
            laxpc_orb,
            BaryEphemeris::analytical(),
        );
        let bc = calc.correct_time(orb_mid).expect("barycenter Crab timestamp");
        assert!(
            bc.roemer_earth_s > 400.0 && bc.roemer_earth_s < 510.0,
            "Crab Jan-2016 Roemer should be ~+460 s, got {} s",
            bc.roemer_earth_s
        );
    }

    #[test]
    fn safety_edge_cases_nan_double_bary_and_primary_dates() {
        // 1. ISO-8601 never panics on NaN / Inf / extreme values.
        assert_eq!(astrosat_utc_to_iso8601(f64::NAN), "0000-00-00T00:00:00.000");
        assert_eq!(astrosat_utc_to_iso8601(f64::INFINITY), "0000-00-00T00:00:00.000");
        assert_eq!(astrosat_utc_to_iso8601(1.0e30), "0000-00-00T00:00:00.000");

        // 2. CztiTct rejects NaN/Inf cleanly.
        let czti = CztiTct::new(
            TctTable::from_pairs(vec![100.0, 200.0], vec![10.0, 110.0]).unwrap(),
        )
        .unwrap();
        assert!(czti.interpolate(f64::NAN, 2.0, None).is_err());
        assert!(czti.interpolate(f64::INFINITY, 2.0, None).is_err());

        // 3. OrbitTable rejects NaN interpolation & empty 0-row ORBIT FITS cleanly.
        let orb = OrbitTable::from_arrays(
            vec![100.0, 110.0],
            vec![[7000.0, 0.0, 0.0], [6999.0, 75.0, 0.0]],
            vec![[0.0, 7.5, 0.0], [-0.08, 7.5, 0.0]],
        )
        .unwrap();
        assert!(orb.interpolate(f64::NAN).is_err());

        let dir = tempdir().unwrap();
        let empty_orb = dir.path().join("empty.orb");
        let meta = OgipMeta::default();
        let p0 = Hdu::empty_primary();
        let empty_hdu = build_orbit_hdu(&[], &[], &[], &meta).unwrap();
        write_fits(&empty_orb, &mut [p0.clone(), empty_hdu], false).unwrap();
        assert!(OrbitTable::from_fits(&empty_orb).is_err());

        // 4. barycenter_fits updates Primary DATE-OBS/DATE-END from global TSTART/TSTOP
        //    and rejects double-barycentering a file already in TIMESYS='TDB'.
        let t0 = mjd_utc_to_astrosat_utc(57394.0);
        let orb_path = dir.path().join("ok.orb");
        let times = vec![t0, t0 + 10.0, t0 + 20.0, t0 + 30.0];
        let pos = vec![[7000.0, 0.0, 0.0]; 4];
        let vel = vec![[0.0, 7.5, 0.0]; 4];
        let orb_hdu = build_orbit_hdu(&times, &pos, &vel, &meta).unwrap();
        write_fits(&orb_path, &mut [p0.clone(), orb_hdu], false).unwrap();

        let in_evt = dir.path().join("in.evt");
        let out_bary = dir.path().join("bary.evt");
        let out_bary2 = dir.path().join("bary2.evt");
        let mut p_in = Hdu::empty_primary();
        p_in.header.set_f64("TSTART", t0, "");
        p_in.header.set_f64("TSTOP", t0 + 29.0, "");
        let evt = BinTableBuilder::new("EVENTS")
            .add_f64("TIME", Some("s"), vec![t0 + 5.0, t0 + 15.0])
            .header_str("TIMESYS", "UTC", "")
            .header_i64("MJDREFI", 55197, "")
            .header_f64("MJDREFF", 0.0, "")
            .build()
            .unwrap();
        write_fits(&in_evt, &mut [p_in, evt], false).unwrap();

        let s = barycenter_fits(
            &in_evt,
            &out_bary,
            &orb_path,
            83.6331,
            22.0145,
            &BaryFitsOptions::default(),
        )
        .unwrap();
        let f = FitsFile::open(&out_bary).unwrap();
        let p_hdr = f.hdu(0).unwrap().header();
        // When Primary already had its own TSTART/TSTOP, Primary keeps its own barycentered bounds
        // and DATE-OBS/DATE-END match them consistently (never inflated by aux HDUs).
        let p_tstart = p_hdr.get_f64("TSTART").unwrap();
        let p_tstop = p_hdr.get_f64("TSTOP").unwrap();
        assert!(p_tstart < s.tstart_tdb.unwrap() && p_tstop > s.tstop_tdb.unwrap());
        assert_eq!(
            p_hdr.get_str("DATE-OBS"),
            Some(astrosat_utc_to_iso8601(p_tstart).as_str())
        );
        assert_eq!(
            p_hdr.get_str("DATE-END"),
            Some(astrosat_utc_to_iso8601(p_tstop).as_str())
        );

        // When Primary lacks TSTART/TSTOP, it inherits global_tstart/global_tstop + matching dates.
        let in_no_ts = dir.path().join("in_nots.evt");
        let out_no_ts = dir.path().join("bary_nots.evt");
        let evt2 = BinTableBuilder::new("EVENTS")
            .add_f64("TIME", Some("s"), vec![t0 + 5.0, t0 + 15.0])
            .header_str("TIMESYS", "UTC", "")
            .header_i64("MJDREFI", 55197, "")
            .header_f64("MJDREFF", 0.0, "")
            .build()
            .unwrap();
        write_fits(&in_no_ts, &mut [Hdu::empty_primary(), evt2], false).unwrap();
        let s2 = barycenter_fits(
            &in_no_ts,
            &out_no_ts,
            &orb_path,
            83.6331,
            22.0145,
            &BaryFitsOptions::default(),
        )
        .unwrap();
        let f2 = FitsFile::open(&out_no_ts).unwrap();
        let p2 = f2.hdu(0).unwrap().header();
        assert_eq!(
            p2.get_str("DATE-OBS"),
            Some(astrosat_utc_to_iso8601(s2.tstart_tdb.unwrap()).as_str())
        );
        assert_eq!(
            p2.get_str("DATE-END"),
            Some(astrosat_utc_to_iso8601(s2.tstop_tdb.unwrap()).as_str())
        );

        // Re-running barycenter_fits on an already-TDB file MUST error (no silent double-shift).
        assert!(barycenter_fits(
            &out_bary,
            &out_bary2,
            &orb_path,
            83.6331,
            22.0145,
            &BaryFitsOptions::default(),
        )
        .is_err());

        // 4b. 0-row EVENTS extension carrying UTC TSTART/TSTOP MUST still barycenter its header
        //     TSTART/TSTOP & DATE-OBS/DATE-END (and unsorted GTI-only FITS uses min/max bounds).
        let in_empty_evt = dir.path().join("in_empty_evt.fits");
        let out_empty_evt = dir.path().join("out_empty_evt.fits");
        let empty_evt_hdu = BinTableBuilder::new("Q3")
            .add_f64("TIME", Some("s"), vec![])
            .header_str("TIMESYS", "UTC", "")
            .header_i64("MJDREFI", 55197, "")
            .header_f64("MJDREFF", 0.0, "")
            .header_f64("TSTART", t0 + 2.0, "")
            .header_f64("TSTOP", t0 + 25.0, "")
            .build()
            .unwrap();
        let unsorted_gti = BinTableBuilder::new("GTI")
            .add_f64("START", Some("s"), vec![t0 + 15.0, t0 + 3.0])
            .add_f64("STOP", Some("s"), vec![t0 + 28.0, t0 + 10.0])
            .header_str("TIMESYS", "UTC", "")
            .header_i64("MJDREFI", 55197, "")
            .header_f64("MJDREFF", 0.0, "")
            .build()
            .unwrap();
        write_fits(
            &in_empty_evt,
            &mut [Hdu::empty_primary(), empty_evt_hdu, unsorted_gti],
            false,
        )
        .unwrap();
        let s_empty = barycenter_fits(
            &in_empty_evt,
            &out_empty_evt,
            &orb_path,
            83.6331,
            22.0145,
            &BaryFitsOptions::default(),
        )
        .unwrap();
        let f_empty = FitsFile::open(&out_empty_evt).unwrap();
        let q3_hdr = f_empty.hdu_by_name("Q3").unwrap().header();
        let q3_ts = q3_hdr.get_f64("TSTART").unwrap();
        let q3_te = q3_hdr.get_f64("TSTOP").unwrap();
        // Crab Jan-2016 shift is ~+530 s (TT-UTC ~68 s + Roemer ~460 s); must NOT remain UTC!
        assert!(q3_ts > t0 + 400.0 && q3_te > t0 + 400.0);
        assert_eq!(
            q3_hdr.get_str("DATE-OBS"),
            Some(astrosat_utc_to_iso8601(q3_ts).as_str())
        );
        // Unsorted GTI min/max propagated cleanly to GTI header and Primary fallback.
        let gti_hdr = f_empty.hdu_by_name("GTI").unwrap().header();
        assert!(gti_hdr.get_f64("TSTART").unwrap() < gti_hdr.get_f64("TSTOP").unwrap());
        assert!((s_empty.tstart_tdb.unwrap() - gti_hdr.get_f64("TSTART").unwrap()).abs() < 1e-6);
        assert!((s_empty.tstop_tdb.unwrap() - gti_hdr.get_f64("TSTOP").unwrap()).abs() < 1e-6);

        // 5. Real NASA HEASARC FITS JPLEPH.405 & JPLEPH.200 support (when test data is present).
        let ephem_dir = test_data_root().join("ephem");
        let jpl405 = ephem_dir.join("JPLEPH.405");
        let jpl200 = ephem_dir.join("JPLEPH.200");
        if jpl405.exists() && jpl200.exists() {
            let e405 = BaryEphemeris::from_jpl_file(&jpl405).expect("load HEASARC FITS JPLEPH.405");
            assert_eq!(e405.plephem_name(), "JPL-DE405");
            let e200 = BaryEphemeris::from_jpl_file(&jpl200).expect("load HEASARC FITS JPLEPH.200");
            assert_eq!(e200.plephem_name(), "JPL-DE200");

            let orb_ok = OrbitTable::from_fits(&orb_path).unwrap();
            let calc405 = BaryCalculator::new(83.6331, 22.0145, CoordFrame::Icrs, orb_ok.clone(), e405);
            let calc200 = BaryCalculator::new(83.6331, 22.0145, CoordFrame::Fk5, orb_ok, e200);
            // DE200 + FK5 skips FK5->ICRS bias rotation (matches raw unit vector).
            assert_eq!(
                calc200.unit_vector_icrs(),
                source_unit_vector_icrs(83.6331, 22.0145, CoordFrame::Icrs)
            );
            let b405 = calc405.correct_time(t0 + 10.0).unwrap();
            let b200 = calc200.correct_time(t0 + 10.0).unwrap();
            // Physical DE200 (1982) vs DE405 (1998) Voyager outer-planet mass shift is ~1.5 ms.
            assert!(
                (b405.t_tdb_s - b200.t_tdb_s).abs() < 3.0e-3,
                "DE405 vs DE200 diff too large: {} s",
                (b405.t_tdb_s - b200.t_tdb_s).abs()
            );

            // HEASARC FITS JPLEPH.405 must match raw JPL binary lnxp1600p2200.405 to < 1 ps!
            let raw405 = ephem_dir.join("lnxp1600p2200.405");
            if raw405.exists() {
                let eraw = BaryEphemeris::from_jpl_file(&raw405).expect("load raw binary DE405");
                let orb_ok2 = OrbitTable::from_fits(&orb_path).unwrap();
                let calc_raw = BaryCalculator::new(83.6331, 22.0145, CoordFrame::Icrs, orb_ok2, eraw);
                let braw = calc_raw.correct_time(t0 + 10.0).unwrap();
                assert!(
                    (b405.t_tdb_s - braw.t_tdb_s).abs() < 1e-11,
                    "HEASARC FITS JPLEPH.405 vs raw binary DE405 mismatch: {} s",
                    (b405.t_tdb_s - braw.t_tdb_s).abs()
                );
            }
        }
    }
}
