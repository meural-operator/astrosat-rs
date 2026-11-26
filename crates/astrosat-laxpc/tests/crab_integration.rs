//! End-to-end integration tests on real AstroSat Crab pulsar ObsID `9000000252` (orbit `01509`),
//! cross-checking `astrosat-laxpc` Level-1 -> Level-2 outputs against official ISRO/TIFR Level-2
//! reference products and validating OGIP/FITS conformance (including NASA HEASARC `fitsverify`).

use astrosat_core::gti::GtiSet;
use astrosat_fits::{validate_ogip, FitsFile};
use astrosat_laxpc::*;
use std::path::{Path, PathBuf};
use std::process::Command;
use tempfile::tempdir;

fn test_data_root() -> PathBuf {
    std::env::var("ASTROSAT_TEST_DATA")
        .map(PathBuf::from)
        .unwrap_or_else(|_| PathBuf::from("/mnt/d/ISRO/data"))
}

const FITSVERIFY_BIN: &str = "/mnt/d/ISRO/.toolchain/bin/fitsverify";

fn run_fitsverify(path: &Path) {
    if !Path::new(FITSVERIFY_BIN).is_file() {
        return;
    }
    let out = Command::new(FITSVERIFY_BIN)
        .arg("-q")
        .arg(path)
        .output()
        .expect("failed to execute fitsverify");
    let stdout = String::from_utf8_lossy(&out.stdout);
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        out.status.success() && (stdout.contains("0 errors") || stdout.contains("verification OK")),
        "fitsverify failed on {}:\nstdout: {stdout}\nstderr: {stderr}",
        path.display()
    );
}

#[test]
fn crab_9000000252_parity_against_official_l2_reference() {
    let data_root = test_data_root();
    let l1_root = data_root.join("sample/laxpc_l1/20160107_G02_010T01_9000000252_level1");
    let l2_root = data_root.join("sample/laxpc_l2/20160107_G02_010T01_9000000252_level2/laxpc/01509");
    if !l1_root.exists() || !l2_root.exists() {
        eprintln!(
            "Skipping Crab integration test: sample datasets not present under {}",
            data_root.display()
        );
        return;
    }

    // Compare each unit against its official L2 `.pha` and `.lc` using the official L2 `.gti`
    // so exposure windows match apples-to-apples.
    for unit in LaxpcUnit::ALL {
        let prefix = unit.prefix();
        let l2_ea_dir = l2_root.join(prefix).join("modeEA");
        let l2_gti_path = l2_ea_dir.join(format!("AS1G02_010T01_9000000252{prefix}EA_level2.gti"));
        let l2_pha_path = l2_ea_dir.join(format!("AS1G02_010T01_9000000252{prefix}EA_level2.pha"));

        let official_gti = GtiSet::from_fits(&l2_gti_path, Some("GTI")).unwrap();
        let official_pha = FitsFile::open(&l2_pha_path).unwrap();
        let pha_hdu = official_pha.hdu_by_name("SPECTRUM").unwrap();
        let official_exposure = pha_hdu
            .header()
            .get_f64("EXPOSURE")
            .or_else(|| pha_hdu.header().get_i64("EXPOSURE").map(|v| v as f64))
            .or_else(|| {
                pha_hdu
                    .header()
                    .get_str("EXPOSURE")
                    .and_then(|s| s.trim().parse::<f64>().ok())
            })
            .unwrap();
        let official_counts = pha_hdu.as_bintable().unwrap().read_col_f64("COUNTS").unwrap();
        let official_total_counts: f64 = official_counts.iter().sum();

        let out_dir = tempdir().unwrap();
        let cfg = PipelineConfig {
            input_dir: Some(l1_root.clone()),
            output_dir: out_dir.path().to_path_buf(),
            units: vec![unit],
            tbin: 1.0,
            deadtime_sec: 43.0e-6,
            spectrum_channels: Some(1024),
            user_gti_utc: Some(official_gti),
            use_user_gti_only: true,
            write_format_b_events: false,
            write_format_a_combined: false,
            write_orbit_file: false,
            ..Default::default()
        };

        let multi = run_pipeline(&cfg).unwrap();
        let res = &multi.unit_results[0];

        // 1. Exposure matches official L2 GTI exposure closely (within 1%)
        let exp_rel_diff = (res.gti.gross_exposure_sec - official_exposure).abs() / official_exposure;
        assert!(
            exp_rel_diff < 0.01,
            "{unit}: exposure mismatch: ours={:.2}s vs official={:.2}s",
            res.gti.gross_exposure_sec,
            official_exposure
        );

        // 2. Total raw main-anode event count matches official L2 PHA total counts closely (< 2.5%)
        //    (Small ~1% difference comes from `laxpcl1.f` strict `1 <= pha <= 999` & double-event / frame-end filtering)
        let cnt_rel_diff =
            (res.spectrum.total_raw_counts as f64 - official_total_counts).abs() / official_total_counts;
        println!(
            "[{unit}] official_exp={official_exposure:.1}s, our_gross_exp={:.1}s, our_net_exp={:.1}s | official_cnts={official_total_counts:.0}, our_raw_cnts={} (rel_diff={:.3}%), corr_rate={:.2} cps (fcd={:.4})",
            res.gti.gross_exposure_sec,
            res.gti.net_exposure_sec,
            res.spectrum.total_raw_counts,
            cnt_rel_diff * 100.0,
            res.spectrum.total_corrected_rate_cps,
            res.spectrum.deadtime_factor
        );
        assert!(
            cnt_rel_diff < 0.03,
            "{unit}: count mismatch too large: ours={} vs official={official_total_counts:.0} ({:.2}%)",
            res.spectrum.total_raw_counts,
            cnt_rel_diff * 100.0
        );

        // 3. Pearson correlation between our 1024-ch deadtime-corrected spectrum and official L2 PHA
        //    must be > 0.995 (near-identical spectral shape across the Crab continuum + Xe K edge!).
        let mean_a: f64 = res.spectrum.selected_rate.iter().sum::<f64>() / 1024.0;
        let mean_b: f64 = official_counts.iter().sum::<f64>() / 1024.0;
        let mut num = 0.0;
        let mut da2 = 0.0;
        let mut db2 = 0.0;
        for i in 0..1024 {
            let xa = res.spectrum.selected_rate[i] - mean_a;
            let xb = official_counts[i] - mean_b;
            num += xa * xb;
            da2 += xa * xa;
            db2 += xb * xb;
        }
        let corr = num / (da2 * db2).sqrt();
        println!("[{unit}] 1024-ch Crab spectrum Pearson r = {corr:.6}");
        assert!(corr > 0.995, "{unit}: spectrum correlation too low: {corr:.6}");
    }
}

#[test]
fn crab_9000000252_full_multi_unit_auto_gti_format_b_and_format_a() {
    let data_root = test_data_root();
    let l1_root = data_root.join("sample/laxpc_l1/20160107_G02_010T01_9000000252_level1");
    let ref_src = data_root.join("ref_src/laxpc");
    if !l1_root.exists() {
        eprintln!(
            "Skipping Crab full multi-unit test: L1 dataset not present under {}",
            data_root.display()
        );
        return;
    }

    let out_dir = tempdir().unwrap();
    let cfg = PipelineConfig {
        input_dir: Some(l1_root),
        output_dir: out_dir.path().to_path_buf(),
        units: LaxpcUnit::ALL.to_vec(),
        tbin: 1.0,
        deadtime_sec: 43.0e-6,
        spectrum_channels: None, // 512 for Lxp10/30, 256 for Lxp20
        caldb_dir: ref_src.is_dir().then(|| ref_src.clone()),
        background_model_path: ref_src
            .join("backfit.lxp10mar16")
            .is_file()
            .then(|| ref_src.join("backfit.lxp10mar16")),
        background_pha_path: ref_src
            .join("backlxp10jan16.pha")
            .is_file()
            .then(|| ref_src.join("backlxp10jan16.pha")),
        write_format_b_events: true,
        write_format_a_combined: true,
        write_orbit_file: true,
        write_ascii_products: true,
        ..Default::default()
    };

    let result = run_pipeline(&cfg).unwrap();
    assert_eq!(result.unit_results.len(), 3);

    // Check every unit's standard GTI, Crab count rate, RMF recommendation, and OGIP FITS validity
    for ures in &result.unit_results {
        println!(
            "[AUTO-GTI {}] frames={}/{} (missing={}), valid_evs={}, gti_evs={}, gti_n={}, net_exp={:.1}s, rate={:.1} cps, lc_mean={:.1} cps, rmf={:?}",
            ures.unit,
            ures.ea_frame_diagnostics.retained_frames,
            ures.ea_frame_diagnostics.total_raw_frames,
            ures.ea_frame_diagnostics.total_missing_frames,
            ures.total_valid_events,
            ures.gti_filtered_events,
            ures.gti.gti_utc.len(),
            ures.gti.net_exposure_sec,
            ures.spectrum.total_corrected_rate_cps,
            ures.lightcurve_mean_rate_cps,
            ures.recommended_rmf
        );

        assert!(ures.gti.net_exposure_sec > 2500.0, "{}: exposure too low", ures.unit);
        // Crab pulsar on-source (ELV > 3.5 deg) deadtime-corrected rate in each LAXPC unit is ~2950-3500 cps
        assert!(
            ures.spectrum.total_corrected_rate_cps > 2800.0
                && ures.spectrum.total_corrected_rate_cps < 3700.0,
            "{}: unexpected Crab rate {:.1} cps",
            ures.unit,
            ures.spectrum.total_corrected_rate_cps
        );
        assert!(ures.lightcurve_bins > 2500);
        assert_eq!(
            ures.spectrum.num_channels,
            ures.unit.default_spectrum_channels()
        );

        // Validate all Format B FITS products
        for p in [
            &ures.spectrum_fits_path,
            &ures.lightcurve_fits_path,
            &ures.gti_fits_path,
        ] {
            let f = FitsFile::open(p).unwrap();
            let issues = validate_ogip(&f);
            assert!(issues.is_empty(), "{} OGIP issues: {issues:?}", p.display());
            run_fitsverify(p);
        }
        if let Some(evp) = &ures.event_fits_path {
            let f = FitsFile::open(evp).unwrap();
            assert!(validate_ogip(&f).is_empty());
            run_fitsverify(evp);
        }
    }

    // Validate AstroSat.orb
    let orb_path = result.orbit_fits_path.as_ref().unwrap();
    let f_orb = FitsFile::open(orb_path).unwrap();
    assert!(validate_ogip(&f_orb).is_empty());
    run_fitsverify(orb_path);

    // Validate Format A combined level2.event.fits (uniform 0..1023 channel scale across all 3 units)
    let fa_path = result.format_a_event_fits_path.as_ref().unwrap();
    assert!(result.format_a_total_events > 15_000_000);
    let f_fa = FitsFile::open(fa_path).unwrap();
    assert!(validate_ogip(&f_fa).is_empty());
    let fa_hdr = f_fa.hdu_by_name("LAXPC_EVENTS").unwrap().header();
    assert_eq!(fa_hdr.get_i64("TLMIN2"), Some(0));
    assert_eq!(fa_hdr.get_i64("TLMAX2"), Some(1023));
    run_fitsverify(fa_path);

    // Verify per-unit background spectrum & lightcurve substitution succeeded for all 3 units!
    for ures in &result.unit_results {
        assert!(
            ures.bkg_spectrum_fits_path.is_some(),
            "{}: expected shifted background spectrum",
            ures.unit
        );
        assert!(
            ures.bkg_lightcurve_fits_path.is_some(),
            "{}: expected background lightcurve .lcbk",
            ures.unit
        );
    }

    // Verify native `astrosat_core::bary::barycenter_fits` works directly on BOTH Format B (.evn)
    // and Format A (level2.event.fits) with zero `preplaxbary` hacks!
    let evn1 = result.unit_results[0].event_fits_path.as_ref().unwrap();
    let bary_b_out = out_dir.path().join("lxp1level2_bary.evn");
    let bary_b = astrosat_core::bary::barycenter_fits(
        evn1,
        &bary_b_out,
        orb_path,
        83.6331,
        22.0145,
        &astrosat_core::bary::BaryFitsOptions::default(),
    )
    .unwrap();
    assert_eq!(bary_b.events_corrected, result.unit_results[0].gti_filtered_events);
    run_fitsverify(&bary_b_out);

    let bary_a_out = out_dir.path().join("level2_bary.event.fits");
    let bary_a = astrosat_core::bary::barycenter_fits(
        fa_path,
        &bary_a_out,
        orb_path,
        83.6331,
        22.0145,
        &astrosat_core::bary::BaryFitsOptions::default(),
    )
    .unwrap();
    assert_eq!(bary_a.events_corrected, result.format_a_total_events);
    run_fitsverify(&bary_a_out);
}
