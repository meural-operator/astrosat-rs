//! Unit and integration tests for `astrosat-cli` (CLI tree, config parser, inspect,
//! validate, timing periodicity, and help output across all subcommands).

use std::process::Command;

use astrosat_core::time::ASTROSAT_MJDREFI;
use astrosat_fits::{build_gti_hdu, write_fits, BinTableBuilder, Hdu, OgipMeta};
use tempfile::tempdir;

fn astrosat_bin() -> &'static str {
    env!("CARGO_BIN_EXE_astrosat")
}

#[test]
fn help_succeeds_for_every_subcommand() {
    let subcommands: &[&[&str]] = &[
        &["--help"],
        &["czti", "--help"],
        &["czti", "run", "--help"],
        &["czti", "decode", "--help"],
        &["czti", "bunchclean", "--help"],
        &["czti", "pha2energy", "--help"],
        &["czti", "gtigen", "--help"],
        &["czti", "gaas", "--help"],
        &["czti", "datasel", "--help"],
        &["czti", "clean", "--help"],
        &["czti", "evtencor", "--help"],
        &["czti", "evtclean", "--help"],
        &["czti", "dpigen", "--help"],
        &["czti", "image", "--help"],
        &["czti", "bindata", "--help"],
        &["czti", "addspec", "--help"],
        &["laxpc", "--help"],
        &["laxpc", "run", "--help"],
        &["laxpc", "decode", "--help"],
        &["laxpc", "gtigen", "--help"],
        &["laxpc", "lightcurve", "--help"],
        &["laxpc", "spectrum", "--help"],
        &["laxpc", "backshift", "--help"],
        &["laxpc", "orbit", "--help"],
        &["bary", "--help"],
        &["inspect", "--help"],
        &["validate", "--help"],
        &["timing", "--help"],
        &["ui", "--help"],
    ];

    for args in subcommands {
        let out = Command::new(astrosat_bin())
            .args(*args)
            .output()
            .expect("failed to execute astrosat binary");
        assert!(
            out.status.success(),
            "`astrosat {}` failed: {}",
            args.join(" "),
            String::from_utf8_lossy(&out.stderr)
        );
        let stdout = String::from_utf8_lossy(&out.stdout);
        assert!(
            stdout.contains("Usage:"),
            "`astrosat {}` missing Usage block",
            args.join(" ")
        );
    }
}

#[test]
fn inspect_validate_and_timing_on_synthetic_pulsar_fits() {
    let dir = tempdir().unwrap();
    let fits_path = dir.path().join("pulsar_test.fits");

    // Generate 2000 events modulated at 29.66 Hz over 200 seconds
    let freq_true = 29.66_f64;
    let t0 = 190_000_000.0_f64;
    let mut times = Vec::with_capacity(2000);
    let mut chans = Vec::with_capacity(2000);
    for i in 0..2000 {
        let dt_coarse = (i as f64) * 0.1;
        // Concentrate photons near phase 0.25 (pulsed)
        let cycle = (dt_coarse * freq_true).floor();
        let phase_offset = if i % 4 == 0 {
            ((i as f64) * 0.37) % 1.0
        } else {
            0.25 + 0.04 * (((i % 11) as f64) - 5.0) / 5.0
        };
        times.push(t0 + (cycle + phase_offset) / freq_true);
        chans.push((i % 256) as i16);
    }
    times.sort_by(f64::total_cmp);

    let meta = OgipMeta {
        instrume: "LAXPC".into(),
        detnam: Some("LX10".into()),
        object: Some("CRAB_SYNTH".into()),
        ra_obj: Some(83.6331),
        dec_obj: Some(22.0145),
        tstart: t0,
        tstop: t0 + 200.0,
        exposure: 200.0,
        mjdrefi: ASTROSAT_MJDREFI,
        mjdreff: 0.0,
        timesys: "TDB".into(),
        timeref: "SOLARSYSTEM".into(),
        ..Default::default()
    };

    let mut primary = Hdu::empty_primary();
    astrosat_fits::apply_ogip_meta(&mut primary.header, &meta);
    let mut ev = BinTableBuilder::new("EVENTS")
        .add_f64("TIME", Some("s"), times)
        .add_i16("CHANNEL", Some("chan"), chans)
        .build()
        .unwrap();
    astrosat_fits::apply_ogip_meta(&mut ev.header, &meta);
    let gti = build_gti_hdu("GTI", &[t0], &[t0 + 200.0], &meta).unwrap();

    write_fits(&fits_path, &mut [primary, ev, gti], true).unwrap();

    // 1. `astrosat --json inspect --full`
    let insp = Command::new(astrosat_bin())
        .args(["--json", "inspect", "--full", fits_path.to_str().unwrap()])
        .output()
        .unwrap();
    assert!(
        insp.status.success(),
        "inspect failed: {}",
        String::from_utf8_lossy(&insp.stderr)
    );
    let insp_json: serde_json::Value = serde_json::from_slice(&insp.stdout).unwrap();
    assert_eq!(insp_json["files"][0]["health"]["all_checksums_ok"], true);
    assert_eq!(insp_json["files"][0]["health"]["total_events"], 2000);

    // 2. `astrosat --json validate --strict`
    let val = Command::new(astrosat_bin())
        .args(["--json", "validate", "--strict", fits_path.to_str().unwrap()])
        .output()
        .unwrap();
    assert!(
        val.status.success(),
        "validate --strict failed: {} / {}",
        String::from_utf8_lossy(&val.stdout),
        String::from_utf8_lossy(&val.stderr)
    );
    let val_json: serde_json::Value = serde_json::from_slice(&val.stdout).unwrap();
    assert_eq!(val_json["passed"], true);
    assert_eq!(val_json["total_errors"], 0);

    // 3. `astrosat --threads 2 --json timing`
    let tim = Command::new(astrosat_bin())
        .args([
            "--threads",
            "2",
            "--json",
            "timing",
            fits_path.to_str().unwrap(),
            "--freq",
            "29.66",
            "--search-span",
            "0.02",
            "--harmonics",
            "2",
            "--nbins",
            "32",
        ])
        .output()
        .unwrap();
    assert!(
        tim.status.success(),
        "timing failed: {}",
        String::from_utf8_lossy(&tim.stderr)
    );
    let tim_json: serde_json::Value = serde_json::from_slice(&tim.stdout).unwrap();
    let best_f = tim_json["best_freq_hz"].as_f64().unwrap();
    let z2 = tim_json["max_z2_n"].as_f64().unwrap();
    assert!((best_f - freq_true).abs() < 0.005, "best_f={best_f}");
    assert!(z2 > 100.0, "expected strong pulsar Z^2_2, got {z2}");
}

#[test]
fn config_toml_and_json_work_with_cli() {
    let dir = tempdir().unwrap();
    let toml_path = dir.path().join("cfg.toml");
    std::fs::write(
        &toml_path,
        r#"
# AstroSat CLI config test
threads = 2

[timing]
freq_hz = 29.66
harmonics = 2
nbins = 16
"#,
    )
    .unwrap();

    // Verify --config loads cleanly even with inspect --help or validate
    let fits_path = dir.path().join("tiny.fits");
    let p = Hdu::empty_primary();
    let g = build_gti_hdu("GTI", &[100.0], &[200.0], &OgipMeta::default()).unwrap();
    write_fits(&fits_path, &mut [p, g], true).unwrap();

    let out = Command::new(astrosat_bin())
        .args([
            "--config",
            toml_path.to_str().unwrap(),
            "--json",
            "validate",
            fits_path.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(out.status.success());
}

#[test]
fn bary_command_corrects_events_and_passes_validate() {
    use astrosat_core::time::mjd_utc_to_astrosat_utc;
    use astrosat_fits::build_orbit_hdu;

    let dir = tempdir().unwrap();
    let orb_path = dir.path().join("test.orb");
    let in_evt = dir.path().join("evt.fits");
    let out_evt = dir.path().join("evt_bary.fits");

    let t0 = mjd_utc_to_astrosat_utc(57394.0);
    let r_km = 7028.0_f64;
    let v_kms = 7.53_f64;
    let omega = v_kms / r_km;
    let mut times = Vec::new();
    let mut pos = Vec::new();
    let mut vel = Vec::new();
    for i in 0..30 {
        let dt = (i as f64) * 10.0;
        let th = omega * dt;
        times.push(t0 + dt);
        pos.push([r_km * th.cos(), r_km * th.sin(), 0.0]);
        vel.push([-v_kms * th.sin(), v_kms * th.cos(), 0.0]);
    }
    let meta = OgipMeta {
        instrume: "LAXPC".into(),
        ra_obj: Some(83.6331),
        dec_obj: Some(22.0145),
        tstart: t0,
        tstop: t0 + 200.0,
        exposure: 200.0,
        ..Default::default()
    };
    let orb_hdu = build_orbit_hdu(&times, &pos, &vel, &meta).unwrap();
    write_fits(&orb_path, &mut [Hdu::empty_primary(), orb_hdu], true).unwrap();

    let mut ev = BinTableBuilder::new("EVENTS")
        .add_f64("TIME", Some("s"), vec![t0 + 10.0, t0 + 50.0, t0 + 100.0])
        .add_i16("CHANNEL", None, vec![10, 20, 30])
        .build()
        .unwrap();
    astrosat_fits::apply_ogip_meta(&mut ev.header, &meta);
    let gti = build_gti_hdu("GTI", &[t0 + 10.0], &[t0 + 100.0], &meta).unwrap();
    write_fits(&in_evt, &mut [Hdu::empty_primary(), ev, gti], true).unwrap();

    let out = Command::new(astrosat_bin())
        .args([
            "--json",
            "bary",
            in_evt.to_str().unwrap(),
            "-o",
            out_evt.to_str().unwrap(),
            "--orbit",
            orb_path.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "bary failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let j: serde_json::Value = serde_json::from_slice(&out.stdout).unwrap();
    assert_eq!(j["events_corrected"], 3);
    assert_eq!(j["gti_rows_corrected"], 1);
    assert_eq!(j["plephem"], "ASTROSAT-ANALYTICAL");

    let val = Command::new(astrosat_bin())
        .args(["validate", "--strict", out_evt.to_str().unwrap()])
        .output()
        .unwrap();
    assert!(val.status.success());
}

#[test]
fn validate_returns_nonzero_on_invalid_file() {
    let dir = tempdir().unwrap();
    let bad = dir.path().join("bad.fits");
    std::fs::write(&bad, b"NOT_A_FITS_FILE").unwrap();

    let out = Command::new(astrosat_bin())
        .args(["--json", "validate", bad.to_str().unwrap()])
        .output()
        .unwrap();
    assert!(!out.status.success());
    let j: serde_json::Value = serde_json::from_slice(&out.stdout).unwrap();
    assert_eq!(j["passed"], false);
    assert!(j["total_errors"].as_u64().unwrap() >= 1);
}

#[test]
fn cli_defaults_symmetrically_override_config_and_energy_range_requires_caldb() {
    use astrosat_core::time::mjd_utc_to_astrosat_utc;
    use astrosat_fits::build_orbit_hdu;

    let dir = tempdir().unwrap();
    let cfg_path = dir.path().join("override.toml");
    std::fs::write(
        &cfg_path,
        r#"
[bary]
frame = "fk5"
ephemeris = "/nonexistent/de405"
write_checksums = false

[timing]
fdot = 1.0e-5
harmonics = 4
nbins = 16
"#,
    )
    .unwrap();

    let orb_path = dir.path().join("test.orb");
    let in_evt = dir.path().join("evt.fits");
    let out_evt = dir.path().join("evt_bary.fits");

    let t0 = mjd_utc_to_astrosat_utc(57394.0);
    let mut times = Vec::new();
    let mut pos = Vec::new();
    let mut vel = Vec::new();
    for i in 0..10 {
        times.push(t0 + (i as f64) * 10.0);
        pos.push([7028.0, 0.0, 0.0]);
        vel.push([0.0, 7.53, 0.0]);
    }
    let meta = OgipMeta {
        instrume: "LAXPC".into(),
        ra_obj: Some(83.6331),
        dec_obj: Some(22.0145),
        tstart: t0,
        tstop: t0 + 90.0,
        exposure: 90.0,
        ..Default::default()
    };
    let orb_hdu = build_orbit_hdu(&times, &pos, &vel, &meta).unwrap();
    write_fits(&orb_path, &mut [Hdu::empty_primary(), orb_hdu], true).unwrap();

    let ev_times: Vec<f64> = (0..16).map(|i| t0 + (i as f64) * 0.1).collect();
    let chans: Vec<i16> = vec![10; 16];
    let mut ev = BinTableBuilder::new("EVENTS")
        .add_f64("TIME", Some("s"), ev_times)
        .add_i16("CHANNEL", None, chans)
        .build()
        .unwrap();
    astrosat_fits::apply_ogip_meta(&mut ev.header, &meta);
    let gti = build_gti_hdu("GTI", &[t0], &[t0 + 90.0], &meta).unwrap();
    write_fits(&in_evt, &mut [Hdu::empty_primary(), ev, gti], true).unwrap();

    // Explicitly passing CLI defaults (--frame icrs, --ephem analytical, --checksum) must override config!
    let bary_out = Command::new(astrosat_bin())
        .args([
            "--config",
            cfg_path.to_str().unwrap(),
            "--json",
            "bary",
            in_evt.to_str().unwrap(),
            "-o",
            out_evt.to_str().unwrap(),
            "--orbit",
            orb_path.to_str().unwrap(),
            "--frame",
            "icrs",
            "--ephem",
            "analytical",
            "--checksum",
        ])
        .output()
        .unwrap();
    assert!(
        bary_out.status.success(),
        "bary override failed: {}",
        String::from_utf8_lossy(&bary_out.stderr)
    );
    let bj: serde_json::Value = serde_json::from_slice(&bary_out.stdout).unwrap();
    assert_eq!(bj["frame"], "ICRS");
    assert_eq!(bj["plephem"], "ASTROSAT-ANALYTICAL");

    // Explicitly passing CLI defaults (--fdot 0, --harmonics 2, --nbins 32) must override config!
    let tim_out = Command::new(astrosat_bin())
        .args([
            "--config",
            cfg_path.to_str().unwrap(),
            "--json",
            "timing",
            out_evt.to_str().unwrap(),
            "--freq",
            "10.0",
            "--fdot",
            "0",
            "--harmonics",
            "2",
            "--nbins",
            "32",
        ])
        .output()
        .unwrap();
    assert!(
        tim_out.status.success(),
        "timing override failed: {}",
        String::from_utf8_lossy(&tim_out.stderr)
    );
    let tj: serde_json::Value = serde_json::from_slice(&tim_out.stdout).unwrap();
    assert_eq!(tj["fdot"], 0.0);
    assert_eq!(tj["harmonics"], 2);
    assert_eq!(tj["nbins"], 32);

    // `--energy-range` without `--caldb` must error explicitly instead of silently ignoring!
    let er_out = Command::new(astrosat_bin())
        .args([
            "laxpc",
            "run",
            "--input",
            dir.path().to_str().unwrap(),
            "--energy-range",
            "3:80",
        ])
        .output()
        .unwrap();
    assert!(!er_out.status.success());
    let stderr = String::from_utf8_lossy(&er_out.stderr);
    assert!(
        stderr.contains("--caldb") || stderr.contains("channel"),
        "expected actionable --caldb error, got: {stderr}"
    );

    // `laxpc spectrum --energy-range 3:20` without `--caldb` must also fail explicitly!
    let spec_er = Command::new(astrosat_bin())
        .args([
            "laxpc",
            "spectrum",
            "--ea",
            in_evt.to_str().unwrap(),
            "--mkf",
            in_evt.to_str().unwrap(),
            "--tct",
            in_evt.to_str().unwrap(),
            "--unit",
            "10",
            "--energy-range",
            "3:20",
            "-o",
            dir.path().join("out.spec").to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(!spec_er.status.success());
}

#[test]
fn neon_review_regressions_negative_args_mkf_precedence_and_timing_guards() {
    use astrosat_core::time::mjd_utc_to_astrosat_utc;
    use astrosat_fits::build_orbit_hdu;

    let dir = tempdir().unwrap();
    let t0 = mjd_utc_to_astrosat_utc(57394.0);

    // 1. Negative --dec and --fdot (hyphen values) + explicit --period overriding config freq_hz
    let orb_path = dir.path().join("good.orb");
    let in_evt = dir.path().join("evt.fits");
    let out_evt = dir.path().join("evt_south_bary.fits");
    let times: Vec<f64> = (0..10).map(|i| t0 + (i as f64) * 10.0).collect();
    let pos = vec![[7028.0, 0.0, 0.0]; 10];
    let vel = vec![[0.0, 7.53, 0.0]; 10];
    let meta = OgipMeta {
        instrume: "LAXPC".into(),
        tstart: t0,
        tstop: t0 + 90.0,
        exposure: 90.0,
        ..Default::default()
    };
    let orb_hdu = build_orbit_hdu(&times, &pos, &vel, &meta).unwrap();
    write_fits(&orb_path, &mut [Hdu::empty_primary(), orb_hdu], true).unwrap();

    let ev_times: Vec<f64> = (0..16).map(|i| t0 + (i as f64) * 0.1).collect();
    let chans: Vec<i16> = vec![10; 16];
    let mut ev = BinTableBuilder::new("EVENTS")
        .add_f64("TIME", Some("s"), ev_times)
        .add_i16("CHANNEL", None, chans)
        .build()
        .unwrap();
    astrosat_fits::apply_ogip_meta(&mut ev.header, &meta);
    let gti = build_gti_hdu("GTI", &[t0], &[t0 + 90.0], &meta).unwrap();
    write_fits(&in_evt, &mut [Hdu::empty_primary(), ev, gti], true).unwrap();

    let bary_south = Command::new(astrosat_bin())
        .args([
            "--json",
            "bary",
            in_evt.to_str().unwrap(),
            "-o",
            out_evt.to_str().unwrap(),
            "--orbit",
            orb_path.to_str().unwrap(),
            "--ra",
            "266.4168",
            "--dec",
            "-29.0078",
        ])
        .output()
        .unwrap();
    assert!(
        bary_south.status.success(),
        "negative --dec failed: {}",
        String::from_utf8_lossy(&bary_south.stderr)
    );
    let bj: serde_json::Value = serde_json::from_slice(&bary_south.stdout).unwrap();
    assert!((bj["dec_deg"].as_f64().unwrap() - (-29.0078)).abs() < 1e-6);

    // 2. Explicit --mkf must take precedence over a sibling AstroSat.orb (even if corrupt)
    let corrupt_sibling = dir.path().join("AstroSat.orb");
    std::fs::write(&corrupt_sibling, b"NOT_A_VALID_ORBIT_FITS").unwrap();
    let mkf_path = dir.path().join("custom.mkf");
    let mkf_t: Vec<f64> = (0..10).map(|i| t0 + (i as f64) * 10.0).collect();
    let mkf_hdu = BinTableBuilder::new("MKF")
        .add_f64("TIME", Some("s"), mkf_t)
        .add_f64("POSX", Some("km"), vec![7028.0; 10])
        .add_f64("POSY", Some("km"), vec![0.0; 10])
        .add_f64("POSZ", Some("km"), vec![0.0; 10])
        .add_f64("VELX", Some("km/s"), vec![0.0; 10])
        .add_f64("VELY", Some("km/s"), vec![7.53; 10])
        .add_f64("VELZ", Some("km/s"), vec![0.0; 10])
        .build()
        .unwrap();
    write_fits(&mkf_path, &mut [Hdu::empty_primary(), mkf_hdu], true).unwrap();

    let out_mkf_bary = dir.path().join("evt_mkf_bary.fits");
    let bary_mkf = Command::new(astrosat_bin())
        .args([
            "--json",
            "bary",
            in_evt.to_str().unwrap(),
            "-o",
            out_mkf_bary.to_str().unwrap(),
            "--mkf",
            mkf_path.to_str().unwrap(),
            "--ra",
            "83.6331",
            "--dec",
            "22.0145",
        ])
        .output()
        .unwrap();
    assert!(
        bary_mkf.status.success(),
        "--mkf precedence over sibling AstroSat.orb failed: {}",
        String::from_utf8_lossy(&bary_mkf.stderr)
    );
    let mj: serde_json::Value = serde_json::from_slice(&bary_mkf.stdout).unwrap();
    assert_eq!(mj["orbit"].as_str().unwrap(), mkf_path.to_str().unwrap());

    // 3. Negative --fdot and CLI --period overriding [timing] freq_hz in config
    let cfg_path = dir.path().join("timing_override.toml");
    std::fs::write(&cfg_path, "[timing]\nfreq_hz = 29.66\n").unwrap();
    let tim_out = Command::new(astrosat_bin())
        .args([
            "--config",
            cfg_path.to_str().unwrap(),
            "--json",
            "timing",
            out_evt.to_str().unwrap(),
            "--period",
            "0.2",
            "--fdot",
            "-3.68e-10",
        ])
        .output()
        .unwrap();
    assert!(
        tim_out.status.success(),
        "negative --fdot / --period override failed: {}",
        String::from_utf8_lossy(&tim_out.stderr)
    );
    let tj: serde_json::Value = serde_json::from_slice(&tim_out.stdout).unwrap();
    assert!((tj["best_freq_hz"].as_f64().unwrap() - 5.0).abs() < 1e-9);
    assert!((tj["fdot"].as_f64().unwrap() - (-3.68e-10)).abs() < 1e-18);

    // 4. Non-finite / OOM timing inputs rejected cleanly
    for bad_args in [
        vec!["timing", out_evt.to_str().unwrap(), "--freq", "NaN"],
        vec!["timing", out_evt.to_str().unwrap(), "--freq", "10.0", "--fdot", "inf"],
        vec!["timing", out_evt.to_str().unwrap(), "--freq", "10.0", "--n-trials", "0"],
        vec!["timing", out_evt.to_str().unwrap(), "--freq", "10.0", "--n-trials", "5000000"],
    ] {
        let res = Command::new(astrosat_bin()).args(&bad_args).output().unwrap();
        assert!(
            !res.status.success(),
            "expected failure for {:?}, got success",
            bad_args
        );
    }

    // 5. Invalid config bary.frame cleanly rejected
    let bad_cfg = dir.path().join("bad_frame.toml");
    std::fs::write(&bad_cfg, "[bary]\nframe = \"galactic\"\n").unwrap();
    let bad_bary = Command::new(astrosat_bin())
        .args([
            "--config",
            bad_cfg.to_str().unwrap(),
            "bary",
            in_evt.to_str().unwrap(),
            "-o",
            dir.path().join("should_fail.fits").to_str().unwrap(),
            "--orbit",
            orb_path.to_str().unwrap(),
            "--ra",
            "83.6",
            "--dec",
            "22.0",
        ])
        .output()
        .unwrap();
    assert!(!bad_bary.status.success());
}
