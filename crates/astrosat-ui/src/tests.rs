//! Comprehensive unit & live HTTP integration tests for `astrosat-ui`.

use std::f64::consts::TAU;
use std::io::{Read, Write};
use std::net::{SocketAddr, TcpStream};
use std::path::Path;
use std::time::Duration;

use astrosat_core::time::mjd_utc_to_astrosat_utc;
use astrosat_czti::bindata::{CztiResponse, CztiSpectrum};
use astrosat_czti::fits_io::{write_response_file, write_spectrum_file, CztiObsMeta};
use astrosat_czti::types::{GtiType, PI_CHANNELS};
use astrosat_fits::{
    build_gti_hdu, build_lightcurve_hdu, build_orbit_hdu, build_response_hdus, build_spectrum_hdu,
    write_fits, BinTableBuilder, FitsFile, Hdu, ImageBuilder, OgipMeta, ResponseMatrixData,
    SpectrumData, SpectrumMeta, WcsParams,
};
use serde_json::{json, Value};
use tempfile::tempdir;

use crate::jobs::JobManager;
use crate::pulse::{circular_xcorr_phase_lag, handle_pulse_multiband, BandSpec, MultibandInputSpec, MultibandRequest};
use crate::science::rebin_lightcurve;
use crate::{UiConfig, UiServer};

// ---------------------------------------------------------------------------
// Minimal zero-dep HTTP/1.1 test client over raw TcpStream
// ---------------------------------------------------------------------------

struct HttpResponse {
    status: u16,
    content_type: String,
    body: String,
}

fn http_request(addr: SocketAddr, method: &str, path: &str, body: Option<&str>) -> HttpResponse {
    http_request_with_ct(addr, method, path, "application/json", body)
}

fn http_request_with_ct(
    addr: SocketAddr,
    method: &str,
    path: &str,
    content_type: &str,
    body: Option<&str>,
) -> HttpResponse {
    let mut stream = TcpStream::connect(addr).expect("connect to UiServer");
    stream
        .set_read_timeout(Some(Duration::from_secs(15)))
        .unwrap();

    let payload = body.unwrap_or("");
    let req = if body.is_some() {
        format!(
            "{method} {path} HTTP/1.1\r\nHost: 127.0.0.1\r\nContent-Type: {content_type}\r\nContent-Length: {}\r\nConnection: close\r\n\r\n{payload}",
            payload.len()
        )
    } else {
        format!("{method} {path} HTTP/1.1\r\nHost: 127.0.0.1\r\nConnection: close\r\n\r\n")
    };
    stream.write_all(req.as_bytes()).unwrap();
    stream.flush().unwrap();

    let mut raw = Vec::new();
    stream.read_to_end(&mut raw).unwrap();
    let text = String::from_utf8_lossy(&raw).into_owned();
    let (head, rest) = text
        .split_once("\r\n\r\n")
        .expect("HTTP response should contain header/body separator");

    let mut lines = head.lines();
    let status_line = lines.next().expect("status line");
    let status: u16 = status_line
        .split_whitespace()
        .nth(1)
        .expect("status code")
        .parse()
        .expect("numeric status");

    let mut content_type = String::new();
    let mut is_chunked = false;
    for line in lines {
        if let Some((k, v)) = line.split_once(':') {
            if k.eq_ignore_ascii_case("content-type") {
                content_type = v.trim().to_string();
            } else if k.eq_ignore_ascii_case("transfer-encoding")
                && v.to_ascii_lowercase().contains("chunked")
            {
                is_chunked = true;
            }
        }
    }

    let decoded_body = if is_chunked {
        decode_chunked(rest)
    } else {
        rest.to_string()
    };

    HttpResponse {
        status,
        content_type,
        body: decoded_body,
    }
}

fn decode_chunked(s: &str) -> String {
    let mut out = String::new();
    let mut cur = s;
    while let Some((len_hex, after)) = cur.split_once("\r\n") {
        let Ok(len) = usize::from_str_radix(len_hex.trim(), 16) else {
            break;
        };
        if len == 0 {
            break;
        }
        if after.len() < len {
            break;
        }
        out.push_str(&after[..len]);
        cur = after[len..].strip_prefix("\r\n").unwrap_or(&after[len..]);
    }
    out
}

fn url_encode_path(p: &Path) -> String {
    let s = p.to_string_lossy();
    let mut out = String::with_capacity(s.len());
    for b in s.bytes() {
        match b {
            b'A'..=b'Z' | b'a'..=b'z' | b'0'..=b'9' | b'-' | b'_' | b'.' | b'~' | b'/' => {
                out.push(b as char);
            }
            _ => out.push_str(&format!("%{:02X}", b)),
        }
    }
    out
}

// ---------------------------------------------------------------------------
// Synthetic test FITS fixture builders
// ---------------------------------------------------------------------------

fn write_test_lightcurve(path: &Path) {
    let meta = OgipMeta {
        instrume: "CZTI".into(),
        detnam: Some("Q0".into()),
        object: Some("CRAB".into()),
        ra_obj: Some(83.6331),
        dec_obj: Some(22.0145),
        tstart: 190_000_000.0,
        tstop: 190_000_008.0,
        exposure: 8.0,
        ..Default::default()
    };
    let p0 = Hdu::empty_primary();
    let times = vec![
        190_000_000.5,
        190_000_001.5,
        190_000_002.5,
        190_000_003.5,
        190_000_004.5,
        190_000_005.5,
        190_000_006.5,
        190_000_007.5,
    ];
    let rates = vec![10.0, 14.0, 18.0, 22.0, 20.0, 24.0, 28.0, 32.0];
    let errs = vec![1.0, 2.0, 3.0, 4.0, 2.0, 2.0, 2.0, 2.0];
    let frac = vec![1.0; 8];
    let lc = build_lightcurve_hdu("RATE", &times, &rates, &errs, Some(&frac), 1.0, &meta).unwrap();
    let gti = build_gti_hdu("GTI", &[190_000_000.0], &[190_000_008.0], &meta).unwrap();
    write_fits(path, &mut [p0, lc, gti], true).unwrap();
}

fn write_test_spectrum(path: &Path) {
    let meta = OgipMeta {
        instrume: "CZTI".into(),
        detnam: Some("Q0".into()),
        object: Some("CRAB".into()),
        tstart: 190_000_000.0,
        tstop: 190_001_000.0,
        exposure: 1000.0,
        ..Default::default()
    };
    let p0 = Hdu::empty_primary();
    let spec_meta = SpectrumMeta {
        detchans: 6,
        tlmin1: 0,
        tlmax1: 5,
        ..Default::default()
    };
    let spec = build_spectrum_hdu(
        "SPECTRUM",
        &[0, 1, 2, 3, 4, 5],
        SpectrumData::Rate(vec![1.0, 2.5, 4.0, 3.0, 1.5, 0.5]),
        Some(&[0.1, 0.2, 0.3, 0.2, 0.1, 0.05]),
        Some(&[5, 0, 0, 0, 0, 1]), // channels 0 and 5 are bad (QUALITY != 0)
        None,
        &meta,
        &spec_meta,
    )
    .unwrap();
    let rmf = ResponseMatrixData {
        chantype: "PI".into(),
        detchans: 6,
        tlmin4: 0,
        energ_lo: vec![20.0, 40.0],
        energ_hi: vec![40.0, 60.0],
        n_grp: vec![1, 1],
        f_chan: vec![vec![0], vec![3]],
        n_chan: vec![vec![3], vec![3]],
        matrix: vec![vec![0.5, 0.3, 0.2], vec![0.4, 0.4, 0.2]],
        channels: vec![0, 1, 2, 3, 4, 5],
        e_min: vec![20.0, 25.0, 30.0, 35.0, 40.0, 45.0],
        e_max: vec![25.0, 30.0, 35.0, 40.0, 45.0, 50.0],
        lo_thresh: 1e-6,
    };
    let (_mat, eb) = build_response_hdus(&rmf, &meta).unwrap();
    write_fits(path, &mut [p0, spec, eb], true).unwrap();
}

fn write_test_image(path: &Path) {
    let mut pixels = vec![10.0f32; 16 * 12];
    // Plant sharp peak at (x = 9, y = 4) -> idx = 4 * 16 + 9 = 73
    pixels[4 * 16 + 9] = 250.0;
    let wcs = WcsParams {
        crpix1: 8.5,
        crpix2: 6.5,
        crval1: 83.6331,
        crval2: 22.0145,
        cdelt1: -0.1,
        cdelt2: 0.1,
        ..Default::default()
    };
    let h0 = ImageBuilder::from_f32(None, vec![16, 12], &pixels)
        .unwrap()
        .with_wcs(&wcs)
        .header_str("INSTRUME", "CZTI", "instrument")
        .header_str("OBJECT", "CRAB", "target")
        .build_primary();
    write_fits(path, &mut [h0], true).unwrap();
}

/// Write a deterministic synthetic event FITS with two energy bands:
/// - Band 1 (`ENERGY = 15.0 keV`): double-peak Crab-like pulse centered at phase `0.25`
/// - Band 2 (`ENERGY = 65.0 keV`): identical pulse shifted LATER by `+shift_cycles`!
fn write_planted_multiband_events(path: &Path, freq_hz: f64, epoch: f64, shift_cycles: f64) {
    let mut times = Vec::new();
    let mut energies = Vec::new();
    let mut pis = Vec::new();
    let mut weights = Vec::new();

    // Deterministic synthetic event generator over 1000 cycles
    let n_cycles = 200usize;
    let period = 1.0 / freq_hz;
    let mut push_band_phase = |cycle: usize, phi: f64, energy: f32, pi: i16, weight: f64| {
        let t = epoch + ((cycle as f64) + phi.rem_euclid(1.0)) * period;
        times.push(t);
        energies.push(energy);
        pis.push(pi);
        weights.push(weight);
    };

    // Exact planted double-peak template (main peak at 0.25, interpulse at 0.65)
    let template_phases: &[(f64, usize)] = &[
        (0.05, 1),
        (0.15, 2),
        (0.24, 8),
        (0.25, 14),
        (0.26, 8),
        (0.35, 2),
        (0.45, 1),
        (0.55, 1),
        (0.64, 5),
        (0.65, 8),
        (0.66, 5),
        (0.75, 1),
        (0.85, 1),
        (0.95, 1),
    ];

    for c in 0..n_cycles {
        for &(phi, mult) in template_phases {
            for k in 0..mult {
                let jitter = ((k as f64) - (mult as f64) * 0.5) * 1e-4;
                push_band_phase(c, phi + jitter, 15.0, 20, 1.0);
                push_band_phase(c, phi + shift_cycles + jitter, 65.0, 120, 0.5);
            }
        }
    }

    let p0 = Hdu::empty_primary();
    let evt = BinTableBuilder::new("EVENTS")
        .add_f64("TIME", Some("s"), times)
        .add_f32("ENERGY", Some("keV"), energies)
        .add_i16("PI", Some("chan"), pis)
        .add_f64("WEIGHT", None, weights)
        .header_str("INSTRUME", "CZTI", "instrument")
        .header_str("TIMESYS", "TDB", "barycentric time")
        .header_f64("EXPOSURE", (n_cycles as f64) * period, "exposure")
        .build()
        .unwrap();
    write_fits(path, &mut [p0, evt], true).unwrap();
}

// ---------------------------------------------------------------------------
// Unit tests: exact mathematics
// ---------------------------------------------------------------------------

#[test]
fn rebin_lightcurve_exact_mean_and_quadrature_errors() {
    let t = vec![0.5, 1.5, 2.5, 3.5, 4.5, 5.5];
    let r = vec![10.0, 14.0, 18.0, 22.0, 30.0, 50.0];
    let e = vec![1.0, 2.0, 3.0, 4.0, 3.0, 4.0];
    let f = vec![1.0; 6];

    // Uniform rebin by 4: group 0 = [0..4] (4 bins), group 1 = [4..6] (2 bins)
    let (rt, rr, re, rf) = rebin_lightcurve(&t, &r, &e, &f, 4, true);
    assert_eq!(rt.len(), 2);
    assert!((rt[0] - 2.0).abs() < 1e-12);
    assert!((rr[0] - 16.0).abs() < 1e-12);
    // Exact quadrature error: sqrt(1^2 + 2^2 + 3^2 + 4^2) / 4 = sqrt(30) / 4
    assert!((re[0] - (30.0f64.sqrt() / 4.0)).abs() < 1e-12);
    assert!((rf[0] - 1.0).abs() < 1e-12);

    assert!((rt[1] - 5.0).abs() < 1e-12);
    assert!((rr[1] - 40.0).abs() < 1e-12);
    // sqrt(3^2 + 4^2) / 2 = 2.5
    assert!((re[1] - 2.5).abs() < 1e-12);

    // Non-uniform fracexp weighting: [10.0 (w=0.25), 20.0 (w=0.75)]
    let (wt, wr, we, wf) = rebin_lightcurve(
        &[0.5, 1.5],
        &[10.0, 20.0],
        &[4.0, 2.0],
        &[0.25, 0.75],
        2,
        true,
    );
    assert_eq!(wt.len(), 1);
    assert!((wr[0] - 17.5).abs() < 1e-12);
    // sqrt((0.25*4)^2 + (0.75*2)^2) / 1.0 = sqrt(1 + 2.25) = sqrt(3.25)
    assert!((we[0] - 3.25f64.sqrt()).abs() < 1e-12);
    assert!((wf[0] - 0.5).abs() < 1e-12);
}

#[test]
fn circular_xcorr_recovers_subbin_planted_phase_lag() {
    let nbins = 64usize;
    let true_lag = 0.1234_f64;
    let make_profile = |shift: f64| -> Vec<f64> {
        (0..nbins)
            .map(|i| {
                let phi = ((i as f64) + 0.5) / (nbins as f64) - shift;
                100.0 + 35.0 * (TAU * phi).cos() + 18.0 * (2.0 * TAU * (phi - 0.1)).cos()
            })
            .collect()
    };
    let a = make_profile(0.0);
    let b = make_profile(true_lag);
    let (lag, corr) = circular_xcorr_phase_lag(&a, &b);
    assert!(
        (lag - true_lag).abs() < 5e-4,
        "expected lag {true_lag}, got {lag} (corr={corr})"
    );
    assert!(corr > 0.99);

    // Negative lag test (-0.215 cycles)
    let c = make_profile(-0.215);
    let (neg_lag, neg_corr) = circular_xcorr_phase_lag(&a, &c);
    assert!(
        (neg_lag - (-0.215)).abs() < 5e-4,
        "expected lag -0.215, got {neg_lag} (corr={neg_corr})"
    );
}

#[test]
fn multiband_folding_recovers_exact_counts_weights_and_phase_lag() {
    let dir = tempdir().unwrap();
    let evt_path = dir.path().join("planted.evt");
    let freq_hz = 29.658_f64;
    let epoch = 190_000_000.0_f64;
    let planted_shift = 0.125_f64; // exactly +8 bins out of 64!
    write_planted_multiband_events(&evt_path, freq_hz, epoch, planted_shift);

    // 1. Unweighted folding
    let req = MultibandRequest {
        inputs: vec![MultibandInputSpec {
            label: Some("Synth".into()),
            path: evt_path.to_string_lossy().into_owned(),
            hdu: None,
            caldb: None,
            use_weights: Some(false),
            bands: vec![
                BandSpec {
                    label: "Low 10-30 keV".into(),
                    emin_kev: Some(10.0),
                    emax_kev: Some(30.0),
                    pi_min: None,
                    pi_max: None,
                    chan_min: None,
                    chan_max: None,
                    unit: None,
                    layer: None,
                },
                BandSpec {
                    label: "High 40-100 keV".into(),
                    emin_kev: Some(40.0),
                    emax_kev: Some(100.0),
                    pi_min: None,
                    pi_max: None,
                    chan_min: None,
                    chan_max: None,
                    unit: None,
                    layer: None,
                },
            ],
        }],
        freq_hz,
        fdot: Some(0.0),
        fddot: Some(0.0),
        epoch: Some(epoch),
        nbins: Some(64),
        harmonics: Some(2),
    };

    let res = handle_pulse_multiband(&req, dir.path()).unwrap();
    assert_eq!(res.bands.len(), 2);
    let b0 = &res.bands[0];
    let b1 = &res.bands[1];
    assert_eq!(b0.num_events, b1.num_events);
    assert!(b0.num_events > 10_000);
    // Exact Poisson error check
    for (&c, &e) in b0.counts.iter().zip(b0.errors.iter()) {
        assert!((e - c.sqrt()).abs() < 1e-12);
    }
    assert!(b0.z2_n > 5_000.0);
    assert!((b0.phase_lag_cycles).abs() < 1e-12);
    assert!(
        (b1.phase_lag_cycles - planted_shift).abs() < 1e-3,
        "expected +{planted_shift} cycles, got {}",
        b1.phase_lag_cycles
    );
    let expected_ms = planted_shift * (1.0 / freq_hz) * 1000.0;
    assert!((b1.phase_lag_ms - expected_ms).abs() < 0.05);

    // 2. Mask-weighted folding (Band 2 has weight = 0.5 -> total_weight = 0.5 * N, error = 0.5 * sqrt(N_b))
    let mut wreq = req;
    wreq.inputs[0].use_weights = Some(true);
    let wres = handle_pulse_multiband(&wreq, dir.path()).unwrap();
    let wb1 = &wres.bands[1];
    assert!(wb1.weighted);
    assert!((wb1.total_weight - 0.5 * (wb1.num_events as f64)).abs() < 1e-9);
    for (bin_idx, &w_err) in wb1.errors.iter().enumerate() {
        let expected_err = 0.5 * b1.counts[bin_idx].sqrt();
        assert!((w_err - expected_err).abs() < 1e-10);
    }
}

// ---------------------------------------------------------------------------
// Full End-to-End Live HTTP Integration Test Suite
// ---------------------------------------------------------------------------

#[test]
fn live_http_endpoints_and_adversarial_error_handling() {
    let dir = tempdir().unwrap();
    let lc_path = dir.path().join("crab_q0.lc");
    let spec_path = dir.path().join("crab_q0.pha");
    let img_path = dir.path().join("crab.img_Q0");
    let evt_path = dir.path().join("crab_bary.evt");
    let bad_fits = dir.path().join("corrupt.fits");

    write_test_lightcurve(&lc_path);
    write_test_spectrum(&spec_path);
    write_test_image(&img_path);
    write_planted_multiband_events(&evt_path, 29.658, 190_000_000.0, 0.125);
    std::fs::write(&bad_fits, b"NOT_A_FITS_FILE_GARBAGE").unwrap();

    // Also create a UTC event + orbit file for real background bary job test
    let t0 = mjd_utc_to_astrosat_utc(57394.0);
    let orb_path = dir.path().join("test.orb");
    let utc_evt = dir.path().join("utc_in.evt");
    let bary_out = dir.path().join("bary_job_out.evt");
    {
        let meta = OgipMeta {
            instrume: "LAXPC".into(),
            tstart: t0,
            tstop: t0 + 100.0,
            exposure: 100.0,
            ..Default::default()
        };
        let p0 = Hdu::empty_primary();
        let orb_t: Vec<f64> = (0..12).map(|i| t0 + (i as f64) * 10.0).collect();
        let orb_pos = vec![[7000.0, 0.0, 0.0]; 12];
        let orb_vel = vec![[0.0, 7.5, 0.0]; 12];
        let orb_hdu = build_orbit_hdu(&orb_t, &orb_pos, &orb_vel, &meta).unwrap();
        write_fits(&orb_path, &mut [p0.clone(), orb_hdu], true).unwrap();

        let evt_hdu = BinTableBuilder::new("EVENTS")
            .add_f64("TIME", Some("s"), vec![t0 + 10.0, t0 + 20.0, t0 + 30.0])
            .add_i16("CHANNEL", None, vec![10, 20, 30])
            .header_str("TIMESYS", "UTC", "")
            .header_i64("MJDREFI", 55197, "")
            .header_f64("MJDREFF", 0.0, "")
            .header_f64("TSTART", t0 + 10.0, "")
            .header_f64("TSTOP", t0 + 30.0, "")
            .build()
            .unwrap();
        write_fits(&utc_evt, &mut [p0, evt_hdu], true).unwrap();
    }

    let server = UiServer::spawn(UiConfig {
        host: "127.0.0.1".into(),
        port: 0,
        data_dir: dir.path().to_path_buf(),
        open_browser: false,
        threads: Some(2),
    })
    .expect("spawn UiServer on ephemeral port");
    let addr = server.addr();

    // 1. Static assets (Studio UI default at `/` and `/studio`; Classic UI at `/classic`)
    for (route, ct) in [
        ("/", "text/html"),
        ("/index.html", "text/html"),
        ("/studio", "text/html"),
        ("/studio/", "text/html"),
        ("/studio/index.html", "text/html"),
        ("/studio/studio.css", "text/css"),
        ("/studio/studio.js", "application/javascript"),
        ("/classic", "text/html"),
        ("/classic/", "text/html"),
        ("/classic/index.html", "text/html"),
        ("/styles.css", "text/css"),
        ("/classic/styles.css", "text/css"),
        ("/charts.js", "application/javascript"),
        ("/classic/charts.js", "application/javascript"),
        ("/app.js", "application/javascript"),
        ("/classic/app.js", "application/javascript"),
    ] {
        let r = http_request(addr, "GET", route, None);
        assert_eq!(r.status, 200, "static route {route} failed");
        assert!(
            r.content_type.contains(ct),
            "route {route} expected {ct}, got {}",
            r.content_type
        );
        let body_no_svg_ns = r.body.replace("http://www.w3.org/2000/svg", "");
        assert!(
            !body_no_svg_ns.contains("https://") && !body_no_svg_ns.contains("http://"),
            "route {route} must remain 100% offline with zero external URLs"
        );
    }
    // Verify `/` serves Studio UI (matching `/studio`) and `/classic` serves Classic UI
    let root_html = http_request(addr, "GET", "/", None).body;
    let studio_html = http_request(addr, "GET", "/studio", None).body;
    let classic_html = http_request(addr, "GET", "/classic", None).body;
    assert!(root_html.contains("/studio/studio.js"));
    assert!(root_html.contains("AstroSat Studio"));
    assert_eq!(root_html, studio_html);
    assert!(classic_html.contains("/app.js"));
    assert!(!classic_html.contains("AstroSat Studio"));

    // 2. Health & Browse
    let r = http_request(addr, "GET", "/api/health", None);
    assert_eq!(r.status, 200);
    let health: Value = serde_json::from_str(&r.body).unwrap();
    assert_eq!(health["status"], "ok");

    let r = http_request(addr, "GET", "/api/browse", None);
    assert_eq!(r.status, 200);
    let browse: Value = serde_json::from_str(&r.body).unwrap();
    let entries = browse["entries"].as_array().unwrap();
    assert!(entries.iter().any(|e| e["name"] == "crab_q0.lc" && e["product"] == "lightcurve"));
    assert!(entries.iter().any(|e| e["name"] == "crab_q0.pha" && e["product"] == "spectrum"));
    assert!(entries.iter().any(|e| e["name"] == "crab.img_Q0" && e["product"] == "image"));

    // 3. Inspect & Validate
    let r = http_request(
        addr,
        "GET",
        &format!("/api/inspect?path={}&full=true&hdu=1&limit=4", url_encode_path(&lc_path)),
        None,
    );
    assert_eq!(r.status, 200);
    let insp: Value = serde_json::from_str(&r.body).unwrap();
    assert_eq!(insp["instrume"], "CZTI");
    assert_eq!(insp["health"]["all_checksums_ok"], true);
    assert_eq!(insp["selected_hdu"]["table"]["rows"].as_array().unwrap().len(), 4);

    let r = http_request(
        addr,
        "GET",
        &format!("/api/validate?path={}&strict=true", url_encode_path(&lc_path)),
        None,
    );
    assert_eq!(r.status, 200);
    let val: Value = serde_json::from_str(&r.body).unwrap();
    assert_eq!(val["valid"], true);
    assert_eq!(val["errors"], 0);
    assert_eq!(val["warnings"], 0);

    // 4. Lightcurve with rebin=4
    let r = http_request(
        addr,
        "GET",
        &format!("/api/lightcurve?path={}&rebin=4", url_encode_path(&lc_path)),
        None,
    );
    assert_eq!(r.status, 200);
    let lc: Value = serde_json::from_str(&r.body).unwrap();
    assert_eq!(lc["effective_rebin"], 4);
    assert_eq!(lc["rate"].as_array().unwrap().len(), 2);
    assert!((lc["rate"][0].as_f64().unwrap() - 16.0).abs() < 1e-9);
    assert!((lc["error"][0].as_f64().unwrap() - (30.0f64.sqrt() / 4.0)).abs() < 1e-9);

    // 5. Spectrum with good_only=true vs false
    let r = http_request(
        addr,
        "GET",
        &format!("/api/spectrum?path={}&good_only=true", url_encode_path(&spec_path)),
        None,
    );
    assert_eq!(r.status, 200);
    let sp: Value = serde_json::from_str(&r.body).unwrap();
    assert_eq!(sp["total_channels"], 6);
    assert_eq!(sp["good_channels"], 4);
    assert_eq!(sp["channel"].as_array().unwrap().len(), 4);
    assert_eq!(sp["energy_mapped"], true);
    assert!((sp["energy_kev"][0].as_f64().unwrap() - 27.5).abs() < 1e-9);

    // 6. 2D Image & WCS
    let r = http_request(
        addr,
        "GET",
        &format!("/api/image?path={}", url_encode_path(&img_path)),
        None,
    );
    assert_eq!(r.status, 200);
    let im: Value = serde_json::from_str(&r.body).unwrap();
    assert_eq!(im["width"], 16);
    assert_eq!(im["height"], 12);
    assert_eq!(im["stats"]["peak"]["x"], 9);
    assert_eq!(im["stats"]["peak"]["y"], 4);
    assert!((im["stats"]["peak"]["val"].as_f64().unwrap() - 250.0).abs() < 1e-6);
    assert_eq!(im["wcs"]["has_wcs"], true);
    assert!((im["wcs"]["crval1"].as_f64().unwrap() - 83.6331).abs() < 1e-6);

    // 7. Pulse Search & Multi-band Folding via HTTP POST
    let search_body = json!({
        "path": evt_path.to_string_lossy(),
        "freq_hz": 29.658,
        "search_span_hz": 0.01,
        "n_trials": 41,
        "harmonics": 2,
        "nbins": 32,
        "pi_min": 0,
        "pi_max": 50
    })
    .to_string();
    let r = http_request(addr, "POST", "/api/pulse/search", Some(&search_body));
    assert_eq!(r.status, 200, "pulse search failed: {}", r.body);
    let ps: Value = serde_json::from_str(&r.body).unwrap();
    assert!((ps["best_freq_hz"].as_f64().unwrap() - 29.658).abs() < 0.002);
    assert!(ps["max_z2_n"].as_f64().unwrap() > 1000.0);

    let mb_body = json!({
        "inputs": [{
            "label": "CZTI",
            "path": evt_path.to_string_lossy(),
            "bands": [
                { "label": "10-30 keV", "emin_kev": 10.0, "emax_kev": 30.0 },
                { "label": "40-100 keV", "emin_kev": 40.0, "emax_kev": 100.0 }
            ]
        }],
        "freq_hz": 29.658,
        "epoch": 190_000_000.0,
        "nbins": 64,
        "harmonics": 2
    })
    .to_string();
    let r = http_request(addr, "POST", "/api/pulse/multiband", Some(&mb_body));
    assert_eq!(r.status, 200, "multiband failed: {}", r.body);
    let mb: Value = serde_json::from_str(&r.body).unwrap();
    assert!((mb["bands"][1]["phase_lag_cycles"].as_f64().unwrap() - 0.125).abs() < 1e-3);

    // 8. CLI Preview & Real Background Bary Job Execution
    let job_body = json!({
        "pipeline": "bary",
        "threads": 2,
        "bary": {
            "input": utc_evt.to_string_lossy(),
            "output": bary_out.to_string_lossy(),
            "orbit": orb_path.to_string_lossy(),
            "ra_deg": 83.6331,
            "dec_deg": 22.0145,
            "frame": "icrs",
            "ephemeris": "analytical"
        }
    })
    .to_string();

    let r = http_request(addr, "POST", "/api/cli-preview", Some(&job_body));
    assert_eq!(r.status, 200);
    let prev: Value = serde_json::from_str(&r.body).unwrap();
    assert!(prev["cli_command"].as_str().unwrap().contains("astrosat --threads 2 bary"));
    assert!(prev["toml_config"].as_str().unwrap().contains("[bary]"));

    // Verify CZTI CLI preview uses valid CztiRunArgs flags (--tbin, --resp / --no-resp)
    let czti_preview_body = json!({
        "pipeline": "czti",
        "czti": {
            "input_dir": "/data/czti_l1",
            "caldb_dir": "/data/caldb/czti",
            "output_dir": "/data/out",
            "time_binsize": 2.0,
            "gen_resp": true,
            "mask_weight": true
        }
    })
    .to_string();
    let r = http_request(addr, "POST", "/api/cli-preview", Some(&czti_preview_body));
    assert_eq!(r.status, 200);
    let cp: Value = serde_json::from_str(&r.body).unwrap();
    let czti_cmd = cp["cli_command"].as_str().unwrap();
    assert!(czti_cmd.contains("--tbin 2"), "got: {czti_cmd}");
    assert!(czti_cmd.contains("--resp"), "got: {czti_cmd}");
    assert!(!czti_cmd.contains("--time-binsize"), "got: {czti_cmd}");
    assert!(!czti_cmd.contains("--gen-resp"), "got: {czti_cmd}");

    let r = http_request(addr, "POST", "/api/jobs/run", Some(&job_body));
    assert_eq!(r.status, 202, "job spawn failed: {}", r.body);
    let spawned: Value = serde_json::from_str(&r.body).unwrap();
    let job_id = spawned["id"].as_str().unwrap().to_string();

    // Poll until succeeded
    let mut final_job = Value::Null;
    for _ in 0..100 {
        let poll = http_request(addr, "GET", &format!("/api/jobs/{job_id}"), None);
        assert_eq!(poll.status, 200);
        let j: Value = serde_json::from_str(&poll.body).unwrap();
        if j["status"] != "running" {
            final_job = j;
            break;
        }
        std::thread::sleep(Duration::from_millis(25));
    }
    assert_eq!(final_job["status"], "succeeded", "job did not succeed: {final_job}");
    assert!(bary_out.is_file());
    let bary_fits = FitsFile::open(&bary_out).unwrap();
    assert_eq!(bary_fits.hdu(0).unwrap().header().get_str("TIMESYS"), Some("TDB"));

    let list_r = http_request(addr, "GET", "/api/jobs", None);
    assert_eq!(list_r.status, 200);
    let list_v: Value = serde_json::from_str(&list_r.body).unwrap();
    assert!(!list_v["jobs"].as_array().unwrap().is_empty());

    // 9. Adversarial & Security Checks (all must return 4xx JSON, zero panics)
    let adversarial_gets = [
        "/api/browse?path=../../etc/passwd",
        "/api/browse?path=%2e%2e/%2e%2e/etc",
        "/api/inspect?path=/etc/passwd",
        "/api/inspect?path=does_not_exist.fits",
        &format!("/api/inspect?path={}", url_encode_path(&bad_fits)),
        &format!("/api/lightcurve?path={}&rebin=0", url_encode_path(&lc_path)),
        "/api/unknown_endpoint",
    ];
    for url in adversarial_gets {
        let res = http_request(addr, "GET", url, None);
        assert!(
            (400..500).contains(&res.status),
            "expected 4xx for {url}, got {} ({})",
            res.status,
            res.body
        );
        let v: Value = serde_json::from_str(&res.body).expect("error must be valid JSON");
        assert!(v["error"].is_string());
    }

    let adversarial_posts = [
        ("/api/pulse/search", "{not valid json"),
        (
            "/api/pulse/search",
            &json!({"path": evt_path.to_string_lossy(), "freq_hz": -5.0, "nbins": 0}).to_string(),
        ),
        (
            "/api/pulse/multiband",
            &json!({"inputs": [], "freq_hz": 0.0, "nbins": 0}).to_string(),
        ),
        ("/api/jobs/run", &json!({"pipeline": "unknown"}).to_string()),
    ];
    for (url, payload) in adversarial_posts {
        let res = http_request(addr, "POST", url, Some(payload));
        assert!(
            (400..500).contains(&res.status),
            "expected 4xx for POST {url}, got {} ({})",
            res.status,
            res.body
        );
        let v: Value = serde_json::from_str(&res.body).expect("error must be valid JSON");
        assert!(v["error"].is_string());
    }

    // Confirm server is still healthy after all adversarial inputs
    let r = http_request(addr, "GET", "/api/health", None);
    assert_eq!(r.status, 200);
    server.stop();
}

fn test_data_root() -> Option<std::path::PathBuf> {
    std::env::var("ASTROSAT_TEST_DATA")
        .ok()
        .map(std::path::PathBuf::from)
        .filter(|p| p.is_dir())
}

#[test]
fn real_astrosat_datasets_ui_smoke() {
    let Some(data_root_buf) = test_data_root() else {
        return;
    };
    let data_root = data_root_buf.as_path();

    let t0 = std::time::Instant::now();
    let health = crate::browse::handle_health(data_root, 4);
    let detect_ms = t0.elapsed().as_millis();
    eprintln!(
        "real_astrosat_datasets_ui_smoke: detect_presets={}ms, lc={}, spec={}, img={}, evt={}",
        detect_ms,
        health.presets.sample_lightcurves.len(),
        health.presets.sample_spectra.len(),
        health.presets.sample_images.len(),
        health.presets.sample_events.len()
    );
    assert_eq!(health.status, "ok");
    assert!(health.presets.czti_l1_dir.is_some());
    assert!(health.presets.laxpc_l1_dir.is_some());
    assert!(health.presets.czti_caldb.is_some());
    assert!(health.presets.jpleph_405.is_some());
    assert!(!health.presets.sample_lightcurves.is_empty());
    assert!(!health.presets.sample_spectra.is_empty());
    assert!(!health.presets.sample_images.is_empty());

    // Verify LAXPC channel map auto-discovery works from data_root
    let cmap = crate::science::find_laxpc_channel_map(
        astrosat_laxpc::types::LaxpcUnit::Lxp10,
        100.0,
        None,
        data_root,
    );
    assert!(cmap.is_some(), "should auto-discover lxp10.channel under ASTROSAT_TEST_DATA");

    // Verify sample CZTI L2 image + spectrum + lightcurve if discovered
    if let Some(img_path) = health.presets.sample_images.first() {
        let mut p = std::collections::HashMap::new();
        p.insert("path".to_string(), img_path.clone());
        let img = crate::science::handle_image(&p, data_root).expect("read sample image");
        assert!(img.width > 0 && img.height > 0);
        assert_eq!(img.pixels.len(), img.width * img.height);
    }
    if let Some(lc_path) = health.presets.sample_lightcurves.first() {
        let mut p = std::collections::HashMap::new();
        p.insert("path".to_string(), lc_path.clone());
        p.insert("rebin".to_string(), "4".to_string());
        let lc = crate::science::handle_lightcurve(&p, data_root).expect("read sample lc");
        assert!(!lc.time.is_empty());
    }
    if let Some(spec_path) = health.presets.sample_spectra.first() {
        let mut p = std::collections::HashMap::new();
        p.insert("path".to_string(), spec_path.clone());
        let sp = crate::science::handle_spectrum(&p, data_root).expect("read sample spectrum");
        assert!(sp.total_channels > 0);
        assert!(sp.energy_mapped);
    }

    // Verify czti-addspec and laxpc-backshift against real ASTROSAT_TEST_DATA files into tempdir()
    let czti_mode_m0 = data_root.join(
        "sample/czti_2018/20181014_C04_007T02_9000002432_level2/czti/modeM0",
    );
    let q0_pha = czti_mode_m0.join("AS1C04_007T02_9000002432cztM0_level2_quad_clean_Q0.pha");
    let q1_pha = czti_mode_m0.join("AS1C04_007T02_9000002432cztM0_level2_quad_clean_Q1.pha");
    let q0_rsp = czti_mode_m0.join("AS1C04_007T02_9000002432cztM0_level2_quad_clean_Q0.rsp");
    let q1_rsp = czti_mode_m0.join("AS1C04_007T02_9000002432cztM0_level2_quad_clean_Q1.rsp");
    let laxpc_ref = data_root.join("ref_src/laxpc");
    let laxpc_bkg_pha = laxpc_ref.join("backlxp10jan16.pha");

    if q0_pha.is_file() && q1_pha.is_file() && q0_rsp.is_file() && q1_rsp.is_file() && laxpc_bkg_pha.is_file() {
        let tmp = tempdir().unwrap();
        let out_pha = tmp.path().join("coadd_q01.pha");
        let out_rsp = tmp.path().join("coadd_q01.rsp");
        let out_lxp_pha = tmp.path().join("shifted_bkg.pha");

        let srv = UiServer::spawn(UiConfig {
            host: "127.0.0.1".into(),
            port: 0,
            data_dir: data_root.to_path_buf(),
            open_browser: false,
            threads: Some(2),
        })
        .unwrap();
        let addr = srv.addr();

        let addspec_req = json!({
            "pipeline": "czti-addspec",
            "czti_addspec": {
                "pha": [q0_pha.to_string_lossy(), q1_pha.to_string_lossy()],
                "rsp": [q0_rsp.to_string_lossy(), q1_rsp.to_string_lossy()],
                "out_pha": out_pha.to_string_lossy(),
                "out_rsp": out_rsp.to_string_lossy(),
                "add_exposure": false
            }
        })
        .to_string();
        let r = http_request(addr, "POST", "/api/jobs/run", Some(&addspec_req));
        assert_eq!(r.status, 202, "real czti-addspec spawn failed: {}", r.body);
        let jid = serde_json::from_str::<Value>(&r.body).unwrap()["id"]
            .as_str()
            .unwrap()
            .to_string();
        let done = poll_job_until_done(addr, &jid);
        assert_eq!(done["status"], "succeeded", "real czti-addspec failed: {done}");
        assert!(out_pha.is_file());
        assert!(out_rsp.is_file());

        let backshift_req = json!({
            "pipeline": "laxpc-backshift",
            "laxpc_backshift": {
                "unit": "10",
                "caldb_dir": laxpc_ref.to_string_lossy(),
                "mission_day": 250.0,
                "bkg_pha": laxpc_bkg_pha.to_string_lossy(),
                "output": out_lxp_pha.to_string_lossy()
            }
        })
        .to_string();
        let r = http_request(addr, "POST", "/api/jobs/run", Some(&backshift_req));
        assert_eq!(r.status, 202, "real laxpc-backshift spawn failed: {}", r.body);
        let jid = serde_json::from_str::<Value>(&r.body).unwrap()["id"]
            .as_str()
            .unwrap()
            .to_string();
        let done = poll_job_until_done(addr, &jid);
        assert_eq!(done["status"], "succeeded", "real laxpc-backshift failed: {done}");
        assert!(out_lxp_pha.is_file());
        assert!(done["summary"]["lcbk_shift_fraction"].is_number());
        srv.stop();
    }
}

fn poll_job_until_done(addr: SocketAddr, job_id: &str) -> Value {
    for _ in 0..200 {
        let poll = http_request(addr, "GET", &format!("/api/jobs/{job_id}"), None);
        assert_eq!(poll.status, 200);
        let j: Value = serde_json::from_str(&poll.body).unwrap();
        if j["status"] != "running" {
            return j;
        }
        std::thread::sleep(Duration::from_millis(25));
    }
    panic!("job {job_id} did not finish within timeout");
}

fn write_test_czti_pha_and_rsp(pha: &Path, rsp: &Path, exposure: f64, counts_scale: f32) {
    let meta = CztiObsMeta {
        obs_id: Some("9000000001".into()),
        object: Some("CRAB".into()),
        ra_obj: Some(83.6331),
        dec_obj: Some(22.0145),
        tstart: 190_000_000.0,
        tstop: 190_000_000.0 + exposure,
        ..Default::default()
    };
    let counts: Vec<f32> = (0..PI_CHANNELS)
        .map(|ch| {
            if (40..=300).contains(&ch) {
                counts_scale * (10.0 + ((ch % 17) as f32))
            } else {
                0.0
            }
        })
        .collect();
    let stat_err: Vec<f32> = counts.iter().map(|&c| c.sqrt()).collect();
    let spec = CztiSpectrum {
        quad_id: Some(0),
        channel: (0..PI_CHANNELS as i32).collect(),
        counts,
        stat_err,
        sys_err: vec![0.0; PI_CHANNELS],
        quality: vec![0; PI_CHANNELS],
        grouping: vec![1; PI_CHANNELS],
        exposure,
        thetax_deg: 0.1,
        thetay_deg: -0.1,
        mask_weighted: true,
        badpix_threshold: 0,
        gti_type: GtiType::Quad,
        respfile: rsp
            .file_name()
            .map(|s| s.to_string_lossy().into_owned())
            .unwrap_or_else(|| "test.rsp".into()),
        tstart: 190_000_000.0,
        tstop: 190_000_000.0 + exposure,
    };
    write_spectrum_file(pha, &meta, &spec).unwrap();

    let n_ebins = 4usize;
    let mut matrix = vec![0.0f32; n_ebins * PI_CHANNELS];
    for eb in 0..n_ebins {
        let ch = 50 + eb * 50;
        matrix[eb * PI_CHANNELS + ch] = 100.0 * counts_scale;
        matrix[eb * PI_CHANNELS + ch + 1] = 50.0 * counts_scale;
    }
    let czti_rsp = CztiResponse {
        quad_id: Some(0),
        n_ebins,
        n_pibins: PI_CHANNELS,
        energ_lo: vec![10.0, 35.0, 60.0, 85.0],
        energ_hi: vec![35.0, 60.0, 85.0, 110.0],
        matrix,
        ebounds_channel: (0..PI_CHANNELS as i16).collect(),
        ebounds_emin: (0..PI_CHANNELS).map(|i| 5.0 + (i as f32) * 0.5).collect(),
        ebounds_emax: (0..PI_CHANNELS)
            .map(|i| 5.5 + (i as f32) * 0.5)
            .collect(),
    };
    write_response_file(rsp, &meta, &czti_rsp).unwrap();
}

fn write_test_laxpc_backshift_caldb(dir: &Path) -> (std::path::PathBuf, std::path::PathBuf) {
    let gaina8 = "\
# day p30_10 f30_10 p60_10 f60_10 p30_20 f30_20 p60_20 f60_20 p30_30 f30_30 p60_30 f60_30
50.0   250.0  10.0  500.0  15.0  220.0  10.0  440.0  15.0  250.0  10.0  500.0  15.0
200.0  230.0  10.0  460.0  15.0  218.0  10.0  436.0  15.0  230.0  10.0  460.0  15.0
";
    std::fs::write(dir.join("gaina8.dat"), gaina8).unwrap();
    std::fs::write(dir.join("lxhv.dat"), "120.0 1\n").unwrap();

    let mut chan_txt = String::new();
    for i in 0..512 {
        let lo = (i as f64) * 0.25;
        let hi = ((i + 1) as f64) * 0.25;
        chan_txt.push_str(&format!("{lo:.4} {hi:.4}\n"));
    }
    std::fs::write(dir.join("lxp10.channel"), chan_txt).unwrap();

    let bkg_pha = dir.join("backlxp10test.pha");
    let mut pha_txt = String::from("# Exposure = 1000.0 lxp1 Channels: 0 511 day: 100.0\n");
    for ch in 0..512 {
        let r = if (10..=200).contains(&ch) { 0.25 } else { 0.01 };
        pha_txt.push_str(&format!(
            "{ch} {r:.4} {r:.4} {r:.4} {r:.4} {r:.4} {r:.4} {r:.4} 0.0 0.0 0.0 {:.4}\n",
            r * 7.0
        ));
    }
    std::fs::write(&bkg_pha, pha_txt).unwrap();

    let fits_path = dir.join("obs_day150.fits");
    let tstart_day150 = (150.0 + 2096.0) * 86_400.0;
    let p0 = Hdu::empty_primary();
    let h1 = BinTableBuilder::new("EVENTS")
        .add_f64("TIME", Some("s"), vec![tstart_day150 + 1.0, tstart_day150 + 2.0])
        .header_f64("TSTART", tstart_day150, "start")
        .header_f64("TSTOP", tstart_day150 + 10.0, "stop")
        .build()
        .unwrap();
    write_fits(&fits_path, &mut [p0, h1], true).unwrap();

    (bkg_pha, fits_path)
}

#[test]
fn m1_config_parse_and_cli_preview_parity() {
    let dir = tempdir().unwrap();
    let srv = UiServer::spawn(UiConfig {
        host: "127.0.0.1".into(),
        port: 0,
        data_dir: dir.path().to_path_buf(),
        open_browser: false,
        threads: Some(2),
    })
    .unwrap();
    let addr = srv.addr();

    // 1. POST /api/config/parse with TOML payload
    let toml_doc = r#"
threads = 4

[czti]
input_dir = "/data/czti_l1"
oversampling = 32
generate_eventfile = true
quads = [0, 2]

[laxpc]
anodes = "layer1"
saa_entry_lon_deg = 260.5
write_ascii_products = true

[bary]
orbit_file = "/data/test.orb"
ra_deg = 83.6331
dec_deg = 22.0145

[timing]
freq_hz = 29.658
fddot = 1.25e-20
nbins = 64
"#;
    let r = http_request(
        addr,
        "POST",
        "/api/config/parse",
        Some(&json!({ "toml": toml_doc }).to_string()),
    );
    assert_eq!(r.status, 200, "config parse failed: {}", r.body);
    let parsed: Value = serde_json::from_str(&r.body).unwrap();
    assert_eq!(parsed["config"]["threads"], 4);
    assert_eq!(parsed["config"]["czti"]["oversampling"], 32);
    assert_eq!(parsed["config"]["czti"]["generate_eventfile"], true);
    assert_eq!(parsed["config"]["czti"]["quads"], json!([0, 2]));
    assert_eq!(parsed["config"]["laxpc"]["anodes"], "layer1");
    assert_eq!(parsed["config"]["laxpc"]["write_ascii_products"], true);
    assert_eq!(parsed["config"]["bary"]["orbit_file"], "/data/test.orb");
    assert!((parsed["config"]["timing"]["fddot"].as_f64().unwrap() - 1.25e-20).abs() < 1e-30);
    assert_eq!(parsed["sections"], json!(["czti", "laxpc", "bary", "timing"]));

    // Re-parse canonical_toml via raw application/toml
    let canon = parsed["canonical_toml"].as_str().unwrap();
    let r2 = http_request_with_ct(addr, "POST", "/api/config/parse", "application/toml", Some(canon));
    assert_eq!(r2.status, 200, "canonical reparse failed: {}", r2.body);
    let reparsed: Value = serde_json::from_str(&r2.body).unwrap();
    assert_eq!(reparsed["config"], parsed["config"]);

    // Content + json / auto format
    let r_json = http_request(
        addr,
        "POST",
        "/api/config/parse",
        Some(
            &json!({
                "content": serde_json::to_string(&json!({"threads": 2, "timing": {"nbins": 32}})).unwrap(),
                "format": "json"
            })
            .to_string(),
        ),
    );
    assert_eq!(r_json.status, 200);
    assert_eq!(
        serde_json::from_str::<Value>(&r_json.body).unwrap()["sections"],
        json!(["timing"])
    );

    // 2. F6 CLI Preview & Round-Trip via /api/config/parse for CZTI, LAXPC, and Bary
    let czti_req = json!({
        "pipeline": "czti",
        "threads": 4,
        "czti": {
            "input_dir": "/data/czti_l1",
            "m0": "/data/m0.fits",
            "ss": "/data/ss.fits",
            "tct": "/data/tct.tct",
            "mkf": "/data/mkf.mkf",
            "caldb_dir": "/data/caldb",
            "output_dir": "/data/out",
            "mkf_threshold_file": "/data/mkfthresh.txt",
            "user_gti_file": "/data/user.gti",
            "output_type": "spec",
            "generate_eventfile": true,
            "oversampling": 64,
            "quads": [0, 2]
        }
    })
    .to_string();
    let r = http_request(addr, "POST", "/api/cli-preview", Some(&czti_req));
    assert_eq!(r.status, 200);
    let cp: Value = serde_json::from_str(&r.body).unwrap();
    let czti_cmd = cp["cli_command"].as_str().unwrap();
    for flag in [
        "--oversampling 64",
        "--mkf-thresholds /data/mkfthresh.txt",
        "--user-gti /data/user.gti",
        "--gen-evt",
        "--output-type spec",
        "--quads 0,2",
    ] {
        assert!(czti_cmd.contains(flag), "missing '{flag}' in {czti_cmd}");
    }
    let czti_round = http_request(
        addr,
        "POST",
        "/api/config/parse",
        Some(&json!({ "toml": cp["toml_config"] }).to_string()),
    );
    assert_eq!(czti_round.status, 200, "CZTI TOML round-trip failed: {}", czti_round.body);
    let czti_cfg = serde_json::from_str::<Value>(&czti_round.body).unwrap();
    assert_eq!(czti_cfg["config"]["czti"]["oversampling"], 64);
    assert_eq!(czti_cfg["config"]["czti"]["generate_eventfile"], true);
    assert_eq!(czti_cfg["config"]["czti"]["l1_m0_file"], "/data/m0.fits");

    let laxpc_req = json!({
        "pipeline": "laxpc",
        "threads": 4,
        "laxpc": {
            "input_dir": "/data/lxp_l1",
            "output_dir": "/data/lxp_out",
            "caldb_dir": "/data/caldb/laxpc",
            "units": ["10", "20"],
            "anodes": "1,2,3",
            "saa_entry_lon_deg": 255.0,
            "saa_exit_lon_deg": 335.0,
            "merge_gap_sec": 5.0,
            "min_duration_sec": 30.0,
            "user_gti_file": "/data/user.gti",
            "user_gti_only": true,
            "write_ascii_products": true,
            "background_model_path": "/data/backfit.lxp10",
            "background_pha_path": "/data/backlxp10.pha",
            "background_peak_ch": 237.3
        }
    })
    .to_string();
    let r = http_request(addr, "POST", "/api/cli-preview", Some(&laxpc_req));
    assert_eq!(r.status, 200);
    let lp: Value = serde_json::from_str(&r.body).unwrap();
    let lxp_cmd = lp["cli_command"].as_str().unwrap();
    for flag in [
        "--anodes '1,2,3'",
        "--saa-entry-lon 255",
        "--saa-exit-lon 335",
        "--merge-gap 5",
        "--min-duration 30",
        "--user-gti /data/user.gti",
        "--user-gti-only",
        "--ascii",
        "--bkg-model /data/backfit.lxp10",
        "--bkg-pha /data/backlxp10.pha",
        "--bkg-peak-ch 237.3",
    ] {
        assert!(lxp_cmd.contains(flag), "missing '{flag}' in {lxp_cmd}");
    }
    let lxp_round = http_request(
        addr,
        "POST",
        "/api/config/parse",
        Some(&json!({ "toml": lp["toml_config"] }).to_string()),
    );
    assert_eq!(lxp_round.status, 200, "LAXPC TOML round-trip failed: {}", lxp_round.body);
    let lxp_cfg = serde_json::from_str::<Value>(&lxp_round.body).unwrap();
    assert_eq!(lxp_cfg["config"]["laxpc"]["anodes"], "1,2,3");
    assert_eq!(lxp_cfg["config"]["laxpc"]["user_gti_only"], true);
    assert_eq!(lxp_cfg["config"]["laxpc"]["write_ascii_products"], true);

    // Bary --mkf only preview produces valid TOML; both orbit + mkf returns 400
    let bary_mkf = json!({
        "pipeline": "bary",
        "bary": {
            "input": "/data/in.evt",
            "output": "/data/out.evt",
            "mkf": "/data/test.mkf",
            "ra_deg": 83.6331,
            "dec_deg": 22.0145
        }
    })
    .to_string();
    let r = http_request(addr, "POST", "/api/cli-preview", Some(&bary_mkf));
    assert_eq!(r.status, 200);
    let bp: Value = serde_json::from_str(&r.body).unwrap();
    assert!(bp["cli_command"].as_str().unwrap().contains("--mkf /data/test.mkf"));
    assert!(!bp["cli_command"].as_str().unwrap().contains("--orbit"));
    let bary_round = http_request(
        addr,
        "POST",
        "/api/config/parse",
        Some(&json!({ "toml": bp["toml_config"] }).to_string()),
    );
    assert_eq!(bary_round.status, 200);

    let bary_both = json!({
        "pipeline": "bary",
        "bary": {
            "input": "/data/in.evt",
            "output": "/data/out.evt",
            "orbit": "/data/test.orb",
            "mkf": "/data/test.mkf"
        }
    })
    .to_string();
    assert_eq!(
        http_request(addr, "POST", "/api/cli-preview", Some(&bary_both)).status,
        400
    );

    // 3. Adversarial /api/config/parse & /api/jobs/run validation
    assert_eq!(http_request(addr, "GET", "/api/config/parse", None).status, 405);
    for bad in [
        json!({ "toml": "[czti]\nbogus_field = 1" }).to_string(),
        json!({ "toml": "[bogus_section]\nx = 1" }).to_string(),
        json!({ "toml": "threads = 0" }).to_string(),
        json!({ "toml": "threads = \"not_a_number\"" }).to_string(),
        json!({ "toml": "[unclosed" }).to_string(),
        json!({ "toml": "a = 1", "content": "b = 2" }).to_string(),
        json!({ "content": "threads = 2", "format": "yaml" }).to_string(),
        "{}".to_string(),
    ] {
        let res = http_request(addr, "POST", "/api/config/parse", Some(&bad));
        assert_eq!(res.status, 400, "expected 400 for {bad}, got {}: {}", res.status, res.body);
    }

    // F6 validation errors on /api/jobs/run
    let dummy_file = dir.path().join("dummy.fits");
    std::fs::write(&dummy_file, b"dummy").unwrap();
    for bad_job in [
        json!({
            "pipeline": "czti",
            "czti": {
                "m0": dummy_file.to_string_lossy(),
                "tct": dummy_file.to_string_lossy(),
                "mkf": dummy_file.to_string_lossy(),
                "caldb_dir": dir.path().to_string_lossy(),
                "output_dir": dir.path().join("out").to_string_lossy(),
                "oversampling": 65
            }
        }),
        json!({
            "pipeline": "czti",
            "czti": {
                "m0": dummy_file.to_string_lossy(),
                "tct": dummy_file.to_string_lossy(),
                "mkf": dummy_file.to_string_lossy(),
                "caldb_dir": dir.path().to_string_lossy(),
                "output_dir": dir.path().join("out").to_string_lossy(),
                "oversampling": 0
            }
        }),
        json!({
            "pipeline": "czti",
            "czti": {
                "m0": dummy_file.to_string_lossy(),
                "tct": dummy_file.to_string_lossy(),
                "mkf": dummy_file.to_string_lossy(),
                "caldb_dir": dir.path().to_string_lossy(),
                "output_dir": dir.path().join("out").to_string_lossy(),
                "ra_deg": 83.6331
            }
        }),
        json!({
            "pipeline": "laxpc",
            "laxpc": {
                "input_dir": dir.path().to_string_lossy(),
                "output_dir": dir.path().join("out").to_string_lossy(),
                "anodes": "99"
            }
        }),
        json!({
            "pipeline": "laxpc",
            "laxpc": {
                "input_dir": dir.path().to_string_lossy(),
                "output_dir": dir.path().join("out").to_string_lossy(),
                "user_gti_only": true
            }
        }),
        json!({
            "pipeline": "laxpc",
            "laxpc": {
                "input_dir": dir.path().to_string_lossy(),
                "output_dir": dir.path().join("out").to_string_lossy(),
                "merge_gap_sec": -1.0
            }
        }),
        json!({
            "pipeline": "laxpc",
            "laxpc": {
                "input_dir": dir.path().to_string_lossy(),
                "output_dir": dir.path().join("out").to_string_lossy(),
                "energy_range": "3:20"
            }
        }),
        json!({
            "pipeline": "bary",
            "bary": {
                "input": dummy_file.to_string_lossy(),
                "output": dir.path().join("bary.evt").to_string_lossy(),
                "orbit": dummy_file.to_string_lossy(),
                "mkf": dummy_file.to_string_lossy()
            }
        }),
    ] {
        let res = http_request(addr, "POST", "/api/jobs/run", Some(&bad_job.to_string()));
        assert_eq!(res.status, 400, "expected 400 for {bad_job}, got {}: {}", res.status, res.body);
    }

    // Confirm JobManager default works without unused warnings
    let _mgr = JobManager::default();
    srv.stop();
}

#[test]
fn m1_czti_addspec_laxpc_backshift_and_pulse_parity() {
    let dir = tempdir().unwrap();

    // Setup CZTI addspec fixtures
    let pha1 = dir.path().join("q0.pha");
    let rsp1 = dir.path().join("q0.rsp");
    let pha2 = dir.path().join("q1.pha");
    let rsp2 = dir.path().join("q1.rsp");
    write_test_czti_pha_and_rsp(&pha1, &rsp1, 100.0, 1.0);
    write_test_czti_pha_and_rsp(&pha2, &rsp2, 200.0, 2.0);

    // Setup LAXPC backshift fixtures
    let caldb_dir = dir.path().join("laxpc_caldb");
    std::fs::create_dir_all(&caldb_dir).unwrap();
    let (bkg_pha, obs_fits) = write_test_laxpc_backshift_caldb(&caldb_dir);

    // Setup Pulse search + GTI fixtures
    let evt_path = dir.path().join("planted.evt");
    let freq_hz = 29.658;
    let epoch = 190_000_000.0;
    write_planted_multiband_events(&evt_path, freq_hz, epoch, 0.125);
    let half_span = 100.0 / freq_hz;
    let gti_half = dir.path().join("half.gti");
    let gti_empty = dir.path().join("narrow.gti");
    {
        let meta = OgipMeta::default();
        let h_half = build_gti_hdu("GTI", &[epoch], &[epoch + half_span], &meta).unwrap();
        write_fits(&gti_half, &mut [Hdu::empty_primary(), h_half], true).unwrap();

        let h_empty = build_gti_hdu("GTI", &[epoch - 10.0], &[epoch - 5.0], &meta).unwrap();
        write_fits(&gti_empty, &mut [Hdu::empty_primary(), h_empty], true).unwrap();
    }

    let srv = UiServer::spawn(UiConfig {
        host: "127.0.0.1".into(),
        port: 0,
        data_dir: dir.path().to_path_buf(),
        open_browser: false,
        threads: Some(2),
    })
    .unwrap();
    let addr = srv.addr();

    // 1. F4: czti-addspec preview + add_exposure=true + add_exposure=false
    let out_pha_sum = dir.path().join("sub/coadd_sum.pha");
    let out_rsp_sum = dir.path().join("sub/coadd_sum.rsp");
    let addspec_sum_req = json!({
        "pipeline": "czti-addspec",
        "threads": 2,
        "czti_addspec": {
            "pha": [pha1.to_string_lossy(), pha2.to_string_lossy()],
            "rsp": [rsp1.to_string_lossy(), rsp2.to_string_lossy()],
            "out_pha": out_pha_sum.to_string_lossy(),
            "out_rsp": out_rsp_sum.to_string_lossy(),
            "add_exposure": true
        }
    })
    .to_string();
    let prev = http_request(addr, "POST", "/api/cli-preview", Some(&addspec_sum_req));
    assert_eq!(prev.status, 200);
    let pv: Value = serde_json::from_str(&prev.body).unwrap();
    let cmd = pv["cli_command"].as_str().unwrap();
    assert!(cmd.contains("czti addspec"));
    assert!(cmd.contains("--pha"));
    assert!(cmd.contains("--rsp"));
    assert!(cmd.contains("--out-rsp"));
    assert!(cmd.contains("--add-exposure"));

    let run1 = http_request(addr, "POST", "/api/jobs/run", Some(&addspec_sum_req));
    assert_eq!(run1.status, 202, "czti-addspec failed: {}", run1.body);
    let jid1 = serde_json::from_str::<Value>(&run1.body).unwrap()["id"]
        .as_str()
        .unwrap()
        .to_string();
    let done1 = poll_job_until_done(addr, &jid1);
    assert_eq!(done1["status"], "succeeded", "addspec job failed: {done1}");
    assert!(out_pha_sum.is_file());
    assert!(out_rsp_sum.is_file());
    assert!((done1["summary"]["exposure_sec"].as_f64().unwrap() - 300.0).abs() < 1e-6);
    assert_eq!(done1["outputs"][0]["product"], "spectrum");
    assert_eq!(done1["outputs"][1]["product"], "response");

    // GET /api/spectrum reads co-added PHA cleanly
    let sp_r = http_request(
        addr,
        "GET",
        &format!("/api/spectrum?path={}", url_encode_path(&out_pha_sum)),
        None,
    );
    assert_eq!(sp_r.status, 200, "spectrum read failed: {}", sp_r.body);

    // add_exposure = false -> averaged exposure (100 + 200)/2 = 150
    let out_pha_avg = dir.path().join("coadd_avg.pha");
    let addspec_avg_req = json!({
        "pipeline": "czti-addspec",
        "czti_addspec": {
            "pha": [pha1.to_string_lossy(), pha2.to_string_lossy()],
            "out_pha": out_pha_avg.to_string_lossy(),
            "add_exposure": false
        }
    })
    .to_string();
    let run2 = http_request(addr, "POST", "/api/jobs/run", Some(&addspec_avg_req));
    assert_eq!(run2.status, 202);
    let jid2 = serde_json::from_str::<Value>(&run2.body).unwrap()["id"]
        .as_str()
        .unwrap()
        .to_string();
    let done2 = poll_job_until_done(addr, &jid2);
    assert_eq!(done2["status"], "succeeded");
    assert!((done2["summary"]["exposure_sec"].as_f64().unwrap() - 150.0).abs() < 1e-6);

    // Adversarial czti-addspec
    for bad in [
        json!({"pipeline": "czti-addspec", "czti_addspec": {"pha": [], "out_pha": "a.pha"}}),
        json!({"pipeline": "czti-addspec", "czti_addspec": {"pha": [pha1.to_string_lossy()], "rsp": [rsp1.to_string_lossy(), rsp2.to_string_lossy()], "out_pha": "a.pha"}}),
        json!({"pipeline": "czti-addspec", "czti_addspec": {"pha": [pha1.to_string_lossy()], "out_pha": "a.pha", "out_rsp": "a.rsp"}}),
        json!({"pipeline": "czti-addspec", "czti_addspec": {"pha": ["../escape.pha"], "out_pha": "a.pha"}}),
    ] {
        assert_eq!(
            http_request(addr, "POST", "/api/jobs/run", Some(&bad.to_string())).status,
            400,
            "expected 400 for {bad}"
        );
    }

    // 2. F5: laxpc-backshift preview + gain-shift-only + fits+bkg_pha+output
    let bs_preview = http_request(
        addr,
        "POST",
        "/api/cli-preview",
        Some(
            &json!({
                "pipeline": "laxpc-backshift",
                "laxpc_backshift": {
                    "unit": "10",
                    "caldb_dir": caldb_dir.to_string_lossy(),
                    "mission_day": 150.0
                }
            })
            .to_string(),
        ),
    );
    assert_eq!(bs_preview.status, 200);
    let bs_cmd = serde_json::from_str::<Value>(&bs_preview.body).unwrap()["cli_command"]
        .as_str()
        .unwrap()
        .to_string();
    assert!(bs_cmd.contains("astrosat laxpc backshift --unit 10"));

    // Gain-shift-only (day 150 vs 150, no bkg_pha)
    let bs_gain_only = http_request(
        addr,
        "POST",
        "/api/jobs/run",
        Some(
            &json!({
                "pipeline": "laxpc-backshift",
                "laxpc_backshift": {
                    "unit": "10",
                    "caldb_dir": caldb_dir.to_string_lossy(),
                    "mission_day": 150.0
                }
            })
            .to_string(),
        ),
    );
    assert_eq!(bs_gain_only.status, 202);
    let jid_bs1 = serde_json::from_str::<Value>(&bs_gain_only.body).unwrap()["id"]
        .as_str()
        .unwrap()
        .to_string();
    let done_bs1 = poll_job_until_done(addr, &jid_bs1);
    assert_eq!(done_bs1["status"], "succeeded", "bs1 failed: {done_bs1}");
    assert!(done_bs1["outputs"].as_array().unwrap().is_empty());
    assert!(done_bs1["summary"]["delta_ch_30kev"].is_number());
    assert!(done_bs1["summary"]["gain_ratio"].is_number());
    assert!(done_bs1["summary"]["secular_scale"].is_number());
    assert!(done_bs1["summary"]["recommended_rmf"].is_string());

    // fits + bkg_pha + output -> crosses HV jump at day 120.0 (bkg_day=100 -> src_day=150)
    let shifted_out = dir.path().join("shifted_out/lxp10_shifted.pha");
    let bs_full = http_request(
        addr,
        "POST",
        "/api/jobs/run",
        Some(
            &json!({
                "pipeline": "laxpc-backshift",
                "laxpc_backshift": {
                    "unit": "10",
                    "caldb_dir": caldb_dir.to_string_lossy(),
                    "fits": obs_fits.to_string_lossy(),
                    "bkg_pha": bkg_pha.to_string_lossy(),
                    "output": shifted_out.to_string_lossy()
                }
            })
            .to_string(),
        ),
    );
    assert_eq!(bs_full.status, 202);
    let jid_bs2 = serde_json::from_str::<Value>(&bs_full.body).unwrap()["id"]
        .as_str()
        .unwrap()
        .to_string();
    let done_bs2 = poll_job_until_done(addr, &jid_bs2);
    assert_eq!(done_bs2["status"], "succeeded", "bs2 failed: {done_bs2}");
    assert!(shifted_out.is_file());
    assert_eq!(done_bs2["outputs"][0]["product"], "spectrum");
    assert_eq!(done_bs2["summary"]["crossed_hv_jumps"], json!([120.0]));
    assert!(done_bs2["summary"]["lcbk_shift_fraction"].is_number());
    let sp_bs = http_request(
        addr,
        "GET",
        &format!("/api/spectrum?path={}", url_encode_path(&shifted_out)),
        None,
    );
    assert_eq!(sp_bs.status, 200, "shifted spectrum read failed: {}", sp_bs.body);

    // Adversarial laxpc-backshift
    for bad in [
        json!({"pipeline": "laxpc-backshift", "laxpc_backshift": {"unit": "10", "caldb_dir": caldb_dir.to_string_lossy(), "mission_day": 150.0, "fits": obs_fits.to_string_lossy()}}),
        json!({"pipeline": "laxpc-backshift", "laxpc_backshift": {"unit": "10", "caldb_dir": caldb_dir.to_string_lossy()}}),
        json!({"pipeline": "laxpc-backshift", "laxpc_backshift": {"unit": "99", "caldb_dir": caldb_dir.to_string_lossy(), "mission_day": 150.0}}),
        json!({"pipeline": "laxpc-backshift", "laxpc_backshift": {"unit": "10", "caldb_dir": caldb_dir.to_string_lossy(), "mission_day": 150.0, "pha_min": 500, "pha_max": 100}}),
        json!({"pipeline": "laxpc-backshift", "laxpc_backshift": {"unit": "10", "caldb_dir": caldb_dir.to_string_lossy(), "mission_day": 150.0, "output": "out.pha"}}),
        json!({"pipeline": "laxpc-backshift", "laxpc_backshift": {"unit": "10", "caldb_dir": "../caldb", "mission_day": 150.0}}),
    ] {
        assert_eq!(
            http_request(addr, "POST", "/api/jobs/run", Some(&bad.to_string())).status,
            400,
            "expected 400 for {bad}"
        );
    }

    // 3. F7: Pulse timing parity (gti + df_hz + fddot)
    let full_search = http_request(
        addr,
        "POST",
        "/api/pulse/search",
        Some(
            &json!({
                "path": evt_path.to_string_lossy(),
                "freq_min": 29.60,
                "freq_max": 29.70,
                "df_hz": 0.005,
                "fddot": 0.0
            })
            .to_string(),
        ),
    );
    assert_eq!(full_search.status, 200);
    let full_v: Value = serde_json::from_str(&full_search.body).unwrap();
    let full_events = full_v["num_events"].as_u64().unwrap();
    assert_eq!(full_v["n_trials"], 20);

    let gti_search = http_request(
        addr,
        "POST",
        "/api/pulse/search",
        Some(
            &json!({
                "path": evt_path.to_string_lossy(),
                "gti": gti_half.to_string_lossy(),
                "freq_min": 29.60,
                "freq_max": 29.70,
                "df_hz": 0.005,
                "fddot": 0.0
            })
            .to_string(),
        ),
    );
    assert_eq!(gti_search.status, 200, "gti pulse search failed: {}", gti_search.body);
    let gti_v: Value = serde_json::from_str(&gti_search.body).unwrap();
    let gti_events = gti_v["num_events"].as_u64().unwrap();
    assert!(
        gti_events < full_events && (gti_events as f64) <= (full_events as f64) * 0.55,
        "expected ~half events after GTI filter, got {gti_events} of {full_events}"
    );
    assert!((gti_v["best_freq_hz"].as_f64().unwrap() - 29.658).abs() < 0.01);

    // Adversarial F7
    assert_eq!(
        http_request(
            addr,
            "POST",
            "/api/pulse/search",
            Some(&json!({"path": evt_path.to_string_lossy(), "freq_hz": 29.658, "df_hz": 0.0}).to_string())
        )
        .status,
        400
    );
    assert_eq!(
        http_request(
            addr,
            "POST",
            "/api/pulse/search",
            Some(&json!({"path": evt_path.to_string_lossy(), "freq_hz": 29.658, "gti": "missing.gti"}).to_string())
        )
        .status,
        404
    );
    assert_eq!(
        http_request(
            addr,
            "POST",
            "/api/pulse/search",
            Some(&json!({"path": evt_path.to_string_lossy(), "freq_hz": 29.658, "gti": gti_empty.to_string_lossy()}).to_string())
        )
        .status,
        400
    );

    srv.stop();
}
