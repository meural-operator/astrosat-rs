//! Opaque-box, requirement-driven 4-tier E2E parity & adversarial test suite
//! for AstroSat Studio v2 (`F1`–`F7` backend, API, routing, TOML & CLI parity).
//!
//! Run baseline harness tests (green during M1/M2):
//!   cargo test -p astrosat-ui --test e2e_parity -- --nocapture
//!
//! Run full M3 E2E parity suite (Tiers 1–4):
//!   cargo test -p astrosat-ui --test e2e_parity -- --include-ignored --nocapture

use std::collections::HashSet;
use std::f64::consts::TAU;
use std::io::{Read, Write};
use std::net::{SocketAddr, TcpStream};
use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::Duration;

use astrosat_core::time::mjd_utc_to_astrosat_utc;
use astrosat_fits::{
    build_gti_hdu, build_lightcurve_hdu, build_orbit_hdu, build_response_hdus, build_spectrum_hdu,
    write_fits, BinTableBuilder, FitsFile, Hdu, ImageBuilder, OgipMeta, ResponseMatrixData,
    SpectrumData, SpectrumMeta, WcsParams,
};
use astrosat_ui::{UiConfig, UiServer};
use serde_json::{json, Value};
use tempfile::tempdir;

// ---------------------------------------------------------------------------
// Minimal zero-dep HTTP/1.1 client over std::net::TcpStream
// ---------------------------------------------------------------------------

struct HttpResponse {
    status: u16,
    content_type: String,
    body: String,
}

fn http_request_with_headers(
    addr: SocketAddr,
    method: &str,
    path: &str,
    content_type: Option<&str>,
    body: Option<&str>,
) -> HttpResponse {
    let mut stream = TcpStream::connect(addr).expect("connect to UiServer");
    stream
        .set_read_timeout(Some(Duration::from_secs(120)))
        .expect("set_read_timeout");
    stream
        .set_write_timeout(Some(Duration::from_secs(30)))
        .expect("set_write_timeout");

    let req = if let Some(payload) = body {
        let ct = content_type.unwrap_or("application/json");
        format!(
            "{method} {path} HTTP/1.1\r\nHost: 127.0.0.1\r\nContent-Type: {ct}\r\nContent-Length: {}\r\nConnection: close\r\n\r\n{payload}",
            payload.len()
        )
    } else {
        format!("{method} {path} HTTP/1.1\r\nHost: 127.0.0.1\r\nConnection: close\r\n\r\n")
    };
    stream.write_all(req.as_bytes()).expect("write HTTP request");
    stream.flush().expect("flush HTTP request");

    read_http_response(&mut stream)
}

fn http_request_allow_early_close(
    addr: SocketAddr,
    method: &str,
    path: &str,
    content_type: Option<&str>,
    body: &str,
) -> Option<HttpResponse> {
    let mut stream = TcpStream::connect(addr).expect("connect to UiServer");
    stream
        .set_read_timeout(Some(Duration::from_secs(120)))
        .expect("set_read_timeout");
    stream
        .set_write_timeout(Some(Duration::from_secs(30)))
        .expect("set_write_timeout");

    let ct = content_type.unwrap_or("application/json");
    let req = format!(
        "{method} {path} HTTP/1.1\r\nHost: 127.0.0.1\r\nContent-Type: {ct}\r\nContent-Length: {}\r\nConnection: close\r\n\r\n{body}",
        body.len()
    );
    if let Err(e) = stream.write_all(req.as_bytes()) {
        if !matches!(
            e.kind(),
            std::io::ErrorKind::BrokenPipe | std::io::ErrorKind::ConnectionReset
        ) {
            panic!("unexpected write error: {e}");
        }
    }
    let _ = stream.flush();

    let mut raw = Vec::new();
    match stream.read_to_end(&mut raw) {
        Ok(_) => {}
        Err(e)
            if matches!(
                e.kind(),
                std::io::ErrorKind::ConnectionReset | std::io::ErrorKind::BrokenPipe
            ) => {}
        Err(e) => panic!("unexpected read error: {e}"),
    }
    if raw.is_empty() {
        return None;
    }
    Some(parse_http_response_bytes(&raw))
}

fn read_http_response(stream: &mut TcpStream) -> HttpResponse {
    let mut raw = Vec::new();
    stream.read_to_end(&mut raw).expect("read HTTP response");
    parse_http_response_bytes(&raw)
}

fn parse_http_response_bytes(raw: &[u8]) -> HttpResponse {
    let text = String::from_utf8_lossy(raw).into_owned();
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
        .expect("numeric HTTP status");

    let mut resp_ct = String::new();
    let mut is_chunked = false;
    let mut content_len: Option<usize> = None;
    for line in lines {
        if let Some((k, v)) = line.split_once(':') {
            if k.eq_ignore_ascii_case("content-type") {
                resp_ct = v.trim().to_string();
            } else if k.eq_ignore_ascii_case("transfer-encoding")
                && v.to_ascii_lowercase().contains("chunked")
            {
                is_chunked = true;
            } else if k.eq_ignore_ascii_case("content-length") {
                content_len = v.trim().parse::<usize>().ok();
            }
        }
    }

    let decoded_body = if is_chunked {
        decode_chunked(rest)
    } else if let Some(n) = content_len {
        let bytes = rest.as_bytes();
        String::from_utf8_lossy(&bytes[..n.min(bytes.len())]).into_owned()
    } else {
        rest.to_string()
    };

    HttpResponse {
        status,
        content_type: resp_ct,
        body: decoded_body,
    }
}

fn http_request(addr: SocketAddr, method: &str, path: &str, body: Option<&str>) -> HttpResponse {
    http_request_with_headers(addr, method, path, Some("application/json"), body)
}

fn decode_chunked(s: &str) -> String {
    let mut out = String::new();
    let mut cur = s;
    while let Some((len_hex, after)) = cur.split_once("\r\n") {
        let len_tok = len_hex.split(';').next().unwrap_or(len_hex).trim();
        let Ok(len) = usize::from_str_radix(len_tok, 16) else {
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

fn url_encode(s: &str) -> String {
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

fn spawn_test_server(data_dir: &Path) -> UiServer {
    UiServer::spawn(UiConfig {
        host: "127.0.0.1".into(),
        port: 0,
        data_dir: data_dir.to_path_buf(),
        open_browser: false,
        threads: Some(2),
    })
    .expect("spawn UiServer on ephemeral port 127.0.0.1:0")
}

fn poll_job_until_done(addr: SocketAddr, job_id: &str, timeout_secs: u64) -> Value {
    let max_iters = (timeout_secs * 1000) / 50;
    for _ in 0..max_iters.max(20) {
        let r = http_request(addr, "GET", &format!("/api/jobs/{job_id}"), None);
        assert_eq!(r.status, 200, "job poll failed: {}", r.body);
        let j: Value = serde_json::from_str(&r.body).expect("valid job JSON");
        if j["status"] != "running" {
            return j;
        }
        std::thread::sleep(Duration::from_millis(50));
    }
    panic!("job {job_id} did not finish within {timeout_secs}s");
}

fn submit_and_wait(addr: SocketAddr, req_body: &Value, timeout_secs: u64) -> Value {
    let payload = serde_json::to_string(req_body).unwrap();
    let r = http_request(addr, "POST", "/api/jobs/run", Some(&payload));
    assert_eq!(r.status, 202, "expected 202 job spawn, got {}: {}", r.status, r.body);
    let spawned: Value = serde_json::from_str(&r.body).unwrap();
    let id = spawned["id"].as_str().expect("job id");
    poll_job_until_done(addr, id, timeout_secs)
}

fn assert_offline_asset(body: &str, route: &str) {
    let cleaned = body.replace("http://www.w3.org/2000/svg", "");
    assert!(
        !cleaned.contains("https://") && !cleaned.contains("http://"),
        "route {route} leaked external URL"
    );
}

fn find_astrosat_binary() -> PathBuf {
    if let Ok(explicit) = std::env::var("ASTROSAT_BIN") {
        let p = PathBuf::from(explicit);
        if p.is_file() {
            return p;
        }
    }
    let exe_name = if cfg!(windows) { "astrosat.exe" } else { "astrosat" };
    let mut candidates = Vec::new();
    if let Ok(tdir) = std::env::var("CARGO_TARGET_DIR") {
        let base = PathBuf::from(tdir);
        candidates.push(base.join("release").join(exe_name));
        candidates.push(base.join("debug").join(exe_name));
    }
    let workspace_root = Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(|p| p.parent())
        .expect("workspace root");
    candidates.push(workspace_root.join("target").join("release").join(exe_name));
    candidates.push(workspace_root.join("target").join("debug").join(exe_name));

    let newest = candidates
        .iter()
        .filter(|c| c.is_file())
        .max_by_key(|c| {
            std::fs::metadata(c)
                .and_then(|m| m.modified())
                .unwrap_or(std::time::SystemTime::UNIX_EPOCH)
        });
    if let Some(best) = newest {
        return best.clone();
    }
    panic!(
        "astrosat CLI binary not found in {:?}; build it first with `cargo build -p astrosat-cli`",
        candidates
    );
}

/// Split a shell command preview line respecting `'...'` and `"..."` quotes.
fn split_shell_words(cmd: &str) -> Vec<String> {
    let mut args = Vec::new();
    let mut cur = String::new();
    let mut in_sq = false;
    let mut in_dq = false;
    let mut esc = false;
    let mut has_token = false;

    for ch in cmd.chars() {
        if esc {
            cur.push(ch);
            esc = false;
            has_token = true;
            continue;
        }
        if in_sq {
            if ch == '\'' {
                in_sq = false;
            } else {
                cur.push(ch);
            }
            has_token = true;
            continue;
        }
        if in_dq {
            if ch == '\\' {
                esc = true;
            } else if ch == '"' {
                in_dq = false;
            } else {
                cur.push(ch);
            }
            has_token = true;
            continue;
        }
        match ch {
            '\'' => {
                in_sq = true;
                has_token = true;
            }
            '"' => {
                in_dq = true;
                has_token = true;
            }
            '\\' => {
                esc = true;
                has_token = true;
            }
            c if c.is_whitespace() => {
                if has_token {
                    args.push(std::mem::take(&mut cur));
                    has_token = false;
                }
            }
            c => {
                cur.push(c);
                has_token = true;
            }
        }
    }
    if has_token {
        args.push(cur);
    }
    args
}

// ---------------------------------------------------------------------------
// Synthetic FITS & CALDB Fixture Builders
// ---------------------------------------------------------------------------

fn write_synthetic_czti_pha_512(
    path: &Path,
    counts: &[f32],
    errs: &[f32],
    exposure: f64,
    respfile: &str,
    quad_id: Option<i64>,
) {
    assert_eq!(counts.len(), 512);
    assert_eq!(errs.len(), 512);
    let mut meta = OgipMeta {
        instrume: "CZTI".into(),
        detnam: quad_id.map(|q| format!("QUADRANT{q}")),
        object: Some("CRAB".into()),
        ra_obj: Some(83.6331),
        dec_obj: Some(22.0145),
        tstart: 190_000_000.0,
        tstop: 190_000_000.0 + exposure,
        exposure,
        ..Default::default()
    };
    meta.timesys = "UTC".into();
    let p0 = Hdu::empty_primary();
    let channels: Vec<i32> = (0..512).collect();
    let spec_meta = SpectrumMeta {
        detchans: 512,
        tlmin1: 0,
        tlmax1: 511,
        respfile: respfile.to_string(),
        ..Default::default()
    };
    let counts_i32: Vec<i32> = counts.iter().map(|&v| v.round() as i32).collect();
    let mut spec = build_spectrum_hdu(
        "SPECTRUM",
        &channels,
        SpectrumData::Counts(counts_i32),
        Some(errs),
        None,
        None,
        &meta,
        &spec_meta,
    )
    .unwrap();
    if let Some(q) = quad_id {
        spec.header.set_i64("QUADID", q, "Quadrant ID");
    }
    spec.header.set_str("MASKWT", "YES", "Mask weighting applied");
    spec.header.set_str("GTITYPE", "QUAD", "GTI type");

    // Also include EBOUNDS so `/api/spectrum` can map energy
    let ch_i16: Vec<i16> = (0..512).map(|i| i as i16).collect();
    let e_min: Vec<f32> = (0..512).map(|i| 5.0 + 0.5 * (i as f32)).collect();
    let e_max: Vec<f32> = (0..512).map(|i| 5.5 + 0.5 * (i as f32)).collect();
    let rmf = ResponseMatrixData {
        chantype: "PI".into(),
        detchans: 512,
        tlmin4: 0,
        energ_lo: vec![10.0, 100.0],
        energ_hi: vec![100.0, 250.0],
        n_grp: vec![1, 1],
        f_chan: vec![vec![0], vec![200]],
        n_chan: vec![vec![10], vec![10]],
        matrix: vec![vec![0.1; 10], vec![0.1; 10]],
        channels: ch_i16,
        e_min,
        e_max,
        lo_thresh: 1e-10,
    };
    let (_mat, eb) = build_response_hdus(&rmf, &meta).unwrap();
    write_fits(path, &mut [p0, spec, eb], true).unwrap();
}

fn write_synthetic_czti_rsp_512(path: &Path, area_val: f32) {
    let meta = OgipMeta {
        instrume: "CZTI".into(),
        detnam: Some("QUADRANT0".into()),
        object: Some("CRAB".into()),
        tstart: 190_000_000.0,
        tstop: 190_001_000.0,
        exposure: 0.0,
        ..Default::default()
    };
    let p0 = Hdu::empty_primary();
    let ch_i16: Vec<i16> = (0..512).map(|i| i as i16).collect();
    let e_min: Vec<f32> = (0..512).map(|i| 5.0 + 0.5 * (i as f32)).collect();
    let e_max: Vec<f32> = (0..512).map(|i| 5.5 + 0.5 * (i as f32)).collect();
    let rmf = ResponseMatrixData {
        chantype: "PI".into(),
        detchans: 512,
        tlmin4: 0,
        energ_lo: vec![10.0, 50.0, 100.0],
        energ_hi: vec![50.0, 100.0, 200.0],
        n_grp: vec![1, 1, 1],
        f_chan: vec![vec![10], vec![90], vec![190]],
        n_chan: vec![vec![5], vec![5], vec![5]],
        matrix: vec![
            vec![area_val; 5],
            vec![area_val * 0.8; 5],
            vec![area_val * 0.5; 5],
        ],
        channels: ch_i16,
        e_min,
        e_max,
        lo_thresh: 1e-10,
    };
    let (mat, eb) = build_response_hdus(&rmf, &meta).unwrap();
    write_fits(path, &mut [p0, eb, mat], true).unwrap();
}

fn write_synthetic_laxpc_caldb(dir: &Path, with_channel_map: bool) {
    std::fs::create_dir_all(dir).unwrap();
    let gaina8 = "\
# day p30_1 f30_1 p60_1 f60_1 p30_2 f30_2 p60_2 f60_2 p30_3 f30_3 p60_3 f60_3
0.0   220.0 25.0 440.0 35.0  215.0 25.0 430.0 35.0  210.0 25.0 420.0 35.0
100.0 225.0 25.5 450.0 35.5  218.0 25.2 436.0 35.2  212.0 25.1 424.0 35.1
500.0 240.0 26.5 480.0 37.0  226.0 26.0 452.0 36.0  220.0 26.0 440.0 36.0
1000.0 260.0 28.0 520.0 39.0 235.0 27.0 470.0 37.0
";
    std::fs::write(dir.join("gaina8.dat"), gaina8).unwrap();
    if with_channel_map {
        let mut ch_lines = String::new();
        for ch in 0..512 {
            let e0 = 0.2 * (ch as f64);
            let e1 = 0.2 * ((ch + 1) as f64);
            ch_lines.push_str(&format!("{e0:.4} {e1:.4}\n"));
        }
        std::fs::write(dir.join("lxp10.channel"), &ch_lines).unwrap();
        std::fs::write(dir.join("lxp20.channel"), &ch_lines).unwrap();
        std::fs::write(dir.join("lxp30.channel"), &ch_lines).unwrap();
    }
}

fn write_synthetic_laxpc_bkg_pha(path: &Path, mission_day: f64) {
    let mut out = format!(
        "#  Exposure =   1.000000E+04  lxp1  Channels:    0  511  day:  {mission_day:.1}\n"
    );
    for ch in 0..512 {
        let base = if (100..220).contains(&ch) { 0.05 } else { 0.01 };
        let sum7 = base * 7.0;
        out.push_str(&format!(
            "{ch:4} {base:.5} {base:.5} {base:.5} {base:.5} {base:.5} {base:.5} {base:.5} 0.00100 0.00100 0.00100 {sum7:.5}\n"
        ));
    }
    std::fs::write(path, out).unwrap();
}

fn write_synthetic_tstart_fits(path: &Path, tstart: f64) {
    let mut p0 = Hdu::empty_primary();
    p0.header.set_str("INSTRUME", "LAXPC", "Instrument");
    p0.header.set_str("OBJECT", "CRAB", "Target");
    p0.header.set_f64("TSTART", tstart, "Observation start MET");
    let evt = BinTableBuilder::new("EVENTS")
        .add_f64("TIME", Some("s"), vec![tstart, tstart + 10.0])
        .header_f64("TSTART", tstart, "Observation start MET")
        .build()
        .unwrap();
    write_fits(path, &mut [p0, evt], true).unwrap();
}

fn write_synthetic_primary_without_tstart(path: &Path) {
    let mut p0 = Hdu::empty_primary();
    p0.header.set_str("OBJECT", "EMPTY", "No TSTART");
    write_fits(path, &mut [p0], true).unwrap();
}

fn write_synthetic_utc_events_and_orbit(evt_path: &Path, orb_path: &Path, timesys: &str) {
    let t0 = mjd_utc_to_astrosat_utc(57400.0);
    let mut times = Vec::new();
    for i in 0..24 {
        times.push(t0 + (i as f64) * 5.0);
    }
    let p0 = Hdu::empty_primary();
    let evt = BinTableBuilder::new("EVENTS")
        .add_f64("TIME", Some("s"), times)
        .add_i16("PI", Some("chan"), vec![40; 24])
        .header_str("INSTRUME", "CZTI", "instrument")
        .header_str("TIMESYS", timesys, "time system")
        .header_f64("RA_OBJ", 83.6331, "Crab RA")
        .header_f64("DEC_OBJ", 22.0145, "Crab Dec")
        .header_f64("TSTART", t0, "start")
        .header_f64("TSTOP", t0 + 120.0, "stop")
        .build()
        .unwrap();
    let mut p0_mut = p0;
    p0_mut.header.set_str("TIMESYS", timesys, "time system");
    write_fits(evt_path, &mut [p0_mut, evt], true).unwrap();

    // Circular LEO orbit around Earth (~7000 km radius)
    let mut ot = Vec::new();
    let mut pos = Vec::new();
    let mut vel = Vec::new();
    let r_m = 7_028_000.0_f64;
    let omega = 0.00107_f64;
    for i in 0..30 {
        let dt = (i as f64) * 10.0 - 30.0;
        let ang = omega * dt;
        ot.push(t0 + dt);
        pos.push([r_m * ang.cos(), r_m * ang.sin(), 0.0]);
        vel.push([-r_m * omega * ang.sin(), r_m * omega * ang.cos(), 0.0]);
    }
    let orb_p0 = Hdu::empty_primary();
    let orb_meta = OgipMeta {
        instrume: "ASTROSAT".into(),
        tstart: t0 - 30.0,
        tstop: t0 + 260.0,
        ..Default::default()
    };
    let orb_hdu = build_orbit_hdu(&ot, &pos, &vel, &orb_meta).unwrap();
    write_fits(orb_path, &mut [orb_p0, orb_hdu], true).unwrap();
}

fn write_synthetic_mkf_for_bary(mkf_path: &Path) {
    let t0 = mjd_utc_to_astrosat_utc(57400.0);
    let n = 30usize;
    let mut time = Vec::with_capacity(n);
    let mut x = Vec::with_capacity(n);
    let mut y = Vec::with_capacity(n);
    let mut z = Vec::with_capacity(n);
    let mut vx = Vec::with_capacity(n);
    let mut vy = Vec::with_capacity(n);
    let mut vz = Vec::with_capacity(n);
    let r_m = 7_028_000.0_f64;
    let omega = 0.00107_f64;
    for i in 0..n {
        let dt = (i as f64) * 10.0 - 30.0;
        let ang = omega * dt;
        time.push(t0 + dt);
        x.push(r_m * ang.cos());
        y.push(r_m * ang.sin());
        z.push(0.0);
        vx.push(-r_m * omega * ang.sin());
        vy.push(r_m * omega * ang.cos());
        vz.push(0.0);
    }
    let p0 = Hdu::empty_primary();
    let mkf = BinTableBuilder::new("MKF")
        .add_f64("TIME", Some("s"), time)
        .add_f64("POSX", Some("m"), x)
        .add_f64("POSY", Some("m"), y)
        .add_f64("POSZ", Some("m"), z)
        .add_f64("VELX", Some("m/s"), vx)
        .add_f64("VELY", Some("m/s"), vy)
        .add_f64("VELZ", Some("m/s"), vz)
        .build()
        .unwrap();
    write_fits(mkf_path, &mut [p0, mkf], true).unwrap();
}

/// Planted pulsar events with optional fdot, fddot, LAXPC_No, Layer, PI columns and 2 energy bands.
fn write_planted_pulsar_events(
    path: &Path,
    epoch: f64,
    freq_hz: f64,
    fdot: f64,
    fddot: f64,
    shift_cycles: f64,
) {
    write_planted_pulsar_events_ext(path, epoch, freq_hz, fdot, fddot, shift_cycles, 250, 1);
}

#[allow(clippy::too_many_arguments)]
fn write_planted_pulsar_events_ext(
    path: &Path,
    epoch: f64,
    freq_hz: f64,
    fdot: f64,
    fddot: f64,
    shift_cycles: f64,
    n_pulses: usize,
    cycle_stride: usize,
) {
    let mut times = Vec::new();
    let mut energies = Vec::new();
    let mut pis = Vec::new();
    let mut units = Vec::new();
    let mut layers = Vec::new();

    let period = 1.0 / freq_hz;
    let total_cycles = n_pulses * cycle_stride.max(1);

    // Solve approx phase(t) = f*(dt) + 0.5*fdot*dt^2 + (1/6)*fddot*dt^3 = c + phi
    let solve_t = |target_cyc: f64| -> f64 {
        let mut dt = target_cyc / freq_hz;
        for _ in 0..8 {
            let ph = freq_hz * dt + 0.5 * fdot * dt * dt + (1.0 / 6.0) * fddot * dt * dt * dt;
            let inst_f = (freq_hz + fdot * dt + 0.5 * fddot * dt * dt).max(1e-6);
            dt -= (ph - target_cyc) / inst_f;
        }
        epoch + dt
    };

    let template: &[(f64, usize)] = &[
        (0.05, 1),
        (0.24, 6),
        (0.25, 12),
        (0.26, 6),
        (0.55, 1),
        (0.65, 5),
        (0.85, 1),
    ];

    for p in 0..n_pulses {
        let c = p * cycle_stride.max(1);
        for &(phi, mult) in template {
            for k in 0..mult {
                let jitter = ((k as f64) - (mult as f64) * 0.5) * 1e-4;
                let t1 = solve_t((c as f64) + phi + jitter);
                times.push(t1);
                energies.push(12.0_f32);
                pis.push(25_i16);
                units.push(1_u8);
                layers.push(1_u8);

                let t2 = solve_t((c as f64) + phi + shift_cycles + jitter);
                times.push(t2);
                energies.push(55.0_f32);
                pis.push(140_i16);
                units.push(2_u8);
                layers.push(2_u8);
            }
        }
    }

    let p0 = Hdu::empty_primary();
    let evt = BinTableBuilder::new("EVENTS")
        .add_f64("TIME", Some("s"), times)
        .add_f32("ENERGY", Some("keV"), energies)
        .add_i16("PI", Some("chan"), pis)
        .add_u8("LAXPC_No", None, units)
        .add_u8("Layer", None, layers)
        .header_str("INSTRUME", "LAXPC", "instrument")
        .header_str("TIMESYS", "TDB", "barycentric")
        .header_f64("TSTART", epoch, "start")
        .header_f64("TSTOP", epoch + (total_cycles as f64) * period, "stop")
        .header_f64("EXPOSURE", (total_cycles as f64) * period, "exposure")
        .build()
        .unwrap();
    write_fits(path, &mut [p0, evt], true).unwrap();
}

fn write_tiny_pulsar_events(path: &Path, epoch: f64, freq_hz: f64) {
    let mut times = Vec::with_capacity(20);
    for i in 0..20usize {
        let cyc = (i * 15) as f64 + 0.25;
        times.push(epoch + cyc / freq_hz);
    }
    let p0 = Hdu::empty_primary();
    let evt = BinTableBuilder::new("EVENTS")
        .add_f64("TIME", Some("s"), times)
        .add_i16("PI", Some("chan"), vec![25_i16; 20])
        .header_str("INSTRUME", "LAXPC", "instrument")
        .header_str("TIMESYS", "TDB", "barycentric")
        .header_f64("TSTART", epoch, "start")
        .header_f64("TSTOP", epoch + 11.0, "stop")
        .build()
        .unwrap();
    write_fits(path, &mut [p0, evt], true).unwrap();
}

fn write_synthetic_gti(path: &Path, starts: &[f64], stops: &[f64]) {
    let meta = OgipMeta {
        instrume: "LAXPC".into(),
        tstart: *starts.first().unwrap_or(&0.0),
        tstop: *stops.last().unwrap_or(&0.0),
        ..Default::default()
    };
    let p0 = Hdu::empty_primary();
    let gti = build_gti_hdu("GTI", starts, stops, &meta).unwrap();
    write_fits(path, &mut [p0, gti], true).unwrap();
}

// ===========================================================================
// ACTIVE BASELINE TESTS (Always run & pass during M1/M2)
// ===========================================================================

#[test]
fn e2e_harness_health_security_and_dataset_paths() {
    let dir = tempdir().unwrap();
    let server = spawn_test_server(dir.path());
    let addr = server.addr();

    // 1. Health & CORS/OPTIONS & Favicon
    let opt = http_request(addr, "OPTIONS", "/api/health", None);
    assert_eq!(opt.status, 204);

    let fav = http_request(addr, "GET", "/favicon.ico", None);
    assert_eq!(fav.status, 204);

    let health = http_request(addr, "GET", "/api/health", None);
    assert_eq!(health.status, 200);
    let hv: Value = serde_json::from_str(&health.body).unwrap();
    assert_eq!(hv["status"], "ok");

    // 2. Currently served static assets remain 100% offline (zero external URLs)
    for route in ["/", "/studio", "/studio/studio.css", "/studio/studio.js", "/styles.css", "/charts.js", "/app.js"] {
        let r = http_request(addr, "GET", route, None);
        assert_eq!(r.status, 200, "static route {route} failed");
        assert_offline_asset(&r.body, route);
    }

    // 3. Path traversal & adversarial HTTP checks (must return structured 4xx JSON, zero panics)
    let bad_gets = [
        "/api/browse?path=../../etc/passwd",
        "/api/browse?path=%2e%2e/%2e%2e/etc",
        "/api/browse?path=..\\..\\etc\\passwd",
        "/api/inspect?path=/etc/passwd",
        "/api/inspect?path=nonexistent_file.fits",
        "/api/unknown_e2e_route",
    ];
    for route in bad_gets {
        let r = http_request(addr, "GET", route, None);
        assert!((400..500).contains(&r.status), "route {route} returned {}", r.status);
        let j: Value = serde_json::from_str(&r.body).expect("4xx must be JSON");
        assert!(j["error"].is_string());
    }

    // 4. Method mismatch returns 405 JSON (GET /api/jobs/run hits /api/jobs/{id} -> 404/405 JSON)
    for post_route in ["/api/cli-preview", "/api/pulse/search", "/api/pulse/multiband"] {
        let r = http_request(addr, "GET", post_route, None);
        assert_eq!(r.status, 405, "GET {post_route} should return 405");
        let j: Value = serde_json::from_str(&r.body).unwrap();
        assert!(j["error"].is_string());
    }
    let r_jobs_run = http_request(addr, "GET", "/api/jobs/run", None);
    assert!(
        matches!(r_jobs_run.status, 404 | 405),
        "GET /api/jobs/run should return 404 or 405, got {}",
        r_jobs_run.status
    );
    let j_jobs_run: Value = serde_json::from_str(&r_jobs_run.body).unwrap();
    assert!(j_jobs_run["error"].is_string());

    // 5. Verify Tier-4 AstroSat dataset paths if ASTROSAT_TEST_DATA is set
    if let Ok(data_root_str) = std::env::var("ASTROSAT_TEST_DATA") {
        let data_root = Path::new(&data_root_str);
        if data_root.is_dir() {
            let required = [
                "sample/czti_2018/20181014_C04_007T02_9000002432_level2/czti/modeM0/AS1C04_007T02_9000002432cztM0_level2_quad_clean_Q0.pha",
                "sample/czti_2018/20181014_C04_007T02_9000002432_level2/czti/modeM0/AS1C04_007T02_9000002432cztM0_level2_quad_clean_Q0.rsp",
                "ref_src/laxpc/backlxp20jan16.pha",
                "sample/laxpc_l2/20160107_G02_010T01_9000000252_level2/laxpc/01509/lxp2/modeEA/AS1G02_010T01_9000000252lxp2EA_level2.fits",
                "sample/czti/czti_sample_data/20160331_T01_112T01_9000000406_level1",
                "sample/czti/czti_sample_data/20160331_T01_112T01_9000000406_level1/AS1T01_112T01_9000000406czt_level2.mkf",
                "sample/czti/czti_sample_data/20160331_T01_112T01_9000000406_level2/modeM0/AS1T01_112T01_9000000406czt_level2_quad_clean.evt",
                "sample/czti/czti_sample_data/20160331_T01_112T01_9000000406_level2/modeM0/AS1T01_112T01_9000000406czt_level2.gti",
                "sample/czti/czti_sample_data/20160331_T01_112T01_9000000406_level2/modeM0/AS1T01_112T01_9000000406cztM0_level2_quad_clean.img_Q0",
                "sample/czti_orb/AS1T01_112T01_9000000406czt_level1.orb",
                "caldb/czti",
                "ephem/JPLEPH.405",
            ];
            for rel in required {
                let p = data_root.join(rel);
                assert!(p.exists(), "missing Tier-4 dataset path: {}", p.display());
            }
        }
    }

    server.stop();
}

#[test]
fn e2e_harness_baseline_science_and_bary_roundtrip() {
    let dir = tempdir().unwrap();
    let lc_path = dir.path().join("test.lc");
    let img_path = dir.path().join("test.img");
    let evt_path = dir.path().join("pulsar.evt");
    let utc_evt = dir.path().join("utc.evt");
    let orb_path = dir.path().join("sat.orb");
    let bary_out = dir.path().join("bary.evt");

    let meta = OgipMeta {
        instrume: "CZTI".into(),
        tstart: 190_000_000.0,
        tstop: 190_000_004.0,
        exposure: 4.0,
        ..Default::default()
    };
    let lc = build_lightcurve_hdu(
        "RATE",
        &[190_000_000.5, 190_000_001.5, 190_000_002.5, 190_000_003.5],
        &[10.0, 20.0, 30.0, 40.0],
        &[1.0, 2.0, 3.0, 4.0],
        Some(&[1.0; 4]),
        1.0,
        &meta,
    )
    .unwrap();
    write_fits(&lc_path, &mut [Hdu::empty_primary(), lc], true).unwrap();

    let mut pix = vec![5.0_f32; 64];
    pix[2 * 8 + 3] = 100.0;
    let wcs = WcsParams {
        crpix1: 4.5,
        crpix2: 4.5,
        crval1: 83.6331,
        crval2: 22.0145,
        cdelt1: -0.05,
        cdelt2: 0.05,
        ..Default::default()
    };
    let ih = ImageBuilder::from_f32(None, vec![8, 8], &pix)
        .unwrap()
        .with_wcs(&wcs)
        .build_primary();
    write_fits(&img_path, &mut [ih], true).unwrap();

    let epoch = 190_000_000.0;
    write_planted_pulsar_events(&evt_path, epoch, 29.658, 0.0, 0.0, 0.25);
    write_synthetic_utc_events_and_orbit(&utc_evt, &orb_path, "UTC");

    let server = spawn_test_server(dir.path());
    let addr = server.addr();

    let lc_r = http_request(
        addr,
        "GET",
        &format!("/api/lightcurve?path={}&rebin=2", url_encode(&lc_path.to_string_lossy())),
        None,
    );
    assert_eq!(lc_r.status, 200);
    let lc_v: Value = serde_json::from_str(&lc_r.body).unwrap();
    assert_eq!(lc_v["rate"].as_array().unwrap().len(), 2);

    let img_r = http_request(
        addr,
        "GET",
        &format!("/api/image?path={}", url_encode(&img_path.to_string_lossy())),
        None,
    );
    assert_eq!(img_r.status, 200);
    let img_v: Value = serde_json::from_str(&img_r.body).unwrap();
    assert_eq!(img_v["stats"]["peak"]["x"], 3);
    assert_eq!(img_v["stats"]["peak"]["y"], 2);
    assert_eq!(img_v["wcs"]["has_wcs"], true);

    let mb_r = http_request(
        addr,
        "POST",
        "/api/pulse/multiband",
        Some(
            &json!({
                "inputs": [{
                    "label": "LAXPC",
                    "path": evt_path.to_string_lossy(),
                    "bands": [
                        { "label": "Soft", "emin_kev": 5.0, "emax_kev": 20.0 },
                        { "label": "Hard", "emin_kev": 30.0, "emax_kev": 80.0 }
                    ]
                }],
                "freq_hz": 29.658,
                "epoch": epoch,
                "nbins": 64,
                "harmonics": 2
            })
            .to_string(),
        ),
    );
    assert_eq!(mb_r.status, 200);
    let mb_v: Value = serde_json::from_str(&mb_r.body).unwrap();
    let lag = mb_v["bands"][1]["phase_lag_cycles"].as_f64().unwrap();
    assert!((lag - 0.25).abs() < 0.03, "expected lag ~0.25, got {lag}");

    let job_done = submit_and_wait(
        addr,
        &json!({
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
        }),
        20,
    );
    assert_eq!(job_done["status"], "succeeded", "bary job failed: {job_done}");
    let f = FitsFile::open(&bary_out).unwrap();
    assert_eq!(f.hdu(0).unwrap().header().get_str("TIMESYS"), Some("TDB"));
    let _ = TAU;

    server.stop();
}

// ===========================================================================
// TIER 1 — FEATURE COVERAGE (F1–F7, M3 E2E Gate)
// ===========================================================================

#[test]
#[ignore = "M3 E2E gate: run with -- --include-ignored"]
fn tier1_f1_f2_toml_parser_and_config_parse_api() {
    let dir = tempdir().unwrap();
    let server = spawn_test_server(dir.path());
    let addr = server.addr();

    // T1-F1-01 & T1-F2-01: Multi-section TOML via {"toml": "..."}
    let multi_toml = r#"
threads = 4

[czti]
input_dir = "/data/czti_l1"
caldb_dir = "/data/caldb/czti"
output_dir = "/data/out_czti"
output_prefix = "run#1" # inline comment after quoted '#'
oversampling = 32
quads = [0, 1, 2, 3]
generate_eventfile = true
mask_weight = true

[laxpc]
input_dir = "/data/laxpc_l1"
output_dir = "/data/out_laxpc"
units = ["10", "20"]
energy_range = '3:20'
deadtime_sec = 0.000043
spectrum_channels = 256

[bary]
orbit_file = "/data/sat.orb"
ra_deg = 83.6331 # Crab RA
dec_deg = 22.0145
frame = "icrs"
write_checksums = true

[timing]
freq_hz = 29.658051471
fdot = -3.68e-10
fddot = 1.2e-20
harmonics = 2
nbins = 64

[ui]
host = "127.0.0.1"
port = 8085
open_browser = false
"#;

    let r1 = http_request(
        addr,
        "POST",
        "/api/config/parse",
        Some(&json!({ "toml": multi_toml }).to_string()),
    );
    assert_eq!(r1.status, 200, "T1-F2-01 failed: {}", r1.body);
    let v1: Value = serde_json::from_str(&r1.body).unwrap();
    assert_eq!(v1["config"]["threads"], 4);
    assert_eq!(v1["config"]["czti"]["output_prefix"], "run#1");
    assert_eq!(v1["config"]["czti"]["oversampling"], 32);
    assert_eq!(v1["config"]["czti"]["quads"], json!([0, 1, 2, 3]));
    assert_eq!(v1["config"]["czti"]["generate_eventfile"], true);
    assert_eq!(v1["config"]["laxpc"]["units"], json!(["10", "20"]));
    assert_eq!(v1["config"]["laxpc"]["energy_range"], "3:20");
    assert_eq!(v1["config"]["laxpc"]["spectrum_channels"], 256);
    assert_eq!(v1["config"]["bary"]["orbit_file"], "/data/sat.orb");
    assert!((v1["config"]["bary"]["ra_deg"].as_f64().unwrap() - 83.6331).abs() < 1e-9);
    assert!((v1["config"]["timing"]["fddot"].as_f64().unwrap() - 1.2e-20).abs() < 1e-26);
    let canonical = v1["canonical_toml"].as_str().expect("canonical_toml string");
    assert!(!canonical.is_empty());

    // T1-F1-03 & T1-F1-04 & T1-F1-05: Escapes, underscores, dotted keys
    let tricky_toml = r#"
threads = 2
czti.caldb_dir = "C:\\data\\\"test\""
czti.time_binsize = 1_000.5
bary.ra_deg = 83.633_1
timing.nbins = 32
"#;
    let r2 = http_request(
        addr,
        "POST",
        "/api/config/parse",
        Some(&json!({ "content": tricky_toml, "format": "toml" }).to_string()),
    );
    assert_eq!(r2.status, 200, "T1-F2-02 failed: {}", r2.body);
    let v2: Value = serde_json::from_str(&r2.body).unwrap();
    assert_eq!(v2["config"]["threads"], 2);
    assert_eq!(v2["config"]["czti"]["caldb_dir"], r#"C:\data\"test""#);
    assert!((v2["config"]["czti"]["time_binsize"].as_f64().unwrap() - 1000.5).abs() < 1e-9);
    assert!((v2["config"]["bary"]["ra_deg"].as_f64().unwrap() - 83.6331).abs() < 1e-9);

    // T1-F2-03: JSON format input via {"content": "...", "format": "json"}
    let json_cfg = r#"{"threads": 3, "bary": {"frame": "fk5", "dec_deg": 22.0145}}"#;
    let r3 = http_request(
        addr,
        "POST",
        "/api/config/parse",
        Some(&json!({ "content": json_cfg, "format": "json" }).to_string()),
    );
    assert_eq!(r3.status, 200, "T1-F2-03 failed: {}", r3.body);
    let v3: Value = serde_json::from_str(&r3.body).unwrap();
    assert_eq!(v3["config"]["threads"], 3);
    assert_eq!(v3["config"]["bary"]["frame"], "fk5");

    // T1-F2-04: Raw text/plain & application/toml body
    for ct in ["text/plain", "application/toml"] {
        let r4 = http_request_with_headers(
            addr,
            "POST",
            "/api/config/parse",
            Some(ct),
            Some("threads = 2\n[timing]\nharmonics = 4\n"),
        );
        assert_eq!(r4.status, 200, "raw {ct} failed: {}", r4.body);
        let v4: Value = serde_json::from_str(&r4.body).unwrap();
        assert_eq!(v4["config"]["threads"], 2);
        assert_eq!(v4["config"]["timing"]["harmonics"], 4);
    }

    // T1-F2-05: Round-trip re-parse of canonical_toml
    let r5 = http_request(
        addr,
        "POST",
        "/api/config/parse",
        Some(&json!({ "toml": canonical }).to_string()),
    );
    assert_eq!(r5.status, 200, "canonical round-trip failed: {}", r5.body);
    let v5: Value = serde_json::from_str(&r5.body).unwrap();
    assert_eq!(v5["config"], v1["config"]);
    assert_eq!(v5["sections"], v1["sections"]);

    server.stop();
}

#[test]
#[ignore = "M3 E2E gate: run with -- --include-ignored"]
fn tier1_f3_studio_default_and_classic_routing() {
    let dir = tempdir().unwrap();
    let server = spawn_test_server(dir.path());
    let addr = server.addr();

    // T1-F3-01: Studio entrypoint aliases
    let mut studio_bodies = Vec::new();
    for route in ["/", "/index.html", "/studio", "/studio/", "/studio/index.html"] {
        let r = http_request(addr, "GET", route, None);
        assert_eq!(r.status, 200, "Studio route {route} failed");
        assert!(r.content_type.contains("text/html"));
        assert!(r.body.contains("/studio/studio.js"), "{route} must serve Studio");
        assert!(r.body.contains("/classic"), "{route} must link to /classic");
        assert_offline_asset(&r.body, route);
        studio_bodies.push(r.body);
    }
    for b in &studio_bodies[1..] {
        assert_eq!(b, &studio_bodies[0], "Studio route alias body mismatch");
    }

    // T1-F3-02: Studio static assets
    for (route, ct) in [
        ("/studio/studio.css", "text/css"),
        ("/studio/studio.js", "application/javascript"),
    ] {
        let r = http_request(addr, "GET", route, None);
        assert_eq!(r.status, 200);
        assert!(r.content_type.contains(ct));
        assert_offline_asset(&r.body, route);
    }

    // T1-F3-03: Classic entrypoint aliases
    let mut classic_bodies = Vec::new();
    for route in ["/classic", "/classic/", "/classic/index.html"] {
        let r = http_request(addr, "GET", route, None);
        assert_eq!(r.status, 200, "Classic route {route} failed");
        assert!(r.content_type.contains("text/html"));
        assert!(r.body.contains("/app.js"), "{route} must serve Classic");
        assert!(r.body.contains("version-pill"));
        assert!(r.body.contains("runner-autofill"));
        assert!(r.body.contains("href=\"/\""), "{route} must link back to Studio /");
        assert_offline_asset(&r.body, route);
        classic_bodies.push(r.body);
    }
    for b in &classic_bodies[1..] {
        assert_eq!(b, &classic_bodies[0], "Classic route alias body mismatch");
    }

    // T1-F3-04: Classic asset pairs (/styles.css vs /classic/styles.css etc.)
    for (root_route, classic_route, ct) in [
        ("/styles.css", "/classic/styles.css", "text/css"),
        ("/charts.js", "/classic/charts.js", "application/javascript"),
        ("/app.js", "/classic/app.js", "application/javascript"),
    ] {
        let a = http_request(addr, "GET", root_route, None);
        let b = http_request(addr, "GET", classic_route, None);
        assert_eq!(a.status, 200);
        assert_eq!(b.status, 200);
        assert!(a.content_type.contains(ct));
        assert!(b.content_type.contains(ct));
        assert_eq!(a.body, b.body, "{root_route} and {classic_route} must match");
        assert_offline_asset(&a.body, root_route);
    }

    // T1-F3-05: Method & unknown route handling (including GET /api/config/parse -> 405)
    let r_405 = http_request(addr, "GET", "/api/config/parse", None);
    assert_eq!(r_405.status, 405);
    let r_404 = http_request(addr, "GET", "/nonexistent-route", None);
    assert_eq!(r_404.status, 404);

    server.stop();
}

#[test]
#[ignore = "M3 E2E gate: run with -- --include-ignored"]
fn tier1_f4_czti_addspec_backend_coadd_and_ogip() {
    let dir = tempdir().unwrap();
    let pha1 = dir.path().join("q0.pha");
    let pha2 = dir.path().join("q1.pha");
    let rsp1 = dir.path().join("q0.rsp");
    let rsp2 = dir.path().join("q1.rsp");
    let out_pha_avg = dir.path().join("comb_avg.pha");
    let out_pha_sum = dir.path().join("comb_sum.pha");
    let out_rsp_sum = dir.path().join("comb_sum.rsp");

    let mut c1 = vec![2.0_f32; 512];
    let mut c2 = vec![3.0_f32; 512];
    c1[20] = 10.0;
    c2[20] = 15.0;
    let e1 = vec![1.0_f32; 512];
    let e2 = vec![2.0_f32; 512];
    write_synthetic_czti_pha_512(&pha1, &c1, &e1, 100.0, "q0.rsp", Some(0));
    write_synthetic_czti_pha_512(&pha2, &c2, &e2, 200.0, "q1.rsp", Some(1));
    write_synthetic_czti_rsp_512(&rsp1, 10.0);
    write_synthetic_czti_rsp_512(&rsp2, 20.0);

    let server = spawn_test_server(dir.path());
    let addr = server.addr();

    // T1-F4-01: CLI preview without RSP / without --add-exposure
    let p1 = http_request(
        addr,
        "POST",
        "/api/cli-preview",
        Some(
            &json!({
                "pipeline": "czti-addspec",
                "czti_addspec": {
                    "pha": [pha1.to_string_lossy(), pha2.to_string_lossy()],
                    "rsp": [],
                    "out_pha": out_pha_avg.to_string_lossy(),
                    "add_exposure": false
                }
            })
            .to_string(),
        ),
    );
    assert_eq!(p1.status, 200, "addspec preview 1 failed: {}", p1.body);
    let pv1: Value = serde_json::from_str(&p1.body).unwrap();
    let cmd1 = pv1["cli_command"].as_str().unwrap();
    assert!(cmd1.contains("astrosat czti addspec"), "{cmd1}");
    assert!(cmd1.contains("--pha"), "{cmd1}");
    assert!(cmd1.contains("-o"), "{cmd1}");
    assert!(!cmd1.contains("--out-rsp"), "{cmd1}");
    assert!(!cmd1.contains("--add-exposure"), "{cmd1}");

    // T1-F4-02: CLI preview with RSP & --add-exposure
    let p2 = http_request(
        addr,
        "POST",
        "/api/cli-preview",
        Some(
            &json!({
                "pipeline": "czti-addspec",
                "czti_addspec": {
                    "pha": [pha1.to_string_lossy(), pha2.to_string_lossy()],
                    "rsp": [rsp1.to_string_lossy(), rsp2.to_string_lossy()],
                    "out_pha": out_pha_sum.to_string_lossy(),
                    "out_rsp": out_rsp_sum.to_string_lossy(),
                    "add_exposure": true
                }
            })
            .to_string(),
        ),
    );
    assert_eq!(p2.status, 200, "addspec preview 2 failed: {}", p2.body);
    let cmd2 = serde_json::from_str::<Value>(&p2.body).unwrap()["cli_command"]
        .as_str()
        .unwrap()
        .to_string();
    assert!(cmd2.contains("--rsp"), "{cmd2}");
    assert!(cmd2.contains("--out-rsp"), "{cmd2}");
    assert!(cmd2.contains("--add-exposure"), "{cmd2}");

    // T1-F4-03: Co-add PHAs only (exposure averaged: (100+200)/2 = 150.0, counts summed)
    let j1 = submit_and_wait(
        addr,
        &json!({
            "pipeline": "czti-addspec",
            "czti_addspec": {
                "pha": [pha1.to_string_lossy(), pha2.to_string_lossy()],
                "rsp": [],
                "out_pha": out_pha_avg.to_string_lossy(),
                "add_exposure": false
            }
        }),
        30,
    );
    assert_eq!(j1["status"], "succeeded", "addspec avg job failed: {j1}");
    assert!(out_pha_avg.is_file());
    let f_avg = FitsFile::open(&out_pha_avg).unwrap();
    let avg_shdu = f_avg.hdu_by_name("SPECTRUM").unwrap();
    let exp_avg = avg_shdu.header().get_f64("EXPOSURE").unwrap();
    assert!((exp_avg - 150.0).abs() < 1e-6, "expected avg exp 150, got {exp_avg}");
    let avg_counts = avg_shdu.as_bintable().unwrap().read_col_f64("COUNTS").unwrap();
    assert_eq!(avg_counts.len(), 512);
    assert!((avg_counts[20] - 25.0).abs() < 1e-6, "expected c1[20]+c2[20]=25, got {}", avg_counts[20]);
    assert!((avg_counts[0] - 5.0).abs() < 1e-6, "expected c1[0]+c2[0]=5, got {}", avg_counts[0]);
    assert!(j1["outputs"]
        .as_array()
        .unwrap()
        .iter()
        .any(|o| o["product"] == "spectrum"));

    // T1-F4-04: Co-add PHAs + RSPs with add_exposure: true (sum = 300.0, RSP exposure-weighted)
    let j2 = submit_and_wait(
        addr,
        &json!({
            "pipeline": "czti-addspec",
            "czti_addspec": {
                "pha": [pha1.to_string_lossy(), pha2.to_string_lossy()],
                "rsp": [rsp1.to_string_lossy(), rsp2.to_string_lossy()],
                "out_pha": out_pha_sum.to_string_lossy(),
                "out_rsp": out_rsp_sum.to_string_lossy(),
                "add_exposure": true
            }
        }),
        30,
    );
    assert_eq!(j2["status"], "succeeded", "addspec sum job failed: {j2}");
    assert!(out_pha_sum.is_file());
    assert!(out_rsp_sum.is_file());
    let f_sum = FitsFile::open(&out_pha_sum).unwrap();
    let shdu = f_sum.hdu_by_name("SPECTRUM").unwrap();
    let exp_sum = shdu.header().get_f64("EXPOSURE").unwrap();
    assert!((exp_sum - 300.0).abs() < 1e-6, "expected sum exp 300, got {exp_sum}");
    assert_eq!(shdu.header().get_str("RESPFILE"), Some("comb_sum.rsp"));
    let sum_counts = shdu.as_bintable().unwrap().read_col_f64("COUNTS").unwrap();
    assert!((sum_counts[20] - 25.0).abs() < 1e-6, "expected c1[20]+c2[20]=25, got {}", sum_counts[20]);
    assert!((sum_counts[0] - 5.0).abs() < 1e-6, "expected c1[0]+c2[0]=5, got {}", sum_counts[0]);

    let (_rmeta, comb_rsp) = astrosat_czti::fits_io::read_response_file(&out_rsp_sum).unwrap();
    let expected_area_bin0 = (100.0 * 10.0 + 200.0 * 20.0) / 300.0;
    for (ch, &val) in comb_rsp.matrix.iter().enumerate().take(15).skip(10) {
        let actual = val as f64;
        assert!(
            (actual - expected_area_bin0).abs() < 1e-3,
            "expected exposure-weighted RSP[{ch}]={expected_area_bin0}, got {actual}"
        );
    }
    for ch in 90..95 {
        let actual = comb_rsp.matrix[512 + ch] as f64;
        let exp1 = expected_area_bin0 * 0.8;
        assert!(
            (actual - exp1).abs() < 1e-3,
            "expected exposure-weighted RSP bin1[{ch}]={exp1}, got {actual}"
        );
    }

    // T1-F4-05: Co-added spectrum readable via /api/spectrum & /api/validate
    let sp_r = http_request(
        addr,
        "GET",
        &format!(
            "/api/spectrum?path={}&good_only=false",
            url_encode(&out_pha_sum.to_string_lossy())
        ),
        None,
    );
    assert_eq!(sp_r.status, 200, "spectrum read failed: {}", sp_r.body);
    let sp_v: Value = serde_json::from_str(&sp_r.body).unwrap();
    assert_eq!(sp_v["total_channels"], 512);
    assert_eq!(sp_v["channel"][20], 20);
    let sp_ch20 = sp_v["rate_or_counts"][20].as_f64().unwrap();
    let expected_sp20 = if sp_v["is_rate"].as_bool().unwrap_or(false) {
        25.0 / 300.0
    } else {
        25.0
    };
    assert!((sp_ch20 - expected_sp20).abs() < 1e-4, "expected ch20={expected_sp20}, got {sp_ch20}");

    let val_r = http_request(
        addr,
        "GET",
        &format!(
            "/api/validate?path={}&strict=true",
            url_encode(&out_pha_sum.to_string_lossy())
        ),
        None,
    );
    assert_eq!(val_r.status, 200);
    let val_v: Value = serde_json::from_str(&val_r.body).unwrap();
    assert_eq!(val_v["valid"], true, "strict validate failed: {val_v}");

    server.stop();
}

#[test]
#[ignore = "M3 E2E gate: run with -- --include-ignored"]
fn tier1_f5_laxpc_backshift_backend_preview_metrics_and_pha() {
    let dir = tempdir().unwrap();
    let caldb = dir.path().join("caldb_laxpc");
    write_synthetic_laxpc_caldb(&caldb, true);
    let bkg_pha = dir.path().join("backlxp10.pha");
    // Use day 120.0 (not in builtin catalog) so background peak interpolates cleanly from gaina8.dat
    write_synthetic_laxpc_bkg_pha(&bkg_pha, 120.0);
    let obs_fits = dir.path().join("lxp1_l2.fits");
    // mission_day = 500 => TSTART = (500 + 2096) * 86400 => non-zero delta_ch vs bkg day 120
    write_synthetic_tstart_fits(&obs_fits, (500.0 + 2096.0) * 86_400.0);
    let out_pha = dir.path().join("shifted_bkg.pha");

    let server = spawn_test_server(dir.path());
    let addr = server.addr();

    // T1-F5-01: CLI preview with mission_day
    let p1 = http_request(
        addr,
        "POST",
        "/api/cli-preview",
        Some(
            &json!({
                "pipeline": "laxpc-backshift",
                "laxpc_backshift": {
                    "unit": "10",
                    "caldb_dir": caldb.to_string_lossy(),
                    "mission_day": 100.0,
                    "layer1_only": true
                }
            })
            .to_string(),
        ),
    );
    assert_eq!(p1.status, 200, "backshift preview 1 failed: {}", p1.body);
    let cmd1 = serde_json::from_str::<Value>(&p1.body).unwrap()["cli_command"]
        .as_str()
        .unwrap()
        .to_string();
    assert!(cmd1.contains("astrosat laxpc backshift"), "{cmd1}");
    assert!(cmd1.contains("--unit 10"), "{cmd1}");
    assert!(cmd1.contains("--mission-day 100"), "{cmd1}");
    assert!(cmd1.contains("--layer1-only"), "{cmd1}");

    // T1-F5-02: CLI preview with fits + bkg_pha + output + filters
    let p2 = http_request(
        addr,
        "POST",
        "/api/cli-preview",
        Some(
            &json!({
                "pipeline": "laxpc-backshift",
                "laxpc_backshift": {
                    "unit": "lxp2",
                    "caldb_dir": caldb.to_string_lossy(),
                    "fits": obs_fits.to_string_lossy(),
                    "bkg_pha": bkg_pha.to_string_lossy(),
                    "bkg_peak_ch": 223.6,
                    "output": out_pha.to_string_lossy(),
                    "single_only": true,
                    "pha_min": 10,
                    "pha_max": 500
                }
            })
            .to_string(),
        ),
    );
    assert_eq!(p2.status, 200, "backshift preview 2 failed: {}", p2.body);
    let cmd2 = serde_json::from_str::<Value>(&p2.body).unwrap()["cli_command"]
        .as_str()
        .unwrap()
        .to_string();
    for flag in [
        "--fits",
        "--bkg-pha",
        "--bkg-peak-ch",
        "-o",
        "--single-only",
        "--pha-min 10",
        "--pha-max 500",
    ] {
        assert!(cmd2.contains(flag), "missing '{flag}' in '{cmd2}'");
    }

    // T1-F5-03: Gain-shift metrics only (mission_day mode, no bkg_pha)
    let j1 = submit_and_wait(
        addr,
        &json!({
            "pipeline": "laxpc-backshift",
            "laxpc_backshift": {
                "unit": "10",
                "caldb_dir": caldb.to_string_lossy(),
                "mission_day": 500.0
            }
        }),
        20,
    );
    assert_eq!(j1["status"], "succeeded", "backshift metrics failed: {j1}");
    let s1 = &j1["summary"];
    assert!((s1["source_peak_ch"].as_f64().unwrap() - 240.0).abs() < 1e-6);
    assert!((s1["delta_ch_30kev"].as_f64().unwrap() - 0.0).abs() < 1e-6);
    assert!((s1["gain_ratio"].as_f64().unwrap() - 1.0).abs() < 1e-6);
    assert!((s1["secular_scale"].as_f64().unwrap() - 1.0).abs() < 1e-6);
    let rmf1 = s1["recommended_rmf"].as_str().unwrap();
    assert!(rmf1.starts_with("lx10") && rmf1.ends_with(".rmf"), "unexpected rmf1={rmf1}");

    // T1-F5-04: FITS TSTART (day 500) vs bkg_pha (day 120) + explicit bkg_peak_ch + shifted .pha output
    let j2 = submit_and_wait(
        addr,
        &json!({
            "pipeline": "laxpc-backshift",
            "laxpc_backshift": {
                "unit": "10",
                "caldb_dir": caldb.to_string_lossy(),
                "fits": obs_fits.to_string_lossy(),
                "bkg_pha": bkg_pha.to_string_lossy(),
                "bkg_peak_ch": 225.0,
                "output": out_pha.to_string_lossy(),
                "pha_min": 0,
                "pha_max": 1023
            }
        }),
        20,
    );
    assert_eq!(j2["status"], "succeeded", "backshift shifted pha failed: {j2}");
    let s2 = &j2["summary"];
    let dch = s2["delta_ch_30kev"].as_f64().unwrap();
    let gr = s2["gain_ratio"].as_f64().unwrap();
    let sec = s2["secular_scale"].as_f64().unwrap();
    let fsh = s2["lcbk_shift_fraction"].as_f64().unwrap();
    assert!((dch - 15.0).abs() < 1e-4, "expected delta_ch_30kev=15.0, got {dch}");
    assert!((gr - (240.0 / 225.0)).abs() < 1e-5, "expected gain_ratio=1.0667, got {gr}");
    assert!(sec.is_finite() && (sec - 1.0).abs() > 0.01, "expected non-unit secular_scale, got {sec}");
    assert!(fsh.is_finite() && fsh.abs() > 1e-4, "expected non-zero lcbk_shift_fraction, got {fsh}");
    let rmf2 = s2["recommended_rmf"].as_str().unwrap();
    assert!(rmf2.starts_with("lx10") && rmf2.ends_with(".rmf"), "unexpected rmf2={rmf2}");
    assert!(out_pha.is_file());
    assert!(j2["outputs"]
        .as_array()
        .unwrap()
        .iter()
        .any(|o| o["product"] == "spectrum"));

    // T1-F5-05: Shifted .pha readable via /api/spectrum & /api/inspect with non-zero shifted counts
    let sp_r = http_request(
        addr,
        "GET",
        &format!("/api/spectrum?path={}", url_encode(&out_pha.to_string_lossy())),
        None,
    );
    assert_eq!(sp_r.status, 200);
    let sp_v: Value = serde_json::from_str(&sp_r.body).unwrap();
    assert_eq!(sp_v["total_channels"], 512);
    let sum_spec: f64 = sp_v["rate_or_counts"]
        .as_array()
        .unwrap()
        .iter()
        .filter_map(|v| v.as_f64())
        .sum();
    assert!(sum_spec > 1.0, "expected non-zero shifted PHA counts/rate, got sum={sum_spec}");
    assert!(sp_v["total_rate"].as_f64().unwrap() > 0.0);

    let ins_r = http_request(
        addr,
        "GET",
        &format!("/api/inspect?path={}", url_encode(&out_pha.to_string_lossy())),
        None,
    );
    assert_eq!(ins_r.status, 200);
    let ins_v: Value = serde_json::from_str(&ins_r.body).unwrap();
    assert_eq!(ins_v["detnam"], "LX10");

    server.stop();
}

#[test]
#[ignore = "M3 E2E gate: run with -- --include-ignored"]
fn tier1_f6_advanced_pipeline_overrides_and_bary_mkf() {
    let dir = tempdir().unwrap();
    let utc_evt = dir.path().join("utc.evt");
    let dummy_orb = dir.path().join("dummy.orb");
    let mkf_path = dir.path().join("sat.mkf");
    let bary_out = dir.path().join("bary_mkf.evt");
    write_synthetic_utc_events_and_orbit(&utc_evt, &dummy_orb, "UTC");
    write_synthetic_mkf_for_bary(&mkf_path);

    let server = spawn_test_server(dir.path());
    let addr = server.addr();

    // T1-F6-01 & T1-F6-04: CZTI extended CLI/TOML preview + oversampling=64
    let czti_req = json!({
        "pipeline": "czti",
        "threads": 4,
        "czti": {
            "input_dir": "/data/czti_l1",
            "m0": "/data/m0.fits",
            "ss": "/data/ss.fits",
            "tct": "/data/tct.tct",
            "mkf": "/data/obs.mkf",
            "caldb_dir": "/data/caldb/czti",
            "output_dir": "/data/out_czti",
            "output_prefix": "crab",
            "mkf_threshold_file": "/data/mkfThresholds.txt",
            "user_gti_file": "/data/user_gti.txt",
            "generate_eventfile": true,
            "oversampling": 64,
            "quads": [0, 2]
        }
    });
    let r1 = http_request(addr, "POST", "/api/cli-preview", Some(&czti_req.to_string()));
    assert_eq!(r1.status, 200, "CZTI preview failed: {}", r1.body);
    let v1: Value = serde_json::from_str(&r1.body).unwrap();
    let cmd1 = v1["cli_command"].as_str().unwrap();
    for flag in [
        "--m0",
        "--ss",
        "--tct",
        "--mkf ",
        "--mkf-thresholds",
        "--user-gti",
        "--gen-evt",
        "--oversampling 64",
        "--quads 0,2",
    ] {
        assert!(cmd1.contains(flag), "CZTI preview missing '{flag}' in: {cmd1}");
    }
    let toml1 = v1["toml_config"].as_str().unwrap();
    let p_toml = http_request(
        addr,
        "POST",
        "/api/config/parse",
        Some(&json!({ "toml": toml1 }).to_string()),
    );
    assert_eq!(p_toml.status, 200, "CZTI toml parse failed: {}", p_toml.body);
    let pt: Value = serde_json::from_str(&p_toml.body).unwrap();
    assert_eq!(pt["config"]["czti"]["oversampling"], 64);
    assert_eq!(pt["config"]["czti"]["quads"], json!([0, 2]));
    assert_eq!(pt["config"]["czti"]["generate_eventfile"], true);
    assert_eq!(pt["config"]["czti"]["mkf_threshold_file"], "/data/mkfThresholds.txt");
    assert_eq!(pt["config"]["czti"]["user_gti_file"], "/data/user_gti.txt");

    // T1-F6-02: LAXPC extended CLI/TOML preview flags
    let laxpc_req = json!({
        "pipeline": "laxpc",
        "threads": 2,
        "laxpc": {
            "input_dir": "/data/laxpc_l1",
            "output_dir": "/data/out_laxpc",
            "caldb_dir": "/data/caldb/laxpc",
            "units": ["10", "20"],
            "deadtime_sec": 4.3e-5,
            "anodes": "1,2,3",
            "energy_range": "3:20",
            "pha_min": 5,
            "pha_max": 600,
            "spectrum_channels": 256,
            "saa_entry_lon_deg": 280.0,
            "saa_exit_lon_deg": 330.0,
            "merge_gap_sec": 2.0,
            "min_duration_sec": 10.0,
            "user_gti_file": "/data/laxpc_gti.txt",
            "user_gti_only": true,
            "write_ascii_products": true,
            "background_model_path": "/data/backfit.lxp10",
            "background_pha_path": "/data/backlxp10.pha",
            "background_peak_ch": 225.0
        }
    });
    let r2 = http_request(addr, "POST", "/api/cli-preview", Some(&laxpc_req.to_string()));
    assert_eq!(r2.status, 200, "LAXPC preview failed: {}", r2.body);
    let v2: Value = serde_json::from_str(&r2.body).unwrap();
    let cmd2 = v2["cli_command"].as_str().unwrap();
    for flag in [
        "--deadtime",
        "--anodes '1,2,3'",
        "--energy-range 3:20",
        "--pha-min 5",
        "--pha-max 600",
        "--channels 256",
        "--saa-entry-lon",
        "--saa-exit-lon",
        "--merge-gap",
        "--min-duration",
        "--user-gti",
        "--user-gti-only",
        "--ascii",
        "--bkg-model",
        "--bkg-pha",
        "--bkg-peak-ch",
    ] {
        assert!(cmd2.contains(flag), "LAXPC preview missing '{flag}' in: {cmd2}");
    }
    let cmd2_tokens = split_shell_words(cmd2);
    assert!(
        cmd2_tokens.windows(2).any(|w| w == ["--anodes", "1,2,3"]),
        "shell tokens missing [--anodes, 1,2,3]: {cmd2_tokens:?}"
    );

    // T1-F6 negative validation checks on POST /api/jobs/run
    let laxpc_l1_dummy = dir.path().join("laxpc_l1_dummy");
    std::fs::create_dir_all(&laxpc_l1_dummy).unwrap();
    let bad_laxpc_overrides = [
        ("invalid_anode_id", json!({
            "pipeline": "laxpc",
            "laxpc": { "input_dir": laxpc_l1_dummy.to_string_lossy(), "output_dir": dir.path().join("o").to_string_lossy(), "anodes": "1,99" }
        })),
        ("non_numeric_anode", json!({
            "pipeline": "laxpc",
            "laxpc": { "input_dir": laxpc_l1_dummy.to_string_lossy(), "output_dir": dir.path().join("o").to_string_lossy(), "anodes": "foo" }
        })),
        ("user_gti_only_without_file", json!({
            "pipeline": "laxpc",
            "laxpc": { "input_dir": laxpc_l1_dummy.to_string_lossy(), "output_dir": dir.path().join("o").to_string_lossy(), "user_gti_only": true }
        })),
        ("energy_range_without_caldb", json!({
            "pipeline": "laxpc",
            "laxpc": { "input_dir": laxpc_l1_dummy.to_string_lossy(), "output_dir": dir.path().join("o").to_string_lossy(), "energy_range": "3:20" }
        })),
        ("background_peak_ch_zero", json!({
            "pipeline": "laxpc",
            "laxpc": { "input_dir": laxpc_l1_dummy.to_string_lossy(), "output_dir": dir.path().join("o").to_string_lossy(), "background_peak_ch": 0.0 }
        })),
        ("background_peak_ch_negative", json!({
            "pipeline": "laxpc",
            "laxpc": { "input_dir": laxpc_l1_dummy.to_string_lossy(), "output_dir": dir.path().join("o").to_string_lossy(), "background_peak_ch": -5.0 }
        })),
        ("merge_gap_sec_negative", json!({
            "pipeline": "laxpc",
            "laxpc": { "input_dir": laxpc_l1_dummy.to_string_lossy(), "output_dir": dir.path().join("o").to_string_lossy(), "merge_gap_sec": -1.0 }
        })),
    ];
    for (lbl, req) in bad_laxpc_overrides {
        let r_bad = http_request(addr, "POST", "/api/jobs/run", Some(&req.to_string()));
        assert_eq!(
            r_bad.status, 400,
            "negative LAXPC override {lbl} should return 400, got {}: {}",
            r_bad.status, r_bad.body
        );
        let j_bad: Value = serde_json::from_str(&r_bad.body).expect("400 JSON");
        assert!(j_bad["error"].is_string(), "missing error string for {lbl}");
    }

    // T1-F6-03: Barycenter --mkf mode (without --orbit)
    let bary_mkf_req = json!({
        "pipeline": "bary",
        "threads": 2,
        "bary": {
            "input": utc_evt.to_string_lossy(),
            "output": bary_out.to_string_lossy(),
            "mkf": mkf_path.to_string_lossy(),
            "ra_deg": 83.6331,
            "dec_deg": 22.0145,
            "frame": "icrs",
            "ephemeris": "analytical"
        }
    });
    let bp = http_request(addr, "POST", "/api/cli-preview", Some(&bary_mkf_req.to_string()));
    assert_eq!(bp.status, 200);
    let bcmd = serde_json::from_str::<Value>(&bp.body).unwrap()["cli_command"]
        .as_str()
        .unwrap()
        .to_string();
    assert!(bcmd.contains("--mkf"), "{bcmd}");
    assert!(!bcmd.contains("--orbit"), "{bcmd}");

    let jb = submit_and_wait(addr, &bary_mkf_req, 20);
    assert_eq!(jb["status"], "succeeded", "bary --mkf failed: {jb}");
    assert!(bary_out.is_file());
    assert!(!bary_out.with_extension("tmp_synth.orb").exists());
    let fb = FitsFile::open(&bary_out).unwrap();
    assert_eq!(fb.hdu(0).unwrap().header().get_str("TIMESYS"), Some("TDB"));

    // T1-F6-05: Job ordering (newest first) & GET /api/jobs/{id}
    let j_id = jb["id"].as_str().unwrap();
    let list_r = http_request(addr, "GET", "/api/jobs", None);
    assert_eq!(list_r.status, 200);
    let lv: Value = serde_json::from_str(&list_r.body).unwrap();
    assert_eq!(lv["jobs"][0]["id"], j_id);

    server.stop();
}

#[test]
#[ignore = "M3 E2E gate: run with -- --include-ignored"]
fn tier1_f7_pulse_timing_cli_parity_gti_df_hz_and_fddot() {
    let dir = tempdir().unwrap();
    let evt_path = dir.path().join("pulsar_fddot.evt");
    let gti_path = dir.path().join("subset.gti");

    let epoch = 190_000_000.0;
    let f0 = 29.658;
    let fdot = -3.68e-10;
    // Choose span ~84.3s (250 pulses * stride 10) and fddot=1.0e-5 so (1/6)*fddot*T^3 ≈ 1.0 cycle >= 0.1
    let fddot = 1.0e-5;
    write_planted_pulsar_events_ext(&evt_path, epoch, f0, fdot, fddot, 0.2, 250, 10);
    // GTI covers only first ~30 seconds (excluding later cycles)
    write_synthetic_gti(&gti_path, &[epoch - 0.1], &[epoch + 30.0]);

    let server = spawn_test_server(dir.path());
    let addr = server.addr();

    // T1-F7-01: period_s search mode (without freq_hz)
    let r1 = http_request(
        addr,
        "POST",
        "/api/pulse/search",
        Some(
            &json!({
                "path": evt_path.to_string_lossy(),
                "period_s": 1.0 / f0,
                "search_span_hz": 0.02,
                "n_trials": 41,
                "fdot": fdot,
                "fddot": fddot,
                "epoch": epoch,
                "harmonics": 2,
                "nbins": 32
            })
            .to_string(),
        ),
    );
    assert_eq!(r1.status, 200, "period_s search failed: {}", r1.body);
    let v1: Value = serde_json::from_str(&r1.body).unwrap();
    assert!((v1["best_freq_hz"].as_f64().unwrap() - f0).abs() < 0.005);
    assert!((v1["best_period_sec"].as_f64().unwrap() - 1.0 / f0).abs() < 1e-4);

    // T1-F7-02: freq_min / freq_max + explicit df_hz
    let r2 = http_request(
        addr,
        "POST",
        "/api/pulse/search",
        Some(
            &json!({
                "path": evt_path.to_string_lossy(),
                "freq_min": 29.648,
                "freq_max": 29.668,
                "df_hz": 0.001,
                "epoch": epoch,
                "harmonics": 2,
                "nbins": 32
            })
            .to_string(),
        ),
    );
    assert_eq!(r2.status, 200, "df_hz search failed: {}", r2.body);
    let v2: Value = serde_json::from_str(&r2.body).unwrap();
    let nt = v2["n_trials"].as_u64().unwrap();
    assert!((20..=22).contains(&nt), "expected ~21 trials for df=0.001, got {nt}");

    // T1-F7-03: gti FITS filtering reduces num_events strictly to GTI window
    let total_events = v1["num_events"].as_u64().unwrap();
    let r3 = http_request(
        addr,
        "POST",
        "/api/pulse/search",
        Some(
            &json!({
                "path": evt_path.to_string_lossy(),
                "gti": gti_path.to_string_lossy(),
                "freq_hz": f0,
                "search_span_hz": 0.02,
                "n_trials": 21,
                "epoch": epoch
            })
            .to_string(),
        ),
    );
    assert_eq!(r3.status, 200, "gti search failed: {}", r3.body);
    let v3: Value = serde_json::from_str(&r3.body).unwrap();
    let gti_events = v3["num_events"].as_u64().unwrap();
    assert!(gti_events >= 4 && gti_events < total_events, "gti_events={gti_events}, total={total_events}");

    // T1-F7-04: fddot + pi_min/pi_max + unit/layer filter + mathematical sensitivity check
    let r4 = http_request(
        addr,
        "POST",
        "/api/pulse/search",
        Some(
            &json!({
                "path": evt_path.to_string_lossy(),
                "freq_hz": f0,
                "search_span_hz": 0.01,
                "n_trials": 21,
                "fdot": fdot,
                "fddot": fddot,
                "epoch": epoch,
                "pi_min": 0,
                "pi_max": 50,
                "unit": 1,
                "layer": 1
            })
            .to_string(),
        ),
    );
    assert_eq!(r4.status, 200, "fddot+unit/layer search failed: {}", r4.body);
    let v4: Value = serde_json::from_str(&r4.body).unwrap();
    assert!((v4["fddot"].as_f64().unwrap() - fddot).abs() < 1e-20);
    assert_eq!(v4["num_events"].as_u64().unwrap(), total_events / 2);
    let t_span = v4["t_span_sec"].as_f64().unwrap();
    let cubic_cycles = (1.0 / 6.0) * fddot.abs() * t_span.powi(3);
    assert!(
        cubic_cycles >= 0.1,
        "expected cubic phase drift >= 0.1 cycles for fddot sensitivity, got {cubic_cycles} (T={t_span})"
    );
    let r4_no_fddot = http_request(
        addr,
        "POST",
        "/api/pulse/search",
        Some(
            &json!({
                "path": evt_path.to_string_lossy(),
                "freq_hz": f0,
                "search_span_hz": 0.01,
                "n_trials": 21,
                "fdot": fdot,
                "fddot": 0.0,
                "epoch": epoch,
                "pi_min": 0,
                "pi_max": 50,
                "unit": 1,
                "layer": 1
            })
            .to_string(),
        ),
    );
    assert_eq!(r4_no_fddot.status, 200);
    let v4_no: Value = serde_json::from_str(&r4_no_fddot.body).unwrap();
    let z2_with = v4["max_z2_n"].as_f64().unwrap();
    let z2_without = v4_no["max_z2_n"].as_f64().unwrap();
    assert!(
        z2_with > z2_without,
        "expected z2_with_fddot ({z2_with}) > z2_without_fddot ({z2_without})"
    );

    // T1-F7-05: multiband with fddot + unit/layer/pi_min/pi_max
    let r5 = http_request(
        addr,
        "POST",
        "/api/pulse/multiband",
        Some(
            &json!({
                "inputs": [{
                    "label": "LAXPC",
                    "path": evt_path.to_string_lossy(),
                    "bands": [
                        { "label": "U1 L1", "unit": 1, "layer": 1, "pi_min": 0, "pi_max": 50 },
                        { "label": "U2 L2", "unit": 2, "layer": 2, "pi_min": 100, "pi_max": 200 }
                    ]
                }],
                "freq_hz": f0,
                "fdot": fdot,
                "fddot": fddot,
                "epoch": epoch,
                "nbins": 64,
                "harmonics": 2
            })
            .to_string(),
        ),
    );
    assert_eq!(r5.status, 200, "multiband fddot failed: {}", r5.body);
    let v5: Value = serde_json::from_str(&r5.body).unwrap();
    assert!((v5["fddot"].as_f64().unwrap() - fddot).abs() < 1e-20);
    let lag = v5["bands"][1]["phase_lag_cycles"].as_f64().unwrap();
    assert!((lag - 0.2).abs() < 0.03, "expected ~0.2 lag, got {lag}");

    server.stop();
}

// ===========================================================================
// TIER 2 — BOUNDARY & CORNER CASES (M3 E2E Gate)
// ===========================================================================

#[test]
#[ignore = "M3 E2E gate: run with -- --include-ignored"]
fn tier2_f1_f2_toml_malformed_and_boundary_rejections() {
    let dir = tempdir().unwrap();
    let server = spawn_test_server(dir.path());
    let addr = server.addr();

    let bad_tomls = [
        // T2-01: array of tables
        ("array_of_tables", "[[czti]]\nthreads = 1\n"),
        ("empty_header", "[]\nthreads = 1\n"),
        // T2-02: unterminated string & bad escape
        ("unterminated_str", "[czti]\ninput_dir = \"/abc\n"),
        ("bad_escape", "[czti]\ninput_dir = \"/abc\\q\"\n"),
        ("unquoted_val", "[czti]\ninput_dir = unquoted_path\n"),
        // T2-03: scalar-table conflict & non-finite floats
        ("scalar_table_conflict", "czti = 5\n[czti]\nra_deg = 10.0\n"),
        ("nan_float", "[czti]\nra_deg = NaN\n"),
        ("inf_float", "[czti]\nra_deg = inf\n"),
        // T2-04: schema type mismatch
        ("schema_mismatch", "threads = \"not_a_number\"\n"),
    ];

    for (label, toml) in bad_tomls {
        let r = http_request(
            addr,
            "POST",
            "/api/config/parse",
            Some(&json!({ "toml": toml }).to_string()),
        );
        assert_eq!(r.status, 400, "{label} should return 400, got {}: {}", r.status, r.body);
        let j: Value = serde_json::from_str(&r.body).expect("error JSON");
        assert!(j["error"].is_string());
    }

    // Empty body -> 400; whitespace-only valid empty TOML -> 200; oversized > 4 MiB body -> 400/413 (or early close)
    let r_empty = http_request_with_headers(addr, "POST", "/api/config/parse", Some("text/plain"), Some(""));
    assert_eq!(r_empty.status, 400, "empty body must return 400, got {}: {}", r_empty.status, r_empty.body);

    let r_ws = http_request_with_headers(addr, "POST", "/api/config/parse", Some("text/plain"), Some("   "));
    assert_eq!(r_ws.status, 200, "whitespace TOML should parse as empty config, got {}: {}", r_ws.status, r_ws.body);

    let big = format!("threads = 1\n# {}\n", "x".repeat(4 * 1024 * 1024 + 512));
    if let Some(r_big) = http_request_allow_early_close(addr, "POST", "/api/config/parse", Some("text/plain"), &big) {
        assert!(
            matches!(r_big.status, 400 | 413),
            "oversized body returned unexpected status {}: {}",
            r_big.status,
            r_big.body
        );
    }

    server.stop();
}

#[test]
#[ignore = "M3 E2E gate: run with -- --include-ignored"]
fn tier2_f4_f5_addspec_and_backshift_boundaries() {
    let dir = tempdir().unwrap();
    let p1 = dir.path().join("a.pha");
    let p2 = dir.path().join("b.pha");
    let r1 = dir.path().join("a.rsp");
    let r2 = dir.path().join("b.rsp");
    write_synthetic_czti_pha_512(&p1, &[1.0; 512], &[0.5; 512], 50.0, "a.rsp", Some(0));
    write_synthetic_czti_pha_512(&p2, &[2.0; 512], &[0.5; 512], 50.0, "b.rsp", Some(1));
    write_synthetic_czti_rsp_512(&r1, 5.0);
    write_synthetic_czti_rsp_512(&r2, 5.0);

    let caldb_no_chan = dir.path().join("caldb_no_chan");
    write_synthetic_laxpc_caldb(&caldb_no_chan, false);
    let bkg_pha = dir.path().join("backlxp10.pha");
    write_synthetic_laxpc_bkg_pha(&bkg_pha, 100.0);
    let no_tstart_fits = dir.path().join("no_tstart.fits");
    write_synthetic_primary_without_tstart(&no_tstart_fits);
    let tstart_fits = dir.path().join("tstart.fits");
    write_synthetic_tstart_fits(&tstart_fits, (100.0 + 2096.0) * 86_400.0);

    let server = spawn_test_server(dir.path());
    let addr = server.addr();

    // T2-06..08, T2-10: czti-addspec invalid requests
    let out_o = dir.path().join("out.pha");
    let bad_addspecs = [
        json!({
            "pipeline": "czti-addspec",
            "czti_addspec": {
                "pha": [p1.to_string_lossy(), p2.to_string_lossy()],
                "rsp": [r1.to_string_lossy()],
                "out_pha": out_o.to_string_lossy()
            }
        }),
        json!({
            "pipeline": "czti-addspec",
            "czti_addspec": {
                "pha": [p1.to_string_lossy()],
                "rsp": [r1.to_string_lossy(), r2.to_string_lossy()],
                "out_pha": out_o.to_string_lossy()
            }
        }),
        json!({
            "pipeline": "czti-addspec",
            "czti_addspec": {
                "pha": [],
                "rsp": [],
                "out_pha": out_o.to_string_lossy()
            }
        }),
        json!({
            "pipeline": "czti-addspec",
            "czti_addspec": {
                "pha": [p1.to_string_lossy()],
                "rsp": [],
                "out_pha": ""
            }
        }),
        json!({
            "pipeline": "czti-addspec",
            "czti_addspec": {
                "pha": ["../../etc/passwd"],
                "rsp": [],
                "out_pha": out_o.to_string_lossy()
            }
        }),
        json!({
            "pipeline": "czti-addspec",
            "czti_addspec": {
                "pha": [dir.path().join("missing.pha").to_string_lossy()],
                "rsp": [],
                "out_pha": out_o.to_string_lossy()
            }
        }),
    ];
    for req in bad_addspecs {
        let res = http_request(addr, "POST", "/api/jobs/run", Some(&req.to_string()));
        assert!((400..500).contains(&res.status), "expected 4xx for {req}, got {}: {}", res.status, res.body);
    }

    // T2-09: out_rsp set when rsp is empty -> either 400 OR succeeds without creating a bogus .rsp
    let rsp_orphan = dir.path().join("orphan.rsp");
    let r_orphan = http_request(
        addr,
        "POST",
        "/api/jobs/run",
        Some(
            &json!({
                "pipeline": "czti-addspec",
                "czti_addspec": {
                    "pha": [p1.to_string_lossy()],
                    "rsp": [],
                    "out_pha": dir.path().join("single_out.pha").to_string_lossy(),
                    "out_rsp": rsp_orphan.to_string_lossy()
                }
            })
            .to_string(),
        ),
    );
    assert!(
        matches!(r_orphan.status, 400 | 202),
        "expected T2-09 status 400 or 202, got {}: {}",
        r_orphan.status,
        r_orphan.body
    );
    if r_orphan.status == 202 {
        let id = serde_json::from_str::<Value>(&r_orphan.body).unwrap()["id"]
            .as_str()
            .unwrap()
            .to_string();
        let _ = poll_job_until_done(addr, &id, 15);
    }
    assert!(!rsp_orphan.exists(), "must not create bogus orphan.rsp when rsp: []");

    // T2-11..15: laxpc-backshift invalid / boundary inputs
    let bad_backshifts = [
        // T2-11: invalid unit
        json!({
            "pipeline": "laxpc-backshift",
            "laxpc_backshift": { "unit": "40", "caldb_dir": caldb_no_chan.to_string_lossy(), "mission_day": 100.0 }
        }),
        // T2-12: missing caldb dir
        json!({
            "pipeline": "laxpc-backshift",
            "laxpc_backshift": { "unit": "10", "caldb_dir": dir.path().join("no_caldb").to_string_lossy(), "mission_day": 100.0 }
        }),
        // T2-13a: BOTH mission_day and fits
        json!({
            "pipeline": "laxpc-backshift",
            "laxpc_backshift": { "unit": "10", "caldb_dir": caldb_no_chan.to_string_lossy(), "mission_day": 100.0, "fits": tstart_fits.to_string_lossy() }
        }),
        // T2-13b: NEITHER mission_day nor fits
        json!({
            "pipeline": "laxpc-backshift",
            "laxpc_backshift": { "unit": "10", "caldb_dir": caldb_no_chan.to_string_lossy() }
        }),
        // T2-15b: invalid PHA range (pha_min > pha_max)
        json!({
            "pipeline": "laxpc-backshift",
            "laxpc_backshift": { "unit": "10", "caldb_dir": caldb_no_chan.to_string_lossy(), "mission_day": 100.0, "pha_min": 500, "pha_max": 100 }
        }),
        // T2-15c: pha_max > 1023
        json!({
            "pipeline": "laxpc-backshift",
            "laxpc_backshift": { "unit": "10", "caldb_dir": caldb_no_chan.to_string_lossy(), "mission_day": 100.0, "pha_max": 1024 }
        }),
    ];
    for req in bad_backshifts {
        let res = http_request(addr, "POST", "/api/jobs/run", Some(&req.to_string()));
        assert!((400..500).contains(&res.status), "expected 4xx for {req}, got {}: {}", res.status, res.body);
    }

    // T2-14a: Lxp30 after mission day 892 -> job fails with explicit message
    let j_lxp30 = submit_and_wait(
        addr,
        &json!({
            "pipeline": "laxpc-backshift",
            "laxpc_backshift": {
                "unit": "30",
                "caldb_dir": caldb_no_chan.to_string_lossy(),
                "mission_day": 900.0
            }
        }),
        15,
    );
    assert_eq!(j_lxp30["status"], "failed");
    assert!(j_lxp30["error"].as_str().unwrap().contains("892"));

    // T2-14b: FITS missing TSTART -> job fails with TSTART message
    let j_notstart = submit_and_wait(
        addr,
        &json!({
            "pipeline": "laxpc-backshift",
            "laxpc_backshift": {
                "unit": "10",
                "caldb_dir": caldb_no_chan.to_string_lossy(),
                "fits": no_tstart_fits.to_string_lossy()
            }
        }),
        15,
    );
    assert_eq!(j_notstart["status"], "failed");
    assert!(j_notstart["error"].as_str().unwrap().contains("TSTART"));

    // T2-15a: bkg_pha with CALDB lacking .channel -> fails cleanly with channel map context
    let j_nochan = submit_and_wait(
        addr,
        &json!({
            "pipeline": "laxpc-backshift",
            "laxpc_backshift": {
                "unit": "10",
                "caldb_dir": caldb_no_chan.to_string_lossy(),
                "mission_day": 100.0,
                "bkg_pha": bkg_pha.to_string_lossy()
            }
        }),
        15,
    );
    assert_eq!(j_nochan["status"], "failed");
    assert!(j_nochan["error"].as_str().unwrap().to_ascii_lowercase().contains("channel"));

    server.stop();
}

#[test]
#[ignore = "M3 E2E gate: run with -- --include-ignored"]
fn tier2_f6_f7_pipeline_and_pulse_boundary_cases() {
    let dir = tempdir().unwrap();
    let evt_path = dir.path().join("pulsar.evt");
    let tiny_evt = dir.path().join("tiny_pulsar.evt");
    let tdb_evt = dir.path().join("already_tdb.evt");
    let orb_path = dir.path().join("sat.orb");
    let mkf_path = dir.path().join("sat.mkf");
    let empty_gti = dir.path().join("empty.gti");
    let corrupt_gti = dir.path().join("corrupt.gti");
    let czti_l1_dummy = dir.path().join("czti_l1_dummy");
    std::fs::create_dir_all(&czti_l1_dummy).unwrap();
    // Create minimal valid modeM0 FITS, .tct, and .mkf inside czti_l1_dummy so stage/energy
    // validation is exercised directly rather than failing early on missing files.
    write_fits(
        czti_l1_dummy.join("AS1TESTcztM0_level1.fits"),
        &mut [Hdu::empty_primary()],
        true,
    )
    .unwrap();
    write_fits(
        czti_l1_dummy.join("AS1TEST.tct"),
        &mut [Hdu::empty_primary()],
        true,
    )
    .unwrap();
    write_fits(
        czti_l1_dummy.join("AS1TEST.mkf"),
        &mut [Hdu::empty_primary()],
        true,
    )
    .unwrap();
    std::fs::write(&corrupt_gti, b"NOT_A_FITS_FILE").unwrap();

    let epoch = 190_000_000.0;
    write_planted_pulsar_events(&evt_path, epoch, 29.658, 0.0, 0.0, 0.25);
    write_tiny_pulsar_events(&tiny_evt, epoch, 29.658);
    write_synthetic_utc_events_and_orbit(&tdb_evt, &orb_path, "TDB");
    write_synthetic_mkf_for_bary(&mkf_path);
    write_synthetic_gti(&empty_gti, &[epoch - 500.0], &[epoch - 400.0]);

    let server = spawn_test_server(dir.path());
    let addr = server.addr();

    // T2-16: CZTI oversampling exact boundaries 0, 1, 64, 65
    for bad_ov in [0_usize, 65] {
        let r = http_request(
            addr,
            "POST",
            "/api/jobs/run",
            Some(
                &json!({
                    "pipeline": "czti",
                    "czti": {
                        "input_dir": czti_l1_dummy.to_string_lossy(),
                        "caldb_dir": czti_l1_dummy.to_string_lossy(),
                        "output_dir": dir.path().join("out").to_string_lossy(),
                        "oversampling": bad_ov
                    }
                })
                .to_string(),
            ),
        );
        assert_eq!(r.status, 400, "oversampling={bad_ov} must be rejected: {}", r.body);
    }
    for ok_ov in [1_usize, 64] {
        let r = http_request(
            addr,
            "POST",
            "/api/cli-preview",
            Some(
                &json!({
                    "pipeline": "czti",
                    "czti": {
                        "input_dir": czti_l1_dummy.to_string_lossy(),
                        "caldb_dir": czti_l1_dummy.to_string_lossy(),
                        "output_dir": dir.path().join("out").to_string_lossy(),
                        "oversampling": ok_ov
                    }
                })
                .to_string(),
            ),
        );
        assert_eq!(r.status, 200, "oversampling={ok_ov} should preview cleanly");
    }

    // T2-17 & T2-18: Invalid CZTI & LAXPC overrides
    // Explicitly verify CZTI stage & energy validation errors against valid czti_l1_dummy
    for (req, expected_sub) in [
        (
            json!({
                "pipeline": "czti",
                "czti": { "input_dir": czti_l1_dummy.to_string_lossy(), "caldb_dir": czti_l1_dummy.to_string_lossy(), "output_dir": dir.path().join("o").to_string_lossy(), "entry_stage": 3, "exit_stage": 1 }
            }),
            "entry_stage",
        ),
        (
            json!({
                "pipeline": "czti",
                "czti": { "input_dir": czti_l1_dummy.to_string_lossy(), "caldb_dir": czti_l1_dummy.to_string_lossy(), "output_dir": dir.path().join("o").to_string_lossy(), "emin_kev": 100.0, "emax_kev": 20.0 }
            }),
            "emin_kev",
        ),
        (
            json!({
                "pipeline": "czti",
                "czti": { "input_dir": czti_l1_dummy.to_string_lossy(), "caldb_dir": czti_l1_dummy.to_string_lossy(), "output_dir": dir.path().join("o").to_string_lossy(), "emin_kev": 50.0, "emax_kev": 50.0 }
            }),
            "emin_kev",
        ),
        (
            json!({
                "pipeline": "czti",
                "czti": { "input_dir": czti_l1_dummy.to_string_lossy(), "caldb_dir": czti_l1_dummy.to_string_lossy(), "output_dir": dir.path().join("o").to_string_lossy(), "emin_kev": -5.0, "emax_kev": 20.0 }
            }),
            "emin_kev",
        ),
    ] {
        let res = http_request(addr, "POST", "/api/jobs/run", Some(&req.to_string()));
        assert_eq!(res.status, 400, "expected 400 for {req}, got {}: {}", res.status, res.body);
        assert!(
            res.body.contains(expected_sub),
            "expected '{expected_sub}' in 400 body, got: {}",
            res.body
        );
    }

    let bad_jobs = [
        // partial RA without Dec
        json!({
            "pipeline": "czti",
            "czti": { "input_dir": czti_l1_dummy.to_string_lossy(), "caldb_dir": czti_l1_dummy.to_string_lossy(), "output_dir": dir.path().join("o").to_string_lossy(), "ra_deg": 83.63 }
        }),
        // invalid quad
        json!({
            "pipeline": "czti",
            "czti": { "input_dir": czti_l1_dummy.to_string_lossy(), "caldb_dir": czti_l1_dummy.to_string_lossy(), "output_dir": dir.path().join("o").to_string_lossy(), "quads": [4] }
        }),
        // LAXPC invalid spectrum_channels, negative deadtime, reversed energy_range, pha_max > 1023
        json!({
            "pipeline": "laxpc",
            "laxpc": { "input_dir": czti_l1_dummy.to_string_lossy(), "output_dir": dir.path().join("o").to_string_lossy(), "spectrum_channels": 128 }
        }),
        json!({
            "pipeline": "laxpc",
            "laxpc": { "input_dir": czti_l1_dummy.to_string_lossy(), "output_dir": dir.path().join("o").to_string_lossy(), "deadtime_sec": -1e-6 }
        }),
        json!({
            "pipeline": "laxpc",
            "laxpc": { "input_dir": czti_l1_dummy.to_string_lossy(), "output_dir": dir.path().join("o").to_string_lossy(), "energy_range": "80:10" }
        }),
        json!({
            "pipeline": "laxpc",
            "laxpc": { "input_dir": czti_l1_dummy.to_string_lossy(), "output_dir": dir.path().join("o").to_string_lossy(), "pha_max": 1024 }
        }),
        // T2-19: bary orbit + mkf mutually exclusive
        json!({
            "pipeline": "bary",
            "bary": {
                "input": evt_path.to_string_lossy(),
                "output": dir.path().join("b1.evt").to_string_lossy(),
                "orbit": orb_path.to_string_lossy(),
                "mkf": mkf_path.to_string_lossy()
            }
        }),
        // T2-20a: dec_deg > 90
        json!({
            "pipeline": "bary",
            "bary": {
                "input": evt_path.to_string_lossy(),
                "output": dir.path().join("b2.evt").to_string_lossy(),
                "orbit": orb_path.to_string_lossy(),
                "ra_deg": 83.6,
                "dec_deg": 91.0
            }
        }),
    ];
    for req in bad_jobs {
        let res = http_request(addr, "POST", "/api/jobs/run", Some(&req.to_string()));
        assert_eq!(res.status, 400, "expected 400 for {req}, got {}: {}", res.status, res.body);
    }

    // T2-20b: double-barycentering already-TDB file fails cleanly
    let j_tdb = submit_and_wait(
        addr,
        &json!({
            "pipeline": "bary",
            "bary": {
                "input": tdb_evt.to_string_lossy(),
                "output": dir.path().join("double_bary.evt").to_string_lossy(),
                "orbit": orb_path.to_string_lossy(),
                "ra_deg": 83.6331,
                "dec_deg": 22.0145
            }
        }),
        15,
    );
    assert_eq!(j_tdb["status"], "failed");
    assert!(j_tdb["error"].as_str().is_some_and(|s| !s.is_empty()));

    // T2-21..23: Pulse search & multiband boundary rejections
    let bad_pulses = [
        // missing gti
        json!({ "path": evt_path.to_string_lossy(), "gti": dir.path().join("missing.gti").to_string_lossy(), "freq_hz": 29.658 }),
        // corrupt gti
        json!({ "path": evt_path.to_string_lossy(), "gti": corrupt_gti.to_string_lossy(), "freq_hz": 29.658 }),
        // gti leaving < 4 events
        json!({ "path": evt_path.to_string_lossy(), "gti": empty_gti.to_string_lossy(), "freq_hz": 29.658 }),
        // df_hz <= 0
        json!({ "path": evt_path.to_string_lossy(), "freq_min": 29.6, "freq_max": 29.7, "df_hz": 0.0 }),
        json!({ "path": evt_path.to_string_lossy(), "freq_min": 29.6, "freq_max": 29.7, "df_hz": -0.001 }),
        // harmonics / nbins / period_s / freq_min >= freq_max
        json!({ "path": evt_path.to_string_lossy(), "freq_hz": 29.658, "harmonics": 0 }),
        json!({ "path": evt_path.to_string_lossy(), "freq_hz": 29.658, "harmonics": 1025 }),
        json!({ "path": evt_path.to_string_lossy(), "freq_hz": 29.658, "nbins": 1 }),
        json!({ "path": evt_path.to_string_lossy(), "period_s": 0.0 }),
        json!({ "path": evt_path.to_string_lossy(), "freq_min": 30.0, "freq_max": 29.0 }),
    ];
    for req in bad_pulses {
        let r = http_request(addr, "POST", "/api/pulse/search", Some(&req.to_string()));
        assert!((400..500).contains(&r.status), "expected 4xx for pulse {req}, got {}: {}", r.status, r.body);
    }

    // Tiny df_hz either rejects cleanly with 400 (trials > MAX) or bounds n_trials <= 200_000
    // Use tiny 20-event FITS fixture so 200,000 trials completes in <1s in debug mode
    let r_huge = http_request(
        addr,
        "POST",
        "/api/pulse/search",
        Some(
            &json!({
                "path": tiny_evt.to_string_lossy(),
                "freq_min": 29.0,
                "freq_max": 30.0,
                "df_hz": 1e-12
            })
            .to_string(),
        ),
    );
    if r_huge.status == 200 {
        let vh: Value = serde_json::from_str(&r_huge.body).unwrap();
        assert!(vh["n_trials"].as_u64().unwrap() <= 200_000);
    } else {
        assert_eq!(r_huge.status, 400);
    }

    server.stop();
}

// ===========================================================================
// TIER 3 — CLI `--config`, `addspec`, `backshift` PARITY
// ===========================================================================

#[test]
#[ignore = "M3 E2E gate: run with -- --include-ignored"]
fn tier3_cli_config_addspec_backshift_and_bary_parity() {
    let bin = find_astrosat_binary();
    let dir = tempdir().unwrap();

    // Prepare fixtures
    let utc_evt = dir.path().join("utc.evt");
    let orb_path = dir.path().join("sat.orb");
    let ui_bary_out = dir.path().join("ui_bary.evt");
    let cli_bary_out = dir.path().join("cli_bary.evt");
    write_synthetic_utc_events_and_orbit(&utc_evt, &orb_path, "UTC");

    let p1 = dir.path().join("q0.pha");
    let p2 = dir.path().join("q1.pha");
    let r1 = dir.path().join("q0.rsp");
    let r2 = dir.path().join("q1.rsp");
    let ui_pha = dir.path().join("ui_comb.pha");
    let ui_rsp = dir.path().join("ui_comb.rsp");
    let cli_pha = dir.path().join("cli_comb.pha");
    let cli_rsp = dir.path().join("cli_comb.rsp");
    write_synthetic_czti_pha_512(&p1, &[4.0; 512], &[1.0; 512], 120.0, "q0.rsp", Some(0));
    write_synthetic_czti_pha_512(&p2, &[6.0; 512], &[1.5; 512], 80.0, "q1.rsp", Some(1));
    write_synthetic_czti_rsp_512(&r1, 10.0);
    write_synthetic_czti_rsp_512(&r2, 20.0);

    let caldb = dir.path().join("caldb_laxpc");
    write_synthetic_laxpc_caldb(&caldb, true);
    let bkg_pha = dir.path().join("backlxp10.pha");
    write_synthetic_laxpc_bkg_pha(&bkg_pha, 100.0);
    let ui_shift_pha = dir.path().join("ui_shift.pha");
    let cli_shift_pha = dir.path().join("cli_shift.pha");

    let server = spawn_test_server(dir.path());
    let addr = server.addr();

    // 1. Bary UI job -> export TOML -> astrosat --config <file.toml> bary ...
    let bary_req = json!({
        "pipeline": "bary",
        "threads": 2,
        "bary": {
            "input": utc_evt.to_string_lossy(),
            "output": ui_bary_out.to_string_lossy(),
            "orbit": orb_path.to_string_lossy(),
            "ra_deg": 83.6331,
            "dec_deg": 22.0145,
            "frame": "icrs",
            "ephemeris": "analytical",
            "write_checksums": true
        }
    });
    let prev = http_request(addr, "POST", "/api/cli-preview", Some(&bary_req.to_string()));
    assert_eq!(prev.status, 200);
    let pv: Value = serde_json::from_str(&prev.body).unwrap();
    let toml_text = pv["toml_config"].as_str().unwrap();

    // Parse TOML via /api/config/parse -> verify round-trip
    let parsed = http_request(
        addr,
        "POST",
        "/api/config/parse",
        Some(&json!({ "toml": toml_text }).to_string()),
    );
    assert_eq!(parsed.status, 200);

    let ui_bary_job = submit_and_wait(addr, &bary_req, 20);
    assert_eq!(ui_bary_job["status"], "succeeded", "UI bary failed: {ui_bary_job}");

    let cfg_path = dir.path().join("roundtrip.toml");
    std::fs::write(&cfg_path, toml_text).unwrap();
    let cli_bary_out_run = Command::new(&bin)
        .args([
            "--config",
            cfg_path.to_str().unwrap(),
            "bary",
            utc_evt.to_str().unwrap(),
            "-o",
            cli_bary_out.to_str().unwrap(),
        ])
        .output()
        .expect("run astrosat --config bary");
    assert!(
        cli_bary_out_run.status.success(),
        "CLI --config bary failed: {}",
        String::from_utf8_lossy(&cli_bary_out_run.stderr)
    );
    let ui_f = FitsFile::open(&ui_bary_out).unwrap();
    let cli_f = FitsFile::open(&cli_bary_out).unwrap();
    let ui_t = ui_f.hdu_by_name("EVENTS").unwrap().as_bintable().unwrap().read_col_f64("TIME").unwrap();
    let cli_t = cli_f.hdu_by_name("EVENTS").unwrap().as_bintable().unwrap().read_col_f64("TIME").unwrap();
    assert_eq!(ui_t.len(), cli_t.len());
    for (a, b) in ui_t.iter().zip(cli_t.iter()) {
        assert!((a - b).abs() < 1e-12, "bary timestamp mismatch: UI={a} CLI={b}");
    }

    // 2. Malformed TOML fails BOTH /api/config/parse (400) AND astrosat --config (non-zero exit)
    let bad_cfg = dir.path().join("bad.toml");
    std::fs::write(&bad_cfg, "[[czti]]\nthreads = 1\n").unwrap();
    let bad_api = http_request(
        addr,
        "POST",
        "/api/config/parse",
        Some(&json!({ "toml": "[[czti]]\nthreads = 1\n" }).to_string()),
    );
    assert_eq!(bad_api.status, 400);
    let bad_cli = Command::new(&bin)
        .args(["--config", bad_cfg.to_str().unwrap(), "inspect", utc_evt.to_str().unwrap()])
        .output()
        .expect("run astrosat --config bad.toml");
    assert!(!bad_cli.status.success(), "CLI must reject malformed TOML");

    // 3. czti addspec preview CLI execution parity vs UI job
    let add_req = json!({
        "pipeline": "czti-addspec",
        "czti_addspec": {
            "pha": [p1.to_string_lossy(), p2.to_string_lossy()],
            "rsp": [r1.to_string_lossy(), r2.to_string_lossy()],
            "out_pha": ui_pha.to_string_lossy(),
            "out_rsp": ui_rsp.to_string_lossy(),
            "add_exposure": true
        }
    });
    let ui_add = submit_and_wait(addr, &add_req, 20);
    assert_eq!(ui_add["status"], "succeeded", "UI addspec failed: {ui_add}");

    let cli_add_req = json!({
        "pipeline": "czti-addspec",
        "czti_addspec": {
            "pha": [p1.to_string_lossy(), p2.to_string_lossy()],
            "rsp": [r1.to_string_lossy(), r2.to_string_lossy()],
            "out_pha": cli_pha.to_string_lossy(),
            "out_rsp": cli_rsp.to_string_lossy(),
            "add_exposure": true
        }
    });
    let add_prev = http_request(addr, "POST", "/api/cli-preview", Some(&cli_add_req.to_string()));
    let add_cmd = serde_json::from_str::<Value>(&add_prev.body).unwrap()["cli_command"]
        .as_str()
        .unwrap()
        .to_string();
    let add_args = split_shell_words(&add_cmd);
    assert_eq!(add_args.first().map(String::as_str), Some("astrosat"));
    let add_run = Command::new(&bin)
        .args(&add_args[1..])
        .output()
        .expect("run previewed czti addspec command");
    assert!(
        add_run.status.success(),
        "CLI czti addspec failed: {}",
        String::from_utf8_lossy(&add_run.stderr)
    );
    let ui_exp = FitsFile::open(&ui_pha).unwrap().hdu_by_name("SPECTRUM").unwrap().header().get_f64("EXPOSURE").unwrap();
    let cli_exp = FitsFile::open(&cli_pha).unwrap().hdu_by_name("SPECTRUM").unwrap().header().get_f64("EXPOSURE").unwrap();
    assert!((ui_exp - cli_exp).abs() < 1e-9);

    // 4. laxpc backshift preview CLI --json execution parity vs UI job
    let bs_req = json!({
        "pipeline": "laxpc-backshift",
        "laxpc_backshift": {
            "unit": "10",
            "caldb_dir": caldb.to_string_lossy(),
            "mission_day": 500.0,
            "bkg_pha": bkg_pha.to_string_lossy(),
            "output": ui_shift_pha.to_string_lossy(),
            "pha_min": 5,
            "pha_max": 500
        }
    });
    let ui_bs = submit_and_wait(addr, &bs_req, 20);
    assert_eq!(ui_bs["status"], "succeeded", "UI backshift failed: {ui_bs}");

    let cli_bs_req = json!({
        "pipeline": "laxpc-backshift",
        "laxpc_backshift": {
            "unit": "10",
            "caldb_dir": caldb.to_string_lossy(),
            "mission_day": 500.0,
            "bkg_pha": bkg_pha.to_string_lossy(),
            "output": cli_shift_pha.to_string_lossy(),
            "pha_min": 5,
            "pha_max": 500
        }
    });
    let bs_prev = http_request(addr, "POST", "/api/cli-preview", Some(&cli_bs_req.to_string()));
    let bs_cmd = serde_json::from_str::<Value>(&bs_prev.body).unwrap()["cli_command"]
        .as_str()
        .unwrap()
        .to_string();
    let mut bs_args = split_shell_words(&bs_cmd);
    bs_args.push("--json".to_string());
    let bs_run = Command::new(&bin)
        .args(&bs_args[1..])
        .output()
        .expect("run previewed laxpc backshift command");
    assert!(
        bs_run.status.success(),
        "CLI laxpc backshift failed: {}",
        String::from_utf8_lossy(&bs_run.stderr)
    );
    let cli_bs_json: Value = serde_json::from_slice(&bs_run.stdout).expect("backshift --json stdout");
    let ui_dch = ui_bs["summary"]["delta_ch_30kev"].as_f64().unwrap();
    let cli_dch = cli_bs_json["gain_shift"]["delta_ch_30kev"].as_f64().unwrap();
    assert!((ui_dch - cli_dch).abs() < 1e-9, "delta_ch_30kev UI={ui_dch} CLI={cli_dch}");
    assert_eq!(
        ui_bs["summary"]["recommended_rmf"],
        cli_bs_json["gain_shift"]["recommended_rmf"]
    );

    server.stop();
}

// ===========================================================================
// TIER 4 — REAL-WORLD ASTROSAT WORKLOADS (set ASTROSAT_TEST_DATA)
// ===========================================================================

fn tier4_data_root() -> Option<std::path::PathBuf> {
    std::env::var("ASTROSAT_TEST_DATA")
        .ok()
        .map(std::path::PathBuf::from)
        .filter(|p| p.is_dir())
}

#[test]
#[ignore = "M3 E2E gate: run with ASTROSAT_TEST_DATA=... -- --include-ignored"]
fn tier4_czti_2018_four_quad_addspec_parity() {
    let Some(data_root) = tier4_data_root() else {
        eprintln!("Skipping tier4_czti_2018_four_quad_addspec_parity (ASTROSAT_TEST_DATA not set)");
        return;
    };
    let base = data_root.join("sample/czti_2018/20181014_C04_007T02_9000002432_level2/czti/modeM0");
    assert!(base.is_dir(), "CZTI 2018 dataset missing: {}", base.display());

    let phas: Vec<String> = (0..4)
        .map(|q| {
            base.join(format!(
                "AS1C04_007T02_9000002432cztM0_level2_quad_clean_Q{q}.pha"
            ))
            .to_string_lossy()
            .into_owned()
        })
        .collect();
    let rsps: Vec<String> = (0..4)
        .map(|q| {
            base.join(format!(
                "AS1C04_007T02_9000002432cztM0_level2_quad_clean_Q{q}.rsp"
            ))
            .to_string_lossy()
            .into_owned()
        })
        .collect();

    let quad_exps: Vec<f64> = phas
        .iter()
        .map(|p| {
            FitsFile::open(p)
                .unwrap()
                .hdu_by_name("SPECTRUM")
                .unwrap()
                .header()
                .get_f64("EXPOSURE")
                .unwrap()
        })
        .collect();
    let expected_mean_exp = quad_exps.iter().sum::<f64>() / 4.0;

    let tmp = tempdir().unwrap();
    let out_pha = tmp.path().join("czti_2018_4q.pha");
    let out_rsp = tmp.path().join("czti_2018_4q.rsp");

    let server = spawn_test_server(&data_root);
    let addr = server.addr();

    let job = submit_and_wait(
        addr,
        &json!({
            "pipeline": "czti-addspec",
            "czti_addspec": {
                "pha": phas,
                "rsp": rsps,
                "out_pha": out_pha.to_string_lossy(),
                "out_rsp": out_rsp.to_string_lossy(),
                "add_exposure": false
            }
        }),
        30,
    );
    assert_eq!(job["status"], "succeeded", "Tier4 addspec failed: {job}");
    assert!(out_pha.is_file());
    assert!(out_rsp.is_file());

    let comb_exp = FitsFile::open(&out_pha)
        .unwrap()
        .hdu_by_name("SPECTRUM")
        .unwrap()
        .header()
        .get_f64("EXPOSURE")
        .unwrap();
    assert!(
        (comb_exp - expected_mean_exp).abs() < 1e-3,
        "expected mean exposure {expected_mean_exp}, got {comb_exp}"
    );

    // Verify co-added spectrum counts equal sum of Q0..Q3 counts, and RSP matrix equals sum of Q0..Q3 RSPs
    let mut expected_sum_counts = vec![0.0_f64; 512];
    for p in &phas {
        let (_m, sp) = astrosat_czti::fits_io::read_spectrum_file(p).unwrap();
        for (dst, &src) in expected_sum_counts.iter_mut().zip(sp.counts.iter()) {
            *dst += src as f64;
        }
    }
    let (_cm, comb_sp) = astrosat_czti::fits_io::read_spectrum_file(&out_pha).unwrap();
    for (ch, (&act, &exp)) in comb_sp.counts.iter().zip(expected_sum_counts.iter()).enumerate() {
        assert!(
            ((act as f64) - exp).abs() < 1e-2,
            "Tier4 addspec channel {ch} mismatch: actual={act}, expected={exp}"
        );
    }

    let mut expected_sum_rsp: Option<Vec<f64>> = None;
    for r in &rsps {
        let (_m, rsp) = astrosat_czti::fits_io::read_response_file(r).unwrap();
        match &mut expected_sum_rsp {
            None => expected_sum_rsp = Some(rsp.matrix.iter().map(|&v| v as f64).collect()),
            Some(acc) => {
                for (dst, &src) in acc.iter_mut().zip(rsp.matrix.iter()) {
                    *dst += src as f64;
                }
            }
        }
    }
    let (_rm, comb_rsp) = astrosat_czti::fits_io::read_response_file(&out_rsp).unwrap();
    let exp_rsp_vec = expected_sum_rsp.unwrap();
    assert_eq!(comb_rsp.matrix.len(), exp_rsp_vec.len());
    let mut max_rsp_diff = 0.0_f64;
    for (&act, &exp) in comb_rsp.matrix.iter().zip(exp_rsp_vec.iter()) {
        max_rsp_diff = max_rsp_diff.max(((act as f64) - exp).abs());
    }
    assert!(
        max_rsp_diff < 1e-2,
        "Tier4 addspec RSP matrix sum mismatch: max_diff={max_rsp_diff}"
    );

    let val_r = http_request(
        addr,
        "GET",
        &format!(
            "/api/validate?path={}&strict=true",
            url_encode(&out_pha.to_string_lossy())
        ),
        None,
    );
    assert_eq!(val_r.status, 200);
    let val_v: Value = serde_json::from_str(&val_r.body).unwrap();
    assert_eq!(val_v["valid"], true, "strict validate failed: {val_v}");

    let sp_r = http_request(
        addr,
        "GET",
        &format!("/api/spectrum?path={}", url_encode(&out_pha.to_string_lossy())),
        None,
    );
    assert_eq!(sp_r.status, 200);
    let sp_v: Value = serde_json::from_str(&sp_r.body).unwrap();
    assert_eq!(sp_v["total_channels"], 512);

    server.stop();
}

#[test]
#[ignore = "M3 E2E gate: run with ASTROSAT_TEST_DATA=... -- --include-ignored"]
fn tier4_laxpc_backshift_crab_parity() {
    let Some(data_root) = tier4_data_root() else {
        eprintln!("Skipping tier4_laxpc_backshift_crab_parity (ASTROSAT_TEST_DATA not set)");
        return;
    };
    let caldb = data_root.join("ref_src/laxpc");
    let bkg_pha = caldb.join("backlxp20jan16.pha");
    let l2_fits = data_root.join("sample/laxpc_l2/20160107_G02_010T01_9000000252_level2/laxpc/01509/lxp2/modeEA/AS1G02_010T01_9000000252lxp2EA_level2.fits");
    assert!(bkg_pha.is_file());
    assert!(l2_fits.is_file());

    let tmp = tempdir().unwrap();
    let out_pha = tmp.path().join("crab_lxp2_shifted_bkg.pha");

    let server = spawn_test_server(&data_root);
    let addr = server.addr();

    let job = submit_and_wait(
        addr,
        &json!({
            "pipeline": "laxpc-backshift",
            "laxpc_backshift": {
                "unit": "20",
                "caldb_dir": caldb.to_string_lossy(),
                "fits": l2_fits.to_string_lossy(),
                "bkg_pha": bkg_pha.to_string_lossy(),
                "output": out_pha.to_string_lossy(),
                "layer1_only": true,
                "single_only": true,
                "pha_min": 0,
                "pha_max": 1023
            }
        }),
        30,
    );
    assert_eq!(job["status"], "succeeded", "Tier4 backshift failed: {job}");
    let s = &job["summary"];
    assert!(s["delta_ch_30kev"].as_f64().unwrap().is_finite());
    assert!(s["lcbk_shift_fraction"].as_f64().unwrap().is_finite());
    let rmf = s["recommended_rmf"].as_str().unwrap().to_string();
    assert!(rmf.starts_with("lx20"), "unexpected rmf for Lxp20: {rmf}");
    assert!(out_pha.is_file());

    let ins = http_request(
        addr,
        "GET",
        &format!("/api/inspect?path={}", url_encode(&out_pha.to_string_lossy())),
        None,
    );
    assert_eq!(ins.status, 200);
    let ins_v: Value = serde_json::from_str(&ins.body).unwrap();
    assert_eq!(ins_v["detnam"], "LX20");

    server.stop();
}

#[test]
#[ignore = "M3 E2E gate: run with ASTROSAT_TEST_DATA=... -- --include-ignored"]
fn tier4_czti_l1_to_l2_and_toml_cli_parity() {
    let Some(data_root) = tier4_data_root() else {
        eprintln!("Skipping tier4_czti_l1_to_l2_and_toml_cli_parity (ASTROSAT_TEST_DATA not set)");
        return;
    };
    let l1 = data_root.join("sample/czti/czti_sample_data/20160331_T01_112T01_9000000406_level1");
    let caldb = data_root.join("caldb/czti");
    assert!(l1.is_dir());
    assert!(caldb.is_dir());

    let tmp = tempdir().unwrap();
    let out_dir = tmp.path().join("czti_l2_out");
    std::fs::create_dir_all(&out_dir).unwrap();

    let server = spawn_test_server(&data_root);
    let addr = server.addr();

    let req = json!({
        "pipeline": "czti",
        "threads": 4,
        "czti": {
            "input_dir": l1.to_string_lossy(),
            "caldb_dir": caldb.to_string_lossy(),
            "output_dir": out_dir.to_string_lossy(),
            "output_prefix": "crab_e2e",
            "entry_stage": 1,
            "exit_stage": 3,
            "gti_type": "quad",
            "algo": "default",
            "image_algo": "balcc",
            "oversampling": 4,
            "mask_weight": true,
            "gen_resp": true
        }
    });
    let prev = http_request(addr, "POST", "/api/cli-preview", Some(&req.to_string()));
    assert_eq!(prev.status, 200);
    let toml_cfg = serde_json::from_str::<Value>(&prev.body).unwrap()["toml_config"]
        .as_str()
        .unwrap()
        .to_string();
    let parsed = http_request(
        addr,
        "POST",
        "/api/config/parse",
        Some(&json!({ "toml": toml_cfg }).to_string()),
    );
    assert_eq!(parsed.status, 200);

    let job = submit_and_wait(addr, &req, 180);
    assert_eq!(job["status"], "succeeded", "Tier4 CZTI L1->L2 failed: {job}");
    let outputs = job["outputs"].as_array().unwrap();
    assert!(
        outputs.len() > 12,
        "expected >12 CZTI L2 outputs, got {}",
        outputs.len()
    );
    let products: HashSet<&str> = outputs.iter().filter_map(|o| o["product"].as_str()).collect();
    assert!(products.contains("image"));
    assert!(products.contains("spectrum"));
    assert!(products.contains("lightcurve"));

    // Verify CLI `--config <toml> czti run` reproduces identical CZTI L2 reduction products
    let bin = find_astrosat_binary();
    let toml_file = tmp.path().join("czti_e2e.toml");
    std::fs::write(&toml_file, &toml_cfg).unwrap();
    let cli_out_dir = tmp.path().join("czti_l2_cli_out");
    std::fs::create_dir_all(&cli_out_dir).unwrap();
    let cli_czti = Command::new(&bin)
        .args([
            "--config",
            toml_file.to_str().unwrap(),
            "czti",
            "run",
            "--output",
            cli_out_dir.to_str().unwrap(),
            "--prefix",
            "crab_cli",
        ])
        .output()
        .expect("astrosat --config czti run");
    assert!(
        cli_czti.status.success(),
        "astrosat --config czti run failed: {}",
        String::from_utf8_lossy(&cli_czti.stderr)
    );
    let ui_q0_pha = out_dir.join("crab_e2e_quad_clean_Q0.pha");
    let cli_q0_pha = cli_out_dir.join("crab_cli_quad_clean_Q0.pha");
    assert!(ui_q0_pha.is_file(), "UI Q0.pha missing: {}", ui_q0_pha.display());
    assert!(cli_q0_pha.is_file(), "CLI Q0.pha missing: {}", cli_q0_pha.display());
    let (_um, ui_sp) = astrosat_czti::fits_io::read_spectrum_file(&ui_q0_pha).unwrap();
    let (_cm, cli_sp) = astrosat_czti::fits_io::read_spectrum_file(&cli_q0_pha).unwrap();
    assert!(
        (ui_sp.exposure - cli_sp.exposure).abs() < 1e-6,
        "CZTI TOML->CLI Q0 exposure mismatch: UI={} CLI={}",
        ui_sp.exposure,
        cli_sp.exposure
    );
    assert_eq!(ui_sp.counts.len(), cli_sp.counts.len());
    for (ch, (&a, &b)) in ui_sp.counts.iter().zip(cli_sp.counts.iter()).enumerate() {
        assert!(
            ((a - b) as f64).abs() < 1e-4,
            "CZTI TOML->CLI Q0 counts mismatch at ch {ch}: UI={a} CLI={b}"
        );
    }

    server.stop();
}

#[test]
#[ignore = "M3 E2E gate: run with ASTROSAT_TEST_DATA=... -- --include-ignored"]
fn tier4_crab_bary_jpleph405_and_pulse_gti_fddot() {
    let Some(data_root) = tier4_data_root() else {
        eprintln!("Skipping tier4_crab_bary_jpleph405_and_pulse_gti_fddot (ASTROSAT_TEST_DATA not set)");
        return;
    };
    let czti_evt = data_root.join("sample/czti/czti_sample_data/20160331_T01_112T01_9000000406_level2/modeM0/AS1T01_112T01_9000000406czt_level2_quad_clean.evt");
    let czti_gti = data_root.join("sample/czti/czti_sample_data/20160331_T01_112T01_9000000406_level2/modeM0/AS1T01_112T01_9000000406czt_level2.gti");
    let czti_orb = data_root.join("sample/czti_orb/AS1T01_112T01_9000000406czt_level1.orb");
    let jpl = data_root.join("ephem/JPLEPH.405");
    assert!(czti_evt.is_file());
    assert!(czti_gti.is_file());
    assert!(czti_orb.is_file());
    assert!(jpl.is_file());

    let tmp = tempdir().unwrap();
    let bary_evt = tmp.path().join("crab_czti_bary.evt");

    let server = spawn_test_server(&data_root);
    let addr = server.addr();

    let bary_job = submit_and_wait(
        addr,
        &json!({
            "pipeline": "bary",
            "threads": 4,
            "bary": {
                "input": czti_evt.to_string_lossy(),
                "output": bary_evt.to_string_lossy(),
                "orbit": czti_orb.to_string_lossy(),
                "ephemeris": jpl.to_string_lossy(),
                "ra_deg": 83.633083,
                "dec_deg": 22.014500,
                "frame": "icrs",
                "write_checksums": true
            }
        }),
        60,
    );
    assert_eq!(bary_job["status"], "succeeded", "Tier4 Crab bary failed: {bary_job}");
    assert_eq!(bary_job["summary"]["plephem"], "JPL-DE405");
    assert_eq!(
        FitsFile::open(&bary_evt).unwrap().hdu(0).unwrap().header().get_str("TIMESYS"),
        Some("TDB")
    );

    // Pulse search with gti + df_hz + fdot + fddot on Crab (~29.655 Hz in March 2016)
    let search_r = http_request(
        addr,
        "POST",
        "/api/pulse/search",
        Some(
            &json!({
                "path": bary_evt.to_string_lossy(),
                "gti": czti_gti.to_string_lossy(),
                "freq_min": 29.648,
                "freq_max": 29.664,
                "df_hz": 0.0002,
                "fdot": -3.68e-10,
                "fddot": 1.0e-20,
                "harmonics": 2,
                "nbins": 32
            })
            .to_string(),
        ),
    );
    assert_eq!(search_r.status, 200, "Tier4 Crab pulse search failed: {}", search_r.body);
    let ps: Value = serde_json::from_str(&search_r.body).unwrap();
    let best_f = ps["best_freq_hz"].as_f64().unwrap();
    let sigma = ps["detection_sigma"].as_f64().unwrap();
    assert!(
        (best_f - 29.655).abs() < 0.01,
        "expected Crab ~29.655 Hz, got {best_f}"
    );
    assert!(sigma > 5.0, "expected Crab detection_sigma > 5, got {sigma}");

    // Multiband fold across CZTI energy bands
    let mb_r = http_request(
        addr,
        "POST",
        "/api/pulse/multiband",
        Some(
            &json!({
                "inputs": [{
                    "label": "CZTI",
                    "path": bary_evt.to_string_lossy(),
                    "bands": [
                        { "label": "20-60 keV", "emin_kev": 20.0, "emax_kev": 60.0 },
                        { "label": "60-150 keV", "emin_kev": 60.0, "emax_kev": 150.0 }
                    ]
                }],
                "freq_hz": best_f,
                "fdot": -3.68e-10,
                "fddot": 1.0e-20,
                "nbins": 32,
                "harmonics": 2
            })
            .to_string(),
        ),
    );
    assert_eq!(mb_r.status, 200, "Tier4 multiband failed: {}", mb_r.body);
    let mb: Value = serde_json::from_str(&mb_r.body).unwrap();
    assert_eq!(mb["bands"].as_array().unwrap().len(), 2);

    server.stop();
}

#[test]
#[ignore = "M3 E2E gate: run with -- --include-ignored"]
fn tier5_phase2_remediation_regressions() {
    let dir = tempdir().unwrap();
    let caldb_dummy = dir.path().join("laxpc_caldb");
    std::fs::create_dir_all(&caldb_dummy).unwrap();
    let l1_dummy = dir.path().join("laxpc_l1");
    std::fs::create_dir_all(&l1_dummy).unwrap();

    let server = spawn_test_server(dir.path());
    let addr = server.addr();

    // 1. laxpc.energy_range hyphenated canonicalization & CLI parity
    let prev = http_request(
        addr,
        "POST",
        "/api/cli-preview",
        Some(
            &json!({
                "pipeline": "laxpc",
                "laxpc": {
                    "input_dir": l1_dummy.to_string_lossy(),
                    "caldb_dir": caldb_dummy.to_string_lossy(),
                    "output_dir": dir.path().join("out").to_string_lossy(),
                    "energy_range": "3-20"
                }
            })
            .to_string(),
        ),
    );
    assert_eq!(prev.status, 200, "hyphenated energy_range preview failed: {}", prev.body);
    let pv: Value = serde_json::from_str(&prev.body).unwrap();
    assert!(
        pv["cli_command"].as_str().unwrap().contains("--energy-range 3:20"),
        "expected canonical --energy-range 3:20, got {}",
        pv["cli_command"]
    );
    assert!(
        pv["toml_config"].as_str().unwrap().contains("energy_range = \"3:20\""),
        "expected canonical energy_range = \"3:20\", got {}",
        pv["toml_config"]
    );

    // Invalid energy_range forms must be rejected with 400 on both preview and run
    for bad_er in ["20-3", "-5:20", "10:10", "abc"] {
        let req = json!({
            "pipeline": "laxpc",
            "laxpc": {
                "input_dir": l1_dummy.to_string_lossy(),
                "caldb_dir": caldb_dummy.to_string_lossy(),
                "output_dir": dir.path().join("out").to_string_lossy(),
                "energy_range": bad_er
            }
        })
        .to_string();
        assert_eq!(
            http_request(addr, "POST", "/api/cli-preview", Some(&req)).status,
            400,
            "expected 400 preview for energy_range={bad_er}"
        );
        assert_eq!(
            http_request(addr, "POST", "/api/jobs/run", Some(&req)).status,
            400,
            "expected 400 run for energy_range={bad_er}"
        );
    }

    // Verify emitted TOML works with real `astrosat --config ... laxpc run` (no --energy-range format error)
    let bin = find_astrosat_binary();
    let laxpc_toml = dir.path().join("laxpc_range.toml");
    std::fs::write(&laxpc_toml, pv["toml_config"].as_str().unwrap()).unwrap();
    let cli_lxp = Command::new(&bin)
        .args(["--config", laxpc_toml.to_str().unwrap(), "laxpc", "run"])
        .output()
        .expect("run astrosat --config laxpc run");
    let lxp_err = String::from_utf8_lossy(&cli_lxp.stderr);
    assert!(
        !lxp_err.contains("must be formatted as EMIN:EMAX"),
        "emitted TOML energy_range failed CLI format check: {lxp_err}"
    );

    // 2 & 4. czti.output_type canonicalization & validation
    for (raw, canon) in [("spectrum", "spec"), ("lightcurve", "lc"), ("all", "both")] {
        let r = http_request(
            addr,
            "POST",
            "/api/cli-preview",
            Some(&json!({ "pipeline": "czti", "czti": { "output_type": raw } }).to_string()),
        );
        assert_eq!(r.status, 200, "output_type={raw} failed: {}", r.body);
        let v: Value = serde_json::from_str(&r.body).unwrap();
        assert!(
            v["cli_command"]
                .as_str()
                .unwrap()
                .contains(&format!("--output-type {canon}")),
            "expected --output-type {canon}, got {}",
            v["cli_command"]
        );
        assert!(
            v["toml_config"]
                .as_str()
                .unwrap()
                .contains(&format!("output_type = \"{canon}\"")),
            "expected output_type = \"{canon}\", got {}",
            v["toml_config"]
        );
    }
    assert_eq!(
        http_request(
            addr,
            "POST",
            "/api/cli-preview",
            Some(&json!({ "pipeline": "czti", "czti": { "output_type": "bogus" } }).to_string()),
        )
        .status,
        400
    );

    // 5. czti-addspec partial/misaligned rows rejected by preview
    let r_addspec_bad = http_request(
        addr,
        "POST",
        "/api/cli-preview",
        Some(
            &json!({
                "pipeline": "czti-addspec",
                "czti_addspec": {
                    "pha": ["/tmp/a.pha", ""],
                    "rsp": ["", "/tmp/b.rsp"]
                }
            })
            .to_string(),
        ),
    );
    assert_eq!(r_addspec_bad.status, 400);

    // 7. Invalid enum validation in /api/config/parse AND astrosat --config
    let dummy_fits = dir.path().join("dummy.fits");
    write_fits(&dummy_fits, &mut [Hdu::empty_primary()], true).unwrap();
    for bad_toml in [
        "[bary]\nframe = \"galactic\"\n",
        "[czti]\nalgo = \"bogus\"\n",
        "[czti]\ngti_type = \"invalid\"\n",
        "[czti]\noutput_type = \"bogus\"\n",
        "[czti]\nimage_algo = \"bogus\"\n",
    ] {
        let r = http_request(
            addr,
            "POST",
            "/api/config/parse",
            Some(&json!({ "toml": bad_toml }).to_string()),
        );
        assert_eq!(r.status, 400, "expected 400 for TOML:\n{bad_toml}");

        let p = dir.path().join("bad_enum.toml");
        std::fs::write(&p, bad_toml).unwrap();
        let cli_res = Command::new(&bin)
            .args(["--config", p.to_str().unwrap(), "inspect", dummy_fits.to_str().unwrap()])
            .output()
            .expect("astrosat --config inspect");
        assert!(
            !cli_res.status.success(),
            "expected CLI failure for bad TOML:\n{bad_toml}"
        );
    }

    // 9. GET /api/jobs/run returns 405 Method Not Allowed
    let r_405 = http_request(addr, "GET", "/api/jobs/run", None);
    assert_eq!(r_405.status, 405, "expected 405 for GET /api/jobs/run, got {}: {}", r_405.status, r_405.body);

    // 10. toml_lite BOM, multiline arrays with comments, and double comma rejection
    let good_bom_multiline = "\u{FEFF}threads = 2\n[czti]\nquads = [\n  # q0\n  0,\n  1, # q1\n]\n";
    let r_ok = http_request(
        addr,
        "POST",
        "/api/config/parse",
        Some(&json!({ "toml": good_bom_multiline }).to_string()),
    );
    assert_eq!(r_ok.status, 200, "BOM+multiline failed: {}", r_ok.body);
    let ok_v: Value = serde_json::from_str(&r_ok.body).unwrap();
    assert_eq!(ok_v["config"]["threads"], 2);
    assert_eq!(ok_v["config"]["czti"]["quads"], json!([0, 1]));

    let bad_comma = "[czti]\nquads = [0, , 1]\n";
    assert_eq!(
        http_request(
            addr,
            "POST",
            "/api/config/parse",
            Some(&json!({ "toml": bad_comma }).to_string()),
        )
        .status,
        400
    );

    server.stop();
}

#[test]
#[ignore = "M3 E2E gate: run with -- --include-ignored"]
fn tier5_iter2_validation_and_preview_parity() {
    let dir = tempdir().unwrap();
    let caldb_dummy = dir.path().join("laxpc_caldb");
    std::fs::create_dir_all(&caldb_dummy).unwrap();
    let l1_dummy = dir.path().join("laxpc_l1");
    std::fs::create_dir_all(&l1_dummy).unwrap();
    let dummy_fits = dir.path().join("dummy.fits");
    write_fits(&dummy_fits, &mut [Hdu::empty_primary()], true).unwrap();
    let pha0 = dir.path().join("q0.pha");
    let rsp0 = dir.path().join("q0.rsp");
    let rsp1 = dir.path().join("q1.rsp");
    write_synthetic_czti_pha_512(&pha0, &[1.0; 512], &[0.5; 512], 50.0, "q0.rsp", Some(0));
    write_synthetic_czti_rsp_512(&rsp0, 5.0);
    write_synthetic_czti_rsp_512(&rsp1, 5.0);

    let bin = find_astrosat_binary();
    let server = spawn_test_server(dir.path());
    let addr = server.addr();

    // Finding #1: AppConfig::validate() parity across POST /api/config/parse and astrosat --config
    let valid_toml = r#"
[czti]
time_binsize = 1.0
ra_deg = 83.6331
dec_deg = 22.0145

[laxpc]
energy_range = "3-20"
units = ["10", "lxp2"]
tbin = 0.5
deadtime_sec = 4.3e-5
anodes = "1,2,3"
merge_gap_sec = 0.0
min_duration_sec = 0.0
background_peak_ch = 213.0

[ui]
host = "127.0.0.1"
port = 8085
"#;
    let ok_api = http_request(
        addr,
        "POST",
        "/api/config/parse",
        Some(&json!({ "toml": valid_toml }).to_string()),
    );
    assert_eq!(ok_api.status, 200, "valid config rejected: {}", ok_api.body);
    let ok_path = dir.path().join("valid_full.toml");
    std::fs::write(&ok_path, valid_toml).unwrap();
    let ok_cli = Command::new(&bin)
        .args(["--config", ok_path.to_str().unwrap(), "inspect", dummy_fits.to_str().unwrap()])
        .output()
        .expect("run astrosat --config valid_full.toml");
    assert!(
        ok_cli.status.success(),
        "valid config failed CLI: {}",
        String::from_utf8_lossy(&ok_cli.stderr)
    );

    for (idx, bad_toml) in [
        "[laxpc]\nenergy_range = \"bogus\"\n",
        "[laxpc]\nenergy_range = \"80:10\"\n",
        "[laxpc]\nenergy_range = \"20-3\"\n",
        "[laxpc]\nunits = [\"99\"]\n",
        "[laxpc]\ntbin = -1.0\n",
        "[laxpc]\ntbin = 0.0\n",
        "[laxpc]\ndeadtime_sec = -0.1\n",
        "[laxpc]\nanodes = \"bogus\"\n",
        "[laxpc]\nmerge_gap_sec = -1.0\n",
        "[laxpc]\nmin_duration_sec = -0.5\n",
        "[laxpc]\nbackground_peak_ch = 0.0\n",
        "[czti]\ntime_binsize = -1.0\n",
        "[czti]\ndec_deg = 120.0\n",
        "[czti]\nra_deg = 400.0\n",
        "[ui]\nport = 0\n",
        "[ui]\nhost = \"\"\n",
    ]
    .iter()
    .enumerate()
    {
        let r = http_request(
            addr,
            "POST",
            "/api/config/parse",
            Some(&json!({ "toml": bad_toml }).to_string()),
        );
        assert_eq!(r.status, 400, "expected 400 for bad TOML #{idx}:\n{bad_toml}");

        let p = dir.path().join(format!("bad_cfg_{idx}.toml"));
        std::fs::write(&p, bad_toml).unwrap();
        let cli_res = Command::new(&bin)
            .args(["--config", p.to_str().unwrap(), "inspect", dummy_fits.to_str().unwrap()])
            .output()
            .expect("run astrosat --config bad_cfg");
        assert!(
            !cli_res.status.success(),
            "expected CLI exit != 0 for bad TOML #{idx}:\n{bad_toml}"
        );
    }

    // Finding #3: /api/cli-preview vs /api/jobs/run 400 parity
    let out_pha = dir.path().join("addspec_out.pha");
    let out_rsp = dir.path().join("addspec_out.rsp");
    let out_dir = dir.path().join("lxp_out");
    let bs_out = dir.path().join("bs_out.pha");
    let invalid_jobs = [
        // czti-addspec
        ("addspec_empty_rsp_entry", json!({
            "pipeline": "czti-addspec",
            "czti_addspec": { "pha": [pha0.to_string_lossy()], "rsp": [""], "out_pha": out_pha.to_string_lossy() }
        })),
        ("addspec_out_rsp_without_rsp", json!({
            "pipeline": "czti-addspec",
            "czti_addspec": { "pha": [pha0.to_string_lossy()], "rsp": [], "out_pha": out_pha.to_string_lossy(), "out_rsp": out_rsp.to_string_lossy() }
        })),
        ("addspec_rsp_len_mismatch", json!({
            "pipeline": "czti-addspec",
            "czti_addspec": { "pha": [pha0.to_string_lossy()], "rsp": [rsp0.to_string_lossy(), rsp1.to_string_lossy()], "out_pha": out_pha.to_string_lossy() }
        })),
        ("addspec_empty_pha_entry", json!({
            "pipeline": "czti-addspec",
            "czti_addspec": { "pha": [""], "out_pha": out_pha.to_string_lossy() }
        })),
        // laxpc-backshift
        ("backshift_invalid_unit", json!({
            "pipeline": "laxpc-backshift",
            "laxpc_backshift": { "unit": "99", "caldb_dir": caldb_dummy.to_string_lossy(), "mission_day": 100.0 }
        })),
        ("backshift_pha_min_gt_max", json!({
            "pipeline": "laxpc-backshift",
            "laxpc_backshift": { "unit": "10", "caldb_dir": caldb_dummy.to_string_lossy(), "mission_day": 100.0, "pha_min": 500, "pha_max": 100 }
        })),
        ("backshift_pha_max_gt_1023", json!({
            "pipeline": "laxpc-backshift",
            "laxpc_backshift": { "unit": "10", "caldb_dir": caldb_dummy.to_string_lossy(), "mission_day": 100.0, "pha_max": 1024 }
        })),
        ("backshift_neg_mission_day", json!({
            "pipeline": "laxpc-backshift",
            "laxpc_backshift": { "unit": "10", "caldb_dir": caldb_dummy.to_string_lossy(), "mission_day": -5.0 }
        })),
        ("backshift_zero_bkg_peak", json!({
            "pipeline": "laxpc-backshift",
            "laxpc_backshift": { "unit": "10", "caldb_dir": caldb_dummy.to_string_lossy(), "mission_day": 100.0, "bkg_peak_ch": 0.0 }
        })),
        ("backshift_output_without_bkg_pha", json!({
            "pipeline": "laxpc-backshift",
            "laxpc_backshift": { "unit": "10", "caldb_dir": caldb_dummy.to_string_lossy(), "mission_day": 100.0, "output": bs_out.to_string_lossy() }
        })),
        // laxpc
        ("laxpc_neg_tbin", json!({
            "pipeline": "laxpc",
            "laxpc": { "input_dir": l1_dummy.to_string_lossy(), "output_dir": out_dir.to_string_lossy(), "tbin": -1.0 }
        })),
        ("laxpc_zero_tbin", json!({
            "pipeline": "laxpc",
            "laxpc": { "input_dir": l1_dummy.to_string_lossy(), "output_dir": out_dir.to_string_lossy(), "tbin": 0.0 }
        })),
        ("laxpc_neg_deadtime", json!({
            "pipeline": "laxpc",
            "laxpc": { "input_dir": l1_dummy.to_string_lossy(), "output_dir": out_dir.to_string_lossy(), "deadtime_sec": -0.1 }
        })),
        ("laxpc_pha_min_gt_max", json!({
            "pipeline": "laxpc",
            "laxpc": { "input_dir": l1_dummy.to_string_lossy(), "output_dir": out_dir.to_string_lossy(), "pha_min": 500, "pha_max": 100 }
        })),
        ("laxpc_pha_max_gt_1023", json!({
            "pipeline": "laxpc",
            "laxpc": { "input_dir": l1_dummy.to_string_lossy(), "output_dir": out_dir.to_string_lossy(), "pha_max": 1024 }
        })),
        ("laxpc_invalid_channels", json!({
            "pipeline": "laxpc",
            "laxpc": { "input_dir": l1_dummy.to_string_lossy(), "output_dir": out_dir.to_string_lossy(), "spectrum_channels": 128 }
        })),
        ("laxpc_invalid_unit", json!({
            "pipeline": "laxpc",
            "laxpc": { "input_dir": l1_dummy.to_string_lossy(), "output_dir": out_dir.to_string_lossy(), "units": ["99"] }
        })),
        ("laxpc_invalid_anodes", json!({
            "pipeline": "laxpc",
            "laxpc": { "input_dir": l1_dummy.to_string_lossy(), "output_dir": out_dir.to_string_lossy(), "anodes": "bogus" }
        })),
        ("laxpc_neg_merge_gap", json!({
            "pipeline": "laxpc",
            "laxpc": { "input_dir": l1_dummy.to_string_lossy(), "output_dir": out_dir.to_string_lossy(), "merge_gap_sec": -1.0 }
        })),
        ("laxpc_neg_min_duration", json!({
            "pipeline": "laxpc",
            "laxpc": { "input_dir": l1_dummy.to_string_lossy(), "output_dir": out_dir.to_string_lossy(), "min_duration_sec": -1.0 }
        })),
    ];

    for (label, req) in invalid_jobs {
        let payload = req.to_string();
        let prev = http_request(addr, "POST", "/api/cli-preview", Some(&payload));
        assert_eq!(
            prev.status, 400,
            "{label}: /api/cli-preview expected 400, got {}: {}",
            prev.status, prev.body
        );
        let run = http_request(addr, "POST", "/api/jobs/run", Some(&payload));
        assert_eq!(
            run.status, 400,
            "{label}: /api/jobs/run expected 400, got {}: {}",
            run.status, run.body
        );
    }

    server.stop();
}

#[test]
#[ignore = "M3 E2E gate: run with -- --include-ignored"]
fn tier5_iter3_validation_and_empty_selection_parity() {
    let dir = tempdir().unwrap();
    let czti_l1_dummy = dir.path().join("czti_l1");
    let czti_caldb_dummy = dir.path().join("czti_caldb");
    let laxpc_l1_dummy = dir.path().join("laxpc_l1");
    std::fs::create_dir_all(&czti_l1_dummy).unwrap();
    std::fs::create_dir_all(&czti_caldb_dummy).unwrap();
    std::fs::create_dir_all(&laxpc_l1_dummy).unwrap();

    let dummy_fits = dir.path().join("dummy.fits");
    write_fits(&dummy_fits, &mut [Hdu::empty_primary()], true).unwrap();
    let m0_dummy = dir.path().join("m0.fits");
    let tct_dummy = dir.path().join("tct.fits");
    let mkf_dummy = dir.path().join("mkf.fits");
    write_fits(&m0_dummy, &mut [Hdu::empty_primary()], true).unwrap();
    write_fits(&tct_dummy, &mut [Hdu::empty_primary()], true).unwrap();
    write_fits(&mkf_dummy, &mut [Hdu::empty_primary()], true).unwrap();

    let bary_in = dir.path().join("bary_in.evt");
    let bary_orb = dir.path().join("bary.orb");
    write_synthetic_utc_events_and_orbit(&bary_in, &bary_orb, "UTC");

    let czti_out = dir.path().join("czti_out");
    let laxpc_out = dir.path().join("laxpc_out");
    let bary_out = dir.path().join("bary_out.evt");

    let bin = find_astrosat_binary();
    let server = spawn_test_server(dir.path());
    let addr = server.addr();

    // Finding #1 & #2: Preview & Job Run 400 parity
    let invalid_requests = [
        ("threads_zero", json!({
            "pipeline": "czti",
            "threads": 0,
            "czti": {
                "m0": m0_dummy.to_string_lossy(),
                "tct": tct_dummy.to_string_lossy(),
                "mkf": mkf_dummy.to_string_lossy(),
                "caldb_dir": czti_caldb_dummy.to_string_lossy(),
                "output_dir": czti_out.to_string_lossy()
            }
        })),
        ("threads_over_512", json!({
            "pipeline": "laxpc",
            "threads": 513,
            "laxpc": {
                "input_dir": laxpc_l1_dummy.to_string_lossy(),
                "output_dir": laxpc_out.to_string_lossy()
            }
        })),
        ("bary_dec_gt_90", json!({
            "pipeline": "bary",
            "bary": {
                "input": bary_in.to_string_lossy(),
                "output": bary_out.to_string_lossy(),
                "orbit": bary_orb.to_string_lossy(),
                "ra_deg": 83.6331,
                "dec_deg": 120.0
            }
        })),
        ("bary_dec_lt_neg90", json!({
            "pipeline": "bary",
            "bary": {
                "input": bary_in.to_string_lossy(),
                "output": bary_out.to_string_lossy(),
                "orbit": bary_orb.to_string_lossy(),
                "ra_deg": 83.6331,
                "dec_deg": -91.0
            }
        })),
        ("bary_ra_gt_360", json!({
            "pipeline": "bary",
            "bary": {
                "input": bary_in.to_string_lossy(),
                "output": bary_out.to_string_lossy(),
                "orbit": bary_orb.to_string_lossy(),
                "ra_deg": 400.0,
                "dec_deg": 22.0145
            }
        })),
        ("bary_ra_lt_neg360", json!({
            "pipeline": "bary",
            "bary": {
                "input": bary_in.to_string_lossy(),
                "output": bary_out.to_string_lossy(),
                "orbit": bary_orb.to_string_lossy(),
                "ra_deg": -400.0,
                "dec_deg": 22.0145
            }
        })),
        ("bary_frame_galactic", json!({
            "pipeline": "bary",
            "bary": {
                "input": bary_in.to_string_lossy(),
                "output": bary_out.to_string_lossy(),
                "orbit": bary_orb.to_string_lossy(),
                "ra_deg": 83.6331,
                "dec_deg": 22.0145,
                "frame": "galactic"
            }
        })),
        ("czti_dec_gt_90", json!({
            "pipeline": "czti",
            "czti": {
                "m0": m0_dummy.to_string_lossy(),
                "tct": tct_dummy.to_string_lossy(),
                "mkf": mkf_dummy.to_string_lossy(),
                "caldb_dir": czti_caldb_dummy.to_string_lossy(),
                "output_dir": czti_out.to_string_lossy(),
                "ra_deg": 83.6331,
                "dec_deg": 120.0
            }
        })),
        ("czti_dec_lt_neg90", json!({
            "pipeline": "czti",
            "czti": {
                "m0": m0_dummy.to_string_lossy(),
                "tct": tct_dummy.to_string_lossy(),
                "mkf": mkf_dummy.to_string_lossy(),
                "caldb_dir": czti_caldb_dummy.to_string_lossy(),
                "output_dir": czti_out.to_string_lossy(),
                "ra_deg": 83.6331,
                "dec_deg": -91.0
            }
        })),
        ("czti_ra_gt_360", json!({
            "pipeline": "czti",
            "czti": {
                "m0": m0_dummy.to_string_lossy(),
                "tct": tct_dummy.to_string_lossy(),
                "mkf": mkf_dummy.to_string_lossy(),
                "caldb_dir": czti_caldb_dummy.to_string_lossy(),
                "output_dir": czti_out.to_string_lossy(),
                "ra_deg": 400.0,
                "dec_deg": 22.0145
            }
        })),
        ("czti_ra_lt_neg360", json!({
            "pipeline": "czti",
            "czti": {
                "m0": m0_dummy.to_string_lossy(),
                "tct": tct_dummy.to_string_lossy(),
                "mkf": mkf_dummy.to_string_lossy(),
                "caldb_dir": czti_caldb_dummy.to_string_lossy(),
                "output_dir": czti_out.to_string_lossy(),
                "ra_deg": -400.0,
                "dec_deg": 22.0145
            }
        })),
        ("czti_empty_quads", json!({
            "pipeline": "czti",
            "czti": {
                "m0": m0_dummy.to_string_lossy(),
                "tct": tct_dummy.to_string_lossy(),
                "mkf": mkf_dummy.to_string_lossy(),
                "caldb_dir": czti_caldb_dummy.to_string_lossy(),
                "output_dir": czti_out.to_string_lossy(),
                "quads": []
            }
        })),
        ("laxpc_empty_units", json!({
            "pipeline": "laxpc",
            "laxpc": {
                "input_dir": laxpc_l1_dummy.to_string_lossy(),
                "output_dir": laxpc_out.to_string_lossy(),
                "units": []
            }
        })),
    ];

    for (label, req) in invalid_requests {
        let payload = req.to_string();
        let prev = http_request(addr, "POST", "/api/cli-preview", Some(&payload));
        assert_eq!(
            prev.status, 400,
            "{label}: /api/cli-preview expected 400, got {}: {}",
            prev.status, prev.body
        );
        let run = http_request(addr, "POST", "/api/jobs/run", Some(&payload));
        assert_eq!(
            run.status, 400,
            "{label}: /api/jobs/run expected 400, got {}: {}",
            run.status, run.body
        );
    }

    // Finding #2: /api/config/parse and astrosat --config parity for empty quads/units
    for (idx, bad_toml) in [
        "[czti]\nquads = []\n",
        "[laxpc]\nunits = []\n",
        "[bary]\ndec_deg = 120.0\n",
        "[bary]\nra_deg = 400.0\n",
        "[bary]\nframe = \"galactic\"\n",
        "threads = 0\n",
    ]
    .iter()
    .enumerate()
    {
        let r = http_request(
            addr,
            "POST",
            "/api/config/parse",
            Some(&json!({ "toml": bad_toml }).to_string()),
        );
        assert_eq!(r.status, 400, "expected 400 for TOML #{idx}:\n{bad_toml}");

        let p = dir.path().join(format!("iter3_bad_{idx}.toml"));
        std::fs::write(&p, bad_toml).unwrap();
        let cli_res = Command::new(&bin)
            .args(["--config", p.to_str().unwrap(), "inspect", dummy_fits.to_str().unwrap()])
            .output()
            .expect("run astrosat --config");
        assert!(
            !cli_res.status.success(),
            "expected CLI exit != 0 for TOML #{idx}:\n{bad_toml}"
        );
    }

    // Omitted quads/units defaults still succeed across TOML, CLI, /api/cli-preview, and /api/jobs/run
    let omitted_toml = "[czti]\nra_deg = 83.6331\ndec_deg = 22.0145\n[laxpc]\ntbin = 1.0\n";
    let omitted_parse = http_request(
        addr,
        "POST",
        "/api/config/parse",
        Some(&json!({ "toml": omitted_toml }).to_string()),
    );
    assert_eq!(omitted_parse.status, 200, "omitted defaults rejected: {}", omitted_parse.body);
    let omitted_path = dir.path().join("iter3_omitted.toml");
    std::fs::write(&omitted_path, omitted_toml).unwrap();
    let omitted_cli = Command::new(&bin)
        .args(["--config", omitted_path.to_str().unwrap(), "inspect", dummy_fits.to_str().unwrap()])
        .output()
        .expect("run astrosat --config omitted defaults");
    assert!(
        omitted_cli.status.success(),
        "omitted defaults CLI failed: {}",
        String::from_utf8_lossy(&omitted_cli.stderr)
    );

    let czti_omitted_req = json!({
        "pipeline": "czti",
        "czti": {
            "m0": m0_dummy.to_string_lossy(),
            "tct": tct_dummy.to_string_lossy(),
            "mkf": mkf_dummy.to_string_lossy(),
            "caldb_dir": czti_caldb_dummy.to_string_lossy(),
            "output_dir": czti_out.to_string_lossy()
        }
    });
    let czti_prev = http_request(addr, "POST", "/api/cli-preview", Some(&czti_omitted_req.to_string()));
    assert_eq!(czti_prev.status, 200, "czti omitted preview failed: {}", czti_prev.body);
    assert!(!czti_prev.body.contains("--quads"));

    let laxpc_omitted_req = json!({
        "pipeline": "laxpc",
        "laxpc": {
            "input_dir": laxpc_l1_dummy.to_string_lossy(),
            "output_dir": laxpc_out.to_string_lossy()
        }
    });
    let laxpc_prev = http_request(addr, "POST", "/api/cli-preview", Some(&laxpc_omitted_req.to_string()));
    assert_eq!(laxpc_prev.status, 200, "laxpc omitted preview failed: {}", laxpc_prev.body);
    assert!(!laxpc_prev.body.contains("--units"));

    let bary_ok_req = json!({
        "pipeline": "bary",
        "bary": {
            "input": bary_in.to_string_lossy(),
            "output": bary_out.to_string_lossy(),
            "orbit": bary_orb.to_string_lossy(),
            "ra_deg": 83.6331,
            "dec_deg": 22.0145
        }
    });
    let bary_prev = http_request(addr, "POST", "/api/cli-preview", Some(&bary_ok_req.to_string()));
    assert_eq!(bary_prev.status, 200, "bary valid preview failed: {}", bary_prev.body);
    let bary_done = submit_and_wait(addr, &bary_ok_req, 30);
    assert_eq!(bary_done["status"], "succeeded", "bary job failed: {bary_done}");

    server.stop();
}
