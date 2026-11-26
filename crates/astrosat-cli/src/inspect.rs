//! `astrosat inspect <fits>` — fast human/JSON summary of HDUs, rows, timing keywords,
//! exposure, and telemetry health.

use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use astrosat_core::time::astrosat_utc_to_iso8601;
use astrosat_fits::{FitsFile, HduKind};
use serde::Serialize;

#[derive(Debug, Clone, Serialize)]
pub struct ColumnSummary {
    pub name: String,
    pub tform: String,
    pub unit: Option<String>,
}

#[derive(Debug, Clone, Serialize)]
pub struct HduSummary {
    pub index: usize,
    pub extname: String,
    pub extver: Option<i64>,
    pub kind: String,
    pub hduclas1: Option<String>,
    pub nrows: Option<usize>,
    pub shape: Option<Vec<usize>>,
    pub columns: Vec<ColumnSummary>,
    pub tstart: Option<f64>,
    pub tstop: Option<f64>,
    pub exposure: Option<f64>,
    pub datasum_ok: bool,
    pub checksum_ok: bool,
    pub has_checksum: bool,
}

#[derive(Debug, Clone, Serialize)]
pub struct TelemetryHealth {
    pub all_checksums_ok: bool,
    pub checksummed_hdus: usize,
    pub total_events: usize,
    pub gti_intervals: usize,
    pub gti_exposure_sec: f64,
    pub non_monotonic_time_hdus: Vec<String>,
    pub non_finite_time_rows: usize,
    pub invalid_gti_rows: usize,
}

#[derive(Debug, Clone, Serialize)]
pub struct FileInspection {
    pub path: PathBuf,
    pub file_size_bytes: u64,
    pub num_hdus: usize,
    pub telescop: Option<String>,
    pub instrume: Option<String>,
    pub detnam: Option<String>,
    pub obs_id: Option<String>,
    pub object: Option<String>,
    pub ra_obj: Option<f64>,
    pub dec_obj: Option<f64>,
    pub timesys: Option<String>,
    pub timeref: Option<String>,
    pub timeunit: Option<String>,
    pub mjdrefi: Option<i64>,
    pub mjdreff: Option<f64>,
    pub tstart: Option<f64>,
    pub tstop: Option<f64>,
    pub date_obs: Option<String>,
    pub date_end: Option<String>,
    pub exposure: Option<f64>,
    pub ontime: Option<f64>,
    pub livetime: Option<f64>,
    pub hdus: Vec<HduSummary>,
    pub health: TelemetryHealth,
}

#[derive(Debug, Clone, Serialize)]
pub struct InspectReport {
    pub files: Vec<FileInspection>,
}

pub fn run_inspect(paths: &[PathBuf], full_scan: bool) -> Result<InspectReport> {
    let mut files = Vec::with_capacity(paths.len());
    for p in paths {
        files.push(inspect_single_file(p, full_scan)?);
    }
    Ok(InspectReport { files })
}

fn inspect_single_file(path: &Path, full_scan: bool) -> Result<FileInspection> {
    let meta = std::fs::metadata(path)
        .with_context(|| format!("failed to stat FITS file {}", path.display()))?;
    let fits = FitsFile::open(path)
        .with_context(|| format!("failed to open FITS file {}", path.display()))?;

    let get_str_any = |key: &str| -> Option<String> {
        for idx in [0usize, 1usize] {
            if let Ok(h) = fits.hdu(idx) {
                if let Some(v) = h.header().get_str(key) {
                    let t = v.trim();
                    if !t.is_empty() {
                        return Some(t.to_string());
                    }
                }
            }
        }
        None
    };

    let get_f64_any = |key: &str| -> Option<f64> {
        for idx in [0usize, 1usize] {
            if let Ok(h) = fits.hdu(idx) {
                if let Some(v) = h.header().get_f64(key) {
                    return Some(v);
                }
            }
        }
        None
    };

    let get_i64_any = |key: &str| -> Option<i64> {
        for idx in [0usize, 1usize] {
            if let Ok(h) = fits.hdu(idx) {
                if let Some(v) = h.header().get_i64(key) {
                    return Some(v);
                }
            }
        }
        None
    };

    let tstart = get_f64_any("TSTART");
    let tstop = get_f64_any("TSTOP");
    let date_obs = get_str_any("DATE-OBS").or_else(|| tstart.map(astrosat_utc_to_iso8601));
    let date_end = get_str_any("DATE-END").or_else(|| tstop.map(astrosat_utc_to_iso8601));

    let mut hdus = Vec::with_capacity(fits.len());
    let mut all_checksums_ok = true;
    let mut checksummed_hdus = 0usize;
    let mut total_events = 0usize;
    let mut gti_intervals = 0usize;
    let mut gti_exposure_sec = 0.0f64;
    let mut non_monotonic_time_hdus = Vec::new();
    let mut non_finite_time_rows = 0usize;
    let mut invalid_gti_rows = 0usize;

    for h in fits.iter_hdus() {
        let hdr = h.header();
        let has_checksum = hdr.get_str("CHECKSUM").is_some() || hdr.get_str("DATASUM").is_some();
        let (datasum_ok, checksum_ok) = h.verify_checksum();
        if has_checksum {
            checksummed_hdus += 1;
            if !datasum_ok || !checksum_ok {
                all_checksums_ok = false;
            }
        }

        let extname = h
            .extname()
            .map(|s| s.trim().to_string())
            .filter(|s| !s.is_empty())
            .unwrap_or_else(|| {
                if h.index() == 0 {
                    "PRIMARY".to_string()
                } else {
                    format!("HDU_{}", h.index())
                }
            });

        let kind_str = match h.kind() {
            HduKind::PrimaryImage => "PrimaryImage".to_string(),
            HduKind::Image => "Image".to_string(),
            HduKind::BinTable => "BinTable".to_string(),
            HduKind::AsciiTable => "AsciiTable".to_string(),
            HduKind::Other(s) => s.clone(),
        };

        let mut nrows = None;
        let mut shape = None;
        let mut columns = Vec::new();

        if *h.kind() == HduKind::BinTable {
            if let Ok(bt) = h.as_bintable() {
                nrows = Some(bt.nrows());
                columns = bt
                    .columns()
                    .iter()
                    .map(|c| ColumnSummary {
                        name: c.name.clone(),
                        tform: c.tform.clone(),
                        unit: c.unit.clone(),
                    })
                    .collect();

                let upper_name = extname.to_ascii_uppercase();
                let is_gti = upper_name.contains("GTI")
                    || hdr
                        .get_str("HDUCLAS1")
                        .map(|s| s.trim().eq_ignore_ascii_case("GTI"))
                        .unwrap_or(false);

                if is_gti {
                    if let (Ok(starts), Ok(stops)) =
                        (bt.read_col_f64("START"), bt.read_col_f64("STOP"))
                    {
                        gti_intervals += starts.len();
                        let mut exp = 0.0;
                        for (&a, &b) in starts.iter().zip(stops.iter()) {
                            if !a.is_finite() || !b.is_finite() || b < a {
                                invalid_gti_rows += 1;
                            } else {
                                exp += b - a;
                            }
                        }
                        if gti_exposure_sec == 0.0 || upper_name == "GTI" {
                            gti_exposure_sec = exp;
                        }
                    }
                } else if bt.column_index("TIME").is_some()
                    && (upper_name.contains("EVENT")
                        || upper_name.contains("MODE_EA")
                        || upper_name.contains("EVT")
                        || (matches!(upper_name.as_str(), "Q0" | "Q1" | "Q2" | "Q3")
                            && (bt.column_index("PHA").is_some()
                                || bt.column_index("PI").is_some()
                                || bt.column_index("ENERGY").is_some())))
                {
                    total_events += bt.nrows();
                }

                if full_scan {
                    if let Ok(times) = bt.read_col_f64("TIME") {
                        let mut prev = f64::NEG_INFINITY;
                        let mut mono = true;
                        for &t in &times {
                            if !t.is_finite() {
                                non_finite_time_rows += 1;
                            } else if t < prev {
                                mono = false;
                            } else {
                                prev = t;
                            }
                        }
                        if !mono {
                            non_monotonic_time_hdus.push(extname.clone());
                        }
                    }
                }
            }
        } else if matches!(h.kind(), HduKind::PrimaryImage | HduKind::Image) {
            if let Ok(img) = h.as_image() {
                shape = Some(img.shape().to_vec());
            }
        }

        hdus.push(HduSummary {
            index: h.index(),
            extname,
            extver: h.extver(),
            kind: kind_str,
            hduclas1: hdr.get_str("HDUCLAS1").map(|s| s.trim().to_string()),
            nrows,
            shape,
            columns,
            tstart: hdr.get_f64("TSTART"),
            tstop: hdr.get_f64("TSTOP"),
            exposure: hdr.get_f64("EXPOSURE"),
            datasum_ok,
            checksum_ok,
            has_checksum,
        });
    }

    Ok(FileInspection {
        path: path.to_path_buf(),
        file_size_bytes: meta.len(),
        num_hdus: fits.len(),
        telescop: get_str_any("TELESCOP"),
        instrume: get_str_any("INSTRUME"),
        detnam: get_str_any("DETNAM"),
        obs_id: get_str_any("OBS_ID"),
        object: get_str_any("OBJECT"),
        ra_obj: get_f64_any("RA_OBJ").or_else(|| get_f64_any("RA_NOM")).or_else(|| get_f64_any("RA_PNT")),
        dec_obj: get_f64_any("DEC_OBJ").or_else(|| get_f64_any("DEC_NOM")).or_else(|| get_f64_any("DEC_PNT")),
        timesys: get_str_any("TIMESYS"),
        timeref: get_str_any("TIMEREF"),
        timeunit: get_str_any("TIMEUNIT"),
        mjdrefi: get_i64_any("MJDREFI"),
        mjdreff: get_f64_any("MJDREFF"),
        tstart,
        tstop,
        date_obs,
        date_end,
        exposure: get_f64_any("EXPOSURE"),
        ontime: get_f64_any("ONTIME"),
        livetime: get_f64_any("LIVETIME"),
        hdus,
        health: TelemetryHealth {
            all_checksums_ok,
            checksummed_hdus,
            total_events,
            gti_intervals,
            gti_exposure_sec,
            non_monotonic_time_hdus,
            non_finite_time_rows,
            invalid_gti_rows,
        },
    })
}

pub fn print_inspect_human(report: &InspectReport) {
    for f in &report.files {
        println!("File: {} ({} bytes, {} HDUs)", f.path.display(), f.file_size_bytes, f.num_hdus);
        println!(
            "  Mission/Instrument : {} / {} (DETNAM: {})",
            f.telescop.as_deref().unwrap_or("-"),
            f.instrume.as_deref().unwrap_or("-"),
            f.detnam.as_deref().unwrap_or("-")
        );
        println!(
            "  Obs/Target         : OBS_ID={}  OBJECT={}  (RA={}, DEC={})",
            f.obs_id.as_deref().unwrap_or("-"),
            f.object.as_deref().unwrap_or("-"),
            f.ra_obj.map(|v| format!("{v:.4}")).unwrap_or_else(|| "-".into()),
            f.dec_obj.map(|v| format!("{v:.4}")).unwrap_or_else(|| "-".into()),
        );
        println!(
            "  Timing             : TIMESYS={}  TIMEREF={}  MJDREF={}+{}",
            f.timesys.as_deref().unwrap_or("-"),
            f.timeref.as_deref().unwrap_or("-"),
            f.mjdrefi.map(|v| v.to_string()).unwrap_or_else(|| "-".into()),
            f.mjdreff.map(|v| format!("{v:.6}")).unwrap_or_else(|| "-".into())
        );
        println!(
            "  Span / Exposure    : TSTART={}  TSTOP={}  EXPOSURE={}",
            f.tstart.map(|v| format!("{v:.3}")).unwrap_or_else(|| "-".into()),
            f.tstop.map(|v| format!("{v:.3}")).unwrap_or_else(|| "-".into()),
            f.exposure.map(|v| format!("{v:.3} s")).unwrap_or_else(|| "-".into()),
        );
        if let (Some(d0), Some(d1)) = (&f.date_obs, &f.date_end) {
            println!("  UTC/ISO Window     : {d0} .. {d1}");
        }
        println!(
            "  Health             : checksums={} ({}/{} HDUs), events={}, GTI={} intervals ({:.3} s)",
            if f.health.all_checksums_ok { "OK" } else { "FAIL" },
            f.health.checksummed_hdus,
            f.num_hdus,
            f.health.total_events,
            f.health.gti_intervals,
            f.health.gti_exposure_sec,
        );
        if !f.health.non_monotonic_time_hdus.is_empty()
            || f.health.non_finite_time_rows > 0
            || f.health.invalid_gti_rows > 0
        {
            println!(
                "  Warnings           : non_monotonic_time={:?}, non_finite_times={}, invalid_gtis={}",
                f.health.non_monotonic_time_hdus,
                f.health.non_finite_time_rows,
                f.health.invalid_gti_rows
            );
        }
        println!("  HDUs:");
        for h in &f.hdus {
            let dim = if let Some(nr) = h.nrows {
                format!("{nr} rows x {} cols", h.columns.len())
            } else if let Some(sh) = &h.shape {
                if sh.is_empty() {
                    "empty".to_string()
                } else {
                    format!(
                        "shape [{}]",
                        sh.iter()
                            .map(|d| d.to_string())
                            .collect::<Vec<_>>()
                            .join("x")
                    )
                }
            } else {
                "-".to_string()
            };
            let csum = if h.has_checksum {
                if h.datasum_ok && h.checksum_ok {
                    "CSUM=OK"
                } else {
                    "CSUM=FAIL"
                }
            } else {
                "CSUM=NONE"
            };
            println!(
                "    [{}] {:<16} {:<12} {:<22} {}",
                h.index, h.extname, h.kind, dim, csum
            );
        }
    }
}
