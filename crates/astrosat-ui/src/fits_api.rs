//! FITS Inspector (`GET /api/inspect`) and OGIP/FITS 4.0 Validator (`GET /api/validate`).

use std::collections::HashMap;
use std::path::Path;

use astrosat_core::time::astrosat_utc_to_iso8601;
use astrosat_fits::{
    validate_ogip, BinTableRef, CardValue, ColumnSpec, ColumnType, FitsFile, HduKind,
    IssueSeverity, VarDescriptor, FITS_BLOCK_SIZE,
};
use serde::Serialize;
use serde_json::Value;

use crate::http_util::{
    query_bool, query_required, query_usize, resolve_existing_file, ApiError, ApiResult,
};
use crate::science::open_fits_clean;

// ---------------------------------------------------------------------------
// FITS Inspect
// ---------------------------------------------------------------------------

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
pub struct GtiIntervalDetail {
    pub start: f64,
    pub stop: f64,
    pub duration: f64,
}

#[derive(Debug, Clone, Serialize)]
pub struct HeaderCardJson {
    pub keyword: String,
    pub value: String,
    pub comment: String,
}

#[derive(Debug, Clone, Serialize)]
pub struct TableSampleJson {
    pub total_rows: usize,
    pub offset: usize,
    pub limit: usize,
    pub columns: Vec<String>,
    pub rows: Vec<Vec<Value>>,
}

#[derive(Debug, Clone, Serialize)]
pub struct SelectedHduDetail {
    pub index: usize,
    pub extname: String,
    pub kind: String,
    pub cards: Vec<HeaderCardJson>,
    pub table: Option<TableSampleJson>,
}

#[derive(Debug, Clone, Serialize)]
pub struct InspectResponse {
    pub path: String,
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
    pub gtis: Vec<GtiIntervalDetail>,
    pub selected_hdu: SelectedHduDetail,
}

/// Handle `GET /api/inspect?path=...&full=false&hdu=0&offset=0&limit=50`.
pub fn handle_inspect(
    params: &HashMap<String, String>,
    data_dir: &Path,
) -> ApiResult<InspectResponse> {
    let raw_path = query_required(params, "path")?;
    let path = resolve_existing_file(raw_path, data_dir)?;
    let full_scan = query_bool(params, "full", false)?;
    let offset = query_usize(params, "offset", 0)?;
    let limit = query_usize(params, "limit", 50)?.min(500);
    let hdu_param = params.get("hdu").map(|s| s.trim()).filter(|s| !s.is_empty());

    let meta = std::fs::metadata(&path)
        .map_err(|e| ApiError::bad_request(format!("failed to stat {}: {e}", path.display())))?;
    let fits = open_fits_clean(&path)?;

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
                    if v.is_finite() {
                        return Some(v);
                    }
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
    let mut gtis = Vec::new();

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

        let extname = hdu_display_name(&h);
        let kind_str = hdu_kind_name(h.kind());

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
                        let mut parsed = Vec::with_capacity(starts.len().min(2000));
                        for (&a, &b) in starts.iter().zip(stops.iter()) {
                            if !a.is_finite() || !b.is_finite() || b < a {
                                invalid_gti_rows += 1;
                            } else {
                                let d = b - a;
                                exp += d;
                                if parsed.len() < 2000 {
                                    parsed.push(GtiIntervalDetail {
                                        start: a,
                                        stop: b,
                                        duration: d,
                                    });
                                }
                            }
                        }
                        if gti_exposure_sec == 0.0 || upper_name == "GTI" {
                            gti_exposure_sec = exp;
                            gtis = parsed;
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
            tstart: hdr.get_f64("TSTART").filter(|v| v.is_finite()),
            tstop: hdr.get_f64("TSTOP").filter(|v| v.is_finite()),
            exposure: hdr.get_f64("EXPOSURE").filter(|v| v.is_finite()),
            datasum_ok,
            checksum_ok,
            has_checksum,
        });
    }

    // Resolve selected HDU for header cards + sample rows
    let sel_idx = match hdu_param {
        None => 0usize,
        Some(s) => {
            if let Ok(i) = s.parse::<usize>() {
                if i >= fits.len() {
                    return Err(ApiError::bad_request(format!(
                        "HDU index {i} out of range (file has {} HDUs)",
                        fits.len()
                    )));
                }
                i
            } else {
                fits.hdu_by_name(s)
                    .map(|h| h.index())
                    .map_err(|_| ApiError::bad_request(format!("HDU '{s}' not found")))?
            }
        }
    };

    let sel_h = fits
        .hdu(sel_idx)
        .map_err(|e| ApiError::bad_request(format!("failed reading HDU {sel_idx}: {e}")))?;
    let sel_cards: Vec<HeaderCardJson> = sel_h
        .header()
        .cards()
        .iter()
        .map(|c| HeaderCardJson {
            keyword: c.keyword.clone(),
            value: format_card_value(&c.value),
            comment: c.comment.clone().unwrap_or_default(),
        })
        .collect();

    let table = if *sel_h.kind() == HduKind::BinTable {
        sel_h
            .as_bintable()
            .ok()
            .map(|bt| extract_table_sample(&bt, offset, limit))
    } else {
        None
    };

    let selected_hdu = SelectedHduDetail {
        index: sel_idx,
        extname: hdu_display_name(&sel_h),
        kind: hdu_kind_name(sel_h.kind()),
        cards: sel_cards,
        table,
    };

    Ok(InspectResponse {
        path: path.to_string_lossy().into_owned(),
        file_size_bytes: meta.len(),
        num_hdus: fits.len(),
        telescop: get_str_any("TELESCOP"),
        instrume: get_str_any("INSTRUME"),
        detnam: get_str_any("DETNAM"),
        obs_id: get_str_any("OBS_ID"),
        object: get_str_any("OBJECT"),
        ra_obj: get_f64_any("RA_OBJ")
            .or_else(|| get_f64_any("RA_NOM"))
            .or_else(|| get_f64_any("RA_PNT")),
        dec_obj: get_f64_any("DEC_OBJ")
            .or_else(|| get_f64_any("DEC_NOM"))
            .or_else(|| get_f64_any("DEC_PNT")),
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
        gtis,
        selected_hdu,
    })
}

fn hdu_display_name(h: &astrosat_fits::HduRef<'_>) -> String {
    h.extname()
        .map(|s| s.trim().to_string())
        .filter(|s| !s.is_empty())
        .unwrap_or_else(|| {
            if h.index() == 0 {
                "PRIMARY".to_string()
            } else {
                format!("HDU_{}", h.index())
            }
        })
}

fn hdu_kind_name(k: &HduKind) -> String {
    match k {
        HduKind::PrimaryImage => "PrimaryImage".to_string(),
        HduKind::Image => "Image".to_string(),
        HduKind::BinTable => "BinTable".to_string(),
        HduKind::AsciiTable => "AsciiTable".to_string(),
        HduKind::Other(s) => s.clone(),
    }
}

fn format_card_value(v: &CardValue) -> String {
    match v {
        CardValue::Logical(b) => if *b { "T".to_string() } else { "F".to_string() },
        CardValue::Integer(i) => i.to_string(),
        CardValue::Float(f) => {
            if f.is_finite() {
                format!("{f}")
            } else {
                "NaN".to_string()
            }
        }
        CardValue::Str(s) => s.clone(),
        CardValue::Empty => String::new(),
    }
}

fn extract_table_sample(bt: &BinTableRef<'_>, offset: usize, limit: usize) -> TableSampleJson {
    let total_rows = bt.nrows();
    let start = offset.min(total_rows);
    let end = start.saturating_add(limit).min(total_rows);
    let cols = bt.columns();
    let col_names: Vec<String> = cols.iter().map(|c| c.name.clone()).collect();
    let mut rows = Vec::with_capacity(end.saturating_sub(start));

    for r in start..end {
        let mut row_vals = Vec::with_capacity(cols.len());
        for (ci, spec) in cols.iter().enumerate() {
            let bytes = bt.cell_bytes(r, ci);
            row_vals.push(decode_cell_json(bytes, spec));
        }
        rows.push(row_vals);
    }

    TableSampleJson {
        total_rows,
        offset: start,
        limit,
        columns: col_names,
        rows,
    }
}

fn decode_cell_json(bytes: &[u8], spec: &ColumnSpec) -> Value {
    if spec.var_desc != VarDescriptor::Fixed {
        return Value::String(format!("[var:{}]", spec.tform));
    }
    if spec.repeat == 0 {
        return Value::Null;
    }
    if spec.col_type == ColumnType::Ascii {
        let s = String::from_utf8_lossy(bytes);
        return Value::String(s.trim_end_matches([' ', '\0']).to_string());
    }
    let tscal = spec.tscal.unwrap_or(1.0);
    let tzero = spec.tzero.unwrap_or(0.0);
    let scaled = tscal != 1.0 || tzero != 0.0;

    if spec.repeat == 1 {
        return decode_scalar_elem(bytes, spec.col_type, scaled, tscal, tzero);
    }

    // Vector column: show compact preview
    let elem_sz = spec.col_type.elem_size();
    let max_show = 6usize;
    let show = spec.repeat.min(max_show);
    let mut parts = Vec::with_capacity(show);
    for i in 0..show {
        let off = i * elem_sz;
        if off + elem_sz <= bytes.len() {
            let v = decode_scalar_elem(&bytes[off..off + elem_sz], spec.col_type, scaled, tscal, tzero);
            parts.push(v.to_string());
        }
    }
    if spec.repeat > max_show {
        Value::String(format!("[{}, ... ({} elems)]", parts.join(", "), spec.repeat))
    } else {
        Value::String(format!("[{}]", parts.join(", ")))
    }
}

fn decode_scalar_elem(
    b: &[u8],
    ty: ColumnType,
    scaled: bool,
    tscal: f64,
    tzero: f64,
) -> Value {
    let f64_to_val = |v: f64| -> Value {
        if v.is_finite() {
            serde_json::Number::from_f64(v)
                .map(Value::Number)
                .unwrap_or_else(|| Value::String(format!("{v}")))
        } else {
            Value::String(format!("{v}"))
        }
    };

    match ty {
        ColumnType::Logical => Value::Bool(b.first().copied() == Some(b'T')),
        ColumnType::Bit => Value::Number((b.first().copied().unwrap_or(0) as u64).into()),
        ColumnType::Byte => {
            let raw = b.first().copied().unwrap_or(0) as i64;
            if scaled {
                f64_to_val(tscal * (raw as f64) + tzero)
            } else {
                Value::Number(raw.into())
            }
        }
        ColumnType::I16 if b.len() >= 2 => {
            let raw = i16::from_be_bytes([b[0], b[1]]) as i64;
            if scaled {
                if (tscal - 1.0).abs() < 1e-12 && (tzero - 32768.0).abs() < 1e-6 {
                    Value::Number((raw + 32768).into())
                } else {
                    f64_to_val(tscal * (raw as f64) + tzero)
                }
            } else {
                Value::Number(raw.into())
            }
        }
        ColumnType::I32 if b.len() >= 4 => {
            let raw = i32::from_be_bytes([b[0], b[1], b[2], b[3]]) as i64;
            if scaled {
                if (tscal - 1.0).abs() < 1e-12 && (tzero - 2_147_483_648.0).abs() < 1e-3 {
                    Value::Number((raw + 2_147_483_648).into())
                } else {
                    f64_to_val(tscal * (raw as f64) + tzero)
                }
            } else {
                Value::Number(raw.into())
            }
        }
        ColumnType::I64 if b.len() >= 8 => {
            let raw = i64::from_be_bytes([b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7]]);
            if scaled {
                f64_to_val(tscal * (raw as f64) + tzero)
            } else {
                Value::Number(raw.into())
            }
        }
        ColumnType::F32 if b.len() >= 4 => {
            let raw = f32::from_be_bytes([b[0], b[1], b[2], b[3]]) as f64;
            f64_to_val(if scaled { tscal * raw + tzero } else { raw })
        }
        ColumnType::F64 if b.len() >= 8 => {
            let raw = f64::from_be_bytes([b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7]]);
            f64_to_val(if scaled { tscal * raw + tzero } else { raw })
        }
        ColumnType::ComplexF32 if b.len() >= 8 => {
            let re = f32::from_be_bytes([b[0], b[1], b[2], b[3]]);
            let im = f32::from_be_bytes([b[4], b[5], b[6], b[7]]);
            Value::String(format!("{re}+{im}i"))
        }
        ColumnType::ComplexF64 if b.len() >= 16 => {
            let re = f64::from_be_bytes([b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7]]);
            let im = f64::from_be_bytes([b[8], b[9], b[10], b[11], b[12], b[13], b[14], b[15]]);
            Value::String(format!("{re}+{im}i"))
        }
        _ => Value::Null,
    }
}

// ---------------------------------------------------------------------------
// FITS Validate
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Serialize)]
pub struct ValidationFinding {
    pub hdu_index: usize,
    pub severity: String,
    pub message: String,
}

#[derive(Debug, Clone, Serialize)]
pub struct ValidateResponse {
    pub path: String,
    pub strict: bool,
    pub valid: bool,
    pub errors: usize,
    pub warnings: usize,
    pub num_hdus: usize,
    pub findings: Vec<ValidationFinding>,
}

/// Handle `GET /api/validate?path=...&strict=true`.
pub fn handle_validate(
    params: &HashMap<String, String>,
    data_dir: &Path,
) -> ApiResult<ValidateResponse> {
    let raw_path = query_required(params, "path")?;
    let path = resolve_existing_file(raw_path, data_dir)?;
    let strict = query_bool(params, "strict", true)?;

    let mut findings = Vec::new();
    let meta = std::fs::metadata(&path)
        .map_err(|e| ApiError::bad_request(format!("cannot stat {}: {e}", path.display())))?;

    if meta.len() == 0 {
        findings.push(ValidationFinding {
            hdu_index: 0,
            severity: "error".into(),
            message: "FITS file is empty (0 bytes)".into(),
        });
    } else if !meta.len().is_multiple_of(FITS_BLOCK_SIZE as u64) {
        findings.push(ValidationFinding {
            hdu_index: 0,
            severity: "error".into(),
            message: format!(
                "FITS file size ({} bytes) is not a multiple of 2880 bytes",
                meta.len()
            ),
        });
    }

    let fits = match FitsFile::open(&path) {
        Ok(f) => f,
        Err(e) => {
            findings.push(ValidationFinding {
                hdu_index: 0,
                severity: "error".into(),
                message: format!("FITS 4.0 parse failure: {e}"),
            });
            return Ok(finalize_validation(&path, 0, findings, strict));
        }
    };

    for issue in validate_ogip(&fits) {
        findings.push(ValidationFinding {
            hdu_index: issue.hdu_index,
            severity: match issue.severity {
                IssueSeverity::Error => "error".into(),
                IssueSeverity::Warning => "warning".into(),
            },
            message: issue.message,
        });
    }

    for hdu in fits.iter_hdus() {
        let idx = hdu.index();
        let hdr = hdu.header();

        if idx == 0 {
            if hdr.get_bool("SIMPLE") != Some(true) {
                findings.push(ValidationFinding {
                    hdu_index: 0,
                    severity: "error".into(),
                    message: "Primary HDU missing SIMPLE = T".into(),
                });
            }
            match hdr.get_i64("BITPIX") {
                Some(8 | 16 | 32 | 64 | -32 | -64) => {}
                other => findings.push(ValidationFinding {
                    hdu_index: 0,
                    severity: "error".into(),
                    message: format!("Primary HDU invalid BITPIX: {other:?}"),
                }),
            }
            if hdr.get_i64("NAXIS").is_none() {
                findings.push(ValidationFinding {
                    hdu_index: 0,
                    severity: "error".into(),
                    message: "Primary HDU missing NAXIS".into(),
                });
            }
            if fits.len() > 1 && hdr.get_bool("EXTEND") != Some(true) {
                findings.push(ValidationFinding {
                    hdu_index: 0,
                    severity: "warning".into(),
                    message: "File has extensions but Primary HDU lacks EXTEND = T".into(),
                });
            }
        } else if hdr.get_str("EXTNAME").is_none() {
            findings.push(ValidationFinding {
                hdu_index: idx,
                severity: "warning".into(),
                message: "Extension HDU missing EXTNAME keyword".into(),
            });
        }

        if hdr.get_str("CHECKSUM").is_none() || hdr.get_str("DATASUM").is_none() {
            findings.push(ValidationFinding {
                hdu_index: idx,
                severity: "warning".into(),
                message: "HDU missing CHECKSUM/DATASUM integrity keywords".into(),
            });
        }

        if let (Some(ts), Some(te)) = (hdr.get_f64("TSTART"), hdr.get_f64("TSTOP")) {
            if !ts.is_finite() || !te.is_finite() || te < ts {
                findings.push(ValidationFinding {
                    hdu_index: idx,
                    severity: "error".into(),
                    message: format!("Invalid time bounds TSTART={ts} > TSTOP={te}"),
                });
            }
        }

        if *hdu.kind() == HduKind::BinTable {
            if let Ok(bt) = hdu.as_bintable() {
                let ext = hdu
                    .extname()
                    .map(|s| s.trim().to_ascii_uppercase())
                    .unwrap_or_default();
                if ext.contains("GTI") {
                    if let (Ok(starts), Ok(stops)) =
                        (bt.read_col_f64("START"), bt.read_col_f64("STOP"))
                    {
                        let mut prev_stop = f64::NEG_INFINITY;
                        for (r, (&a, &b)) in starts.iter().zip(stops.iter()).enumerate() {
                            if !a.is_finite() || !b.is_finite() || b < a {
                                findings.push(ValidationFinding {
                                    hdu_index: idx,
                                    severity: "error".into(),
                                    message: format!(
                                        "GTI row {r} has invalid interval [{a}, {b}]"
                                    ),
                                });
                            } else {
                                if a < prev_stop {
                                    findings.push(ValidationFinding {
                                        hdu_index: idx,
                                        severity: "warning".into(),
                                        message: format!(
                                            "GTI row {r} overlaps or is unsorted ({a} < {prev_stop})"
                                        ),
                                    });
                                }
                                prev_stop = b;
                            }
                        }
                    }
                } else if ext == "ORBIT" {
                    let has_pos =
                        bt.column_index("X").is_some() || bt.column_index("POS_X").is_some();
                    let has_vel =
                        bt.column_index("VX").is_some() || bt.column_index("VEL_X").is_some();
                    if bt.column_index("TIME").is_none() || !has_pos || !has_vel {
                        findings.push(ValidationFinding {
                            hdu_index: idx,
                            severity: "error".into(),
                            message: "ORBIT HDU missing TIME or position/velocity columns".into(),
                        });
                    }
                }
            }
        }
    }

    Ok(finalize_validation(&path, fits.len(), findings, strict))
}

fn finalize_validation(
    path: &Path,
    num_hdus: usize,
    findings: Vec<ValidationFinding>,
    strict: bool,
) -> ValidateResponse {
    let errors = findings.iter().filter(|f| f.severity == "error").count();
    let warnings = findings.iter().filter(|f| f.severity == "warning").count();
    let valid = errors == 0 && (!strict || warnings == 0);
    ValidateResponse {
        path: path.to_string_lossy().into_owned(),
        strict,
        valid,
        errors,
        warnings,
        num_hdus,
        findings,
    }
}
