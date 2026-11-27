//! Fast, exact scientific data endpoints: Light Curves, Energy Spectra, and 2D FITS Images/WCS.

use std::collections::HashMap;
use std::path::{Path, PathBuf};

use astrosat_fits::{FitsFile, HduKind};
use astrosat_laxpc::background::ChannelMap;
use astrosat_laxpc::types::{utc_to_mission_day, LaxpcUnit};
use serde::Serialize;

use crate::http_util::{
    query_bool, query_required, query_usize, resolve_existing_file, ApiError, ApiResult,
};

// ---------------------------------------------------------------------------
// Light Curve API
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Serialize)]
pub struct GtiSpan {
    pub start: f64,
    pub stop: f64,
}

#[derive(Debug, Clone, Serialize)]
pub struct LightCurveResponse {
    pub path: String,
    pub hdu: String,
    pub title: String,
    pub instrume: Option<String>,
    pub detnam: Option<String>,
    pub object: Option<String>,
    pub timesys: Option<String>,
    pub units: String,
    pub tstart: Option<f64>,
    pub tstop: Option<f64>,
    pub timedel: Option<f64>,
    pub effective_timedel: Option<f64>,
    pub rebin: usize,
    pub effective_rebin: usize,
    pub raw_points: usize,
    pub exposure: Option<f64>,
    pub mean_rate: f64,
    pub min_rate: f64,
    pub max_rate: f64,
    pub time: Vec<f64>,
    pub rate: Vec<f64>,
    pub error: Vec<f64>,
    pub fracexp: Vec<f64>,
    pub gtis: Vec<GtiSpan>,
}

/// Handle `GET /api/lightcurve?path=...&hdu=...&rebin=1&max_points=4000`.
pub fn handle_lightcurve(
    params: &HashMap<String, String>,
    data_dir: &Path,
) -> ApiResult<LightCurveResponse> {
    let raw_path = query_required(params, "path")?;
    let path = resolve_existing_file(raw_path, data_dir)?;
    let rebin = query_usize(params, "rebin", 1)?;
    if rebin == 0 {
        return Err(ApiError::bad_request("rebin must be >= 1"));
    }
    let max_points = query_usize(params, "max_points", 4000)?;
    if max_points == 0 {
        return Err(ApiError::bad_request("max_points must be >= 1"));
    }
    let hdu_param = params
        .get("hdu")
        .map(|s| s.trim())
        .filter(|s| !s.is_empty());

    let fits = open_fits_clean(&path)?;
    let (hdu_idx, hdu_name) = select_lightcurve_hdu(&fits, hdu_param)?;
    let hdu = fits
        .hdu(hdu_idx)
        .map_err(|e| ApiError::bad_request(format!("failed to access HDU {hdu_idx}: {e}")))?;
    let hdr = hdu.header();
    let bt = hdu
        .as_bintable()
        .map_err(|e| ApiError::bad_request(format!("HDU '{hdu_name}' is not a BINTABLE: {e}")))?;

    let time_raw = bt
        .read_col_f64("TIME")
        .map_err(|e| ApiError::bad_request(format!("missing TIME column in '{hdu_name}': {e}")))?;

    let timedel = hdr
        .get_f64("TIMEDEL")
        .filter(|v| v.is_finite() && *v > 0.0)
        .or_else(|| infer_timedel(&time_raw));

    let (rate_raw, error_raw, units) = if bt.column_index("RATE").is_some() {
        let r = bt
            .read_col_f64("RATE")
            .map_err(|e| ApiError::bad_request(format!("failed reading RATE column: {e}")))?;
        let e = bt
            .read_col_f64("ERROR")
            .or_else(|_| bt.read_col_f64("STAT_ERR"))
            .unwrap_or_else(|_| {
                let dt = timedel.unwrap_or(1.0);
                r.iter()
                    .map(|&rv| if rv > 0.0 && dt > 0.0 { (rv / dt).sqrt() } else { 0.0 })
                    .collect()
            });
        (r, e, "counts/s".to_string())
    } else if bt.column_index("COUNTS").is_some() {
        let c = bt
            .read_col_f64("COUNTS")
            .map_err(|e| ApiError::bad_request(format!("failed reading COUNTS column: {e}")))?;
        let dt = timedel.unwrap_or(1.0);
        let e = bt
            .read_col_f64("ERROR")
            .or_else(|_| bt.read_col_f64("STAT_ERR"))
            .unwrap_or_else(|_| c.iter().map(|&cv| cv.max(0.0).sqrt()).collect());
        if dt > 0.0 {
            let r: Vec<f64> = c.iter().map(|&cv| cv / dt).collect();
            let err: Vec<f64> = e.iter().map(|&ev| ev / dt).collect();
            (r, err, "counts/s".to_string())
        } else {
            (c, e, "counts".to_string())
        }
    } else {
        return Err(ApiError::bad_request(format!(
            "HDU '{hdu_name}' has neither RATE nor COUNTS column"
        )));
    };

    let has_fracexp = bt.column_index("FRACEXP").is_some();
    let frac_raw = if has_fracexp {
        bt.read_col_f64("FRACEXP")
            .unwrap_or_else(|_| vec![1.0; time_raw.len()])
    } else {
        vec![1.0; time_raw.len()]
    };

    let raw_points = time_raw.len();
    let after_rebin = raw_points.div_ceil(rebin);
    let extra_factor = if after_rebin > max_points {
        after_rebin.div_ceil(max_points)
    } else {
        1
    };
    let effective_rebin = rebin.saturating_mul(extra_factor).max(1);

    let (time, rate, error, fracexp) =
        rebin_lightcurve(&time_raw, &rate_raw, &error_raw, &frac_raw, effective_rebin, has_fracexp);

    let (mean_rate, min_rate, max_rate) = summary_stats(&rate);
    let gtis = extract_gtis(&fits);

    let instrume = header_str_any(&fits, hdu_idx, "INSTRUME");
    let detnam = header_str_any(&fits, hdu_idx, "DETNAM");
    let object = header_str_any(&fits, hdu_idx, "OBJECT");
    let timesys = header_str_any(&fits, hdu_idx, "TIMESYS");
    let tstart = header_f64_any(&fits, hdu_idx, "TSTART")
        .or_else(|| time_raw.iter().copied().find(|v| v.is_finite()));
    let tstop = header_f64_any(&fits, hdu_idx, "TSTOP")
        .or_else(|| time_raw.iter().rev().copied().find(|v| v.is_finite()));
    let exposure = header_f64_any(&fits, hdu_idx, "EXPOSURE").or_else(|| {
        if !gtis.is_empty() {
            Some(gtis.iter().map(|g| (g.stop - g.start).max(0.0)).sum())
        } else {
            None
        }
    });

    let title = build_title(
        object.as_deref(),
        instrume.as_deref(),
        detnam.as_deref(),
        "Light Curve",
        &path,
    );

    Ok(LightCurveResponse {
        path: path.to_string_lossy().into_owned(),
        hdu: hdu_name,
        title,
        instrume,
        detnam,
        object,
        timesys,
        units,
        tstart,
        tstop,
        timedel,
        effective_timedel: timedel.map(|d| d * (effective_rebin as f64)),
        rebin,
        effective_rebin,
        raw_points,
        exposure,
        mean_rate,
        min_rate,
        max_rate,
        time,
        rate,
        error,
        fracexp,
        gtis,
    })
}

/// Exact integer rebinning with optional exposure (`fracexp`) weighting and exact quadrature error propagation.
///
/// When `use_fracexp_weights` is `false` or all `fracexp` in a group are equal, this reduces identically to
/// arithmetic mean rate `sum(r) / m` and quadrature error `sqrt(sum(err^2)) / m`.
pub fn rebin_lightcurve(
    time: &[f64],
    rate: &[f64],
    error: &[f64],
    fracexp: &[f64],
    factor: usize,
    use_fracexp_weights: bool,
) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
    let n = time.len().min(rate.len()).min(error.len()).min(fracexp.len());
    let k = factor.max(1);
    if k == 1 {
        let mut out_t = Vec::with_capacity(n);
        let mut out_r = Vec::with_capacity(n);
        let mut out_e = Vec::with_capacity(n);
        let mut out_f = Vec::with_capacity(n);
        for i in 0..n {
            if time[i].is_finite() && rate[i].is_finite() {
                out_t.push(time[i]);
                out_r.push(rate[i]);
                out_e.push(if error[i].is_finite() { error[i].max(0.0) } else { 0.0 });
                out_f.push(if fracexp[i].is_finite() { fracexp[i] } else { 1.0 });
            }
        }
        return (out_t, out_r, out_e, out_f);
    }

    let groups = n.div_ceil(k);
    let mut out_t = Vec::with_capacity(groups);
    let mut out_r = Vec::with_capacity(groups);
    let mut out_e = Vec::with_capacity(groups);
    let mut out_f = Vec::with_capacity(groups);

    for g in 0..groups {
        let start = g * k;
        let end = (start + k).min(n);

        let mut t_sum = 0.0;
        let mut w_sum = 0.0;
        let mut wr_sum = 0.0;
        let mut w2e2_sum = 0.0;
        let mut f_sum = 0.0;
        let mut valid_cnt = 0usize;

        for i in start..end {
            let ti = time[i];
            let ri = rate[i];
            if !ti.is_finite() || !ri.is_finite() {
                continue;
            }
            let ei = if error[i].is_finite() {
                error[i].max(0.0)
            } else {
                0.0
            };
            let fi = if fracexp[i].is_finite() {
                fracexp[i].max(0.0)
            } else {
                1.0
            };
            let w = if use_fracexp_weights { fi } else { 1.0 };
            t_sum += ti;
            f_sum += fi;
            valid_cnt += 1;
            if w > 0.0 {
                w_sum += w;
                wr_sum += w * ri;
                w2e2_sum += (w * ei) * (w * ei);
            }
        }

        if valid_cnt == 0 {
            continue;
        }
        let m = valid_cnt as f64;
        out_t.push(t_sum / m);
        out_f.push(f_sum / m);
        if w_sum > 0.0 {
            out_r.push(wr_sum / w_sum);
            out_e.push(w2e2_sum.sqrt() / w_sum);
        } else {
            out_r.push(0.0);
            out_e.push(0.0);
        }
    }

    (out_t, out_r, out_e, out_f)
}

fn select_lightcurve_hdu(fits: &FitsFile, requested: Option<&str>) -> ApiResult<(usize, String)> {
    if let Some(req) = requested {
        if let Ok(idx) = req.parse::<usize>() {
            let h = fits
                .hdu(idx)
                .map_err(|e| ApiError::bad_request(format!("HDU index {idx} not found: {e}")))?;
            let name = h
                .extname()
                .map(|s| s.trim().to_string())
                .filter(|s| !s.is_empty())
                .unwrap_or_else(|| format!("HDU_{idx}"));
            return Ok((idx, name));
        }
        let h = fits
            .hdu_by_name(req)
            .map_err(|e| ApiError::bad_request(format!("HDU '{req}' not found: {e}")))?;
        return Ok((h.index(), req.to_string()));
    }

    for candidate in ["RATE", "LIGHTCURVE", "LC"] {
        if let Ok(h) = fits.hdu_by_name(candidate) {
            return Ok((h.index(), candidate.to_string()));
        }
    }
    for h in fits.iter_hdus() {
        if *h.kind() == HduKind::BinTable {
            if let Ok(bt) = h.as_bintable() {
                if bt.column_index("TIME").is_some()
                    && (bt.column_index("RATE").is_some() || bt.column_index("COUNTS").is_some())
                {
                    let name = h
                        .extname()
                        .map(|s| s.trim().to_string())
                        .filter(|s| !s.is_empty())
                        .unwrap_or_else(|| format!("HDU_{}", h.index()));
                    return Ok((h.index(), name));
                }
            }
        }
    }
    Err(ApiError::bad_request(
        "no light curve BINTABLE HDU (with TIME + RATE/COUNTS) found in FITS file",
    ))
}

// ---------------------------------------------------------------------------
// Spectrum API
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Serialize)]
pub struct SpectrumResponse {
    pub path: String,
    pub hdu: String,
    pub title: String,
    pub instrume: Option<String>,
    pub detnam: Option<String>,
    pub object: Option<String>,
    pub chantype: Option<String>,
    pub exposure: Option<f64>,
    pub is_rate: bool,
    pub units: String,
    pub energy_mapped: bool,
    pub energy_source: String,
    pub total_channels: usize,
    pub good_channels: usize,
    pub total_rate: f64,
    pub channel: Vec<i32>,
    pub energy_kev: Vec<f64>,
    pub e_min_kev: Vec<f64>,
    pub e_max_kev: Vec<f64>,
    pub rate_or_counts: Vec<f64>,
    pub stat_err: Vec<f64>,
    pub sys_err: Vec<f64>,
    pub quality: Vec<i32>,
    pub grouping: Vec<i32>,
}

/// Handle `GET /api/spectrum?path=...&hdu=...&good_only=true&caldb=...`.
pub fn handle_spectrum(
    params: &HashMap<String, String>,
    data_dir: &Path,
) -> ApiResult<SpectrumResponse> {
    let raw_path = query_required(params, "path")?;
    let path = resolve_existing_file(raw_path, data_dir)?;
    let good_only = query_bool(params, "good_only", true)?;
    let hdu_param = params
        .get("hdu")
        .map(|s| s.trim())
        .filter(|s| !s.is_empty());
    let caldb_param = params
        .get("caldb")
        .map(|s| s.trim())
        .filter(|s| !s.is_empty())
        .map(PathBuf::from);

    let fits = open_fits_clean(&path)?;
    let (hdu_idx, hdu_name) = select_spectrum_hdu(&fits, hdu_param)?;
    let hdu = fits
        .hdu(hdu_idx)
        .map_err(|e| ApiError::bad_request(format!("failed to access HDU {hdu_idx}: {e}")))?;
    let hdr = hdu.header();
    let bt = hdu
        .as_bintable()
        .map_err(|e| ApiError::bad_request(format!("HDU '{hdu_name}' is not a BINTABLE: {e}")))?;

    let nrows = bt.nrows();
    let channels = bt
        .read_col_i32("CHANNEL")
        .unwrap_or_else(|_| (0..nrows as i32).collect());

    let (is_rate, values, units) = if bt.column_index("RATE").is_some() {
        let r = bt
            .read_col_f64("RATE")
            .map_err(|e| ApiError::bad_request(format!("failed reading RATE column: {e}")))?;
        (true, r, "counts/s".to_string())
    } else if bt.column_index("COUNTS").is_some() {
        let c = bt
            .read_col_f64("COUNTS")
            .map_err(|e| ApiError::bad_request(format!("failed reading COUNTS column: {e}")))?;
        (false, c, "counts".to_string())
    } else {
        return Err(ApiError::bad_request(format!(
            "HDU '{hdu_name}' has neither RATE nor COUNTS column"
        )));
    };

    let exposure = header_f64_any(&fits, hdu_idx, "EXPOSURE");
    let stat_err_all = bt
        .read_col_f64("STAT_ERR")
        .or_else(|_| bt.read_col_f64("ERROR"))
        .unwrap_or_else(|_| {
            let exp = exposure.unwrap_or(1.0);
            values
                .iter()
                .map(|&v| {
                    if is_rate {
                        if v > 0.0 && exp > 0.0 {
                            (v / exp).sqrt()
                        } else {
                            0.0
                        }
                    } else {
                        v.max(0.0).sqrt()
                    }
                })
                .collect()
        });

    let sys_err_all = bt
        .read_col_f64("SYS_ERR")
        .unwrap_or_else(|_| vec![hdr.get_f64("SYS_ERR").unwrap_or(0.0); nrows]);

    let quality_all = bt
        .read_col_i32("QUALITY")
        .unwrap_or_else(|_| vec![hdr.get_i64("QUALITY").unwrap_or(0) as i32; nrows]);

    let grouping_all = bt
        .read_col_i32("GROUPING")
        .unwrap_or_else(|_| vec![hdr.get_i64("GROUPING").unwrap_or(0) as i32; nrows]);

    let instrume = header_str_any(&fits, hdu_idx, "INSTRUME");
    let detnam = header_str_any(&fits, hdu_idx, "DETNAM");
    let object = header_str_any(&fits, hdu_idx, "OBJECT");
    let chantype = hdr
        .get_str("CHANTYPE")
        .map(|s| s.trim().to_string())
        .filter(|s| !s.is_empty());
    let respfile = hdr
        .get_str("RESPFILE")
        .map(|s| s.trim().to_string())
        .filter(|s| !s.is_empty() && !s.eq_ignore_ascii_case("NONE"));
    let tstart = header_f64_any(&fits, hdu_idx, "TSTART");

    let (e_min_all, e_max_all, energy_mapped, energy_source) = resolve_spectrum_energies(
        &fits,
        &path,
        &channels,
        instrume.as_deref(),
        detnam.as_deref(),
        chantype.as_deref(),
        respfile.as_deref(),
        tstart,
        caldb_param.as_deref(),
        data_dir,
    );

    let total_channels = channels.len();
    let good_channels = quality_all.iter().filter(|&&q| q == 0).count();

    let mut out_ch = Vec::with_capacity(total_channels);
    let mut out_ekev = Vec::with_capacity(total_channels);
    let mut out_emin = Vec::with_capacity(total_channels);
    let mut out_emax = Vec::with_capacity(total_channels);
    let mut out_val = Vec::with_capacity(total_channels);
    let mut out_stat = Vec::with_capacity(total_channels);
    let mut out_sys = Vec::with_capacity(total_channels);
    let mut out_qual = Vec::with_capacity(total_channels);
    let mut out_grp = Vec::with_capacity(total_channels);

    let mut total_val_sum = 0.0;
    for i in 0..total_channels {
        let q = quality_all[i];
        if good_only && q != 0 {
            continue;
        }
        let v = sanitize_f64(values[i]);
        let se = sanitize_f64(stat_err_all[i]).max(0.0);
        let sy = sanitize_f64(sys_err_all[i]).max(0.0);
        let elo = sanitize_f64(e_min_all[i]);
        let ehi = sanitize_f64(e_max_all[i]);
        out_ch.push(channels[i]);
        out_emin.push(elo);
        out_emax.push(ehi);
        out_ekev.push(0.5 * (elo + ehi));
        out_val.push(v);
        out_stat.push(se);
        out_sys.push(sy);
        out_qual.push(q);
        out_grp.push(grouping_all[i]);
        total_val_sum += v;
    }

    let total_rate = if is_rate {
        total_val_sum
    } else if let Some(exp) = exposure.filter(|&e| e > 0.0) {
        total_val_sum / exp
    } else {
        total_val_sum
    };

    let title = build_title(
        object.as_deref(),
        instrume.as_deref(),
        detnam.as_deref(),
        "Spectrum",
        &path,
    );

    Ok(SpectrumResponse {
        path: path.to_string_lossy().into_owned(),
        hdu: hdu_name,
        title,
        instrume,
        detnam,
        object,
        chantype,
        exposure,
        is_rate,
        units,
        energy_mapped,
        energy_source,
        total_channels,
        good_channels,
        total_rate,
        channel: out_ch,
        energy_kev: out_ekev,
        e_min_kev: out_emin,
        e_max_kev: out_emax,
        rate_or_counts: out_val,
        stat_err: out_stat,
        sys_err: out_sys,
        quality: out_qual,
        grouping: out_grp,
    })
}

fn select_spectrum_hdu(fits: &FitsFile, requested: Option<&str>) -> ApiResult<(usize, String)> {
    if let Some(req) = requested {
        if let Ok(idx) = req.parse::<usize>() {
            let h = fits
                .hdu(idx)
                .map_err(|e| ApiError::bad_request(format!("HDU index {idx} not found: {e}")))?;
            let name = h
                .extname()
                .map(|s| s.trim().to_string())
                .filter(|s| !s.is_empty())
                .unwrap_or_else(|| format!("HDU_{idx}"));
            return Ok((idx, name));
        }
        let h = fits
            .hdu_by_name(req)
            .map_err(|e| ApiError::bad_request(format!("HDU '{req}' not found: {e}")))?;
        return Ok((h.index(), req.to_string()));
    }
    if let Ok(h) = fits.hdu_by_name("SPECTRUM") {
        return Ok((h.index(), "SPECTRUM".to_string()));
    }
    for h in fits.iter_hdus() {
        if *h.kind() == HduKind::BinTable {
            if let Ok(bt) = h.as_bintable() {
                if bt.column_index("CHANNEL").is_some()
                    && (bt.column_index("RATE").is_some() || bt.column_index("COUNTS").is_some())
                {
                    let name = h
                        .extname()
                        .map(|s| s.trim().to_string())
                        .filter(|s| !s.is_empty())
                        .unwrap_or_else(|| format!("HDU_{}", h.index()));
                    return Ok((h.index(), name));
                }
            }
        }
    }
    Err(ApiError::bad_request(
        "no SPECTRUM BINTABLE HDU found in FITS file",
    ))
}

#[allow(clippy::too_many_arguments)]
fn resolve_spectrum_energies(
    fits: &FitsFile,
    spec_path: &Path,
    channels: &[i32],
    instrume: Option<&str>,
    detnam: Option<&str>,
    chantype: Option<&str>,
    respfile: Option<&str>,
    tstart: Option<f64>,
    caldb_explicit: Option<&Path>,
    data_dir: &Path,
) -> (Vec<f64>, Vec<f64>, bool, String) {
    let n = channels.len();

    // 1. EBOUNDS in the same FITS file
    if let Some((emin, emax)) = read_ebounds_from_fits(fits, channels) {
        return (emin, emax, true, "fits_ebounds".to_string());
    }

    // 2. Sibling or CALDB RESPFILE / .rsp matching stem
    if let Some(rsp_path) = find_response_file(spec_path, respfile, caldb_explicit, data_dir) {
        if let Ok(rfits) = open_fits_clean(&rsp_path) {
            if let Some((emin, emax)) = read_ebounds_from_fits(&rfits, channels) {
                let src = rsp_path
                    .file_name()
                    .map(|s| format!("ebounds:{}", s.to_string_lossy()))
                    .unwrap_or_else(|| "respfile_ebounds".to_string());
                return (emin, emax, true, src);
            }
        }
    }

    let inst_up = instrume.unwrap_or("").to_ascii_uppercase();
    let det_up = detnam.unwrap_or("").to_ascii_uppercase();
    let fname_lo = spec_path
        .file_name()
        .map(|s| s.to_string_lossy().to_ascii_lowercase())
        .unwrap_or_default();

    // 3. LAXPC ChannelMap (`lxp{10,20,30}.channel`)
    let is_laxpc = inst_up.contains("LAXPC")
        || det_up.starts_with("LX")
        || fname_lo.starts_with("lxp")
        || fname_lo.contains("backlxp");
    if is_laxpc {
        let unit = detect_laxpc_unit(detnam, &fname_lo).unwrap_or(LaxpcUnit::Lxp10);
        let ndays = tstart.map(utc_to_mission_day).unwrap_or(100.0);
        if let Some((cmap, src)) = find_laxpc_channel_map(unit, ndays, caldb_explicit, data_dir) {
            let (emin, emax) = map_channels_with_laxpc_map(channels, n, &cmap);
            return (emin, emax, true, src);
        }
    }

    // 4. CZTI standard PI formula: E_min = 5.0 + 0.5 * PI, E_max = 5.5 + 0.5 * PI
    let is_czti = inst_up.contains("CZTI")
        || fname_lo.contains("czt")
        || (chantype.map(|c| c.eq_ignore_ascii_case("PI")).unwrap_or(false) && n == 512);
    if is_czti {
        let emin: Vec<f64> = channels.iter().map(|&c| 5.0 + 0.5 * (c as f64)).collect();
        let emax: Vec<f64> = channels.iter().map(|&c| 5.5 + 0.5 * (c as f64)).collect();
        return (emin, emax, true, "czti_pi_formula".to_string());
    }

    // 5. Fallback: raw channel indices
    let emin: Vec<f64> = channels.iter().map(|&c| c as f64).collect();
    let emax: Vec<f64> = channels.iter().map(|&c| (c as f64) + 1.0).collect();
    (emin, emax, false, "channel_index".to_string())
}

fn read_ebounds_from_fits(fits: &FitsFile, channels: &[i32]) -> Option<(Vec<f64>, Vec<f64>)> {
    let h = fits.hdu_by_name("EBOUNDS").ok()?;
    let bt = h.as_bintable().ok()?;
    let eb_ch = bt.read_col_i32("CHANNEL").ok()?;
    let eb_min = bt.read_col_f64("E_MIN").ok()?;
    let eb_max = bt.read_col_f64("E_MAX").ok()?;
    if eb_ch.is_empty() || eb_ch.len() != eb_min.len() || eb_ch.len() != eb_max.len() {
        return None;
    }
    let mut map: HashMap<i32, (f64, f64)> = HashMap::with_capacity(eb_ch.len());
    for i in 0..eb_ch.len() {
        map.insert(eb_ch[i], (eb_min[i], eb_max[i]));
    }
    let mut out_min = Vec::with_capacity(channels.len());
    let mut out_max = Vec::with_capacity(channels.len());
    let mut matched = 0usize;
    for &c in channels {
        if let Some(&(lo, hi)) = map.get(&c) {
            out_min.push(lo);
            out_max.push(hi);
            matched += 1;
        } else {
            out_min.push(c as f64);
            out_max.push((c as f64) + 1.0);
        }
    }
    if matched > 0 {
        Some((out_min, out_max))
    } else {
        None
    }
}

fn find_response_file(
    spec_path: &Path,
    respfile: Option<&str>,
    caldb_explicit: Option<&Path>,
    data_dir: &Path,
) -> Option<PathBuf> {
    let parent = spec_path.parent();
    if let Some(rf) = respfile {
        let rfp = Path::new(rf);
        if rfp.is_file() {
            return Some(rfp.to_path_buf());
        }
        if let Some(dir) = parent {
            let p = dir.join(rf);
            if p.is_file() {
                return Some(p);
            }
        }
        for base in caldb_search_dirs(caldb_explicit, data_dir) {
            let p = base.join(rf);
            if p.is_file() {
                return Some(p);
            }
        }
    }
    // Check sibling `.rsp` with the same file stem (standard CZTI convention: `foo.pha` -> `foo.rsp`)
    let rsp_sibling = spec_path.with_extension("rsp");
    if rsp_sibling.is_file() {
        return Some(rsp_sibling);
    }
    None
}

pub fn detect_laxpc_unit(detnam: Option<&str>, lower_name: &str) -> Option<LaxpcUnit> {
    let d = detnam.unwrap_or("").to_ascii_uppercase();
    if d.contains("10") || d.ends_with('1') || lower_name.contains("lxp1") || lower_name.contains("lx10") {
        Some(LaxpcUnit::Lxp10)
    } else if d.contains("20") || d.ends_with('2') || lower_name.contains("lxp2") || lower_name.contains("lx20") {
        Some(LaxpcUnit::Lxp20)
    } else if d.contains("30") || d.ends_with('3') || lower_name.contains("lxp3") || lower_name.contains("lx30") {
        Some(LaxpcUnit::Lxp30)
    } else {
        None
    }
}

pub fn caldb_search_dirs(caldb_explicit: Option<&Path>, data_dir: &Path) -> Vec<PathBuf> {
    let mut dirs = Vec::new();
    if let Some(c) = caldb_explicit {
        dirs.push(c.to_path_buf());
    }
    let cwd = std::env::current_dir().unwrap_or_else(|_| data_dir.to_path_buf());
    for root in [data_dir.to_path_buf(), cwd] {
        for sub in [
            "caldb/laxpc",
            "ref_src/laxpc",
            "data/caldb/laxpc",
            "data/ref_src/laxpc",
        ] {
            let p = root.join(sub);
            if p.is_dir() && !dirs.iter().any(|d| d == &p) {
                dirs.push(p);
            }
        }
    }
    if let Some(parent) = data_dir.parent() {
        for sub in ["caldb/laxpc", "ref_src/laxpc"] {
            let p = parent.join(sub);
            if p.is_dir() && !dirs.iter().any(|d| d == &p) {
                dirs.push(p);
            }
        }
    }
    dirs
}

pub fn find_laxpc_channel_map(
    unit: LaxpcUnit,
    ndays: f64,
    caldb_explicit: Option<&Path>,
    data_dir: &Path,
) -> Option<(ChannelMap, String)> {
    for dir in caldb_search_dirs(caldb_explicit, data_dir) {
        if let Some(cm) = ChannelMap::from_caldb_for_day(&dir, unit, ndays) {
            let name = ChannelMap::filename_for_day(unit, ndays);
            return Some((cm, format!("laxpc_channel_map:{name}")));
        }
    }
    None
}

fn map_channels_with_laxpc_map(
    channels: &[i32],
    nch_total: usize,
    cmap: &ChannelMap,
) -> (Vec<f64>, Vec<f64>) {
    let map_len = cmap.len().max(1);
    let mut emin = Vec::with_capacity(channels.len());
    let mut emax = Vec::with_capacity(channels.len());

    if nch_total == map_len {
        for &c in channels {
            let idx = (c.max(0) as usize).min(map_len - 1);
            emin.push(cmap.e_min_kev[idx]);
            emax.push(cmap.e_max_kev[idx]);
        }
    } else if nch_total > map_len && nch_total.is_multiple_of(map_len) {
        // Finer spectrum (e.g. 1024 channels) vs native map (512 or 256): linearly subdivide each native bin.
        let factor = nch_total / map_len;
        let factor_f = factor as f64;
        for &c in channels {
            let cu = c.max(0) as usize;
            let native = (cu / factor).min(map_len - 1);
            let sub = (cu % factor) as f64;
            let lo0 = cmap.e_min_kev[native];
            let hi0 = cmap.e_max_kev[native];
            let de = (hi0 - lo0) / factor_f;
            emin.push(lo0 + sub * de);
            emax.push(lo0 + (sub + 1.0) * de);
        }
    } else if map_len > nch_total && map_len.is_multiple_of(nch_total.max(1)) {
        // Coarser spectrum (e.g. 256 channels) vs 512 native map: group consecutive native bins.
        let factor = map_len / nch_total.max(1);
        for &c in channels {
            let cu = c.max(0) as usize;
            let lo_idx = (cu * factor).min(map_len - 1);
            let hi_idx = ((cu + 1) * factor - 1).min(map_len - 1);
            emin.push(cmap.e_min_kev[lo_idx]);
            emax.push(cmap.e_max_kev[hi_idx]);
        }
    } else {
        // General fractional index mapping
        let ratio = (map_len as f64) / (nch_total.max(1) as f64);
        for &c in channels {
            let idx = (((c.max(0) as f64) * ratio).floor() as usize).min(map_len - 1);
            emin.push(cmap.e_min_kev[idx]);
            emax.push(cmap.e_max_kev[idx]);
        }
    }
    (emin, emax)
}

// ---------------------------------------------------------------------------
// 2D Image & WCS API
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Serialize)]
pub struct ImageHduOption {
    pub index: usize,
    pub extname: String,
    pub width: usize,
    pub height: usize,
}

#[derive(Debug, Clone, Serialize)]
pub struct ImagePeak {
    pub x: usize,
    pub y: usize,
    pub val: f64,
    pub snr: f64,
}

#[derive(Debug, Clone, Serialize)]
pub struct ImageStats {
    pub min: f64,
    pub max: f64,
    pub p1: f64,
    pub p50: f64,
    pub p99: f64,
    pub mean: f64,
    pub std: f64,
    pub peak: ImagePeak,
}

#[derive(Debug, Clone, Serialize)]
pub struct ImageWcs {
    pub has_wcs: bool,
    pub ctype1: Option<String>,
    pub ctype2: Option<String>,
    pub cunit1: Option<String>,
    pub cunit2: Option<String>,
    pub crpix1: Option<f64>,
    pub crpix2: Option<f64>,
    pub crval1: Option<f64>,
    pub crval2: Option<f64>,
    pub cdelt1: Option<f64>,
    pub cdelt2: Option<f64>,
    pub crota2: Option<f64>,
    pub cd: Option<[[f64; 2]; 2]>,
    pub pc: Option<[[f64; 2]; 2]>,
}

#[derive(Debug, Clone, Serialize)]
pub struct ImageResponse {
    pub path: String,
    pub hdu: String,
    pub hdu_index: usize,
    pub title: String,
    pub instrume: Option<String>,
    pub detnam: Option<String>,
    pub object: Option<String>,
    pub width: usize,
    pub height: usize,
    pub available_image_hdus: Vec<ImageHduOption>,
    pub is_czti_detector: bool,
    pub quad_overlay: bool,
    pub non_finite_pixels: usize,
    pub stats: ImageStats,
    pub wcs: ImageWcs,
    pub pixels: Vec<f64>,
}

/// Handle `GET /api/image?path=...&hdu=...`.
pub fn handle_image(
    params: &HashMap<String, String>,
    data_dir: &Path,
) -> ApiResult<ImageResponse> {
    let raw_path = query_required(params, "path")?;
    let path = resolve_existing_file(raw_path, data_dir)?;
    let hdu_param = params
        .get("hdu")
        .map(|s| s.trim())
        .filter(|s| !s.is_empty());

    let fits = open_fits_clean(&path)?;

    let mut available_image_hdus = Vec::new();
    for h in fits.iter_hdus() {
        if matches!(h.kind(), HduKind::PrimaryImage | HduKind::Image) {
            if let Ok(img) = h.as_image() {
                let sh = img.shape();
                if sh.len() == 2 && sh[0] > 0 && sh[1] > 0 {
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
                    available_image_hdus.push(ImageHduOption {
                        index: h.index(),
                        extname,
                        width: sh[0],
                        height: sh[1],
                    });
                }
            }
        }
    }

    if available_image_hdus.is_empty() {
        return Err(ApiError::bad_request(
            "no 2D IMAGE HDU found in FITS file",
        ));
    }

    let selected = if let Some(req) = hdu_param {
        if let Ok(idx) = req.parse::<usize>() {
            available_image_hdus
                .iter()
                .find(|opt| opt.index == idx)
                .cloned()
                .ok_or_else(|| {
                    ApiError::bad_request(format!("HDU {idx} is not a 2D IMAGE HDU"))
                })?
        } else {
            available_image_hdus
                .iter()
                .find(|opt| opt.extname.eq_ignore_ascii_case(req))
                .cloned()
                .ok_or_else(|| {
                    ApiError::bad_request(format!("2D IMAGE HDU '{req}' not found"))
                })?
        }
    } else {
        // Prefer FULL_DPI / FULL_DPH / SKY_IMAGE / IMAGE over single quadrants when available
        available_image_hdus
            .iter()
            .find(|opt| {
                let u = opt.extname.to_ascii_uppercase();
                u == "FULL_DPI" || u == "FULL_DPH" || u == "SKY_IMAGE" || u == "IMAGE"
            })
            .cloned()
            .unwrap_or_else(|| available_image_hdus[0].clone())
    };

    let hdu = fits
        .hdu(selected.index)
        .map_err(|e| ApiError::bad_request(format!("failed to read HDU {}: {e}", selected.index)))?;
    let hdr = hdu.header();
    let img = hdu
        .as_image()
        .map_err(|e| ApiError::bad_request(format!("HDU is not an IMAGE: {e}")))?;

    let width = selected.width;
    let height = selected.height;
    let raw_pixels = img
        .read_f64()
        .map_err(|e| ApiError::bad_request(format!("failed decoding IMAGE pixels: {e}")))?;

    let mut non_finite_pixels = 0usize;
    let mut finite_vals = Vec::with_capacity(raw_pixels.len());
    let mut pixels = Vec::with_capacity(raw_pixels.len());

    let mut sum = 0.0;
    let mut peak_val = f64::NEG_INFINITY;
    let mut peak_idx = 0usize;

    for (idx, &v) in raw_pixels.iter().enumerate() {
        if v.is_finite() {
            pixels.push(v);
            finite_vals.push(v);
            sum += v;
            if v > peak_val {
                peak_val = v;
                peak_idx = idx;
            }
        } else {
            non_finite_pixels += 1;
            pixels.push(0.0);
        }
    }

    let stats = if finite_vals.is_empty() {
        ImageStats {
            min: 0.0,
            max: 0.0,
            p1: 0.0,
            p50: 0.0,
            p99: 0.0,
            mean: 0.0,
            std: 0.0,
            peak: ImagePeak {
                x: 0,
                y: 0,
                val: 0.0,
                snr: 0.0,
            },
        }
    } else {
        let nf = finite_vals.len() as f64;
        let mean = sum / nf;
        let var = finite_vals
            .iter()
            .map(|&v| (v - mean) * (v - mean))
            .sum::<f64>()
            / nf;
        let std = var.sqrt();
        finite_vals.sort_by(|a, b| a.total_cmp(b));
        let min = finite_vals[0];
        let max = finite_vals[finite_vals.len() - 1];
        let p1 = percentile_sorted(&finite_vals, 0.01);
        let p50 = percentile_sorted(&finite_vals, 0.50);
        let p99 = percentile_sorted(&finite_vals, 0.99);
        let snr = if std > 0.0 {
            (peak_val - mean) / std
        } else {
            0.0
        };
        ImageStats {
            min,
            max,
            p1,
            p50,
            p99,
            mean,
            std,
            peak: ImagePeak {
                x: peak_idx % width,
                y: peak_idx / width,
                val: peak_val,
                snr,
            },
        }
    };

    let wcs = extract_wcs(hdr);
    let instrume = header_str_any(&fits, selected.index, "INSTRUME");
    let detnam = header_str_any(&fits, selected.index, "DETNAM");
    let object = header_str_any(&fits, selected.index, "OBJECT");

    let fname_lo = path
        .file_name()
        .map(|s| s.to_string_lossy().to_ascii_lowercase())
        .unwrap_or_default();
    let ext_up = selected.extname.to_ascii_uppercase();
    let is_dph_dpi = fname_lo.ends_with(".dph")
        || fname_lo.ends_with(".dpi")
        || ext_up.contains("DPH")
        || ext_up.contains("DPI");
    let quad_overlay = width == 128 && height == 128 && (is_dph_dpi || !wcs.has_wcs);
    let is_czti_detector = quad_overlay || (width == 64 && height == 64 && is_dph_dpi);

    let title = build_title(
        object.as_deref(),
        instrume.as_deref(),
        Some(&selected.extname),
        "Image",
        &path,
    );

    Ok(ImageResponse {
        path: path.to_string_lossy().into_owned(),
        hdu: selected.extname,
        hdu_index: selected.index,
        title,
        instrume,
        detnam,
        object,
        width,
        height,
        available_image_hdus,
        is_czti_detector,
        quad_overlay,
        non_finite_pixels,
        stats,
        wcs,
        pixels,
    })
}

fn extract_wcs(hdr: &astrosat_fits::Header) -> ImageWcs {
    let ctype1 = hdr.get_str("CTYPE1").map(|s| s.trim().to_string());
    let ctype2 = hdr.get_str("CTYPE2").map(|s| s.trim().to_string());
    let cunit1 = hdr.get_str("CUNIT1").map(|s| s.trim().to_string());
    let cunit2 = hdr.get_str("CUNIT2").map(|s| s.trim().to_string());
    let crpix1 = hdr.get_f64("CRPIX1").filter(|v| v.is_finite());
    let crpix2 = hdr.get_f64("CRPIX2").filter(|v| v.is_finite());
    let crval1 = hdr.get_f64("CRVAL1").filter(|v| v.is_finite());
    let crval2 = hdr.get_f64("CRVAL2").filter(|v| v.is_finite());
    let cdelt1 = hdr.get_f64("CDELT1").filter(|v| v.is_finite());
    let cdelt2 = hdr.get_f64("CDELT2").filter(|v| v.is_finite());
    let crota2 = hdr.get_f64("CROTA2").filter(|v| v.is_finite());

    let cd = match (
        hdr.get_f64("CD1_1"),
        hdr.get_f64("CD1_2"),
        hdr.get_f64("CD2_1"),
        hdr.get_f64("CD2_2"),
    ) {
        (Some(a), b, c, Some(d)) if a.is_finite() && d.is_finite() => Some([
            [a, b.filter(|v| v.is_finite()).unwrap_or(0.0)],
            [c.filter(|v| v.is_finite()).unwrap_or(0.0), d],
        ]),
        _ => None,
    };

    let pc = match (
        hdr.get_f64("PC1_1"),
        hdr.get_f64("PC1_2"),
        hdr.get_f64("PC2_1"),
        hdr.get_f64("PC2_2"),
    ) {
        (Some(a), b, c, Some(d)) if a.is_finite() && d.is_finite() => Some([
            [a, b.filter(|v| v.is_finite()).unwrap_or(0.0)],
            [c.filter(|v| v.is_finite()).unwrap_or(0.0), d],
        ]),
        _ => None,
    };

    let has_wcs = crpix1.is_some()
        && crpix2.is_some()
        && crval1.is_some()
        && crval2.is_some()
        && ((cdelt1.is_some() && cdelt2.is_some()) || cd.is_some());

    ImageWcs {
        has_wcs,
        ctype1,
        ctype2,
        cunit1,
        cunit2,
        crpix1,
        crpix2,
        crval1,
        crval2,
        cdelt1,
        cdelt2,
        crota2,
        cd,
        pc,
    }
}

// ---------------------------------------------------------------------------
// Shared helpers
// ---------------------------------------------------------------------------

pub fn open_fits_clean(path: &Path) -> ApiResult<FitsFile> {
    // Transparently support gzip-compressed FITS (e.g. CALDB `.fits.gz`)
    astrosat_czti::fits_io::open_fits_maybe_gz(path).map_err(|e| {
        ApiError::bad_request(format!("failed to parse FITS file {}: {e}", path.display()))
    })
}

pub fn header_str_any(fits: &FitsFile, primary_hdu: usize, key: &str) -> Option<String> {
    for idx in [primary_hdu, 0, 1] {
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
}

pub fn header_f64_any(fits: &FitsFile, primary_hdu: usize, key: &str) -> Option<f64> {
    for idx in [primary_hdu, 0, 1] {
        if let Ok(h) = fits.hdu(idx) {
            if let Some(v) = h.header().get_f64(key) {
                if v.is_finite() {
                    return Some(v);
                }
            }
        }
    }
    None
}

pub fn extract_gtis(fits: &FitsFile) -> Vec<GtiSpan> {
    for name in ["GTI", "STDGTI", "Q0_GTI"] {
        if let Ok(h) = fits.hdu_by_name(name) {
            if let Ok(bt) = h.as_bintable() {
                if let (Ok(starts), Ok(stops)) = (bt.read_col_f64("START"), bt.read_col_f64("STOP"))
                {
                    return starts
                        .into_iter()
                        .zip(stops)
                        .filter(|&(a, b)| a.is_finite() && b.is_finite() && b >= a)
                        .map(|(start, stop)| GtiSpan { start, stop })
                        .collect();
                }
            }
        }
    }
    Vec::new()
}

fn infer_timedel(times: &[f64]) -> Option<f64> {
    if times.len() < 2 {
        return None;
    }
    let mut best = f64::INFINITY;
    for w in times.windows(2).take(1024) {
        let dt = w[1] - w[0];
        if dt.is_finite() && dt > 0.0 && dt < best {
            best = dt;
        }
    }
    if best.is_finite() {
        Some(best)
    } else {
        None
    }
}

fn summary_stats(vals: &[f64]) -> (f64, f64, f64) {
    let mut sum = 0.0;
    let mut cnt = 0usize;
    let mut min_v = f64::INFINITY;
    let mut max_v = f64::NEG_INFINITY;
    for &v in vals {
        if v.is_finite() {
            sum += v;
            cnt += 1;
            if v < min_v {
                min_v = v;
            }
            if v > max_v {
                max_v = v;
            }
        }
    }
    if cnt == 0 {
        (0.0, 0.0, 0.0)
    } else {
        (sum / (cnt as f64), min_v, max_v)
    }
}

fn percentile_sorted(sorted: &[f64], p: f64) -> f64 {
    if sorted.is_empty() {
        return 0.0;
    }
    if sorted.len() == 1 {
        return sorted[0];
    }
    let pos = p.clamp(0.0, 1.0) * ((sorted.len() - 1) as f64);
    let lo = pos.floor() as usize;
    let hi = (lo + 1).min(sorted.len() - 1);
    let frac = pos - (lo as f64);
    sorted[lo] * (1.0 - frac) + sorted[hi] * frac
}

#[inline]
pub fn sanitize_f64(v: f64) -> f64 {
    if v.is_finite() {
        v
    } else {
        0.0
    }
}

fn build_title(
    object: Option<&str>,
    instrume: Option<&str>,
    detnam: Option<&str>,
    suffix: &str,
    path: &Path,
) -> String {
    let mut parts = Vec::new();
    if let Some(o) = object {
        parts.push(o.to_string());
    }
    match (instrume, detnam) {
        (Some(i), Some(d)) if !d.eq_ignore_ascii_case(i) => parts.push(format!("{i} {d}")),
        (Some(i), _) => parts.push(i.to_string()),
        (None, Some(d)) => parts.push(d.to_string()),
        _ => {}
    }
    if parts.is_empty() {
        let stem = path
            .file_name()
            .map(|s| s.to_string_lossy().into_owned())
            .unwrap_or_else(|| suffix.to_string());
        format!("{stem} — {suffix}")
    } else {
        format!("{} — {suffix}", parts.join(" · "))
    }
}
