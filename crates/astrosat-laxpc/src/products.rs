//! Stage 5: Non-paralyzable deadtime correction and OGIP Level-2 product generation.
//!
//! Supports both **Format B** per-unit TIFR/LaxpcSoft products (`.evn`, `.lcurv`, `.spec`,
//! `.gtis`, `AstroSat.orb`, plus optional ASCII `.pha`/`.lc`) and **Format A** combined
//! multi-unit `level2.event.fits` (`TIME`, `Channel`, `Layer`, `LAXPC_No`) with complete
//! HDU-0, HDU-1, and GTI timing keywords so `astrosat bary` works natively without `preplaxbary`.

use crate::ea::{passes_event_filter, EaEvent, EventKind};
use crate::error::{LaxpcError, Result};
use crate::types::{
    AnodeSelection, EventFilter, LaxpcUnit, DEADTIME_RATE_GUARD_CPS,
};
use astrosat_fits::{
    apply_ogip_meta, build_gti_hdu, build_orbit_hdu, write_fits, BinTableBuilder, Hdu, OgipMeta,
};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

/// Apply non-paralyzable deadtime correction $R_{\text{corr}} = R_{\text{obs}} / (1 - \tau R_{\text{all}})$,
/// guarding against the $R_{\text{all}} \ge 23000\text{ cps}$ singularity (`laxpcl1.f:917,1303,1660`).
#[inline]
pub fn deadtime_correct_rate(obs_rate_cps: f64, all_rate_cps: f64, deadtime_sec: f64) -> f64 {
    if all_rate_cps > 0.0 && all_rate_cps < DEADTIME_RATE_GUARD_CPS {
        let denom = 1.0 - deadtime_sec * all_rate_cps;
        if denom > 0.0 {
            return obs_rate_cps / denom;
        }
    }
    obs_rate_cps
}

/// Compute deadtime correction factor $f_{\text{cd}} = 1 / (1 - \tau R_{\text{all}})$, guarded at 23000 cps.
#[inline]
pub fn deadtime_factor(all_rate_cps: f64, deadtime_sec: f64) -> f64 {
    if all_rate_cps > 0.0 && all_rate_cps < DEADTIME_RATE_GUARD_CPS {
        let denom = 1.0 - deadtime_sec * all_rate_cps;
        if denom > 0.0 {
            return 1.0 / denom;
        }
    }
    1.0
}

/// Compute total overlap duration between a sorted list of `[start, stop]` GTI intervals
/// and a sorted list of `[gap_start, gap_stop]` internal frame gaps (`laxpcl1.f:881-896`).
pub fn compute_gti_frame_gap_loss(gti: &[(f64, f64)], frame_gaps: &[(f64, f64)]) -> f64 {
    let mut tgap = 0.0_f64;
    let mut jpg = 0usize;
    for &(tlo, thi) in gti {
        while jpg < frame_gaps.len() && frame_gaps[jpg].1 <= tlo {
            jpg += 1;
        }
        let mut igp = jpg;
        while igp < frame_gaps.len() {
            let (g0, g1) = frame_gaps[igp];
            if g0 >= thi {
                break;
            }
            if g1 > tlo {
                let lo = tlo.max(g0);
                let hi = thi.min(g1);
                if hi > lo {
                    tgap += hi - lo;
                }
            }
            igp += 1;
        }
    }
    tgap
}

/// Compute fractional frame-gap loss within each uniform light-curve bin `[t0 + i*tbin, t0 + (i+1)*tbin]`
/// (`laxpcl1.f:1250-1263`).
pub fn compute_bin_gap_fractions(
    t0_ins: f64,
    tbin: f64,
    nbins: usize,
    frame_gaps: &[(f64, f64)],
) -> Vec<f64> {
    let mut fracs = vec![0.0_f64; nbins];
    if frame_gaps.is_empty() || tbin <= 0.0 {
        return fracs;
    }
    let mut jpg = 0usize;
    for (i, frac) in fracs.iter_mut().enumerate() {
        let bin_lo = t0_ins + (i as f64) * tbin;
        let bin_hi = bin_lo + tbin;
        while jpg < frame_gaps.len() && frame_gaps[jpg].1 <= bin_lo {
            jpg += 1;
        }
        let mut igp = jpg;
        let mut loss = 0.0_f64;
        while igp < frame_gaps.len() {
            let (g0, g1) = frame_gaps[igp];
            if g0 >= bin_hi {
                break;
            }
            if g1 > bin_lo {
                loss += bin_hi.min(g1) - bin_lo.max(g0);
            }
            igp += 1;
        }
        *frac = (loss / tbin).clamp(0.0, 1.0);
    }
    fracs
}

/// Compute running 1-second mean of `coall` when `tbin < 0.2 s` (`laxpcl1.f:1226,1273-1286,1647-1659`).
/// For `tbin >= 0.2 s`, returns a direct clone of `coall`.
pub fn smooth_all_counts_for_deadtime(coall: &[f64], tbin: f64) -> Vec<f64> {
    let nt = coall.len();
    if nt == 0 || tbin >= 0.2 {
        return coall.to_vec();
    }
    let ntb = (0.5 / tbin + 0.2).floor() as usize;
    if ntb == 0 {
        return coall.to_vec();
    }

    // Prefix sum for O(N) exact running mean instead of O(N * ntb)
    let mut prefix = Vec::with_capacity(nt + 1);
    prefix.push(0.0_f64);
    let mut acc = 0.0_f64;
    for &v in coall {
        acc += v;
        prefix.push(acc);
    }

    let mut out = vec![0.0_f64; nt];
    for (i, slot) in out.iter_mut().enumerate() {
        let mut itbl = i.saturating_sub(ntb);
        if i + ntb >= nt {
            itbl = nt.saturating_sub(2 * ntb + 1);
        }
        let itbu = (itbl + 2 * ntb).min(nt - 1);
        let count = (itbu - itbl + 1) as f64;
        *slot = (prefix[itbu + 1] - prefix[itbl]) / count;
    }
    out
}

/// Level-2 deadtime-corrected LAXPC energy spectrum (all anodes + selected anode spectrum).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LaxpcSpectrum {
    pub unit: LaxpcUnit,
    /// Number of output channels (`256`, `512`, or `1024`).
    pub num_channels: usize,
    /// Net exposure time in seconds (GTI duration minus internal frame gaps).
    pub exposure_sec: f64,
    /// Raw total main-anode counts in GTI (including rejected doubles if `single_only`).
    pub total_raw_counts: u64,
    /// Deadtime-corrected total main-anode count rate (cps).
    pub total_corrected_rate_cps: f64,
    /// Deadtime correction factor applied (`1 / (1 - tau * R_all)`).
    pub deadtime_factor: f64,
    /// Per-channel, per-anode deadtime-corrected count rates `[num_channels][11]`
    /// (columns `0..10` = A1..A10, column `10` = sum A1..A7).
    pub anode_rates: Vec<[f64; 11]>,
    /// Selected-anode spectrum rate (cps) for FITS `SPECTRUM` `RATE` column.
    pub selected_rate: Vec<f64>,
    /// Statistical uncertainty `sqrt(rate / exposure)` for FITS `SPECTRUM` `STAT_ERR` column.
    pub selected_stat_err: Vec<f64>,
}

/// Accumulate and deadtime-correct a LAXPC spectrum from decoded EA events using `rayon` parallel reduction,
/// matching `laxpcl1.f:621-643,905-956`.
pub fn build_spectrum(
    events: &[EaEvent],
    gti_ins: &[(f64, f64)],
    frame_gaps_ins: &[(f64, f64)],
    unit: LaxpcUnit,
    filter: &EventFilter,
    num_channels: usize,
    deadtime_sec: f64,
) -> Result<LaxpcSpectrum> {
    if !matches!(num_channels, 256 | 512 | 1024) {
        return Err(LaxpcError::Config(format!(
            "num_channels must be 256, 512, or 1024 (got {num_channels})"
        )));
    }

    let raw_gti_dur: f64 = gti_ins.iter().map(|&(a, b)| (b - a).max(0.0)).sum();
    let gap_loss = compute_gti_frame_gap_loss(gti_ins, frame_gaps_ins);
    let exposure_sec = (raw_gti_dur - gap_loss).max(0.0);

    // Parallel histogramming over events using rayon fold/reduce
    let (raw_sp, rejected_doubles) = events
        .par_chunks(65_536)
        .map(|chunk| {
            let mut local_sp = vec![[0u64; 10]; 1024];
            let mut local_rej_dbl = 0u64;
            for ev in chunk {
                if !in_sorted_intervals(gti_ins, ev.time_ins) {
                    continue;
                }
                let ch = ev.pha as usize;
                let an = ev.anode as usize;
                if ch >= 1024 || an == 0 || an > 10 {
                    continue;
                }
                if filter.single_only && ev.kind == EventKind::Double {
                    // Track rejected main-anode doubles so deadtime R_all stays exact (`laxpcl1.f:623-625`)
                    if an <= 7 {
                        local_rej_dbl += 1;
                    }
                } else {
                    local_sp[ch][an - 1] += 1;
                }
            }
            (local_sp, local_rej_dbl)
        })
        .reduce(
            || (vec![[0u64; 10]; 1024], 0u64),
            |(mut a, ra), (b, rb)| {
                for ch in 0..1024 {
                    for an in 0..10 {
                        a[ch][an] += b[ch][an];
                    }
                }
                (a, ra + rb)
            },
        );

    let mut total_main_counts: u64 = 0;
    for ch_row in &raw_sp {
        for &c in &ch_row[..7] {
            total_main_counts += c;
        }
    }
    if filter.single_only {
        total_main_counts += rejected_doubles;
    }

    let raw_all_rate = if exposure_sec > 0.0 {
        total_main_counts as f64 / exposure_sec
    } else {
        0.0
    };
    let fcd = deadtime_factor(raw_all_rate, deadtime_sec);
    let total_corrected_rate_cps = raw_all_rate * fcd;

    let ncf = 1024 / num_channels;
    let mut anode_rates = vec![[0.0_f64; 11]; num_channels];
    let mut selected_rate = vec![0.0_f64; num_channels];
    let mut selected_stat_err = vec![0.0_f64; num_channels];

    if exposure_sec > 0.0 {
        let scale = fcd / exposure_sec;
        for out_ch in 0..num_channels {
            let mut sum_main = 0.0_f64;
            for an in 0..10 {
                let mut c = 0u64;
                for sub in 0..ncf {
                    c += raw_sp[out_ch * ncf + sub][an];
                }
                let r = (c as f64) * scale;
                anode_rates[out_ch][an] = r;
                if an < 7 {
                    sum_main += r;
                }
            }
            anode_rates[out_ch][10] = sum_main;

            let sel: f64 = match &filter.anodes {
                AnodeSelection::AllMain => sum_main,
                AnodeSelection::Layer1 => anode_rates[out_ch][0] + anode_rates[out_ch][1],
                AnodeSelection::Layer2 => anode_rates[out_ch][2] + anode_rates[out_ch][3],
                AnodeSelection::SingleAnode(a) => {
                    if (1..=10).contains(a) {
                        anode_rates[out_ch][(*a - 1) as usize]
                    } else {
                        0.0
                    }
                }
                AnodeSelection::Custom(list) => list
                    .iter()
                    .filter(|&&a| (1..=10).contains(&a))
                    .map(|&a| anode_rates[out_ch][(a - 1) as usize])
                    .sum(),
            };
            selected_rate[out_ch] = sel;
            selected_stat_err[out_ch] = (sel / exposure_sec).max(0.0).sqrt();
        }
    }

    Ok(LaxpcSpectrum {
        unit,
        num_channels,
        exposure_sec,
        total_raw_counts: total_main_counts,
        total_corrected_rate_cps,
        deadtime_factor: fcd,
        anode_rates,
        selected_rate,
        selected_stat_err,
    })
}

/// Level-2 deadtime-corrected LAXPC light curve (GTI-filtered bins).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LaxpcLightCurve {
    pub unit: LaxpcUnit,
    /// Bin width in seconds (`TIMEDEL`).
    pub tbin: f64,
    /// Bin center times in AstroSat UTC seconds.
    pub time_utc: Vec<f64>,
    /// Bin center times in instrument seconds.
    pub time_ins: Vec<f64>,
    /// Deadtime-corrected count rate (cps) in each GTI bin.
    pub rate_cps: Vec<f64>,
    /// Statistical uncertainty `sqrt(rate / tbin)` (cps).
    pub error_cps: Vec<f64>,
    /// Fractional bin exposure (`1.0 - frame_gap_frac`).
    pub fracexp: Vec<f64>,
    /// Optional modeled background rate (cps) per bin (same length as `time_utc` when present).
    pub bkg_rate_cps: Option<Vec<f64>>,
    /// Total raw filtered counts across all GTI bins.
    pub total_counts: u64,
    /// Mean deadtime-corrected rate across GTI bins.
    pub mean_rate_cps: f64,
}

/// Raw binned counts across the full observation timeline (before GTI trimming).
#[derive(Debug, Clone)]
pub struct RawBinnedCounts {
    /// Start of bin 0 in instrument time (`t0 = floor(first_t / tbin + 1) * tbin`, matching `laxpcl1.f:560-561`).
    pub t0_ins: f64,
    pub tbin: f64,
    /// Filtered counts per bin (`co`).
    pub selected_counts: Vec<f64>,
    /// All main-anode (A1..A7) counts per bin (`coall`), used for GTI & deadtime.
    pub all_main_counts: Vec<f64>,
}

/// Bin EA events into uniform `tbin` bins in parallel using `rayon`.
pub fn bin_events_raw(
    events: &[EaEvent],
    tbin: f64,
    filter: &EventFilter,
) -> Result<RawBinnedCounts> {
    if tbin <= 0.0 || !tbin.is_finite() {
        return Err(LaxpcError::Config(format!("invalid tbin: {tbin}")));
    }
    if events.is_empty() {
        return Ok(RawBinnedCounts {
            t0_ins: 0.0,
            tbin,
            selected_counts: Vec::new(),
            all_main_counts: Vec::new(),
        });
    }

    let first_t = events[0].time_ins;
    let t0_ins = ((first_t / tbin).floor() + 1.0) * tbin;
    let max_t = events
        .par_iter()
        .map(|e| e.time_ins)
        .reduce(|| first_t, f64::max);

    let max_idx = (((max_t - t0_ins + 1.0e-6) / tbin).floor() as isize).max(0) as usize;
    let nbins = max_idx + 1;

    // Parallel chunked binning (fallback to single-threaded if nbins is huge to bound memory)
    let (selected_counts, all_main_counts) = if nbins <= 2_000_000 {
        events
            .par_chunks(131_072)
            .map(|chunk| {
                let mut sel = vec![0.0_f64; nbins];
                let mut all = vec![0.0_f64; nbins];
                for ev in chunk {
                    let idx_i = ((ev.time_ins - t0_ins + 1.0e-6) / tbin).floor() as isize;
                    if idx_i < 0 || (idx_i as usize) >= nbins {
                        continue;
                    }
                    let idx = idx_i as usize;
                    if (1..=7).contains(&ev.anode) {
                        all[idx] += 1.0;
                    }
                    if passes_event_filter(ev, filter) {
                        sel[idx] += 1.0;
                    }
                }
                (sel, all)
            })
            .reduce(
                || (vec![0.0_f64; nbins], vec![0.0_f64; nbins]),
                |(mut s1, mut a1), (s2, a2)| {
                    for i in 0..nbins {
                        s1[i] += s2[i];
                        a1[i] += a2[i];
                    }
                    (s1, a1)
                },
            )
    } else {
        let mut sel = vec![0.0_f64; nbins];
        let mut all = vec![0.0_f64; nbins];
        for ev in events {
            let idx_i = ((ev.time_ins - t0_ins + 1.0e-6) / tbin).floor() as isize;
            if idx_i < 0 || (idx_i as usize) >= nbins {
                continue;
            }
            let idx = idx_i as usize;
            if (1..=7).contains(&ev.anode) {
                all[idx] += 1.0;
            }
            if passes_event_filter(ev, filter) {
                sel[idx] += 1.0;
            }
        }
        (sel, all)
    };

    Ok(RawBinnedCounts {
        t0_ins,
        tbin,
        selected_counts,
        all_main_counts,
    })
}

/// Build a GTI-filtered, frame-gap-corrected, deadtime-corrected light curve (`laxpcl1.f:1250-1286,1623-1664`).
pub fn build_lightcurve<F>(
    raw: &RawBinnedCounts,
    gti_ins: &[(f64, f64)],
    frame_gaps_ins: &[(f64, f64)],
    unit: LaxpcUnit,
    deadtime_sec: f64,
    ins_to_utc: F,
    raw_bkg_per_bin: Option<&[f64]>,
) -> LaxpcLightCurve
where
    F: Fn(f64) -> f64,
{
    let nt = raw.selected_counts.len();
    let tbin = raw.tbin;
    let gap_fracs = compute_bin_gap_fractions(raw.t0_ins, tbin, nt, frame_gaps_ins);

    // Apply per-bin frame gap compensation (`laxpcl1.f:1265-1271`)
    let mut co_adj = vec![0.0_f64; nt];
    let mut coall_adj = vec![0.0_f64; nt];
    for i in 0..nt {
        if gap_fracs[i] < 0.90 {
            let live = 1.0 - gap_fracs[i];
            co_adj[i] = raw.selected_counts[i] / live;
            coall_adj[i] = raw.all_main_counts[i] / live;
        }
    }

    // Running 1-s mean of all-anode counts when tbin < 0.2 s
    let coall_smooth = smooth_all_counts_for_deadtime(&coall_adj, tbin);

    let mut time_utc = Vec::new();
    let mut time_ins = Vec::new();
    let mut rate_cps = Vec::new();
    let mut error_cps = Vec::new();
    let mut fracexp = Vec::new();
    let mut bkg_out = raw_bkg_per_bin.map(|_| Vec::new());
    let mut total_counts = 0u64;

    for i in 0..nt {
        let tc = raw.t0_ins + (i as f64 + 0.5) * tbin;
        let t_prev = tc - tbin;
        if !bin_inside_gti(gti_ins, t_prev, tc) {
            continue;
        }
        let obs_rate = co_adj[i] / tbin;
        let all_rate = coall_smooth[i] / tbin;
        let corr_rate = deadtime_correct_rate(obs_rate, all_rate, deadtime_sec);
        let err = (corr_rate / tbin).max(0.0).sqrt();

        time_ins.push(tc);
        time_utc.push(ins_to_utc(tc));
        rate_cps.push(corr_rate);
        error_cps.push(err);
        fracexp.push((1.0 - gap_fracs[i]).clamp(0.0, 1.0));
        if let (Some(src_bk), Some(dst_bk)) = (raw_bkg_per_bin, bkg_out.as_mut()) {
            dst_bk.push(src_bk.get(i).copied().unwrap_or(0.0));
        }
        total_counts += co_adj[i].round() as u64;
    }

    let mean_rate_cps = if rate_cps.is_empty() {
        0.0
    } else {
        rate_cps.iter().sum::<f64>() / (rate_cps.len() as f64)
    };

    LaxpcLightCurve {
        unit,
        tbin,
        time_utc,
        time_ins,
        rate_cps,
        error_cps,
        fracexp,
        bkg_rate_cps: bkg_out,
        total_counts,
        mean_rate_cps,
    }
}

/// Level-2 calibrated event in AstroSat UTC time ready for Format B or Format A FITS export.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CalibratedEvent {
    /// AstroSat UTC time (seconds since MJD 55197).
    pub time_utc: f64,
    /// Scaled output channel (`pha / (1024 / num_channels)`).
    pub channel: i16,
    /// Anode ID (`1..=10`).
    pub anode_id: u8,
    /// Physical LAXPC layer (`1..=5` for A1..A7, `6..=8` for A8..A10).
    pub layer: u8,
    /// Detector unit ID (`1`, `2`, or `3`).
    pub detector_id: u8,
}

/// Filter decoded events by GTI and `EventFilter`, convert timestamps via `ins_to_utc`,
/// and scale PHA channels to `num_channels` in parallel using `rayon`.
pub fn build_calibrated_events<F>(
    events: &[EaEvent],
    gti_ins: &[(f64, f64)],
    filter: &EventFilter,
    num_channels: usize,
    ins_to_utc: F,
) -> Vec<CalibratedEvent>
where
    F: Fn(f64) -> f64 + Sync,
{
    let ncf = (1024 / num_channels.max(1)).max(1) as u16;
    events
        .par_iter()
        .filter(|ev| in_sorted_intervals(gti_ins, ev.time_ins) && passes_event_filter(ev, filter))
        .map(|ev| CalibratedEvent {
            time_utc: ins_to_utc(ev.time_ins),
            channel: (ev.pha / ncf) as i16,
            anode_id: ev.anode,
            layer: AnodeSelection::anode_to_layer(ev.anode),
            detector_id: ev.unit_id,
        })
        .collect()
}

/// Write a **Format B** per-unit Level-2 event file (`lxp?level2.evn`: `LAXPC_MODE_EA` + `GTI`),
/// with full OGIP & barycentric-ready timing keywords and FITS `CHECKSUM`/`DATASUM`.
pub fn write_format_b_event_fits(
    path: &Path,
    events: &[CalibratedEvent],
    gti_utc: &[(f64, f64)],
    unit: LaxpcUnit,
    filter: &EventFilter,
    num_channels: usize,
    meta: &OgipMeta,
) -> Result<()> {
    let mut primary = Hdu::empty_primary();
    apply_ogip_meta(&mut primary.header, meta);
    primary.header.set_f64("MJDREF", 55197.0, "2010.0 UTC");
    primary.header.set_f64("TIMEDEL", 1.0e-5, "[s] Time resolution");
    primary.header.set_str("CLOCKCOR", "YES", "Time corrected to UTC");

    let times: Vec<f64> = events.iter().map(|e| e.time_utc).collect();
    let phas: Vec<i16> = events.iter().map(|e| e.channel).collect();
    let anodes: Vec<u8> = events.iter().map(|e| e.anode_id).collect();
    let dets: Vec<u8> = events.iter().map(|e| e.detector_id).collect();

    let mut ev_hdu = BinTableBuilder::new("LAXPC_MODE_EA")
        .add_f64("TIME", Some("s"), times)
        .add_i16("PHA", Some("chan"), phas)
        .with_last_col_range(0, (num_channels as i64) - 1)
        .add_u8("ANODE_ID", None, anodes)
        .add_u8("DETECTOR_ID", None, dets)
        .header_str("HDUCLASS", "OGIP", "format conforms to OGIP")
        .header_str("HDUCLAS1", "EVENTS", "event table extension")
        .header_str("HDUCLAS2", "ALL", "all selected events")
        .header_str("HDUCLAS3", "COUNT", "event counts")
        .header_str("HDUVERS", "1.2.1", "OGIP version")
        .header_str("DATAMODE", "EA", "LAXPC event analysis mode")
        .header_f64("TIMEDEL", 1.0e-5, "[s] Time resolution")
        .header_f64("TIMEPIXR", 0.0, "timestamp refers to start of pixel")
        .header_f64("MJDREF", 55197.0, "2010.0 UTC")
        .header_str("CLOCKCOR", "YES", "Time corrected to UTC")
        .header_str("FILTER", filter.anodes.filter_keyword(filter.single_only), "event filter")
        .build()
        .map_err(|e| LaxpcError::Fits(e.to_string()))?;
    apply_ogip_meta(&mut ev_hdu.header, meta);
    ev_hdu.header.set_str("DETNAM", unit.detnam(), "LAXPC detector unit");

    let starts: Vec<f64> = gti_utc.iter().map(|g| g.0).collect();
    let stops: Vec<f64> = gti_utc.iter().map(|g| g.1).collect();
    let mut gti_hdu = build_gti_hdu("GTI", &starts, &stops, meta)
        .map_err(|e| LaxpcError::Fits(e.to_string()))?;
    gti_hdu.header.set_f64("MJDREF", 55197.0, "2010.0 UTC");

    let mut hdus = vec![primary, ev_hdu, gti_hdu];
    write_fits(path, &mut hdus, true).map_err(|e| LaxpcError::Fits(e.to_string()))
}

/// Write a **Format A** combined multi-unit Level-2 event file (`level2.event.fits`) with columns
/// `TIME`, `Channel`, `Layer`, `LAXPC_No` plus `GTI` extension and complete HDU-0/HDU-1/GTI timing
/// keywords so `astrosat bary` runs natively without `preplaxbary`.
pub fn write_format_a_event_fits(
    path: &Path,
    combined_events: &[CalibratedEvent],
    gti_utc: &[(f64, f64)],
    meta: &OgipMeta,
) -> Result<()> {
    write_format_a_event_fits_with_channels(path, combined_events, gti_utc, 1024, meta)
}

/// Write a **Format A** combined multi-unit Level-2 event file with an explicit uniform channel
/// count (`TLMIN2 = 0`, `TLMAX2 = num_channels - 1`).
pub fn write_format_a_event_fits_with_channels(
    path: &Path,
    combined_events: &[CalibratedEvent],
    gti_utc: &[(f64, f64)],
    num_channels: usize,
    meta: &OgipMeta,
) -> Result<()> {
    let mut primary = Hdu::empty_primary();
    apply_ogip_meta(&mut primary.header, meta);
    primary.header.set_f64("MJDREF", 55197.0, "2010.0 UTC");
    primary.header.set_f64("TIMEDEL", 1.0e-5, "[s] Time resolution");
    primary.header.set_str("CLOCKCOR", "YES", "Time corrected to UTC");
    if let Some(ra) = meta.ra_obj {
        primary.header.set_f64("RA_NOM", ra, "[deg] Nominal RA");
    }
    if let Some(dec) = meta.dec_obj {
        primary.header.set_f64("DEC_NOM", dec, "[deg] Nominal Dec");
    }

    let times: Vec<f64> = combined_events.iter().map(|e| e.time_utc).collect();
    let chans: Vec<i16> = combined_events.iter().map(|e| e.channel).collect();
    let layers: Vec<u8> = combined_events.iter().map(|e| e.layer).collect();
    let units: Vec<u8> = combined_events.iter().map(|e| e.detector_id).collect();
    let tlmax = (num_channels.max(1) as i64) - 1;

    let mut ev_hdu = BinTableBuilder::new("LAXPC_EVENTS")
        .add_f64("TIME", Some("s"), times)
        .add_i16("Channel", Some("chan"), chans)
        .with_last_col_range(0, tlmax)
        .add_u8("Layer", None, layers)
        .add_u8("LAXPC_No", None, units)
        .header_str("HDUCLASS", "OGIP", "format conforms to OGIP")
        .header_str("HDUCLAS1", "EVENTS", "combined multi-unit LAXPC event list")
        .header_str("HDUCLAS2", "ALL", "all selected events")
        .header_str("DATAMODE", "EA", "LAXPC event analysis mode")
        .header_f64("TIMEDEL", 1.0e-5, "[s] Time resolution")
        .header_f64("MJDREF", 55197.0, "2010.0 UTC")
        .header_str("CLOCKCOR", "YES", "Time corrected to UTC")
        .build()
        .map_err(|e| LaxpcError::Fits(e.to_string()))?;
    apply_ogip_meta(&mut ev_hdu.header, meta);
    if let Some(ra) = meta.ra_obj {
        ev_hdu.header.set_f64("RA_NOM", ra, "[deg] Nominal RA");
    }
    if let Some(dec) = meta.dec_obj {
        ev_hdu.header.set_f64("DEC_NOM", dec, "[deg] Nominal Dec");
    }

    let starts: Vec<f64> = gti_utc.iter().map(|g| g.0).collect();
    let stops: Vec<f64> = gti_utc.iter().map(|g| g.1).collect();
    let mut gti_hdu = build_gti_hdu("GTI", &starts, &stops, meta)
        .map_err(|e| LaxpcError::Fits(e.to_string()))?;
    gti_hdu.header.set_f64("MJDREF", 55197.0, "2010.0 UTC");

    let mut hdus = vec![primary, ev_hdu, gti_hdu];
    write_fits(path, &mut hdus, true).map_err(|e| LaxpcError::Fits(e.to_string()))
}

/// Write an OGIP 1.2.0 `SPECTRUM` FITS file (`lxp?level2.spec` or background `.spec`).
#[allow(clippy::too_many_arguments)]
pub fn write_spectrum_fits(
    path: &Path,
    unit: LaxpcUnit,
    rates: &[f64],
    stat_err: &[f64],
    exposure_sec: f64,
    deadtime_sec: f64,
    filter: &EventFilter,
    respfile: &str,
    backfile: &str,
    meta: &OgipMeta,
) -> Result<()> {
    let nch = rates.len();
    let channels: Vec<i16> = (0..nch).map(|c| c as i16).collect();

    let mut primary = Hdu::empty_primary();
    apply_ogip_meta(&mut primary.header, meta);
    primary.header.set_f64("MJDREF", 55197.0, "2010.0 UTC");

    let mut spec_hdu = BinTableBuilder::new("SPECTRUM")
        .add_i16("CHANNEL", None, channels)
        .with_last_col_range(0, (nch as i64) - 1)
        .add_f64("RATE", Some("count/s"), rates.to_vec())
        .add_f64("STAT_ERR", Some("count/s"), stat_err.to_vec())
        .header_str("HDUCLASS", "OGIP", "format conforms to OGIP")
        .header_str("HDUCLAS1", "SPECTRUM", "PHA dataset (OGIP-92-007)")
        .header_str("HDUCLAS2", "TOTAL", "total source+background rate")
        .header_str("HDUCLAS3", "RATE", "data stored as count rate")
        .header_str("HDUVERS", "1.2.0", "OGIP spectrum version")
        .header_str("CHANTYPE", "PHA", "channel type")
        .header_i64("DETCHANS", nch as i64, "total number of detector channels")
        .header_bool("POISSERR", false, "errors given in STAT_ERR column")
        .header_f64("AREASCAL", 1.0, "area scaling factor")
        .header_f64("BACKSCAL", 1.0, "background scaling factor")
        .header_f64("CORRSCAL", 0.0, "correction scaling factor")
        .header_i64("SYS_ERR", 0, "no systematic error specified")
        .header_i64("GROUPING", 0, "no grouping data specified")
        .header_i64("QUALITY", 0, "no data quality specified")
        .header_str("BACKFILE", backfile, "associated background spectrum")
        .header_str("CORRFILE", "NONE", "associated correction file")
        .header_str("RESPFILE", respfile, "associated response matrix file")
        .header_str("ANCRFILE", "NONE", "associated ancillary response file")
        .header_str("XFLT0001", "NONE", "XSPEC selection filter")
        .header_f64("DEADTIME", deadtime_sec, "[s] deadtime constant applied")
        .header_str("FILTER", filter.anodes.filter_keyword(filter.single_only), "anode/event filter")
        .build()
        .map_err(|e| LaxpcError::Fits(e.to_string()))?;

    let mut spec_meta = meta.clone();
    spec_meta.exposure = exposure_sec;
    spec_meta.detnam = Some(unit.detnam().to_string());
    apply_ogip_meta(&mut spec_hdu.header, &spec_meta);
    spec_hdu.header.set_f64("MJDREF", 55197.0, "2010.0 UTC");

    let mut hdus = vec![primary, spec_hdu];
    write_fits(path, &mut hdus, true).map_err(|e| LaxpcError::Fits(e.to_string()))
}

/// Write an OGIP `RATE` light curve FITS file (`lxp?level2.lcurv` or `.lcbk`) with attached `GTI` HDU.
#[allow(clippy::too_many_arguments)]
pub fn write_lightcurve_fits(
    path: &Path,
    lc: &LaxpcLightCurve,
    rates: &[f64],
    errors: &[f64],
    gti_utc: &[(f64, f64)],
    filter: &EventFilter,
    deadtime_sec: f64,
    backapp: bool,
    meta: &OgipMeta,
) -> Result<()> {
    let mut primary = Hdu::empty_primary();
    apply_ogip_meta(&mut primary.header, meta);
    primary.header.set_f64("MJDREF", 55197.0, "2010.0 UTC");

    let mjd_obs = 55197.0 + meta.tstart / 86_400.0;
    let chanrang = format!("{} - {}", filter.pha_min, filter.pha_max);

    let mut lc_hdu = BinTableBuilder::new("RATE")
        .add_f64("TIME", Some("s"), lc.time_utc.clone())
        .add_f64("RATE", Some("count/s"), rates.to_vec())
        .add_f64("ERROR", Some("count/s"), errors.to_vec())
        .add_f64("FRACEXP", None, lc.fracexp.clone())
        .header_str("HDUCLASS", "OGIP", "format conforms to OGIP")
        .header_str("HDUCLAS1", "LIGHTCURVE", "dataset contains a light curve")
        .header_str("HDUCLAS2", "TOTAL", "total counts/rate")
        .header_str("HDUCLAS3", "RATE", "data stored as count rate")
        .header_str("HDUVERS", "1.2.1", "OGIP version")
        .header_f64("TIMEDEL", lc.tbin, "[s] time bin size")
        .header_f64("TIMEPIXR", 0.5, "time stamp refers to bin center")
        .header_f64("MJDREF", 55197.0, "2010.0 UTC")
        .header_f64("MJD-OBS", mjd_obs, "MJD start time")
        .header_f64("DEADTIME", deadtime_sec, "[s] deadtime constant applied")
        .header_str("CLOCKCOR", "YES", "Time corrected to UTC")
        .header_bool("BACKAPP", backapp, "background correction applied")
        .header_bool("DEADAPP", true, "deadtime correction applied")
        .header_bool("VIGNAPP", true, "vignetting correction applied")
        .header_str("CHANRANG", &chanrang, "1024-scale channel range")
        .header_str("FILTER", filter.anodes.filter_keyword(filter.single_only), "anode/event filter")
        .build()
        .map_err(|e| LaxpcError::Fits(e.to_string()))?;

    let mut lc_meta = meta.clone();
    lc_meta.detnam = Some(lc.unit.detnam().to_string());
    apply_ogip_meta(&mut lc_hdu.header, &lc_meta);

    let starts: Vec<f64> = gti_utc.iter().map(|g| g.0).collect();
    let stops: Vec<f64> = gti_utc.iter().map(|g| g.1).collect();
    let mut gti_hdu = build_gti_hdu("GTI", &starts, &stops, &lc_meta)
        .map_err(|e| LaxpcError::Fits(e.to_string()))?;
    gti_hdu.header.set_f64("MJDREF", 55197.0, "2010.0 UTC");

    let mut hdus = vec![primary, lc_hdu, gti_hdu];
    write_fits(path, &mut hdus, true).map_err(|e| LaxpcError::Fits(e.to_string()))
}

/// Write a standalone Level-2 GTI FITS file (`lxp?level2.gtis`).
pub fn write_gti_fits(path: &Path, gti_utc: &[(f64, f64)], meta: &OgipMeta) -> Result<()> {
    let mut primary = Hdu::empty_primary();
    apply_ogip_meta(&mut primary.header, meta);
    primary.header.set_f64("MJDREF", 55197.0, "2010.0 UTC");

    let starts: Vec<f64> = gti_utc.iter().map(|g| g.0).collect();
    let stops: Vec<f64> = gti_utc.iter().map(|g| g.1).collect();
    let mut gti_hdu = build_gti_hdu("GTI", &starts, &stops, meta)
        .map_err(|e| LaxpcError::Fits(e.to_string()))?;
    gti_hdu.header.set_f64("MJDREF", 55197.0, "2010.0 UTC");

    let mut hdus = vec![primary, gti_hdu];
    write_fits(path, &mut hdus, true).map_err(|e| LaxpcError::Fits(e.to_string()))
}

/// Write an AstroSat `ORBIT` FITS file (`AstroSat.orb`) at 1-second cadence interpolated from MKF,
/// matching `laxpcl1.f:orbit` (and using double precision columns compatible with `astrosat-core::orbit`).
pub fn write_orbit_fits(
    path: &Path,
    time_utc: &[f64],
    pos_km: &[[f64; 3]],
    vel_kms: &[[f64; 3]],
    meta: &OgipMeta,
) -> Result<()> {
    let mut primary = Hdu::empty_primary();
    apply_ogip_meta(&mut primary.header, meta);
    primary.header.set_f64("MJDREF", 55197.0, "2010.0 UTC");

    let mut orb_hdu = build_orbit_hdu(time_utc, pos_km, vel_kms, meta)
        .map_err(|e| LaxpcError::Fits(e.to_string()))?;
    orb_hdu.header.set_f64("MJDREF", 55197.0, "2010.0 UTC");
    orb_hdu.header.set_str("MISSION", "ASTROSAT", "mission name");
    orb_hdu.header.set_f64("TIMEDEL", 1.0, "[s] orbit cadence");

    let mut hdus = vec![primary, orb_hdu];
    write_fits(path, &mut hdus, true).map_err(|e| LaxpcError::Fits(e.to_string()))
}

/// Write ASCII TIFR-compatible `.pha` spectrum table (channel + A1..A10 + A1-A7 sum).
pub fn write_ascii_pha(path: &Path, spec: &LaxpcSpectrum, mission_day: i32, single_only: bool) -> Result<()> {
    let file = File::create(path)?;
    let mut w = BufWriter::new(file);
    let iev = if single_only { 1 } else { 2 };
    writeln!(
        w,
        "#  Exposure ={:15.6E}  {}  Channels:{:5}{:5}  day:{:5}  rate:{:15.6E}  event:{:2}",
        spec.exposure_sec,
        spec.unit.prefix(),
        0,
        spec.num_channels.saturating_sub(1),
        mission_day,
        spec.total_corrected_rate_cps,
        iev
    )?;
    for (ch, row) in spec.anode_rates.iter().enumerate() {
        write!(w, "{ch:5}")?;
        for &val in row {
            write!(w, "{val:12.4E}")?;
        }
        writeln!(w)?;
    }
    w.flush()?;
    Ok(())
}

/// Write ASCII `.lc` light curve table (`time_ins, time_utc, rate_cps, error_cps, fracexp`).
pub fn write_ascii_lc(path: &Path, lc: &LaxpcLightCurve) -> Result<()> {
    let file = File::create(path)?;
    let mut w = BufWriter::new(file);
    writeln!(
        w,
        "# {} tbin={:.6} nbins={} mean_rate={:.4} cps",
        lc.unit.prefix(),
        lc.tbin,
        lc.time_utc.len(),
        lc.mean_rate_cps
    )?;
    for i in 0..lc.time_utc.len() {
        let bk = lc
            .bkg_rate_cps
            .as_ref()
            .and_then(|v| v.get(i))
            .copied()
            .unwrap_or(0.0);
        writeln!(
            w,
            "{:20.12E} {:20.12E} {:13.5E} {:13.5E} {:8.5} {:13.5E}",
            lc.time_ins[i], lc.time_utc[i], lc.rate_cps[i], lc.error_cps[i], lc.fracexp[i], bk
        )?;
    }
    w.flush()?;
    Ok(())
}

/// Fast binary-search membership check in sorted non-overlapping `[start, stop]` intervals.
#[inline]
pub fn in_sorted_intervals(intervals: &[(f64, f64)], t: f64) -> bool {
    if intervals.is_empty() {
        return false;
    }
    let idx = intervals.partition_point(|&(_, stop)| stop <= t);
    idx < intervals.len() && t > intervals[idx].0 && t < intervals[idx].1
}

/// Check whether `[t_lo, t_hi]` lies completely within any single GTI interval (`laxpcl1.f:1632-1636`).
#[inline]
fn bin_inside_gti(intervals: &[(f64, f64)], t_lo: f64, t_hi: f64) -> bool {
    if intervals.is_empty() {
        return false;
    }
    let idx = intervals.partition_point(|&(_, stop)| stop < t_hi);
    idx < intervals.len() && t_lo >= intervals[idx].0 && t_hi <= intervals[idx].1
}
