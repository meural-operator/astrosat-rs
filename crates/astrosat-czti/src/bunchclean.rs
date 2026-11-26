//! Stage 1b: `cztbunchclean` — Cosmic-ray bunch identification, post-bunch deadtime,
//! livetime binning, and per-pixel fractional exposure computation.

use serde::{Deserialize, Serialize};

use crate::types::{BunchTable, QuadEvents, QuadLivetime, PIXELS_PER_QUAD};

/// Configuration parameters for `cztbunchclean` (matching CZTI Pipeline v3.0.1 defaults).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BunchCleanConfig {
    /// Maximum inter-event time (microseconds) to belong to a bunch (default: `30` us).
    pub bunch_def_time_us: i64,
    /// Threshold bunch length separating small (module-level) vs heavy (quadrant-level) post-bunch
    /// deadtime (default: `20`).
    pub bunch_length_thresh: usize,
    /// Post-bunch deadtime (seconds) applied to all modules for small bunches (default: `0.0`).
    pub skip_t1: f64,
    /// Post-bunch deadtime (seconds) applied to dominant module for small bunches (default: `1e-3`).
    pub skip_t2: f64,
    /// Post-bunch deadtime (seconds) applied quadrant-wide for heavy bunches (default: `1e-3`).
    pub skip_t3: f64,
    /// Time bin size (seconds) for output livetime curves (default: `1.0`).
    pub livetime_binsize: f64,
}

impl Default for BunchCleanConfig {
    fn default() -> Self {
        Self {
            bunch_def_time_us: 30,
            bunch_length_thresh: 20,
            skip_t1: 0.0,
            skip_t2: 1e-3,
            skip_t3: 1e-3,
            livetime_binsize: 1.0,
        }
    }
}

impl BunchCleanConfig {
    /// Preset for Generalized v3.0 noise cleaning (`skipT1=skipT2=skipT3=0` -> `NCOMPAT=1`).
    pub fn generalized_compat() -> Self {
        Self {
            skip_t1: 0.0,
            skip_t2: 0.0,
            skip_t3: 0.0,
            ..Self::default()
        }
    }

    /// True if `skipT1 == 0 && skipT2 == 0 && skipT3 == 0` (`NCOMPAT=1` in FITS header).
    pub fn is_generalized_compat(&self) -> bool {
        self.skip_t1 == 0.0 && self.skip_t2 == 0.0 && self.skip_t3 == 0.0
    }
}

/// Result of running `bunchclean` on a single quadrant.
#[derive(Debug, Clone)]
pub struct QuadBunchCleanResult {
    pub events: QuadEvents,
    pub out_bunches: BunchTable,
    pub livetime: QuadLivetime,
    /// Fractional exposure (`0.0..=1.0`) for each of the 4096 pixels (`det_id*256 + pix_id`).
    pub pixel_exposure: Vec<f64>,
    pub ncompat: i32,
    pub events_before: usize,
    pub events_after_stage1: usize,
    pub events_after_stage2: usize,
}

/// Subtract or add deadtime `fraction * overlap` across livetime bins (`utils.cpp::updatelivetime`).
#[allow(clippy::too_many_arguments)]
pub fn update_livetime(
    gtstart: f64,
    gtstop: f64,
    livetime: &mut [f64],
    tstarti: i64,
    fraction: f64,
    add: bool,
    binsize: f64,
) {
    if gtstart == 0.0 || gtstop == 0.0 || gtstop <= gtstart || binsize <= 0.0 || livetime.is_empty()
    {
        return;
    }
    let start_bin = ((gtstart - tstarti as f64) / binsize).floor() as isize;
    let stop_bin = ((gtstop - tstarti as f64) / binsize).floor() as isize;
    let n = livetime.len() as isize;
    if start_bin < 0 || stop_bin < 0 || start_bin >= n || stop_bin >= n {
        return;
    }

    let sign = if add { 1.0 } else { -1.0 };
    if start_bin == stop_bin {
        livetime[start_bin as usize] += sign * (gtstop - gtstart) * fraction;
    } else {
        let start_edge =
            (((gtstart - tstarti as f64) / binsize).ceil() * binsize) - (gtstart - tstarti as f64);
        let stop_edge =
            (gtstop - tstarti as f64) - (((gtstop - tstarti as f64) / binsize).floor() * binsize);
        livetime[start_bin as usize] += sign * start_edge * fraction;
        livetime[stop_bin as usize] += sign * stop_edge * fraction;
        for bin in (start_bin + 1)..stop_bin {
            livetime[bin as usize] += sign * binsize * fraction;
        }
    }
}

/// Update livetime weighted by prior `in_livetime` (`utils.cpp::re_updatelivetime`, used by `pixclean`).
#[allow(clippy::too_many_arguments)]
pub fn re_update_livetime(
    gtstart: f64,
    gtstop: f64,
    in_livetime: &[f64],
    out_livetime: &mut [f64],
    tstarti: i64,
    fraction: f64,
    add: bool,
    binsize: f64,
) {
    if gtstart == 0.0
        || gtstop == 0.0
        || gtstop <= gtstart
        || binsize <= 0.0
        || out_livetime.is_empty()
    {
        return;
    }
    let start_bin = ((gtstart - tstarti as f64) / binsize).floor() as isize;
    let stop_bin = ((gtstop - tstarti as f64) / binsize).floor() as isize;
    let n = out_livetime.len() as isize;
    if start_bin < 0 || stop_bin < 0 || start_bin >= n || stop_bin >= n {
        return;
    }
    let sign = if add { 1.0 } else { -1.0 };
    let sb = start_bin as usize;
    let eb = stop_bin as usize;
    if sb == eb {
        out_livetime[sb] += sign * (gtstop - gtstart) * fraction * (in_livetime[sb] / binsize);
    } else {
        let start_edge =
            (((gtstart - tstarti as f64) / binsize).ceil() * binsize) - (gtstart - tstarti as f64);
        let stop_edge =
            (gtstop - tstarti as f64) - (((gtstop - tstarti as f64) / binsize).floor() * binsize);
        out_livetime[sb] += sign * start_edge * fraction * (in_livetime[sb] / binsize);
        out_livetime[eb] += sign * stop_edge * fraction * (in_livetime[eb] / binsize);
        for bin in (sb + 1)..eb {
            out_livetime[bin] += sign * binsize * fraction * (in_livetime[bin] / binsize);
            if out_livetime[bin] < 0.0 {
                out_livetime[bin] = 0.0;
            }
        }
    }
    if out_livetime[sb] < 0.0 {
        out_livetime[sb] = 0.0;
    }
    if out_livetime[eb] < 0.0 {
        out_livetime[eb] = 0.0;
    }
}

/// Apply GTI intervals `(start, stop)` to a livetime array (`utils.cpp::updatelivetime_gti`).
pub fn apply_gti_to_livetime(
    gti_intervals: &[(f64, f64)],
    time_array: &[f64],
    livetime: &mut [f64],
    binsize: f64,
) {
    if time_array.is_empty() || livetime.is_empty() || binsize <= 0.0 {
        return;
    }
    let ntbins = livetime.len();
    let tstart = time_array[0].trunc();
    let tstop = time_array[ntbins - 1].trunc() + 1.0;

    let mut btstart = tstart;
    let nrows = gti_intervals.len();
    #[allow(clippy::needless_range_loop)]
    for i in 0..=nrows {
        let btstop = if i < nrows {
            gti_intervals[i].0
        } else {
            tstop
        };

        let btstartbin = ((btstart - tstart) / binsize).floor() as isize;
        let btstopbin = ((btstop - tstart) / binsize).floor() as isize;

        if btstartbin >= 0
            && btstopbin >= 0
            && (btstartbin as usize) < ntbins
            && (btstopbin as usize) < ntbins
            && btstop > btstart
        {
            let sb = btstartbin as usize;
            let eb = btstopbin as usize;
            if sb == eb {
                livetime[sb] -= btstop - btstart;
            } else {
                for lv in livetime.iter_mut().take(eb).skip(sb + 1) {
                    *lv = 0.0;
                }
                livetime[sb] -= ((sb as f64 + 1.0) * binsize) - (btstart - tstart);
                livetime[eb] -= (btstop - tstart) - ((eb as f64) * binsize);
            }
            if livetime[sb] < 0.0 {
                livetime[sb] = 0.0;
            }
            if livetime[eb] < 0.0 {
                livetime[eb] = 0.0;
            }
        }
        if i < nrows {
            btstart = gti_intervals[i].1;
        }
    }
}

/// Clean cosmic-ray bunches for one quadrant (exact port of `cztbunchcleanv2.cpp`).
///
/// `gti_intervals` is the quadrant's M0 GTI (`Q{qid}_GTI`). If empty, falls back to `[t_first, t_last]`.
pub fn bunchclean_quad(
    events: &QuadEvents,
    in_bunches: &BunchTable,
    gti_intervals: &[(f64, f64)],
    cfg: &BunchCleanConfig,
) -> QuadBunchCleanResult {
    let qid = events.quad_id;
    let nrows = events.len();
    let ncompat = if cfg.is_generalized_compat() { 1 } else { 0 };

    if nrows == 0 {
        return QuadBunchCleanResult {
            events: QuadEvents::new(qid),
            out_bunches: BunchTable::new(qid),
            livetime: QuadLivetime {
                quad_id: qid,
                binsize: cfg.livetime_binsize,
                time: Vec::new(),
                fracexp: Vec::new(),
            },
            pixel_exposure: vec![0.0; PIXELS_PER_QUAD],
            ncompat,
            events_before: 0,
            events_after_stage1: 0,
            events_after_stage2: 0,
        };
    }

    let (gti_first, gti_last, allquad_exp_time) = if !gti_intervals.is_empty() {
        let exp: f64 = gti_intervals.iter().map(|&(a, b)| (b - a).max(0.0)).sum();
        (gti_intervals[0].0, gti_intervals[gti_intervals.len() - 1].1, exp.max(1e-9))
    } else {
        let a = events.time[0];
        let b = events.time[nrows - 1];
        (a, b, (b - a).max(1.0))
    };

    let lv_binsize = if cfg.livetime_binsize.is_finite() && cfg.livetime_binsize > 0.0 {
        cfg.livetime_binsize
    } else {
        1.0
    };
    let tstarti = if gti_first.is_finite() { gti_first as i64 } else { 0 };
    let tstopi = if gti_last.is_finite() { gti_last as i64 } else { tstarti };
    const MAX_LV_BINS: usize = 10_000_000;
    let raw_bins = (((tstopi.saturating_sub(tstarti).max(0)) as f64) / lv_binsize) + 1.0;
    let num_timebins = if raw_bins.is_finite() && raw_bins >= 1.0 {
        (raw_bins as usize).clamp(1, MAX_LV_BINS)
    } else {
        1
    };
    let mut time_array = Vec::with_capacity(num_timebins);
    let mut livetime = vec![lv_binsize; num_timebins];
    for j in 0..num_timebins {
        time_array.push(tstarti as f64 + (j as f64) * lv_binsize + 0.5 * lv_binsize);
    }

    let mut pixel_loss = vec![0.0f64; PIXELS_PER_QUAD];
    let mut evt_index = vec![0usize; nrows];
    let mut evt_flag = vec![0u8; nrows];

    let has_onboard_bunches = !in_bunches.is_empty();
    let nrows_bunch = in_bunches.len();
    let mut bunch_flag = vec![0u8; nrows_bunch];
    let mut detected_bunches = BunchTable::new(qid);

    // Stage 1: Identify bunches & exclude first 400 ms of every 100 s CZT clock tick
    let mut nbunch_evt: usize = 1;
    let mut bunch_tstart = 0.0f64;
    let mut bunch_time_dfs: u32 = 0;
    let mut row_index: usize = 0;

    let mut i = 0usize;
    while i < nrows.saturating_sub(1) {
        let cztt = events.cztseccnt[i];
        let czt_int = cztt as i64;
        let tmod100 = czt_int.rem_euclid(100);
        let in_100s_tick = tmod100 == 0 && (cztt - (czt_int as f64)) <= 0.400;

        if !in_100s_tick {
            let dt = events.time[i + 1] - events.time[i];
            let delta_us: i64 = if dt < 10.0 {
                (dt * 1e6) as i64
            } else {
                10_000_000
            };

            if !has_onboard_bunches {
                if delta_us <= cfg.bunch_def_time_us && i + 2 != nrows {
                    nbunch_evt += 1;
                    if nbunch_evt == 2 {
                        bunch_tstart = events.time[i];
                        bunch_time_dfs = ((delta_us as f64) / 20.0).round().max(0.0) as u32;
                    }
                    if nbunch_evt == 3 {
                        detected_bunches.time.push(events.time[i]);
                        detected_bunches.row_num.push((i - 1) as i32);
                        detected_bunches.time_dfs.push(bunch_time_dfs.min(255) as u8);
                        detected_bunches.detid_fevt.push(events.detid[i - 1]);
                        detected_bunches.pixid_fevt.push(events.pixid[i - 1]);
                        detected_bunches.detid_sevt.push(events.detid[i]);
                        detected_bunches.pixid_sevt.push(events.pixid[i]);
                    }
                    if nbunch_evt == 4 {
                        detected_bunches.detid1.push(events.detid[i]);
                    }
                    if nbunch_evt == 5 {
                        detected_bunches.detid2.push(events.detid[i]);
                    }
                    if nbunch_evt == 6 {
                        detected_bunches.detid3.push(events.detid[i]);
                    }
                    if nbunch_evt == 7 {
                        detected_bunches.detid4.push(events.detid[i]);
                    }
                } else {
                    let bunch_tstop = events.time[i];
                    if nbunch_evt >= 3 {
                        let dsl = (((bunch_tstop - bunch_tstart) * 1e6 / 20.0).round() as i64
                            - bunch_time_dfs as i64)
                            .clamp(0, 255) as u8;
                        detected_bunches.time_dsl.push(dsl);
                        detected_bunches
                            .num_event
                            .push((nbunch_evt - 1).min(255) as u8);
                        detected_bunches.detid_levt.push(events.detid[i]);
                        detected_bunches.pixid_levt.push(events.pixid[i]);
                        if nbunch_evt < 7 {
                            detected_bunches.detid4.push(99);
                        }
                        if nbunch_evt < 6 {
                            detected_bunches.detid3.push(99);
                        }
                        if nbunch_evt < 5 {
                            detected_bunches.detid2.push(99);
                        }
                        if nbunch_evt < 4 {
                            detected_bunches.detid1.push(99);
                        }
                    }

                    let start_j = i + 1 - nbunch_evt;
                    for slot in evt_index.iter_mut().take(i + 1).skip(start_j) {
                        *slot = nbunch_evt;
                    }
                    nbunch_evt = 1;
                }
            } else if row_index < nrows_bunch && (i as i32) == in_bunches.row_num[row_index] {
                let nb = (in_bunches.num_event[row_index] as usize) + 1;
                row_index += 1;
                let end_j = (i + 3).min(nrows);
                for slot in evt_index.iter_mut().take(end_j).skip(i) {
                    *slot = nb;
                }
                i += 2; // loop adds +1 below -> advances by 3
                nbunch_evt = 1;
            } else if delta_us <= cfg.bunch_def_time_us && i + 2 != nrows {
                nbunch_evt += 1;
            } else {
                let start_j = i + 1 - nbunch_evt;
                for slot in evt_index.iter_mut().take(i + 1).skip(start_j) {
                    *slot = nbunch_evt;
                }
                nbunch_evt = 1;
            }
        } else if has_onboard_bunches
            && row_index < nrows_bunch
            && (i as i32) == in_bunches.row_num[row_index]
        {
            bunch_flag[row_index] = 1;
            row_index += 1;
            let end_j = (i + 3).min(nrows);
            for slot in evt_index.iter_mut().take(end_j).skip(i) {
                *slot = 0;
            }
            i += 2;
        } else {
            evt_index[i] = 0;
        }
        i += 1;
    }

    // Handle single-event edge case
    if nrows == 1 {
        let cztt = events.cztseccnt[0];
        let czt_int = cztt as i64;
        if !(czt_int.rem_euclid(100) == 0 && (cztt - czt_int as f64) <= 0.400) {
            evt_index[0] = 1;
        }
    }

    let events_after_stage1 = evt_index
        .iter()
        .filter(|&&idx| idx == 1 || idx == 2)
        .count();

    // Stage 2: Post-bunch deadtime exclusion, livetime & pixel exposure loss update
    let mut prev_bunch_tstop = 0.0f64;
    let mut prev_bunch_module: i32 = -1;
    let mut prev_bunch_length: usize = 0;

    let mut prev2_btstop = 0.0f64;
    let mut prev2_btstartmod = 0.0f64;
    let mut prev2_btstopmod = 0.0f64;
    let mut prev2_bunch_module: i32 = -1;
    let mut bunch_count: usize = 0;

    let limit = if nrows == 1 { 1 } else { nrows - 1 };
    let mut idx = 0usize;
    while idx < limit {
        let eidx = evt_index[idx];
        if eidx == 1 {
            let delta_tstop = events.time[idx] - prev_bunch_tstop;
            let keep = ((delta_tstop > cfg.skip_t1)
                && (prev_bunch_length < cfg.bunch_length_thresh)
                && ((events.detid[idx] as i32 != prev_bunch_module)
                    || (delta_tstop > cfg.skip_t2)))
                || (delta_tstop > cfg.skip_t3);
            if keep {
                evt_flag[idx] = 1;
            }
            idx += 1;
        } else if eidx == 2 {
            let delta_tstop = events.time[idx] - prev_bunch_tstop;
            let keep = ((delta_tstop > cfg.skip_t1)
                && (prev_bunch_length < cfg.bunch_length_thresh)
                && ((events.detid[idx] as i32 != prev_bunch_module)
                    || (delta_tstop > cfg.skip_t2)))
                || (delta_tstop > cfg.skip_t3);
            if keep {
                evt_flag[idx] = 2;
                if idx + 1 < nrows {
                    evt_flag[idx + 1] = 2;
                }
            }
            idx += 2;
        } else if eidx > 2 {
            evt_flag[idx] = 10;
            let prev_bunch_tstart = events.time[idx];
            prev_bunch_tstop = if !has_onboard_bunches {
                events.time[(idx + eidx - 1).min(nrows - 1)]
            } else {
                events.time[(idx + 2).min(nrows - 1)]
            };
            prev_bunch_length = eidx;

            // Dominant module histogram
            let mut modhist = [0usize; 16];
            if !has_onboard_bunches {
                for &det in events.detid.iter().take((idx + eidx).min(nrows)).skip(idx) {
                    let d = det as usize;
                    if d < 16 {
                        modhist[d] += 1;
                    }
                }
            } else {
                for &det in events.detid.iter().take((idx + 3).min(nrows)).skip(idx) {
                    let d = det as usize;
                    if d < 16 {
                        modhist[d] += 1;
                    }
                }
                if prev_bunch_length > 3 && bunch_count < nrows_bunch {
                    let extra_dets = [
                        in_bunches.detid1[bunch_count],
                        in_bunches.detid2[bunch_count],
                        in_bunches.detid3[bunch_count],
                        in_bunches.detid4[bunch_count],
                    ];
                    for (k, &d) in extra_dets.iter().enumerate() {
                        if prev_bunch_length - 3 >= k && (d as usize) < 16 {
                            modhist[d as usize] += 1;
                        }
                    }
                }
                bunch_count += 1;
            }

            // Find module with maximum counts
            prev_bunch_module = -1;
            let mut maxhist = 0usize;
            for (k, &cnt) in modhist.iter().enumerate() {
                if cnt > maxhist {
                    maxhist = cnt;
                    prev_bunch_module = k as i32;
                }
            }

            idx += if !has_onboard_bunches { eidx } else { 3 };

            // Compute deadtime intervals (quadrant-wide, dominant module, other modules)
            let mut btstart: f64;
            let mut btstop: f64;
            let mut btstartmod: f64;
            let mut btstopmod: f64;
            let mut btstartomod = 0.0f64;
            let mut btstopomod = 0.0f64;

            if prev_bunch_length < cfg.bunch_length_thresh {
                btstart = prev_bunch_tstart;
                btstop = prev_bunch_tstop + cfg.skip_t1;
                if cfg.skip_t2 > cfg.skip_t1 {
                    btstartmod = prev_bunch_tstop + cfg.skip_t1;
                    btstopmod = prev_bunch_tstop + cfg.skip_t2;
                } else {
                    btstartmod = 0.0;
                    btstopmod = 0.0;
                }
            } else {
                btstart = prev_bunch_tstart;
                btstop = prev_bunch_tstop + cfg.skip_t1.max(cfg.skip_t3);
                btstartmod = 0.0;
                btstopmod = 0.0;
            }

            if btstopmod == 0.0 {
                if btstart < prev2_btstop {
                    btstart = prev2_btstop;
                }
                if prev2_btstopmod != 0.0 {
                    if btstart < prev2_btstop {
                        btstart = prev2_btstopmod;
                        btstartomod = prev2_btstop;
                        if btstop < prev2_btstopmod {
                            btstopomod = btstop;
                            btstop = prev2_btstopmod;
                        } else {
                            btstopomod = prev2_btstopmod;
                        }
                    } else if btstart < prev2_btstopmod {
                        btstart = prev2_btstopmod;
                        btstartomod = prev2_btstartmod;
                        if btstop < prev2_btstopmod {
                            btstopomod = btstop;
                            btstop = prev2_btstopmod;
                        } else {
                            btstopomod = prev2_btstopmod;
                        }
                    }
                }

                update_livetime(
                    btstart,
                    btstop,
                    &mut livetime,
                    tstarti,
                    1.0,
                    false,
                    lv_binsize,
                );
                update_livetime(
                    btstartomod,
                    btstopomod,
                    &mut livetime,
                    tstarti,
                    15.0 / 16.0,
                    false,
                    lv_binsize,
                );

                if btstop > btstart {
                    let dt = btstop - btstart;
                    for p in &mut pixel_loss {
                        *p += dt;
                    }
                }
                if btstopomod != 0.0 && btstopomod > btstartomod {
                    let dt = btstopomod - btstartomod;
                    for m in 0..16i32 {
                        if m != prev2_bunch_module {
                            let base = (m as usize) * 256;
                            for p in &mut pixel_loss[base..base + 256] {
                                *p += dt;
                            }
                        }
                    }
                }
            } else {
                if prev2_btstopmod == 0.0 {
                    if btstart < prev2_btstop {
                        btstart = prev2_btstop;
                        if btstop < prev2_btstop {
                            btstop = prev2_btstop;
                        }
                    }
                    if btstartmod < prev2_btstop {
                        btstartmod = prev2_btstop;
                        if btstopmod < prev2_btstop {
                            btstopmod = prev2_btstop;
                        }
                    }
                } else if btstart < prev2_btstop {
                    btstartomod = prev2_btstop;
                    btstopomod = btstop;
                    btstart = btstop;
                    if prev_bunch_module == prev2_bunch_module {
                        btstartmod = prev2_btstopmod;
                    }
                } else if btstart < prev2_btstopmod {
                    if btstartmod < prev2_btstopmod {
                        btstartomod = btstart;
                        btstopomod = btstop;
                        if prev_bunch_module == prev2_bunch_module {
                            btstartmod = prev2_btstopmod;
                        }
                        btstart = btstop;
                    } else {
                        btstartomod = btstart;
                        btstopomod = prev2_btstopmod;
                        btstart = prev2_btstopmod;
                    }
                }

                update_livetime(
                    btstart,
                    btstop,
                    &mut livetime,
                    tstarti,
                    1.0,
                    false,
                    lv_binsize,
                );
                update_livetime(
                    btstartmod,
                    btstopmod,
                    &mut livetime,
                    tstarti,
                    1.0 / 16.0,
                    false,
                    lv_binsize,
                );
                update_livetime(
                    btstartomod,
                    btstopomod,
                    &mut livetime,
                    tstarti,
                    15.0 / 16.0,
                    false,
                    lv_binsize,
                );

                if btstop > btstart {
                    let dt = btstop - btstart;
                    for p in &mut pixel_loss {
                        *p += dt;
                    }
                }
                if btstopmod > btstartmod && prev_bunch_module >= 0 {
                    let dt = btstopmod - btstartmod;
                    let base = (prev_bunch_module as usize) * 256;
                    for p in &mut pixel_loss[base..base + 256] {
                        *p += dt;
                    }
                }
                if btstopomod != 0.0 && btstopomod > btstartomod {
                    let dt = btstopomod - btstartomod;
                    for m in 0..16i32 {
                        if m != prev2_bunch_module {
                            let base = (m as usize) * 256;
                            for p in &mut pixel_loss[base..base + 256] {
                                *p += dt;
                            }
                        }
                    }
                }
            }

            prev2_btstop = btstop;
            prev2_btstartmod = btstartmod;
            prev2_btstopmod = btstopmod;
            prev2_bunch_module = prev_bunch_module;
        } else {
            evt_flag[idx] = 0;
            idx += 1;
        }
    }

    let keep_mask: Vec<bool> = evt_flag.iter().map(|&f| f == 1 || f == 2).collect();
    let events_after_stage2 = keep_mask.iter().filter(|&&b| b).count();
    let cleaned_events = events.filter_by_mask(&keep_mask);

    // Clamp livetime bins to [0, binsize]
    for lv in &mut livetime {
        *lv = lv.clamp(0.0, lv_binsize);
    }

    // Fractional pixel exposure
    let mut pixel_exposure = vec![0.0f64; PIXELS_PER_QUAD];
    for j in 0..PIXELS_PER_QUAD {
        pixel_exposure[j] = ((allquad_exp_time - pixel_loss[j]) / allquad_exp_time).clamp(0.0, 1.0);
    }

    // Build output bunch table
    let out_bunches = if !has_onboard_bunches {
        detected_bunches
    } else {
        let mut filtered = BunchTable::new(qid);
        for k in 0..nrows_bunch {
            if bunch_flag[k] == 0 {
                filtered.time.push(in_bunches.time[k]);
                filtered.time_dfs.push(in_bunches.time_dfs[k]);
                filtered.time_dsl.push(in_bunches.time_dsl[k]);
                filtered.num_event.push(in_bunches.num_event[k]);
                filtered.detid1.push(in_bunches.detid1[k]);
                filtered.detid2.push(in_bunches.detid2[k]);
                filtered.detid3.push(in_bunches.detid3[k]);
                filtered.detid4.push(in_bunches.detid4[k]);
                filtered.row_num.push(in_bunches.row_num[k]);
                if k < in_bunches.detid_fevt.len() {
                    filtered.detid_fevt.push(in_bunches.detid_fevt[k]);
                    filtered.pixid_fevt.push(in_bunches.pixid_fevt[k]);
                    filtered.detid_sevt.push(in_bunches.detid_sevt[k]);
                    filtered.pixid_sevt.push(in_bunches.pixid_sevt[k]);
                    filtered.detid_levt.push(in_bunches.detid_levt[k]);
                    filtered.pixid_levt.push(in_bunches.pixid_levt[k]);
                }
            }
        }
        filtered
    };

    QuadBunchCleanResult {
        events: cleaned_events,
        out_bunches,
        livetime: QuadLivetime {
            quad_id: qid,
            binsize: lv_binsize,
            time: time_array,
            fracexp: livetime,
        },
        pixel_exposure,
        ncompat,
        events_before: nrows,
        events_after_stage1,
        events_after_stage2,
    }
}
