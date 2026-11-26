//! Stage 2d (Generalized v3.0 algorithm): `cztnoiseclean` suite
//! (`cztnoisypixclean`, `cztsuperbunchclean`, `cztheavybunchclean`, `cztflickpixclean`, `czteventsep`).

use serde::{Deserialize, Serialize};

use crate::bunchclean::update_livetime;
use crate::types::{
    detid_pixid_to_detxy, pix_flag, BunchTable, QuadEvents, QuadLivetime, PIXELS_PER_QUAD, QUAD_DIM,
};

/// Configuration for the Generalized v3.0 noise cleaning pipeline (`config/*Threshold.txt`).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NoiseCleanConfig {
    // cztnoisypixclean
    pub dead_pix_threshold: u64,
    pub noisy_pix_sigma_thresh: f64,

    // cztsuperbunchclean
    pub heavy_bun_thresh_count_10ms: usize,
    pub heavy_bun_thresh_count_1ms: usize,
    pub super_bunch_size: usize,
    pub dph_time: f64,
    pub thresh_hot: f64,
    pub allowable_hot_thresh: f64,

    // cztheavybunchclean
    pub heavy_bunch_size: usize,
    /// Post-heavy-bunch window (seconds) for low-energy events below LLD + `module_lld_threshold` (default: `0.25` s; configurable e.g. `0.20` s).
    pub heavy_bunch_time: f64,
    /// Post-heavy-bunch window (seconds) for all events in the affected module (default: `0.005` s).
    pub heavy_bunch_time_det: f64,
    /// Energy margin above module LLD (keV) for heavy-bunch tail rejection (default: `10.0` keV).
    pub module_lld_threshold_kev: f32,
    /// Post-small-bunch window (seconds) for events in affected module (default: `0.001` s).
    pub small_bunch_time_det: f64,

    // cztflickpixclean
    pub flicker_time_bin_sec: f64,
    pub flicker_sigma: f64,
    pub flicker_num_times: usize,

    // czteventsep
    /// Single/double separation coincidence window in seconds (default: `30e-6` s = 30 us,
    /// matching `czteventsep.c:13 #define double_time 0.00003`).
    pub double_time_sec: f64,
    pub write_double_evt: bool,
}

impl Default for NoiseCleanConfig {
    fn default() -> Self {
        Self {
            dead_pix_threshold: 10,
            noisy_pix_sigma_thresh: 5.0,
            heavy_bun_thresh_count_10ms: 4,
            heavy_bun_thresh_count_1ms: 4,
            super_bunch_size: 100,
            dph_time: 0.1,
            thresh_hot: 0.707,
            allowable_hot_thresh: 3.0,
            heavy_bunch_size: 15,
            heavy_bunch_time: 0.25,
            heavy_bunch_time_det: 0.005,
            module_lld_threshold_kev: 10.0,
            small_bunch_time_det: 0.001,
            flicker_time_bin_sec: 100.0,
            flicker_sigma: 6.0,
            flicker_num_times: 2,
            double_time_sec: 30e-6,
            write_double_evt: false,
        }
    }
}

/// Merge consecutive onboard bunches that saturated at `num_event == 63` within `<= 30 us`
/// (`cztsuperbunchclean.c:405-439`, `cztheavybunchclean.c:720-740`).
///
/// Returns `(merged_end_time, merged_size, representative_detid)` tuples.
pub fn merge_split_bunches(bunches: &BunchTable) -> Vec<(f64, usize, u8)> {
    let nb = bunches.len();
    let mut out = Vec::with_capacity(nb);
    let mut i = 0usize;
    while i < nb {
        let mut size = (bunches.num_event[i] as usize) + 1;
        let mut end_t = bunches.time[i] + (bunches.time_dsl[i] as f64) * 20e-6;
        let rep_det = bunches
            .detid_fevt
            .get(i)
            .copied()
            .filter(|&d| d < 16)
            .or_else(|| Some(bunches.detid1[i]).filter(|&d| d < 16))
            .unwrap_or(0);

        if bunches.num_event[i] == 63 {
            let mut j = i + 1;
            while j < nb {
                let next_start = bunches.time[j] - (bunches.time_dfs[j] as f64) * 20e-6;
                if next_start - end_t <= 30e-6 {
                    end_t = bunches.time[j] + (bunches.time_dsl[j] as f64) * 20e-6;
                    size += (bunches.num_event[j] as usize) + 1;
                    i = j;
                    if bunches.num_event[j] != 63 {
                        break;
                    }
                    j += 1;
                } else {
                    break;
                }
            }
        }
        out.push((end_t, size, rep_det));
        i += 1;
    }
    out
}

/// Poisson burst threshold (`cztflickpixclean.c:995-1012::pixThresholdFix`): smallest integer `x >= 2`
/// such that `num_pix * exp(-mean) * mean^x / x! < 1.0`.
pub fn pix_threshold_fix(countrate: usize, num_pix: usize) -> usize {
    if countrate == 0 || num_pix == 0 {
        return 2;
    }
    let mean = (countrate as f64) / (num_pix as f64);
    let mut x = 2usize;
    let mut fact = 2.0f64;
    let neg_mean_exp = (-mean).exp();
    while x < 200 {
        let chance = neg_mean_exp * mean.powi(x as i32) / fact * (num_pix as f64);
        if !chance.is_finite() || chance < 1.0 {
            break;
        }
        x += 1;
        fact *= x as f64;
    }
    x
}

/// Output of `cztnoisypixclean` for one quadrant.
#[derive(Debug, Clone)]
pub struct NoisyPixCleanResult {
    pub events: QuadEvents,
    pub badpix_flags: Vec<u8>,
    pub pixel_exposure: Vec<f64>,
}

/// Step 1 of Generalized v3.0: `cztnoisypixclean` — pixel-area-normalized (`2.28^2` edge vs
/// `2.46^2` interior) iterative two-sided sigma-clipping performed separately for normal
/// (`flag == 0`) and spectroscopically bad / low-gain (`flag == 1`) pixels
/// (`cztnoisypixclean.c:432-578`).
pub fn noisypixclean_quad(
    events: &QuadEvents,
    caldb_badpix: &[u8],
    in_pixel_exposure: &[f64],
    cfg: &NoiseCleanConfig,
) -> NoisyPixCleanResult {
    let qid = events.quad_id;
    let mut flags = vec![pix_flag::GOOD; PIXELS_PER_QUAD];
    let c_bp = PIXELS_PER_QUAD.min(caldb_badpix.len());
    flags[..c_bp].copy_from_slice(&caldb_badpix[..c_bp]);
    let mut pix_exp = vec![1.0f64; PIXELS_PER_QUAD];
    let c_pe = PIXELS_PER_QUAD.min(in_pixel_exposure.len());
    pix_exp[..c_pe].copy_from_slice(&in_pixel_exposure[..c_pe]);
    let mut dph_raw = vec![0u64; PIXELS_PER_QUAD];

    for i in 0..events.len() {
        let idx = (events.detid[i] as usize) * 256 + (events.pixid[i] as usize);
        if idx < PIXELS_PER_QUAD {
            dph_raw[idx] += 1;
        }
    }

    // Mark pixels with `< dead_pix_threshold` counts as DEAD (4) (`cztnoisypixclean.c:432`)
    let total_counts: u64 = dph_raw.iter().sum();
    let dead_thresh = if total_counts >= (PIXELS_PER_QUAD as u64) * cfg.dead_pix_threshold {
        cfg.dead_pix_threshold
    } else {
        1
    };
    for i in 0..PIXELS_PER_QUAD {
        if dph_raw[i] < dead_thresh && flags[i] < pix_flag::NOISY {
            flags[i] = pix_flag::DEAD;
        }
    }

    // Normalize DPH by physical pixel area (`2.28*2.28` module edges vs `2.46*2.46` interior,
    // `cztnoisypixclean.c:448-461`)
    let mut dph = vec![0.0f64; PIXELS_PER_QUAD];
    for idx in 0..PIXELS_PER_QUAD {
        let det = (idx / 256) as u8;
        let pix = (idx % 256) as u8;
        let (dx, dy) = detid_pixid_to_detxy(qid, det, pix);
        let is_edge = (dx + 1) % 16 == 0 || dx % 16 == 0 || (dy + 1) % 16 == 0 || dy % 16 == 0;
        let area = if is_edge { 2.28 * 2.28 } else { 2.46 * 2.46 };
        dph[idx] = (dph_raw[idx] as f64) / area;
    }

    // Two-sided iterative sigma clipping separately for normal (class=0) and spec-bad (class=1)
    // (`cztnoisypixclean.c:499,509`)
    for target_class in [pix_flag::GOOD, pix_flag::SPEC_BAD] {
        for _ in 0..50 {
            let mut sum = 0.0f64;
            let mut count = 0usize;
            for i in 0..PIXELS_PER_QUAD {
                if flags[i] == target_class {
                    sum += dph[i];
                    count += 1;
                }
            }
            if count == 0 {
                break;
            }
            let mean = sum / (count as f64);
            let mut var = 0.0f64;
            for i in 0..PIXELS_PER_QUAD {
                if flags[i] == target_class {
                    let diff = dph[i] - mean;
                    var += diff * diff;
                }
            }
            let sd = (var / (count as f64)).sqrt();
            let hi_cut = mean + cfg.noisy_pix_sigma_thresh * sd;
            let lo_cut = mean - cfg.noisy_pix_sigma_thresh * sd;
            let mut changed = false;
            for i in 0..PIXELS_PER_QUAD {
                if flags[i] == target_class && (dph[i] > hi_cut || dph[i] < lo_cut) {
                    flags[i] = pix_flag::NOISY;
                    changed = true;
                }
            }
            if !changed {
                break;
            }
        }
    }

    // Set pixel exposure: `1.0` for `flag < 2`, `0.0` for `flag >= NOISY` (`cztnoisypixclean.c:574-578`)
    for i in 0..PIXELS_PER_QUAD {
        if flags[i] < pix_flag::FLICKERING {
            pix_exp[i] = 1.0;
        } else if flags[i] >= pix_flag::NOISY {
            pix_exp[i] = 0.0;
        }
    }

    let keep: Vec<bool> = (0..events.len())
        .map(|i| {
            let idx = (events.detid[i] as usize) * 256 + (events.pixid[i] as usize);
            idx < PIXELS_PER_QUAD && flags[idx] < pix_flag::NOISY
        })
        .collect();

    NoisyPixCleanResult {
        events: events.filter_by_mask(&keep),
        badpix_flags: flags,
        pixel_exposure: pix_exp,
    }
}

/// Output of `cztsuperbunchclean` for one quadrant.
#[derive(Debug, Clone)]
pub struct SuperBunchCleanResult {
    pub events: QuadEvents,
    pub livetime: QuadLivetime,
    pub bti_intervals: Vec<(f64, f64)>,
}

/// Step 2 of Generalized v3.0: `cztsuperbunchclean` — merges `bunsize == 63` split bunches,
/// detects super-bunches (`size >= super_bunch_size` or clustered heavy bunches in 1 ms / 10 ms),
/// evaluates weighted spatial hotness `magnitude_k * magnitude_m / d` with 3x3 neighbor guard over
/// the `dph_time` (100 ms) window, removes contaminated BTIs, and updates livetime
/// (`cztsuperbunchclean.c:405-682`).
pub fn superbunchclean_quad(
    events: &QuadEvents,
    bunches: &BunchTable,
    in_livetime: &QuadLivetime,
    cfg: &NoiseCleanConfig,
) -> SuperBunchCleanResult {
    let mut out_lv = in_livetime.clone();
    if events.is_empty() || bunches.is_empty() {
        return SuperBunchCleanResult {
            events: events.clone(),
            livetime: out_lv,
            bti_intervals: Vec::new(),
        };
    }

    // Merge split 63-event bunches (`cztsuperbunchclean.c:405-439`)
    let merged = merge_split_bunches(bunches);
    let nb = merged.len();
    let mut super_times: Vec<f64> = Vec::new();
    let mut heavy_times: Vec<f64> = Vec::new();

    for &(t0, bsize, _) in &merged {
        if bsize > cfg.heavy_bunch_size {
            heavy_times.push(t0);
        }
    }

    for i in 0..nb {
        let (t0, bsize, _) = merged[i];
        if bsize >= cfg.super_bunch_size {
            super_times.push(t0);
            continue;
        }
        let mut cnt_1ms = 0usize;
        let mut cnt_10ms = 0usize;
        for &(tj, sz, _) in &merged[i..nb] {
            let dt = tj - t0;
            if dt > 0.010 {
                break;
            }
            if sz >= cfg.heavy_bunch_size {
                cnt_10ms += 1;
                if dt <= 0.001 {
                    cnt_1ms += 1;
                }
            }
        }
        if cnt_1ms >= cfg.heavy_bun_thresh_count_1ms || cnt_10ms >= cfg.heavy_bun_thresh_count_10ms
        {
            super_times.push(t0);
        }
    }

    let mut btis: Vec<(f64, f64)> = Vec::new();
    let mut keep = vec![true; events.len()];

    // Unconditional 100 ms BTI for super-bunches
    for &t_sb in &super_times {
        let win_stop = t_sb + cfg.dph_time;
        let lo = events.time.partition_point(|&t| t < t_sb);
        let hi = events.time.partition_point(|&t| t <= win_stop);
        for flag in keep.iter_mut().take(hi).skip(lo) {
            *flag = false;
        }
        btis.push((t_sb, win_stop));
    }

    // Weighted spatial hotness check (`magnitude_k * magnitude_m / d`, `allowable_hot = total_hotness / n_hot_pix`,
    // `cztsuperbunchclean.c:600-665`) over single `dph_time` window for heavy bunches
    let mut dph = vec![0u32; QUAD_DIM * QUAD_DIM];
    let mut dph_hot = vec![0u8; QUAD_DIM * QUAD_DIM];
    for &t_hb in &heavy_times {
        let win_stop = t_hb + cfg.dph_time;
        let lo = events.time.partition_point(|&t| t < t_hb);
        let hi = events.time.partition_point(|&t| t <= win_stop);
        if hi <= lo + 2 {
            continue;
        }
        dph.fill(0);
        dph_hot.fill(0);
        for k in lo..hi {
            if !keep[k] {
                continue;
            }
            let x = (events.detx[k] as usize).min(QUAD_DIM - 1);
            let y = (events.dety[k] as usize).min(QUAD_DIM - 1);
            dph[x * QUAD_DIM + y] += 1;
        }

        // Collect all non-zero DPH pixels with their precomputed 3x3 non-zero neighbor count,
        // matching `cztsuperbunchclean.c:536-618` pairwise `counter4 > 2` guard.
        let mut active: Vec<(usize, usize, u32, usize)> = Vec::new();
        for x in 0..QUAD_DIM {
            for y in 0..QUAD_DIM {
                let mag = dph[x * QUAD_DIM + y];
                if mag == 0 {
                    continue;
                }
                let mut nbrs = 0usize;
                let x0 = x.saturating_sub(1);
                let x1 = (x + 1).min(QUAD_DIM - 1);
                let y0 = y.saturating_sub(1);
                let y1 = (y + 1).min(QUAD_DIM - 1);
                for nx in x0..=x1 {
                    for ny in y0..=y1 {
                        if dph[nx * QUAD_DIM + ny] > 0 {
                            nbrs += 1;
                        }
                    }
                }
                active.push((x, y, mag, nbrs));
            }
        }

        let mut total_hotness = 0.0f64;
        for a in 0..active.len() {
            let (xa, ya, ma, nbrs_a) = active[a];
            for &(xb, yb, mb, nbrs_b) in &active[(a + 1)..] {
                let dx = (xa as f64) - (xb as f64);
                let dy = (ya as f64) - (yb as f64);
                let d = (dx * dx + dy * dy).sqrt();
                if d > 0.0 {
                    // Legacy `cztsuperbunchclean.c:567-597`: only apply 3x3 neighborhood guard
                    // pairwise when BOTH pixels have `dph == 1` and `d > 1.5`.
                    let counter4 = if ma == 1 && mb == 1 && d > 1.5 {
                        nbrs_a + nbrs_b
                    } else {
                        10
                    };
                    if counter4 > 2 {
                        let cij = ((ma as f64) * (mb as f64)) / d;
                        if cij > cfg.thresh_hot {
                            total_hotness += cij;
                            dph_hot[xa * QUAD_DIM + ya] = 1;
                            dph_hot[xb * QUAD_DIM + yb] = 1;
                        }
                    }
                }
            }
        }
        let n_hot_pix: usize = dph_hot.iter().map(|&v| v as usize).sum();
        let allowable_hot = if n_hot_pix > 0 {
            total_hotness / (n_hot_pix as f64)
        } else {
            0.0
        };
        if allowable_hot > cfg.allowable_hot_thresh {
            for flag in keep.iter_mut().take(hi).skip(lo) {
                *flag = false;
            }
            btis.push((t_hb, win_stop));
        }
    }

    // Merge overlapping BTIs and subtract from livetime
    let btis = crate::gtigen::normalize_intervals(&btis);
    if !out_lv.time.is_empty() {
        let binsize = if out_lv.time.len() >= 2 {
            (out_lv.time[1] - out_lv.time[0]).max(1e-6)
        } else {
            out_lv.binsize.max(1e-6)
        };
        let tstarti = (out_lv.time[0] - 0.5 * binsize) as i64;
        for &(a, b) in &btis {
            update_livetime(a, b, &mut out_lv.fracexp, tstarti, 1.0, false, binsize);
        }
        for lv in &mut out_lv.fracexp {
            *lv = lv.max(0.0);
        }
    }

    SuperBunchCleanResult {
        events: events.filter_by_mask(&keep),
        livetime: out_lv,
        bti_intervals: btis,
    }
}

/// Step 3 of Generalized v3.0: `cztheavybunchclean` (with per-pixel exposure update) —
/// merges `bunsize == 63` split bunches, computes the `25x26` detector box around pixel `119`
/// (`cztheavybunchclean.c:802-840`), rejects pairs within `<= 2 ms` below `LLD + 10 keV` in
/// `heavy_bunch_time` (`cztheavybunchclean.c:754-778`), rejects all events in `heavy_bunch_time_det`
/// / `small_bunch_time_det` inside the 25x26 box, and subtracts deadtime from `pixel_exposure`.
pub fn heavybunchclean_quad_with_exposure(
    events: &QuadEvents,
    bunches: &BunchTable,
    lld_pi: &[i16],
    in_pixel_exposure: &[f64],
    cfg: &NoiseCleanConfig,
) -> (QuadEvents, Vec<f64>) {
    let mut pix_exp = vec![1.0f64; PIXELS_PER_QUAD];
    let c_pe = PIXELS_PER_QUAD.min(in_pixel_exposure.len());
    pix_exp[..c_pe].copy_from_slice(&in_pixel_exposure[..c_pe]);

    if events.is_empty() || bunches.is_empty() {
        return (events.clone(), pix_exp);
    }

    let qid = events.quad_id;
    let merged = merge_split_bunches(bunches);
    let mut keep = vec![true; events.len()];
    let has_energy = events.energy.len() == events.len();

    // Track per-pixel deadtime in detx,dety coordinates (`cztheavybunchclean.c:842-860`)
    let mut pix_loss = vec![0.0f64; QUAD_DIM * QUAD_DIM];
    let mut pix_btistop_prev = vec![0.0f64; QUAD_DIM * QUAD_DIM];

    for (t_end, bsize, rep_det) in merged {
        let (cx, cy) = detid_pixid_to_detxy(qid, rep_det, 119);
        let detx_min = cx.saturating_sub(12) as usize;
        let detx_max = ((cx as usize) + 13).min(QUAD_DIM - 1);
        let dety_min = cy.saturating_sub(12) as usize;
        let dety_max = ((cy as usize) + 13).min(QUAD_DIM - 1);

        let is_heavy = bsize >= cfg.heavy_bunch_size;
        let box_dt = if is_heavy {
            cfg.heavy_bunch_time_det
        } else {
            cfg.small_bunch_time_det
        };

        // Update per-pixel exposure loss inside the 25x26 box
        for x in detx_min..=detx_max {
            for y in dety_min..=dety_max {
                let p = x * QUAD_DIM + y;
                let new_stop = t_end + box_dt;
                if pix_btistop_prev[p] < t_end {
                    pix_loss[p] += box_dt;
                    pix_btistop_prev[p] = new_stop;
                } else if new_stop > pix_btistop_prev[p] {
                    pix_loss[p] += new_stop - pix_btistop_prev[p];
                    pix_btistop_prev[p] = new_stop;
                }
            }
        }

        let max_win = if is_heavy {
            cfg.heavy_bunch_time.max(cfg.heavy_bunch_time_det)
        } else {
            cfg.small_bunch_time_det
        };
        let lo = events.time.partition_point(|&t| t < t_end);
        let hi = events.time.partition_point(|&t| t <= t_end + max_win);

        for k in lo..hi {
            if !keep[k] {
                continue;
            }
            let x = events.detx[k] as usize;
            let y = events.dety[k] as usize;
            if x < detx_min || x > detx_max || y < dety_min || y > dety_max {
                continue;
            }
            let dt = events.time[k] - t_end;
            if dt <= box_dt {
                keep[k] = false;
            } else if is_heavy && dt <= cfg.heavy_bunch_time && has_energy {
                let dpix = (events.detid[k] as usize) * 256 + (events.pixid[k] as usize);
                let lld_ch = lld_pi.get(dpix).copied().unwrap_or(0);
                let lld_kev = 5.0 + 0.5 * (lld_ch as f32);
                if events.energy[k] < lld_kev + cfg.module_lld_threshold_kev {
                    // Reject pairs of low-energy events within `<= 2 ms` (`cztheavybunchclean.c:754-778`)
                    for m in (k + 1)..hi {
                        if events.time[m] - events.time[k] > 0.002 {
                            break;
                        }
                        let xm = events.detx[m] as usize;
                        let ym = events.dety[m] as usize;
                        if xm >= detx_min && xm <= detx_max && ym >= dety_min && ym <= dety_max {
                            let dpix_m =
                                (events.detid[m] as usize) * 256 + (events.pixid[m] as usize);
                            let lld_m = 5.0 + 0.5 * (lld_pi.get(dpix_m).copied().unwrap_or(0) as f32);
                            if events.energy[m] < lld_m + cfg.module_lld_threshold_kev {
                                keep[k] = false;
                                keep[m] = false;
                            }
                        }
                    }
                }
            }
        }
    }

    let tot_exp = if events.exposure > 0.0 {
        events.exposure
    } else {
        (events.tstop - events.tstart).max(1.0)
    };
    for det in 0..16u8 {
        for pix in 0..=255u8 {
            let row = (det as usize) * 256 + (pix as usize);
            let (dx, dy) = detid_pixid_to_detxy(qid, det, pix);
            let loss = pix_loss[(dx as usize) * QUAD_DIM + (dy as usize)];
            if loss > 0.0 {
                pix_exp[row] = (pix_exp[row] - loss / tot_exp).clamp(0.0, 1.0);
            }
        }
    }

    (events.filter_by_mask(&keep), pix_exp)
}

/// Step 3 of Generalized v3.0: `cztheavybunchclean` (events-only wrapper preserving signature).
pub fn heavybunchclean_quad(
    events: &QuadEvents,
    bunches: &BunchTable,
    lld_pi: &[i16],
    cfg: &NoiseCleanConfig,
) -> QuadEvents {
    let dummy_exp = vec![1.0f64; PIXELS_PER_QUAD];
    heavybunchclean_quad_with_exposure(events, bunches, lld_pi, &dummy_exp, cfg).0
}

/// Output of `cztflickpixclean` for one quadrant.
#[derive(Debug, Clone)]
pub struct FlickPixCleanResult {
    pub events: QuadEvents,
    pub badpix_flags: Vec<u8>,
    pub pixel_exposure: Vec<f64>,
    pub newly_flagged_pixels: usize,
}

/// Compute overlap duration (seconds) between sorted normalized GTI intervals and `[a, b]`.
fn gti_overlap_sec(ivs: &[(f64, f64)], a: f64, b: f64) -> f64 {
    let mut total = 0.0f64;
    let start_idx = ivs.partition_point(|&(_, stop)| stop <= a);
    for &(s, e) in &ivs[start_idx..] {
        if s >= b {
            break;
        }
        let lo = s.max(a);
        let hi = e.min(b);
        if hi > lo {
            total += hi - lo;
        }
    }
    total
}

/// Step 4 of Generalized v3.0: `cztflickpixclean` — splits events at `50 keV`, normalizes pixel
/// 100 s light curves by quadrant count rate (`pix_lc / lc`), computes per-pixel **temporal**
/// mean and SD across valid bins with exposure-adaptive `(sigma_thresh_2, flick_count_thresh_2)`
/// (`(5,1) / (6,2) / (8,2)` at `<20k / <120k / >=120k s`), plus 1 s and 10 s `pixThresholdFix`
/// Poisson burst filtering and `pixel_badtime` subtraction (`cztflickpixclean.c:445-845`).
pub fn flickpixclean_quad(
    events: &QuadEvents,
    in_badpix: &[u8],
    in_pixel_exposure: &[f64],
    cfg: &NoiseCleanConfig,
) -> FlickPixCleanResult {
    flickpixclean_quad_with_gti(events, in_badpix, in_pixel_exposure, None, cfg)
}

/// `cztflickpixclean` with optional quadrant GTI intervals for the `flivetime >= 20 s` valid-bin
/// filter (`cztflickpixclean.c:445-528`).
pub fn flickpixclean_quad_with_gti(
    events: &QuadEvents,
    in_badpix: &[u8],
    in_pixel_exposure: &[f64],
    gti: Option<&[(f64, f64)]>,
    cfg: &NoiseCleanConfig,
) -> FlickPixCleanResult {
    let mut flags = vec![pix_flag::GOOD; PIXELS_PER_QUAD];
    let c_bp = PIXELS_PER_QUAD.min(in_badpix.len());
    flags[..c_bp].copy_from_slice(&in_badpix[..c_bp]);
    let mut pix_exp = vec![1.0f64; PIXELS_PER_QUAD];
    let c_pe = PIXELS_PER_QUAD.min(in_pixel_exposure.len());
    pix_exp[..c_pe].copy_from_slice(&in_pixel_exposure[..c_pe]);

    if events.len() < 10
        || !cfg.flicker_time_bin_sec.is_finite()
        || cfg.flicker_time_bin_sec <= 0.0
    {
        return FlickPixCleanResult {
            events: events.clone(),
            badpix_flags: flags,
            pixel_exposure: pix_exp,
            newly_flagged_pixels: 0,
        };
    }

    let t0 = events.time[0];
    let t1 = *events.time.last().unwrap_or(&t0);
    if !t0.is_finite() || !t1.is_finite() || t1 < t0 {
        return FlickPixCleanResult {
            events: events.clone(),
            badpix_flags: flags,
            pixel_exposure: pix_exp,
            newly_flagged_pixels: 0,
        };
    }

    let tot_exposure = if events.exposure > 0.0 {
        events.exposure
    } else {
        (t1 - t0).max(1.0)
    };
    // Exposure-adaptive thresholds (`cztflickpixclean.c:500-514`):
    // use legacy `(5,1) / (6,2) / (8,2)` directly unless caller explicitly overrode defaults.
    let (sigma_thresh_2, flick_count_thresh_2): (f64, usize) = if tot_exposure < 20_000.0 {
        (5.0, 1)
    } else if tot_exposure < 120_000.0 {
        (6.0, 2)
    } else {
        (8.0, 2)
    };
    let default_cfg = NoiseCleanConfig::default();
    let sigma_use = if (cfg.flicker_sigma - default_cfg.flicker_sigma).abs() > 1e-9 {
        cfg.flicker_sigma.max(1.0)
    } else {
        sigma_thresh_2
    };
    let count_use = if cfg.flicker_num_times != default_cfg.flicker_num_times {
        cfg.flicker_num_times.max(1)
    } else {
        flick_count_thresh_2
    };

    let n_active: usize = flags.iter().filter(|&&f| f < pix_flag::NOISY).count();
    if n_active == 0 {
        return FlickPixCleanResult {
            events: events.clone(),
            badpix_flags: flags,
            pixel_exposure: pix_exp,
            newly_flagged_pixels: 0,
        };
    }

    let norm_gti: Option<Vec<(f64, f64)>> = gti
        .filter(|ivs| !ivs.is_empty())
        .map(crate::gtigen::normalize_intervals);
    let min_bin_live = (0.2 * cfg.flicker_time_bin_sec).min(20.0);

    let has_energy = events.energy.len() == events.len();
    let mut bin_lt = vec![0u32; PIXELS_PER_QUAD];
    let mut bin_gt = vec![0u32; PIXELS_PER_QUAD];
    let mut touched_lt: Vec<usize> = Vec::new();
    let mut touched_gt: Vec<usize> = Vec::new();
    // Sparse per-pixel non-zero normalized lightcurve values (`pix_lc[i][j] / lc[j]`),
    // keeping total memory O(N) even across arbitrarily large timestamp gaps.
    let mut pix_vals_lt: Vec<Vec<f64>> = vec![Vec::new(); PIXELS_PER_QUAD];
    let mut pix_vals_gt: Vec<Vec<f64>> = vec![Vec::new(); PIXELS_PER_QUAD];
    let mut n_valid_lt = 0usize;
    let mut n_valid_gt = 0usize;

    // 1. 100 s binned `<50 keV` & `>=50 keV` quadrant-rate normalized temporal lightcurve collection
    let mut cursor = 0usize;
    while cursor < events.len() {
        let tc = events.time[cursor];
        if !tc.is_finite() || tc < t0 {
            cursor += 1;
            continue;
        }
        let bin_idx = ((tc - t0) / cfg.flicker_time_bin_sec).floor();
        let bin_start = t0 + bin_idx * cfg.flicker_time_bin_sec;
        let bin_end = t0 + (bin_idx + 1.0) * cfg.flicker_time_bin_sec;
        let start_cursor = cursor;

        touched_lt.clear();
        touched_gt.clear();
        let mut lc_lt = 0u64;
        let mut lc_gt = 0u64;
        while cursor < events.len() && events.time[cursor] < bin_end {
            let idx = (events.detid[cursor] as usize) * 256 + (events.pixid[cursor] as usize);
            if idx < PIXELS_PER_QUAD && flags[idx] < pix_flag::NOISY {
                let high_e = has_energy && events.energy[cursor] >= 50.0;
                if high_e {
                    if bin_gt[idx] == 0 {
                        touched_gt.push(idx);
                    }
                    bin_gt[idx] += 1;
                    lc_gt += 1;
                } else {
                    if bin_lt[idx] == 0 {
                        touched_lt.push(idx);
                    }
                    bin_lt[idx] += 1;
                    lc_lt += 1;
                }
            }
            cursor += 1;
        }
        if cursor == start_cursor {
            cursor += 1;
        }

        let bin_live_ok = norm_gti
            .as_deref()
            .is_none_or(|ivs| gti_overlap_sec(ivs, bin_start, bin_end) >= min_bin_live);

        if bin_live_ok {
            if lc_lt > 0 {
                let inv = 1.0 / (lc_lt as f64);
                for &idx in &touched_lt {
                    pix_vals_lt[idx].push((bin_lt[idx] as f64) * inv);
                }
                n_valid_lt += 1;
            }
            if lc_gt > 0 {
                let inv = 1.0 / (lc_gt as f64);
                for &idx in &touched_gt {
                    pix_vals_gt[idx].push((bin_gt[idx] as f64) * inv);
                }
                n_valid_gt += 1;
            }
        }
        for &idx in &touched_lt {
            bin_lt[idx] = 0;
        }
        for &idx in &touched_gt {
            bin_gt[idx] = 0;
        }
    }

    // Per-pixel TEMPORAL mean and SD across valid 100 s bins (`cztflickpixclean.c:531-648`)
    let exceeds_temporal = |vals: &[f64], n_valid: usize| -> bool {
        if n_valid < 2 || vals.is_empty() {
            return false;
        }
        let m = n_valid as f64;
        let mean = vals.iter().sum::<f64>() / m;
        let nz_sq: f64 = vals.iter().map(|&v| (v - mean) * (v - mean)).sum();
        let z_sq = ((n_valid - vals.len()) as f64) * mean * mean;
        let sd = ((nz_sq + z_sq) / m).sqrt();
        if sd <= 1e-15 {
            return false;
        }
        let thresh = mean + sigma_use * sd;
        vals.iter().filter(|&&v| v > thresh).count() >= count_use
    };

    let mut newly_flagged = 0usize;
    for i in 0..PIXELS_PER_QUAD {
        if flags[i] < pix_flag::FLICKERING
            && (exceeds_temporal(&pix_vals_lt[i], n_valid_lt)
                || exceeds_temporal(&pix_vals_gt[i], n_valid_gt))
        {
            flags[i] = pix_flag::FLICKERING;
            pix_exp[i] = 0.0;
            newly_flagged += 1;
        }
    }

    let mut keep: Vec<bool> = (0..events.len())
        .map(|i| {
            let idx = (events.detid[i] as usize) * 256 + (events.pixid[i] as usize);
            idx < PIXELS_PER_QUAD && flags[idx] < pix_flag::FLICKERING
        })
        .collect();

    // 2. 1 s and 10 s `pixThresholdFix` Poisson burst filtering + `pixel_badtime` subtraction
    // (`cztflickpixclean.c:685-845`)
    let valid_pix_counter = flags.iter().filter(|&&f| f < pix_flag::FLICKERING).count();
    let mut pixel_badtime = vec![0.0f64; PIXELS_PER_QUAD];
    let mut win_counts = vec![0usize; PIXELS_PER_QUAD];

    for win_sec in [1.0f64, 10.0f64] {
        let mut c = 0usize;
        while c < events.len() {
            let tc = events.time[c];
            if !tc.is_finite() || tc < t0 {
                c += 1;
                continue;
            }
            let b = ((tc - t0) / win_sec).floor();
            let w_end = t0 + (b + 1.0) * win_sec;
            let w_start = c;
            win_counts.fill(0);
            let mut total_w = 0usize;
            while c < events.len() && events.time[c] < w_end {
                if keep[c] {
                    let idx = (events.detid[c] as usize) * 256 + (events.pixid[c] as usize);
                    if idx < PIXELS_PER_QUAD {
                        win_counts[idx] += 1;
                        total_w += 1;
                    }
                }
                c += 1;
            }
            if total_w == 0 {
                continue;
            }
            let thresh = pix_threshold_fix(total_w, valid_pix_counter);
            let mut any_burst = false;
            for p in 0..PIXELS_PER_QUAD {
                if win_counts[p] >= thresh {
                    pixel_badtime[p] += win_sec;
                    any_burst = true;
                }
            }
            if any_burst {
                for j in w_start..c {
                    if keep[j] {
                        let idx = (events.detid[j] as usize) * 256 + (events.pixid[j] as usize);
                        if idx < PIXELS_PER_QUAD && win_counts[idx] >= thresh {
                            keep[j] = false;
                        }
                    }
                }
            }
        }
    }

    for p in 0..PIXELS_PER_QUAD {
        if pixel_badtime[p] > 0.0 && flags[p] < pix_flag::FLICKERING {
            pix_exp[p] = (pix_exp[p] - pixel_badtime[p] / tot_exposure).clamp(0.0, 1.0);
        }
    }

    FlickPixCleanResult {
        events: events.filter_by_mask(&keep),
        badpix_flags: flags,
        pixel_exposure: pix_exp,
        newly_flagged_pixels: newly_flagged,
    }
}

/// Step 5 of Generalized v3.0: `czteventsep` — separates single events and Compton-scattered
/// double events occurring within `double_time_sec` (default `40 us`, pipeline supports `30-40 us`).
/// Residual 3+ event bursts within `double_time_sec` are suppressed from both singles and doubles.
pub fn eventsep_quad(events: &QuadEvents, double_time_sec: f64) -> (QuadEvents, QuadEvents) {
    let n = events.len();
    if n == 0 {
        return (QuadEvents::new(events.quad_id), QuadEvents::new(events.quad_id));
    }
    // 0 = rejected (3+ event residual burst), 1 = single, 2 = double
    let mut evt_class = vec![1u8; n];
    let mut i = 0usize;
    while i < n {
        // Find contiguous cluster of events within `double_time_sec` of `events.time[i]`
        let mut end = i + 1;
        while end < n && (events.time[end] - events.time[i]) <= double_time_sec {
            end += 1;
        }
        let cluster_len = end - i;
        if cluster_len == 2 {
            // Physical double / Compton coincidence within window
            let _ = detid_pixid_to_detxy(events.quad_id, events.detid[i], events.pixid[i]);
            let _ = detid_pixid_to_detxy(events.quad_id, events.detid[i + 1], events.pixid[i + 1]);
            evt_class[i] = 2;
            evt_class[i + 1] = 2;
        } else if cluster_len >= 3 {
            // Suppress all events in a 3+ event residual burst from both singles and doubles
            for slot in evt_class.iter_mut().take(end).skip(i) {
                *slot = 0;
            }
        }
        i = end;
    }

    let is_single: Vec<bool> = evt_class.iter().map(|&c| c == 1).collect();
    let is_double: Vec<bool> = evt_class.iter().map(|&c| c == 2).collect();
    (
        events.filter_by_mask(&is_single),
        events.filter_by_mask(&is_double),
    )
}
