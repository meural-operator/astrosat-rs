//! Stage 2d (Default algorithm): `cztpixclean` — Iterative 5-sigma DPH bad-pixel detection,
//! 1-second module & pixel burst cuts, double-event separation, and livetime/exposure updates.

use serde::{Deserialize, Serialize};

use crate::bunchclean::{apply_gti_to_livetime, re_update_livetime};
use crate::types::{pix_flag, QuadEvents, QuadLivetime, PIXELS_PER_QUAD};

/// Configuration for `cztpixclean` (matching CZTI Pipeline v3.0.1 defaults).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PixCleanConfig {
    /// Iterative sigma-clipping threshold for noisy pixel identification on DPH (default: `5.0`).
    pub nsigma: f64,
    /// Time bin size (seconds) for module-level burst detection (default: `1.0`).
    pub det_tbinsize: f64,
    /// Time bin size (seconds) for pixel-level burst detection (default: `1.0`).
    pub pix_tbinsize: f64,
    /// Maximum counts per module in `det_tbinsize` before rejecting the bin (default: `25`).
    pub det_count_thresh: usize,
    /// Maximum counts per pixel in `pix_tbinsize` before rejecting the bin (default: `2`).
    pub pix_count_thresh: usize,
    /// Whether to output separated double events (default: `false`).
    pub write_double_evt: bool,
}

impl Default for PixCleanConfig {
    fn default() -> Self {
        Self {
            nsigma: 5.0,
            det_tbinsize: 1.0,
            pix_tbinsize: 1.0,
            det_count_thresh: 25,
            pix_count_thresh: 2,
            write_double_evt: false,
        }
    }
}

/// Output of `pixclean` for a single quadrant.
#[derive(Debug, Clone)]
pub struct QuadPixCleanResult {
    pub single_events: QuadEvents,
    pub double_events: QuadEvents,
    pub livetime: QuadLivetime,
    pub badpix_flags: Vec<u8>,
    pub pixel_exposure: Vec<f64>,
}

/// Iterative `nsigma` clipping on a 4096-pixel quadrant DPH (`cztpixclean.cpp::findbadpix`).
///
/// Marks zero-count pixels as `DEAD (4)` and iteratively flags pixels with `dph > mean + nsig*sd`
/// as `NOISY (3)` until convergence.
pub fn find_bad_pixels(dph: &[u64], pix_flags: &mut [u8], nsig: f64) {
    let n_pix = PIXELS_PER_QUAD.min(dph.len()).min(pix_flags.len());
    for i in 0..n_pix {
        if dph[i] == 0 {
            pix_flags[i] = pix_flag::DEAD;
        }
    }

    let compute_mean = |flags: &[u8]| -> f64 {
        let mut sum = 0.0f64;
        let mut n = 0usize;
        for i in 0..n_pix {
            if flags[i] < pix_flag::NOISY {
                sum += dph[i] as f64;
                n += 1;
            }
        }
        if n > 0 {
            sum / (n as f64)
        } else {
            0.0
        }
    };

    let compute_stdev = |flags: &[u8], avg: f64| -> f64 {
        let mut acc = 0.0f64;
        let mut n = 0usize;
        for i in 0..n_pix {
            if flags[i] < pix_flag::NOISY {
                let diff = (dph[i] as f64) - avg;
                acc += diff * diff;
                n += 1;
            }
        }
        if n > 0 {
            (acc / (n as f64)).sqrt()
        } else {
            0.0
        }
    };

    let mut avg = f64::NAN;
    let mut prev_avg = -1.0f64;
    for _ in 0..100 {
        avg = compute_mean(pix_flags);
        if (avg - prev_avg).abs() == 0.0 {
            break;
        }
        prev_avg = avg;
        let stdev = compute_stdev(pix_flags, avg);
        let thresh = avg + nsig * stdev;
        let mut changed = false;
        for i in 0..n_pix {
            if pix_flags[i] < pix_flag::NOISY && (dph[i] as f64) > thresh {
                pix_flags[i] = pix_flag::NOISY;
                changed = true;
            }
        }
        if !changed {
            break;
        }
    }
    let _ = avg;
}

/// Run `cztpixclean` for one quadrant (exact port of `cztpixclean.cpp`).
#[allow(clippy::too_many_arguments)]
pub fn pixclean_quad(
    events: &QuadEvents,
    in_livetime: &QuadLivetime,
    gti_intervals: &[(f64, f64)],
    caldb_badpix: &[u8],
    in_pixel_exposure_frac: &[f64],
    cfg: &PixCleanConfig,
) -> QuadPixCleanResult {
    let qid = events.quad_id;
    let nrows = events.len();
    let mut pix_flags = vec![pix_flag::GOOD; PIXELS_PER_QUAD];
    let copy_bp = PIXELS_PER_QUAD.min(caldb_badpix.len());
    pix_flags[..copy_bp].copy_from_slice(&caldb_badpix[..copy_bp]);
    let mut out_lv = in_livetime.clone();

    // Update livetime with GTI
    let lv_binsize = if out_lv.time.len() >= 2 {
        (out_lv.time[1] - out_lv.time[0]).max(1e-6)
    } else {
        out_lv.binsize.max(1e-6)
    };
    apply_gti_to_livetime(gti_intervals, &out_lv.time, &mut out_lv.fracexp, lv_binsize);

    let exposure_time = if events.exposure > 0.0 {
        events.exposure
    } else {
        let gti_exp: f64 = gti_intervals.iter().map(|&(a, b)| (b - a).max(0.0)).sum();
        gti_exp.max(1.0)
    };

    let det_tbinsize = if cfg.det_tbinsize.is_finite() && cfg.det_tbinsize > 0.0 {
        cfg.det_tbinsize
    } else {
        1.0
    };
    let pix_tbinsize = if cfg.pix_tbinsize.is_finite() && cfg.pix_tbinsize > 0.0 {
        cfg.pix_tbinsize
    } else {
        1.0
    };

    let mut pixel_exposure = vec![0.0f64; PIXELS_PER_QUAD];
    for i in 0..PIXELS_PER_QUAD.min(in_pixel_exposure_frac.len()) {
        pixel_exposure[i] = if pix_flags[i] >= pix_flag::NOISY {
            0.0
        } else {
            in_pixel_exposure_frac[i] * exposure_time
        };
    }

    if nrows == 0 {
        let mut zero_frac = vec![0.0f64; PIXELS_PER_QUAD];
        for i in 0..PIXELS_PER_QUAD.min(in_pixel_exposure_frac.len()) {
            zero_frac[i] = if pix_flags[i] >= pix_flag::NOISY {
                0.0
            } else {
                in_pixel_exposure_frac[i].clamp(0.0, 1.0)
            };
        }
        return QuadPixCleanResult {
            single_events: QuadEvents::new(qid),
            double_events: QuadEvents::new(qid),
            livetime: out_lv,
            badpix_flags: pix_flags,
            pixel_exposure: zero_frac,
        };
    }

    // 1. Flag single (1) vs double (2) events (`dt <= 25 us`)
    let mut evt_flag = vec![0u8; nrows];
    let mut i = 0usize;
    while i + 1 < nrows {
        if events.time[i + 1] - events.time[i] > 2.5e-5 {
            evt_flag[i] = 1;
            i += 1;
        } else {
            evt_flag[i] = 2;
            evt_flag[i + 1] = 2;
            i += 2;
        }
    }
    if nrows > 0 && evt_flag[nrows - 1] != 2 {
        evt_flag[nrows - 1] = 1;
    }

    // 2. Build 4096-pixel DPH and run iterative n-sigma clipping
    let mut dph = vec![0u64; PIXELS_PER_QUAD];
    for row in 0..nrows {
        let idx = (events.detid[row] as usize) * 256 + (events.pixid[row] as usize);
        if idx < PIXELS_PER_QUAD {
            dph[idx] += 1;
        }
    }
    find_bad_pixels(&dph, &mut pix_flags, cfg.nsigma);

    let mut ngoodpix = 0usize;
    for p in 0..PIXELS_PER_QUAD {
        if pix_flags[p] >= pix_flag::NOISY {
            pixel_exposure[p] = 0.0;
        } else {
            ngoodpix += 1;
        }
    }

    // 3. Flag noisy/dead pixels & perform module-level burst cut
    let tstart = events.time[0].trunc();
    let lv_tstarti = if !out_lv.time.is_empty() {
        (out_lv.time[0] - lv_binsize * 0.5) as i64
    } else {
        tstart as i64
    };
    let mut exposure_loss = vec![0.0f64; PIXELS_PER_QUAD];

    let mut detlc = [0usize; 16];
    let mut prev_det_tbin: i64 = -1;
    let mut start_evt: usize = 0;

    let flush_det_bin = |prev_bin: i64,
                         s_evt: usize,
                         e_evt: usize,
                         detlc: &[usize; 16],
                         evt_flag: &mut [u8],
                         exposure_loss: &mut [f64],
                         out_lv: &mut QuadLivetime| {
        if prev_bin < 0 {
            return;
        }
        for j in s_evt..=e_evt {
            let d = events.detid[j] as usize;
            if d < 16 && detlc[d] > cfg.det_count_thresh {
                evt_flag[j] = 0;
            }
        }
        for (d, &cnt) in detlc.iter().enumerate() {
            if cnt > cfg.det_count_thresh {
                for k in 0..256 {
                    if pix_flags[d * 256 + k] < pix_flag::NOISY {
                        exposure_loss[d * 256 + k] += det_tbinsize;
                    }
                }
                let btstart = tstart + (prev_bin as f64) * det_tbinsize;
                let btstop = tstart + ((prev_bin + 1) as f64) * det_tbinsize;
                re_update_livetime(
                    btstart,
                    btstop,
                    &in_livetime.fracexp,
                    &mut out_lv.fracexp,
                    lv_tstarti,
                    1.0 / 16.0,
                    false,
                    lv_binsize,
                );
            }
        }
    };

    for row in 0..nrows {
        let dpix = (events.detid[row] as usize) * 256 + (events.pixid[row] as usize);
        if dpix >= PIXELS_PER_QUAD || pix_flags[dpix] >= pix_flag::NOISY {
            evt_flag[row] = 0;
            continue;
        }
        if (evt_flag[row] == 1 || evt_flag[row] == 2) && events.time[row].is_finite() {
            let det_tbin = ((events.time[row] - tstart) / det_tbinsize).floor() as i64;
            if det_tbin == prev_det_tbin {
                detlc[events.detid[row] as usize] += 1;
            } else {
                if row > 0 {
                    flush_det_bin(
                        prev_det_tbin,
                        start_evt,
                        row - 1,
                        &detlc,
                        &mut evt_flag,
                        &mut exposure_loss,
                        &mut out_lv,
                    );
                }
                detlc = [0usize; 16];
                start_evt = row;
                detlc[events.detid[row] as usize] += 1;
                prev_det_tbin = det_tbin;
            }
        }
    }
    if nrows > 0 && prev_det_tbin >= 0 {
        flush_det_bin(
            prev_det_tbin,
            start_evt,
            nrows - 1,
            &detlc,
            &mut evt_flag,
            &mut exposure_loss,
            &mut out_lv,
        );
    }

    // 4. Pixel-level burst cut (`pix_count_thresh` in `pix_tbinsize`)
    let mut pixlc = vec![0usize; PIXELS_PER_QUAD];
    let mut prev_pix_tbin: i64 = -1;
    start_evt = 0;

    let flush_pix_bin = |prev_bin: i64,
                         s_evt: usize,
                         e_evt: usize,
                         pixlc: &[usize],
                         evt_flag: &mut [u8],
                         exposure_loss: &mut [f64],
                         out_lv: &mut QuadLivetime| {
        if prev_bin < 0 {
            return;
        }
        for j in s_evt..=e_evt {
            let idx = (events.detid[j] as usize) * 256 + (events.pixid[j] as usize);
            if idx < PIXELS_PER_QUAD && pixlc[idx] > cfg.pix_count_thresh {
                evt_flag[j] = 0;
            }
        }
        let pix_frac = if ngoodpix > 0 {
            1.0 / (ngoodpix as f64)
        } else {
            0.0
        };
        for (p, &cnt) in pixlc.iter().enumerate() {
            if cnt > cfg.pix_count_thresh {
                exposure_loss[p] += pix_tbinsize;
                let btstart = tstart + (prev_bin as f64) * pix_tbinsize;
                let btstop = tstart + ((prev_bin + 1) as f64) * pix_tbinsize;
                re_update_livetime(
                    btstart,
                    btstop,
                    &in_livetime.fracexp,
                    &mut out_lv.fracexp,
                    lv_tstarti,
                    pix_frac,
                    false,
                    lv_binsize,
                );
            }
        }
    };

    for row in 0..nrows {
        if (evt_flag[row] == 1 || evt_flag[row] == 2) && events.time[row].is_finite() {
            let pix_tbin = ((events.time[row] - tstart) / pix_tbinsize).floor() as i64;
            let idx = (events.detid[row] as usize) * 256 + (events.pixid[row] as usize);
            if pix_tbin == prev_pix_tbin {
                if idx < PIXELS_PER_QUAD {
                    pixlc[idx] += 1;
                }
            } else {
                if row > 0 {
                    flush_pix_bin(
                        prev_pix_tbin,
                        start_evt,
                        row - 1,
                        &pixlc,
                        &mut evt_flag,
                        &mut exposure_loss,
                        &mut out_lv,
                    );
                }
                pixlc.fill(0);
                start_evt = row;
                if idx < PIXELS_PER_QUAD {
                    pixlc[idx] += 1;
                }
                prev_pix_tbin = pix_tbin;
            }
        }
    }
    if nrows > 0 && prev_pix_tbin >= 0 {
        flush_pix_bin(
            prev_pix_tbin,
            start_evt,
            nrows - 1,
            &pixlc,
            &mut evt_flag,
            &mut exposure_loss,
            &mut out_lv,
        );
    }

    // 5. Companion check for double events (`30 us` window)
    for row in 1..nrows.saturating_sub(1) {
        if evt_flag[row] == 2 {
            let bad_next = evt_flag[row + 1] != 2 && (events.time[row + 1] - events.time[row]) < 3.0e-5;
            let bad_prev = evt_flag[row - 1] != 2 && (events.time[row] - events.time[row - 1]) < 3.0e-5;
            if bad_next || bad_prev {
                evt_flag[row] = 0;
            }
        }
    }
    if nrows >= 2 {
        if evt_flag[0] == 2 && evt_flag[1] != 2 {
            evt_flag[0] = 0;
        }
        if evt_flag[nrows - 1] == 2 && evt_flag[nrows - 2] != 2 {
            evt_flag[nrows - 1] = 0;
        }
    } else if nrows == 1 && evt_flag[0] == 2 {
        evt_flag[0] = 0;
    }

    // 6. Compute updated fractional pixel exposures
    for p in 0..PIXELS_PER_QUAD {
        pixel_exposure[p] = ((pixel_exposure[p] - exposure_loss[p]) / exposure_time).clamp(0.0, 1.0);
    }

    let single_mask: Vec<bool> = evt_flag.iter().map(|&f| f == 1).collect();
    let double_mask: Vec<bool> = evt_flag.iter().map(|&f| f == 2).collect();

    let single_events = events.filter_by_mask(&single_mask);
    let double_events = if cfg.write_double_evt {
        events.filter_by_mask(&double_mask)
    } else {
        QuadEvents::new(qid)
    };

    QuadPixCleanResult {
        single_events,
        double_events,
        livetime: out_lv,
        badpix_flags: pix_flags,
        pixel_exposure,
    }
}
