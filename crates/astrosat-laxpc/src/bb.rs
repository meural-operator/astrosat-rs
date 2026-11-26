//! Stage 3: Broad Band Counting (`modeBB` / BBC) subframe decoding, binning, and ULD quartile rates.
//!
//! Decodes 31 x 64-byte subframes per 2242-byte BB frame (`laxpcl1.f:converb`), extracting
//! `LLD`, `JAN` (valid events), `LDU` (ULD counts), and `IGE` counters, accumulating into
//! 32.768 s (`0.128 * 256 s`) bins, and computing deadtime-corrected ULD rates and quartile fractions.

use crate::frames::OrderedFrame;
use crate::types::{
    LaxpcUnit, BBC_BIN_SEC, DEADTIME_RATE_GUARD_CPS, L1_FRAME_SIZE, STBG_CLOCK_HZ,
};
use serde::{Deserialize, Serialize};

/// A single decoded BBC record (emitted for every odd subframe `j = 1, 3, ..., 31`, 16 per frame).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct BbcRecord {
    /// Instrument timestamp (seconds since STBG switch-on).
    pub time_ins: f64,
    /// Lower Level Discriminator counter (`lld`, 12-bit masked).
    pub lld: u16,
    /// Valid-event counter (`jan`, 12-bit masked).
    pub jan: u16,
    /// Upper Level Discriminator counter (`ldu` / ULD, 11-bit masked).
    pub ldu: u16,
    /// Guard/veto counter (`ige`, 11-bit masked).
    pub ige: u16,
    /// Total rejected/non-science sum `jan + ldu + ige` (`iad` in `laxpcl1.f:2621`).
    pub total_rate_counts: u16,
}

/// Decode one 2242-byte Level-1 BB frame into 16 `BbcRecord`s, bit-exact with `laxpcl1.f:converb`.
pub fn decode_bb_frame(raw: &[u8]) -> Option<Vec<BbcRecord>> {
    if raw.len() < L1_FRAME_SIZE {
        return None;
    }

    let mut records: Vec<BbcRecord> = Vec::with_capacity(16);
    let mut lld: i32 = 0;
    let mut jan: i32 = 0;
    let mut ldu: i32 = 0;
    let mut ige: i32 = 0;

    // 31 subframes (Fortran j = 1..=31), 64 bytes each
    for j in 1..=31usize {
        let base = j * 64;

        // 56-bit big-endian STBG clock at Fortran `j*64 + 48 ..= j*64 + 54` -> 0-based `base+47 .. base+54`
        let mut tbuf = [0u8; 8];
        tbuf[1..8].copy_from_slice(&raw[base + 47..base + 54]);
        let itim = u64::from_be_bytes(tbuf);
        let mut tim = (itim as f64) / STBG_CLOCK_HZ;

        // 16-bit big-endian signed counters matching `laxpcl1.f:2594-2599`
        let i1 = i16::from_be_bytes([raw[base + 56], raw[base + 57]]) as i32;
        let i2 = i16::from_be_bytes([raw[base + 58], raw[base + 59]]) as i32;
        let i3 = i16::from_be_bytes([raw[base + 60], raw[base + 61]]) as i32;
        let i4 = i16::from_be_bytes([raw[base + 64], raw[base + 65]]) as i32;

        lld = (lld + i1).rem_euclid(4096);
        jan = (jan + i2).rem_euclid(4096);
        ldu = (ldu + i3).rem_euclid(2048);
        ige = (ige + i4).rem_euclid(2048);

        if j % 2 == 1 {
            // Glitch guard on timestamp jumps > 2.06 s (`laxpcl1.f:2614-2617`)
            if j > 1 {
                let prev = records[records.len() - 1].time_ins;
                if tim - prev > 2.06 {
                    if j > 3 {
                        let prev2 = records[records.len() - 2].time_ins;
                        tim = prev + (prev - prev2);
                    } else {
                        tim = prev + 1.28e-3;
                    }
                }
            }

            records.push(BbcRecord {
                time_ins: tim,
                lld: lld as u16,
                jan: jan as u16,
                ldu: ldu as u16,
                ige: ige as u16,
                total_rate_counts: (jan + ldu + ige) as u16,
            });

            lld = 0;
            jan = 0;
            ldu = 0;
            ige = 0;
        }
    }

    Some(records)
}

/// Binned BBC time series (32.768 s bins) used for ULD rate interpolation and GTI screening.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BbcSeries {
    /// Reference start instrument time `t0b` (first valid BBC record timestamp).
    pub t0_ins: f64,
    /// Bin width in seconds (`32.768 s` by default).
    pub bin_sec: f64,
    /// Bin center instrument times `tbs[i] = t0_ins + (i + 0.5) * bin_sec`.
    pub bin_times_ins: Vec<f64>,
    /// Accumulated ULD counts per bin (`cdu`).
    pub uld_counts: Vec<f64>,
    /// Accumulated rejected/total counter (`jan + ldu + ige`) per bin (`cdr`).
    pub rejected_counts: Vec<f64>,
    /// Total valid individual BBC records accumulated across all orbits.
    pub num_records: usize,
}

impl BbcSeries {
    /// Build a binned `BbcSeries` from multi-orbit ordered BB frame batches,
    /// matching `laxpcl1.f:830-858,1118-1122`.
    pub fn from_ordered_orbits(orbit_frames: &[Vec<OrderedFrame>]) -> Self {
        let mut records: Vec<BbcRecord> = Vec::new();
        let mut t03 = 0.0_f64;
        let mut t0b: Option<f64> = None;

        for frames in orbit_frames {
            let mut orbit_last_t = 0.0_f64;
            let mut had_any = false;
            for f in frames {
                let decoded = match decode_bb_frame(&f.data) {
                    Some(r) => r,
                    None => continue,
                };
                for rec in decoded {
                    let tim = rec.time_ins;
                    if tim <= t03 || tim <= orbit_last_t || tim <= 1.0e5 || tim > 1.0e9 {
                        continue;
                    }
                    if t0b.is_none() {
                        t0b = Some(tim);
                    }
                    orbit_last_t = tim;
                    had_any = true;
                    records.push(rec);
                }
            }
            if had_any {
                t03 = orbit_last_t;
            }
        }

        Self::from_records(&records, t0b.unwrap_or(0.0), BBC_BIN_SEC)
    }

    /// Accumulate a flat slice of `BbcRecord`s into uniform `bin_sec` bins.
    pub fn from_records(records: &[BbcRecord], t0_ins: f64, bin_sec: f64) -> Self {
        const MAX_BBC_BINS: usize = 5_000_000;
        if records.is_empty()
            || !bin_sec.is_finite()
            || bin_sec <= 0.0
            || !t0_ins.is_finite()
        {
            return Self {
                t0_ins,
                bin_sec,
                bin_times_ins: Vec::new(),
                uld_counts: Vec::new(),
                rejected_counts: Vec::new(),
                num_records: 0,
            };
        }

        let mut max_bin: Option<usize> = None;
        for rec in records {
            if !rec.time_ins.is_finite() {
                continue;
            }
            let kf = ((rec.time_ins - t0_ins + 1.0e-3) / bin_sec).floor();
            if !kf.is_finite() || kf < 0.0 || kf >= (MAX_BBC_BINS as f64) {
                continue;
            }
            let k = kf as usize;
            max_bin = Some(match max_bin {
                Some(prev) => prev.max(k),
                None => k,
            });
        }
        let Some(max_b) = max_bin else {
            return Self {
                t0_ins,
                bin_sec,
                bin_times_ins: Vec::new(),
                uld_counts: Vec::new(),
                rejected_counts: Vec::new(),
                num_records: 0,
            };
        };
        let nbins = max_b + 1;
        let mut uld_counts = vec![0.0_f64; nbins];
        let mut rejected_counts = vec![0.0_f64; nbins];
        let mut accepted = 0usize;

        for rec in records {
            if !rec.time_ins.is_finite() {
                continue;
            }
            let kf = ((rec.time_ins - t0_ins + 1.0e-3) / bin_sec).floor();
            if !kf.is_finite() || kf < 0.0 || kf >= (nbins as f64) {
                continue;
            }
            let idx = kf as usize;
            uld_counts[idx] += rec.ldu as f64;
            rejected_counts[idx] += rec.total_rate_counts as f64;
            accepted += 1;
        }

        let bin_times_ins: Vec<f64> = (0..nbins)
            .map(|i| t0_ins + (i as f64 + 0.5) * bin_sec)
            .collect();

        Self {
            t0_ins,
            bin_sec,
            bin_times_ins,
            uld_counts,
            rejected_counts,
            num_records: accepted,
        }
    }

    /// Interpolate raw ULD count rate (cps) at `t_ins` using linear (`NUSE=2` divided-difference)
    /// interpolation, then apply non-paralyzable deadtime correction (`laxpcl1.f:1306-1309`).
    pub fn uld_rate_cps(&self, t_ins: f64, all_rate_cps: f64, deadtime_sec: f64) -> f64 {
        if self.bin_times_ins.is_empty() {
            return 0.0;
        }
        let raw_counts = interp_linear(&self.bin_times_ins, &self.uld_counts, t_ins);
        let mut rate = raw_counts / self.bin_sec;
        if all_rate_cps < DEADTIME_RATE_GUARD_CPS {
            let denom = 1.0 - deadtime_sec * all_rate_cps;
            if denom > 0.0 {
                rate /= denom;
            }
        }
        rate.max(0.0)
    }

    /// Interpolate deadtime-corrected rejected event rate (cps) at `t_ins` (`cur` in `laxpcl1.f:1301-1304`).
    pub fn rejected_rate_cps(&self, t_ins: f64, all_rate_cps: f64, deadtime_sec: f64) -> f64 {
        if self.bin_times_ins.is_empty() {
            return 0.0;
        }
        let raw_counts = interp_linear(&self.bin_times_ins, &self.rejected_counts, t_ins);
        let mut rate = raw_counts / self.bin_sec;
        if all_rate_cps < DEADTIME_RATE_GUARD_CPS {
            let denom = 1.0 - deadtime_sec * all_rate_cps;
            if denom > 0.0 {
                rate /= denom;
            }
        }
        rate.max(0.0)
    }

    /// Classify a ULD rate (cps) into the 1..=4 quartile bin (0 = below lower bound, 5 = above upper bound)
    /// matching `laxpcl1.f:366-368,1310,1399-1400`.
    pub fn uld_quartile_index(unit: LaxpcUnit, uld_rate_cps: f64) -> i32 {
        let (ull, ulu) = unit.uld_quartile_bounds();
        let dul = (ulu - ull) / 4.0;
        let idx = ((uld_rate_cps - ull) / dul).floor() as i32 + 1;
        idx.clamp(0, 5)
    }
}

/// Summary of ULD quartile exposure fractions across a GTI (`frb[0..=5]` in `laxpcl1.f:1497-1500`).
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct UldQuartileSummary {
    /// Raw bin counts in buckets `0` (below low), `1..=4` (quartiles), `5` (above high).
    pub counts: [u64; 6],
    /// Normalized fractions relative to sum of quartiles `1..=4`.
    pub fractions: [f64; 6],
}

impl UldQuartileSummary {
    pub fn record(&mut self, unit: LaxpcUnit, uld_rate_cps: f64) {
        let idx = BbcSeries::uld_quartile_index(unit, uld_rate_cps) as usize;
        self.counts[idx] += 1;
    }

    pub fn finalize(&mut self) {
        let total: u64 = self.counts[1..=4].iter().sum();
        if total > 0 {
            let denom = total as f64;
            for i in 0..6 {
                self.fractions[i] = self.counts[i] as f64 / denom;
            }
        }
    }
}

/// Linear interpolation matching `laxpcl1.f:DIVDIF` with `NUSE = 2`.
pub(crate) fn interp_linear(xs: &[f64], ys: &[f64], x: f64) -> f64 {
    let n = xs.len();
    if n == 0 {
        return 0.0;
    }
    if n == 1 {
        return ys[0];
    }
    if x <= xs[0] {
        let dx = xs[1] - xs[0];
        return if dx.abs() > 0.0 {
            ys[0] + (x - xs[0]) * (ys[1] - ys[0]) / dx
        } else {
            ys[0]
        };
    }
    if x >= xs[n - 1] {
        let dx = xs[n - 1] - xs[n - 2];
        return if dx.abs() > 0.0 {
            ys[n - 1] + (x - xs[n - 1]) * (ys[n - 1] - ys[n - 2]) / dx
        } else {
            ys[n - 1]
        };
    }
    let idx = xs.partition_point(|&v| v < x);
    let i1 = idx.min(n - 1);
    let i0 = i1.saturating_sub(1);
    let dx = xs[i1] - xs[i0];
    if dx.abs() == 0.0 {
        ys[i0]
    } else {
        ys[i0] + (x - xs[i0]) * (ys[i1] - ys[i0]) / dx
    }
}
