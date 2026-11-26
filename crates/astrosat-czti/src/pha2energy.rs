//! Stage 1c: `cztpha2energy` — Temperature-interpolated PHA -> Energy (keV) & 512-ch PI mapping.

use rayon::prelude::*;

use crate::types::{QuadEvents, TempTable, NUM_DET_PER_QUAD, NUM_TEMP_BINS, PI_CHANNELS, PIXELS_PER_QUAD};

/// CALDB 5-temperature gain/offset calibration table for a single quadrant.
///
/// Temperature bins:
/// - bin 0: `T < 2.5 C` (0 C)
/// - bin 1: `2.5 <= T < 7.5 C` (5 C)
/// - bin 2: `7.5 <= T < 12.5 C` (10 C)
/// - bin 3: `12.5 <= T < 17.5 C` (15 C)
/// - bin 4: `T >= 17.5 C` (20 C)
#[derive(Debug, Clone)]
pub struct QuadGainOffset {
    /// `gain[temp_bin][det_id * 256 + pix_id]` (5 temperature bins x 4096 pixels).
    pub gain: Vec<Vec<f32>>,
    /// `offset[temp_bin][det_id * 256 + pix_id]` (5 temperature bins x 4096 pixels).
    pub offset: Vec<Vec<f32>>,
}

impl Default for QuadGainOffset {
    fn default() -> Self {
        Self {
            gain: vec![vec![0.5; PIXELS_PER_QUAD]; NUM_TEMP_BINS],
            offset: vec![vec![0.0; PIXELS_PER_QUAD]; NUM_TEMP_BINS],
        }
    }
}

/// Select the CALDB temperature bin index (`0..5`) for temperature `temp_c` in Celsius.
#[inline]
pub fn temp_to_bin(temp_c: f32) -> usize {
    if temp_c < 2.5 {
        0
    } else if temp_c < 7.5 {
        1
    } else if temp_c < 12.5 {
        2
    } else if temp_c < 17.5 {
        3
    } else {
        4
    }
}

/// CALDB `EBOUNDS` table mapping energy (keV) to 512 PI channels (`0..511`).
#[derive(Debug, Clone)]
pub struct EboundsTable {
    pub channel: Vec<i16>,
    pub e_min: Vec<f32>,
    pub e_max: Vec<f32>,
}

impl Default for EboundsTable {
    /// Standard CZTI 512-channel EBOUNDS grid: `5.0..=261.0` keV in `0.5` keV bins.
    fn default() -> Self {
        let mut channel = Vec::with_capacity(PI_CHANNELS);
        let mut e_min = Vec::with_capacity(PI_CHANNELS);
        let mut e_max = Vec::with_capacity(PI_CHANNELS);
        for ch in 0..PI_CHANNELS {
            channel.push(ch as i16);
            e_min.push(5.0 + 0.5 * (ch as f32));
            e_max.push(5.0 + 0.5 * ((ch + 1) as f32));
        }
        Self {
            channel,
            e_min,
            e_max,
        }
    }
}

impl EboundsTable {
    /// Map energy (keV) to PI channel via binary search, clamping out-of-range energies to
    /// `0` or `len - 1` (matching CZTI `Ebounds::find_channel_from_energy`).
    #[inline]
    pub fn find_pi(&self, energy_kev: f32) -> i16 {
        if self.e_min.is_empty() {
            return 0;
        }
        if energy_kev <= self.e_min[0] {
            return self.channel[0];
        }
        let last = self.e_max.len() - 1;
        if energy_kev >= self.e_max[last] {
            return self.channel[last];
        }
        // Binary search for first channel with e_max > energy_kev
        let idx = self.e_max.partition_point(|&hi| hi <= energy_kev);
        let clamped = idx.min(last);
        self.channel[clamped]
    }
}

/// Time-sorted temperature profile for a single quadrant extracted from `modeSS` `TEMP` HDU.
#[derive(Debug, Clone, Default)]
pub struct QuadTempSeries {
    pub time: Vec<f64>,
    pub temp: Vec<[f32; NUM_DET_PER_QUAD]>,
}

impl QuadTempSeries {
    /// Extract and time-sort the temperature rows belonging to `quad_id` from a `TempTable`.
    pub fn from_temp_table(table: &TempTable, quad_id: u8) -> Self {
        let mut rows: Vec<(f64, [f32; NUM_DET_PER_QUAD])> = table
            .time
            .iter()
            .zip(table.temperature.iter())
            .zip(table.quad_flag.iter())
            .filter_map(|((&t, &temps), &q)| if q == quad_id { Some((t, temps)) } else { None })
            .collect();
        rows.sort_by(|a, b| a.0.total_cmp(&b.0));
        let (time, temp) = rows.into_iter().unzip();
        Self { time, temp }
    }

    /// Linearly interpolate module `det_id` temperature at `ut_time` (clamped at endpoints;
    /// returns default `10.0 C` if no temperature rows exist).
    #[inline]
    pub fn interpolate_module(&self, ut_time: f64, det_id: usize) -> f32 {
        let n = self.time.len();
        if n == 0 || det_id >= NUM_DET_PER_QUAD {
            return 10.0;
        }
        if n == 1 || ut_time <= self.time[0] {
            return self.temp[0][det_id];
        }
        if ut_time >= self.time[n - 1] {
            return self.temp[n - 1][det_id];
        }
        let hi = self.time.partition_point(|&t| t < ut_time);
        if hi == 0 {
            return self.temp[0][det_id];
        }
        let lo = hi - 1;
        let t0 = self.time[lo];
        let t1 = self.time[hi];
        let dt = t1 - t0;
        if dt <= 0.0 {
            return self.temp[lo][det_id];
        }
        let frac = ((ut_time - t0) / dt) as f32;
        self.temp[lo][det_id] + frac * (self.temp[hi][det_id] - self.temp[lo][det_id])
    }
}

/// Populate `pi` and `energy` columns for a quadrant's events in parallel using `rayon`.
///
/// Implements exact CZTI `cztpha2energyv2.cpp` physics:
/// - Interpolate detector module temperature from `modeSS` `TEMP`.
/// - Select nearest of 5 CALDB temperature bins (`<2.5, [2.5,7.5), [7.5,12.5), [12.5,17.5), >=17.5 C`).
/// - Compute `E = gain * PHA + offset` (`0.0` if `PHA == 0`).
/// - Map `E` to 512-channel `PI` via CALDB `EBOUNDS`.
pub fn pha2energy_quad(
    events: &mut QuadEvents,
    temp_series: &QuadTempSeries,
    gain_offset: &QuadGainOffset,
    ebounds: &EboundsTable,
) {
    let n = events.len();
    let mut energy = vec![0.0f32; n];
    let mut pi = vec![0i16; n];

    energy
        .par_iter_mut()
        .zip(pi.par_iter_mut())
        .enumerate()
        .for_each(|(i, (e_out, pi_out))| {
            let pha = events.pha[i];
            let det = events.detid[i] as usize;
            let pix = events.pixid[i] as usize;
            let e = if pha == 0 || det >= NUM_DET_PER_QUAD {
                0.0
            } else {
                let temp_c = temp_series.interpolate_module(events.time[i], det);
                let tbin = temp_to_bin(temp_c);
                let idx = det * 256 + pix;
                let g = gain_offset.gain[tbin][idx];
                let off = gain_offset.offset[tbin][idx];
                g * (pha as f32) + off
            };
            *e_out = e;
            *pi_out = ebounds.find_pi(e);
        });

    events.energy = energy;
    events.pi = pi;
}
