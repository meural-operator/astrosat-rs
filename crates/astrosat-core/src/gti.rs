//! Good Time Interval (GTI) set algebra and fast event filtering.

use std::path::Path;

use astrosat_fits::{build_gti_hdu, FitsFile, Hdu, HduRef, OgipMeta};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use crate::error::{CoreError, Result};

/// A single closed Good Time Interval `[start, stop]` in seconds.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct GtiInterval {
    pub start: f64,
    pub stop: f64,
}

impl GtiInterval {
    pub fn new(start: f64, stop: f64) -> Self {
        Self { start, stop }
    }

    #[inline]
    pub fn duration(&self) -> f64 {
        (self.stop - self.start).max(0.0)
    }

    #[inline]
    pub fn contains(&self, t: f64) -> bool {
        t >= self.start && t <= self.stop
    }
}

/// A normalized, sorted, disjoint collection of Good Time Intervals.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct GtiSet {
    intervals: Vec<GtiInterval>,
}

impl GtiSet {
    /// Construct and normalize a `GtiSet`: drops invalid (`stop <= start`), sorts, and merges overlaps.
    pub fn new(intervals: Vec<GtiInterval>) -> Self {
        let mut valid: Vec<GtiInterval> = intervals
            .into_iter()
            .filter(|iv| iv.start.is_finite() && iv.stop.is_finite() && iv.stop > iv.start)
            .collect();
        valid.sort_by(|a, b| a.start.total_cmp(&b.start).then(a.stop.total_cmp(&b.stop)));

        let mut merged: Vec<GtiInterval> = Vec::with_capacity(valid.len());
        for iv in valid {
            if let Some(last) = merged.last_mut() {
                if iv.start <= last.stop {
                    if iv.stop > last.stop {
                        last.stop = iv.stop;
                    }
                    continue;
                }
            }
            merged.push(iv);
        }
        Self { intervals: merged }
    }

    pub fn from_pairs(pairs: &[(f64, f64)]) -> Self {
        Self::new(
            pairs
                .iter()
                .map(|&(s, e)| GtiInterval::new(s, e))
                .collect(),
        )
    }

    pub fn from_starts_stops(starts: &[f64], stops: &[f64]) -> Self {
        Self::new(
            starts
                .iter()
                .zip(stops.iter())
                .map(|(&s, &e)| GtiInterval::new(s, e))
                .collect(),
        )
    }

    /// Read a `GtiSet` from an open FITS HDU (`START` and `STOP` columns).
    pub fn from_hdu(hdu: &HduRef<'_>) -> Result<Self> {
        let bt = hdu.as_bintable()?;
        let starts = bt.read_col_f64("START")?;
        let stops = bt.read_col_f64("STOP")?;
        Ok(Self::from_starts_stops(&starts, &stops))
    }

    /// Read a `GtiSet` from a FITS file on disk (defaults to `"GTI"` or first HDU with `START`/`STOP`).
    pub fn from_fits(path: impl AsRef<Path>, extname: Option<&str>) -> Result<Self> {
        let fits = FitsFile::open(path)?;
        if let Some(name) = extname {
            let hdu = fits.hdu_by_name(name)?;
            return Self::from_hdu(&hdu);
        }
        if let Ok(hdu) = fits.hdu_by_name("GTI") {
            return Self::from_hdu(&hdu);
        }
        for hdu in fits.iter_hdus() {
            if let Ok(bt) = hdu.as_bintable() {
                if bt.column_index("START").is_some() && bt.column_index("STOP").is_some() {
                    return Self::from_hdu(&hdu);
                }
            }
        }
        Err(CoreError::InvalidInput(
            "no GTI HDU with START/STOP columns found".into(),
        ))
    }

    /// Build an OGIP-compliant FITS `GTI` HDU from this set.
    pub fn to_fits_hdu(&self, extname: &str, meta: &OgipMeta) -> Result<Hdu> {
        let starts = self.starts();
        let stops = self.stops();
        Ok(build_gti_hdu(extname, &starts, &stops, meta)?)
    }

    pub fn intervals(&self) -> &[GtiInterval] {
        &self.intervals
    }

    pub fn starts(&self) -> Vec<f64> {
        self.intervals.iter().map(|iv| iv.start).collect()
    }

    pub fn stops(&self) -> Vec<f64> {
        self.intervals.iter().map(|iv| iv.stop).collect()
    }

    pub fn len(&self) -> usize {
        self.intervals.len()
    }

    pub fn is_empty(&self) -> bool {
        self.intervals.is_empty()
    }

    pub fn tstart(&self) -> Option<f64> {
        self.intervals.first().map(|iv| iv.start)
    }

    pub fn tstop(&self) -> Option<f64> {
        self.intervals.last().map(|iv| iv.stop)
    }

    /// Total exposure (sum of all GTI interval durations) in seconds.
    pub fn exposure(&self) -> f64 {
        self.intervals.iter().map(|iv| iv.duration()).sum()
    }

    /// Merge adjacent intervals whose separating gap is strictly less than `max_gap` seconds.
    pub fn merge_gaps(&self, max_gap: f64) -> Self {
        if self.intervals.len() <= 1 || max_gap <= 0.0 {
            return self.clone();
        }
        let mut out: Vec<GtiInterval> = Vec::with_capacity(self.intervals.len());
        for &iv in &self.intervals {
            if let Some(last) = out.last_mut() {
                if iv.start - last.stop < max_gap {
                    last.stop = last.stop.max(iv.stop);
                    continue;
                }
            }
            out.push(iv);
        }
        Self { intervals: out }
    }

    /// Drop any interval shorter than `min_exposure` seconds.
    pub fn drop_short(&self, min_exposure: f64) -> Self {
        let out = self
            .intervals
            .iter()
            .copied()
            .filter(|iv| iv.duration() >= min_exposure)
            .collect();
        Self { intervals: out }
    }

    /// Exact set intersection of two `GtiSet`s in `O(N + M)`.
    pub fn intersect(&self, other: &GtiSet) -> Self {
        let a = &self.intervals;
        let b = &other.intervals;
        let mut out = Vec::with_capacity(a.len().min(b.len()));
        let (mut i, mut j) = (0usize, 0usize);
        while i < a.len() && j < b.len() {
            let lo = a[i].start.max(b[j].start);
            let hi = a[i].stop.min(b[j].stop);
            if hi > lo {
                out.push(GtiInterval::new(lo, hi));
            }
            if a[i].stop < b[j].stop {
                i += 1;
            } else {
                j += 1;
            }
        }
        Self { intervals: out }
    }

    /// Exact set union of two `GtiSet`s.
    pub fn union(&self, other: &GtiSet) -> Self {
        let mut all = Vec::with_capacity(self.intervals.len() + other.intervals.len());
        all.extend_from_slice(&self.intervals);
        all.extend_from_slice(&other.intervals);
        Self::new(all)
    }

    /// Exact set difference `self \ other` (e.g., removing Bad Time Intervals or frame gaps).
    pub fn subtract(&self, other: &GtiSet) -> Self {
        if self.is_empty() || other.is_empty() {
            return self.clone();
        }
        let b = &other.intervals;
        let mut out = Vec::with_capacity(self.intervals.len());
        let mut j = 0usize;

        for &iv in &self.intervals {
            let mut cur_start = iv.start;
            let cur_stop = iv.stop;

            while j < b.len() && b[j].stop <= cur_start {
                j += 1;
            }

            let mut k = j;
            while k < b.len() && b[k].start < cur_stop {
                if b[k].start > cur_start {
                    out.push(GtiInterval::new(cur_start, b[k].start));
                }
                cur_start = cur_start.max(b[k].stop);
                if cur_start >= cur_stop {
                    break;
                }
                k += 1;
            }
            if cur_start < cur_stop {
                out.push(GtiInterval::new(cur_start, cur_stop));
            }
        }

        Self { intervals: out }
    }

    /// O(log N) binary-search membership test (`start <= t <= stop`).
    #[inline]
    pub fn contains(&self, t: f64) -> bool {
        self.interval_index(t).is_some()
    }

    /// Return the index of the GTI interval containing `t`, if any.
    #[inline]
    pub fn interval_index(&self, t: f64) -> Option<usize> {
        let idx = self.intervals.partition_point(|iv| iv.stop < t);
        if idx < self.intervals.len() && t >= self.intervals[idx].start {
            Some(idx)
        } else {
            None
        }
    }

    /// Fast boolean mask for event filtering (`rayon`-parallelized for large event lists).
    pub fn filter_mask(&self, times: &[f64]) -> Vec<bool> {
        if self.is_empty() {
            return vec![false; times.len()];
        }
        if times.len() >= 16_384 {
            times.par_iter().map(|&t| self.contains(t)).collect()
        } else {
            times.iter().map(|&t| self.contains(t)).collect()
        }
    }

    /// Return the indices of timestamps falling inside this `GtiSet`.
    pub fn filter_indices(&self, times: &[f64]) -> Vec<usize> {
        if self.is_empty() {
            return Vec::new();
        }
        if times.len() >= 16_384 {
            times
                .par_iter()
                .enumerate()
                .filter_map(|(i, &t)| self.contains(t).then_some(i))
                .collect()
        } else {
            times
                .iter()
                .enumerate()
                .filter_map(|(i, &t)| self.contains(t).then_some(i))
                .collect()
        }
    }

    /// Compute the exact overlap duration (in seconds) between `[start, stop]` and this `GtiSet`.
    /// Useful for light curve fractional exposure (`FRACEXP`).
    pub fn overlap_duration(&self, start: f64, stop: f64) -> f64 {
        if stop <= start || self.is_empty() {
            return 0.0;
        }
        let mut idx = self.intervals.partition_point(|iv| iv.stop <= start);
        let mut total = 0.0f64;
        while idx < self.intervals.len() && self.intervals[idx].start < stop {
            let lo = self.intervals[idx].start.max(start);
            let hi = self.intervals[idx].stop.min(stop);
            if hi > lo {
                total += hi - lo;
            }
            idx += 1;
        }
        total
    }
}
