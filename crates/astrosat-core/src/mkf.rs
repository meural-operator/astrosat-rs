//! Flexible AstroSat MKF (Makefilter) housekeeping table reader.

use std::collections::HashMap;
use std::path::Path;

use astrosat_fits::{FitsFile, HduKind};

use crate::error::{CoreError, Result};
use crate::gti::{GtiInterval, GtiSet};
use crate::tct::linear_interp;

/// In-memory columnar representation of an AstroSat CZTI/LAXPC `.mkf` FITS file.
#[derive(Debug, Clone)]
pub struct MkfTable {
    times: Vec<f64>,
    columns: HashMap<String, Vec<f64>>,
    column_names: Vec<String>,
}

impl MkfTable {
    /// Load all scalar numeric columns from an AstroSat `.mkf` FITS file.
    pub fn from_fits(path: impl AsRef<Path>) -> Result<Self> {
        let fits = FitsFile::open(path)?;
        let hdu = fits
            .hdu_by_name("MKF")
            .or_else(|_| fits.hdu_by_name("MKF_DATA"))
            .or_else(|_| {
                fits.iter_hdus()
                    .find(|h| *h.kind() == HduKind::BinTable)
                    .ok_or_else(|| astrosat_fits::FitsError::HduNotFound("MKF".into()))
            })?;
        let bt = hdu.as_bintable()?;

        let mut columns = HashMap::new();
        let mut column_names = Vec::new();
        for spec in bt.columns() {
            if spec.repeat == 1 {
                if let Ok(vals) = bt.read_col_f64(&spec.name) {
                    let key = spec.name.trim().to_ascii_uppercase();
                    column_names.push(spec.name.clone());
                    columns.insert(key, vals);
                }
            }
        }

        let times = columns
            .get("TIME")
            .or_else(|| columns.get("SPS_TIME"))
            .or_else(|| columns.get("UT"))
            .cloned()
            .ok_or_else(|| CoreError::InvalidMkf("MKF missing TIME column".into()))?;

        Ok(Self {
            times,
            columns,
            column_names,
        })
    }

    /// Construct an `MkfTable` in memory (useful for unit tests).
    pub fn from_columns(times: Vec<f64>, cols: HashMap<String, Vec<f64>>) -> Result<Self> {
        let n = times.len();
        let mut columns = HashMap::with_capacity(cols.len() + 1);
        let mut column_names = Vec::with_capacity(cols.len() + 1);
        for (k, v) in cols {
            if v.len() != n {
                return Err(CoreError::InvalidMkf(format!(
                    "column '{k}' has len {}, expected {n}",
                    v.len()
                )));
            }
            column_names.push(k.clone());
            columns.insert(k.trim().to_ascii_uppercase(), v);
        }
        columns
            .entry("TIME".to_string())
            .or_insert_with(|| times.clone());
        Ok(Self {
            times,
            columns,
            column_names,
        })
    }

    pub fn len(&self) -> usize {
        self.times.len()
    }

    pub fn is_empty(&self) -> bool {
        self.times.is_empty()
    }

    pub fn times(&self) -> &[f64] {
        &self.times
    }

    pub fn column_names(&self) -> &[String] {
        &self.column_names
    }

    pub fn has_column(&self, name: &str) -> bool {
        self.columns
            .contains_key(&name.trim().to_ascii_uppercase())
    }

    /// Case-insensitive column lookup.
    pub fn column(&self, name: &str) -> Option<&[f64]> {
        self.columns
            .get(&name.trim().to_ascii_uppercase())
            .map(|v| v.as_slice())
    }

    /// Lookup the first matching column among a list of candidate aliases.
    pub fn column_any(&self, aliases: &[&str]) -> Option<&[f64]> {
        for &a in aliases {
            if let Some(c) = self.column(a) {
                return Some(c);
            }
        }
        None
    }

    /// Linearly interpolate column `name` at timestamp `t`.
    pub fn interpolate(&self, name: &str, t: f64) -> Option<f64> {
        if self.times.len() < 2 {
            return self.column(name)?.first().copied();
        }
        let col = self.column(name)?;
        Some(linear_interp(&self.times, col, t))
    }

    /// Build a [`GtiSet`] from all MKF rows satisfying `predicate(row_index)`.
    ///
    /// Each accepted row `i` spans `[t_i - half_step, t_i + half_step]` (or uses the median
    /// cadence if `dt_step` is `None`), then adjacent good intervals are merged automatically.
    pub fn gti_where(&self, predicate: impl Fn(usize) -> bool, dt_step: Option<f64>) -> GtiSet {
        let n = self.times.len();
        if n == 0 {
            return GtiSet::default();
        }
        let step = dt_step.unwrap_or_else(|| self.estimate_cadence());
        let half = 0.5 * step;

        let mut intervals = Vec::new();
        let mut run_start: Option<f64> = None;
        let mut run_end: f64 = 0.0;

        for (i, &t) in self.times.iter().enumerate() {
            if predicate(i) {
                let lo = t - half;
                let hi = t + half;
                match run_start {
                    None => {
                        run_start = Some(lo);
                        run_end = hi;
                    }
                    Some(_) => {
                        if lo <= run_end + 1.5 * step {
                            run_end = hi;
                        } else {
                            intervals.push(GtiInterval::new(run_start.unwrap(), run_end));
                            run_start = Some(lo);
                            run_end = hi;
                        }
                    }
                }
            } else if let Some(s) = run_start.take() {
                intervals.push(GtiInterval::new(s, run_end));
            }
        }
        if let Some(s) = run_start {
            intervals.push(GtiInterval::new(s, run_end));
        }
        GtiSet::new(intervals)
    }

    fn estimate_cadence(&self) -> f64 {
        if self.times.len() < 2 {
            return 1.0;
        }
        let mut diffs: Vec<f64> = self
            .times
            .windows(2)
            .take(256)
            .map(|w| w[1] - w[0])
            .filter(|&d| d > 0.0 && d.is_finite())
            .collect();
        if diffs.is_empty() {
            return 1.0;
        }
        diffs.sort_by(|a, b| a.total_cmp(b));
        diffs[diffs.len() / 2]
    }
}
