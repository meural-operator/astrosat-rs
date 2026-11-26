//! OGIP CALDB Calibration Index File (`caldb.indx`) reader and query engine.

use std::path::{Path, PathBuf};

use astrosat_fits::{FitsFile, HduKind};

use crate::error::{CoreError, Result};
use crate::time::{astrosat_utc_to_mjd_utc, parse_date_to_mjd};

/// A single row from an OGIP CALDB CIF (`caldb.indx`) table.
#[derive(Debug, Clone)]
pub struct CaldbEntry {
    pub telescop: String,
    pub instrume: String,
    pub detnam: String,
    pub filter: String,
    pub codename: String,
    pub cal_dir: String,
    pub cal_file: String,
    pub cal_extno: i64,
    pub cal_vsd: String,
    pub cal_vst: String,
    pub ref_time_mjd: f64,
    pub cal_cbd: String,
    pub cal_qual: i64,
    pub full_path: PathBuf,
}

/// Reader and query engine over an AstroSat/OGIP `caldb.indx` Calibration Index File.
#[derive(Debug, Clone)]
pub struct CaldbIndex {
    caldb_root: PathBuf,
    entries: Vec<CaldbEntry>,
}

impl CaldbIndex {
    /// Open a CALDB index file.
    ///
    /// If `index_rel_path` is `None`, searches common AstroSat CALDB locations under `caldb_root`:
    /// `caldb.indx`, `data/astrosat/czti/caldb.indx`, `data/astrosat/laxpc/caldb.indx`,
    /// or `bcf/caldb.indx`.
    pub fn open(caldb_root: impl AsRef<Path>, index_rel_path: Option<&str>) -> Result<Self> {
        let root = caldb_root.as_ref().to_path_buf();
        let index_path = if let Some(rel) = index_rel_path {
            let p = Path::new(rel);
            if p.is_absolute() {
                p.to_path_buf()
            } else {
                root.join(p)
            }
        } else {
            let candidates = [
                root.join("caldb.indx"),
                root.join("data/as1/czti/caldb.indx"),
                root.join("data/as1/laxpc/caldb.indx"),
                root.join("data/astrosat/czti/caldb.indx"),
                root.join("data/astrosat/laxpc/caldb.indx"),
                root.join("bcf/caldb.indx"),
            ];
            candidates
                .into_iter()
                .find(|p| p.exists())
                .ok_or_else(|| {
                    CoreError::CaldbError(format!(
                        "could not locate caldb.indx under {}",
                        root.display()
                    ))
                })?
        };

        let fits = FitsFile::open(&index_path)?;
        let hdu = fits
            .hdu_by_name("CIF")
            .or_else(|_| {
                fits.iter_hdus()
                    .find(|h| *h.kind() == HduKind::BinTable)
                    .ok_or_else(|| astrosat_fits::FitsError::HduNotFound("CIF".into()))
            })?;
        let bt = hdu.as_bintable()?;
        let nrows = bt.nrows();

        let telescop = read_str_or_default(&bt, "TELESCOP", nrows);
        let instrume = read_str_or_default(&bt, "INSTRUME", nrows);
        let detnam = read_str_or_default(&bt, "DETNAM", nrows);
        let filter = read_str_or_default(&bt, "FILTER", nrows);
        let codename = read_str_or_default(&bt, "CAL_CNAM", nrows);
        let cal_dir = read_str_or_default(&bt, "CAL_DIR", nrows);
        let cal_file = read_str_or_default(&bt, "CAL_FILE", nrows);
        let cal_vsd = read_str_or_default(&bt, "CAL_VSD", nrows);
        let cal_vst = read_str_or_default(&bt, "CAL_VST", nrows);
        let cal_cbd = read_str_or_default(&bt, "CAL_CBD", nrows);
        let cal_extno = bt
            .read_col_i64("CAL_XNO")
            .unwrap_or_else(|_| vec![1; nrows]);
        let cal_qual = bt
            .read_col_i64("CAL_QUAL")
            .unwrap_or_else(|_| vec![0; nrows]);
        let ref_time = bt
            .read_col_f64("REF_TIME")
            .unwrap_or_else(|_| vec![0.0; nrows]);

        let mut entries = Vec::with_capacity(nrows);
        for i in 0..nrows {
            let dir_clean = cal_dir[i].trim().trim_start_matches('/');
            let file_clean = cal_file[i].trim();
            let mut full_path = root.join(dir_clean).join(file_clean);
            if !full_path.exists() {
                // Fallback relative to caldb.indx parent directory.
                if let Some(parent) = index_path.parent() {
                    let alt = parent.join(dir_clean).join(file_clean);
                    if alt.exists() {
                        full_path = alt;
                    } else {
                        // Also check stripping leading "data/astrosat/<inst>/" prefix.
                        let alt2 = parent.join("bcf").join(file_clean);
                        if alt2.exists() {
                            full_path = alt2;
                        }
                    }
                }
            }

            let mjd = if ref_time[i] > 10_000.0 {
                ref_time[i]
            } else {
                let date_str = if !cal_vst[i].trim().is_empty() {
                    format!("{}T{}", cal_vsd[i].trim(), cal_vst[i].trim())
                } else {
                    cal_vsd[i].trim().to_string()
                };
                parse_date_to_mjd(&date_str).unwrap_or(0.0)
            };

            entries.push(CaldbEntry {
                telescop: telescop[i].trim().to_string(),
                instrume: instrume[i].trim().to_string(),
                detnam: detnam[i].trim().to_string(),
                filter: filter[i].trim().to_string(),
                codename: codename[i].trim().to_string(),
                cal_dir: dir_clean.to_string(),
                cal_file: file_clean.to_string(),
                cal_extno: cal_extno[i],
                cal_vsd: cal_vsd[i].trim().to_string(),
                cal_vst: cal_vst[i].trim().to_string(),
                ref_time_mjd: mjd,
                cal_cbd: cal_cbd[i].trim().to_string(),
                cal_qual: cal_qual[i],
                full_path,
            });
        }

        Ok(Self {
            caldb_root: root,
            entries,
        })
    }

    pub fn caldb_root(&self) -> &Path {
        &self.caldb_root
    }

    pub fn entries(&self) -> &[CaldbEntry] {
        &self.entries
    }

    /// Query matching CALDB entries sorted with the newest valid calibration first.
    ///
    /// `tstart` / `tstop` are AstroSat UTC seconds since MJD 55197 (if `< 100_000`, treated directly as MJD).
    pub fn query(
        &self,
        mission: &str,
        instrument: &str,
        detnam: Option<&str>,
        codename: &str,
        tstart: Option<f64>,
        _tstop: Option<f64>,
    ) -> Vec<&CaldbEntry> {
        let obs_mjd = tstart.map(|t| {
            if t < 100_000.0 {
                t
            } else {
                astrosat_utc_to_mjd_utc(t)
            }
        });

        let mut matches: Vec<&CaldbEntry> = self
            .entries
            .iter()
            .filter(|e| e.cal_qual == 0)
            .filter(|e| mission.is_empty() || e.telescop.eq_ignore_ascii_case(mission))
            .filter(|e| instrument.is_empty() || e.instrume.eq_ignore_ascii_case(instrument))
            .filter(|e| e.codename.eq_ignore_ascii_case(codename))
            .filter(|e| match detnam {
                None => true,
                Some(d) => {
                    e.detnam.eq_ignore_ascii_case(d)
                        || e.detnam.eq_ignore_ascii_case("NONE")
                        || e.detnam == "-"
                        || e.detnam.is_empty()
                }
            })
            .filter(|e| match obs_mjd {
                Some(mjd) => e.ref_time_mjd <= mjd + 1e-6,
                None => true,
            })
            .collect();

        // Sort newest validity start date first (prefer exact DETNAM match on ties).
        matches.sort_by(|a, b| {
            b.ref_time_mjd
                .total_cmp(&a.ref_time_mjd)
                .then_with(|| {
                    let a_exact = detnam.map(|d| a.detnam.eq_ignore_ascii_case(d)).unwrap_or(false);
                    let b_exact = detnam.map(|d| b.detnam.eq_ignore_ascii_case(d)).unwrap_or(false);
                    b_exact.cmp(&a_exact)
                })
                .then_with(|| b.cal_file.cmp(&a.cal_file))
        });

        matches
    }

    /// Return the single best matching CALDB entry.
    pub fn query_best(
        &self,
        mission: &str,
        instrument: &str,
        detnam: Option<&str>,
        codename: &str,
        tstart: Option<f64>,
        tstop: Option<f64>,
    ) -> Option<&CaldbEntry> {
        let mut res = self.query(mission, instrument, detnam, codename, tstart, tstop);
        if res.is_empty() && tstart.is_some() {
            // Fallback to earliest available entry if observation predates first CAL_VSD.
            res = self.query(mission, instrument, detnam, codename, None, None);
            return res.into_iter().last();
        }
        res.into_iter().next()
    }
}

fn read_str_or_default(bt: &astrosat_fits::BinTableRef<'_>, name: &str, nrows: usize) -> Vec<String> {
    bt.read_col_strings(name)
        .unwrap_or_else(|_| vec![String::new(); nrows])
}
