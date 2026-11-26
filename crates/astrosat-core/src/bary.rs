//! Solar System Barycentric (SSB) timing correction engine (`as1bary` physics).
//!
//! Includes:
//! - Exact relativistic Einstein (`TDB - TT`), Roemer (Earth + spacecraft), spacecraft Lorentz
//!   `(r_sc . v_earth) / c^2`, and solar Shapiro delay terms.
//! - FK5 <-> ICRS frame-bias rotation.
//! - Built-in zero-dependency analytical Earth/Moon/Sun/8-planet barycentric ephemeris
//!   (Standish 1992 Table 2a + 34-term frequency-domain Hill EMB perturbations + ELP2000/Meeus
//!   lunar reflex, ~4–9 ms (< 10 ms) accuracy vs JPL DE405 with zero external files).
//! - Optional JPL DE200/DE405/DE421/DE430 Chebyshev ephemeris reader supporting BOTH raw JPL
//!   binary (`lnxp*.405`) and NASA HEASARC FITS (`JPLEPH.405` / `JPLEPH.200`) formats (< 30 ns parity).
//! - Universal `barycenter_fits(...)` supporting CZTI, LAXPC Format B, and LAXPC Format A natively.

use std::f64::consts::PI;
use std::path::{Path, PathBuf};

use astrosat_fits::{write_fits, FitsFile, Hdu, HduKind};
use rayon::prelude::*;

use crate::error::{CoreError, Result};
use crate::orbit::OrbitTable;
use crate::time::{
    astrosat_utc_to_jd_tt, astrosat_utc_to_tt_sec, tdb_minus_tt, ASTROSAT_MJDREFI, JD_J2000,
    SECS_PER_DAY,
};

/// Exact speed of light in vacuum (m/s).
pub const C_LIGHT: f64 = 299_792_458.0;
/// Astronomical Unit in meters (JPL DE405).
pub const AU_M: f64 = 149_597_870_691.0;
/// Solar gravitational radius divided by `c^3`: `GM_sun / c^3` (seconds).
pub const GM_SUN_OVER_C3: f64 = 4.925_490_947e-6;
/// Earth-to-Moon mass ratio (DE405).
const EMRAT_DE405: f64 = 81.300_56;
/// J2000 mean obliquity of the ecliptic (radians): 84381.448 arcsec.
const EPS_J2000_RAD: f64 = (84_381.448 / 3600.0) * PI / 180.0;
const DEG2RAD: f64 = PI / 180.0;
const ARCSEC2RAD: f64 = PI / (180.0 * 3600.0);

/// Celestial coordinate reference frame of the target `(RA, DEC)`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CoordFrame {
    Icrs,
    Fk5,
}

/// Solar System Barycentric state of Earth and Sun in ICRS/J2000 equatorial meters & m/s.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SolarSystemState {
    pub earth_pos_m: [f64; 3],
    pub earth_vel_ms: [f64; 3],
    pub sun_pos_m: [f64; 3],
}

/// Ephemeris backend selector.
#[derive(Debug, Clone)]
pub enum EphemerisKind {
    /// Built-in zero-dependency analytical ephemeris (~4–9 ms accuracy; use `JplDe` for < 1 us).
    Analytical,
    /// External JPL DE200/DE405/DE421/DE430 ephemeris file (raw binary or HEASARC FITS `JPLEPH.*`, < 30 ns parity).
    JplDe(PathBuf),
}

/// Loaded Solar System ephemeris (built-in analytical or JPL DE binary).
#[derive(Debug, Clone)]
pub enum BaryEphemeris {
    Analytical,
    JplDe(Box<JplDeEphemeris>),
}

impl BaryEphemeris {
    pub fn analytical() -> Self {
        Self::Analytical
    }

    pub fn from_jpl_file(path: impl AsRef<Path>) -> Result<Self> {
        Ok(Self::JplDe(Box::new(JplDeEphemeris::open(path)?)))
    }

    pub fn from_kind(kind: &EphemerisKind) -> Result<Self> {
        match kind {
            EphemerisKind::Analytical => Ok(Self::Analytical),
            EphemerisKind::JplDe(p) => Self::from_jpl_file(p),
        }
    }

    pub fn plephem_name(&self) -> String {
        match self {
            Self::Analytical => "ASTROSAT-ANALYTICAL".to_string(),
            Self::JplDe(de) => format!("JPL-DE{}", de.numde),
        }
    }

    /// Evaluate Earth SSB position/velocity and Sun SSB position at split `JD (TDB)`.
    pub fn earth_and_sun_state(&self, jd_tdb_hi: f64, jd_tdb_lo: f64) -> Result<SolarSystemState> {
        match self {
            Self::Analytical => Ok(analytical_earth_and_sun(jd_tdb_hi, jd_tdb_lo)),
            Self::JplDe(de) => de.earth_and_sun_state(jd_tdb_hi, jd_tdb_lo),
        }
    }
}

/// Convert `(ra_deg, dec_deg)` in `frame` to a unit vector in ICRS.
pub fn source_unit_vector_icrs(ra_deg: f64, dec_deg: f64, frame: CoordFrame) -> [f64; 3] {
    let ra = ra_deg * DEG2RAD;
    let dec = dec_deg * DEG2RAD;
    let cd = dec.cos();
    let v = [cd * ra.cos(), cd * ra.sin(), dec.sin()];
    match frame {
        CoordFrame::Icrs => v,
        CoordFrame::Fk5 => fk5_to_icrs(v),
    }
}

/// Apply Standish's DE200/FK5 -> DE405/ICRS frame rotation (`as1bary/bary.c:c200to405`).
///
/// `dir_icrs = dir_fk5 + eps x dir_fk5`, where `eps = (0.002, 0.012, 0.006) / 206265.0`.
pub fn fk5_to_icrs(v_fk5: [f64; 3]) -> [f64; 3] {
    let rx = 0.002 / 206_265.0;
    let ry = 0.012 / 206_265.0;
    let rz = 0.006 / 206_265.0;

    let x = v_fk5[0] + ry * v_fk5[2] - rz * v_fk5[1];
    let y = v_fk5[1] + rz * v_fk5[0] - rx * v_fk5[2];
    let z = v_fk5[2] + rx * v_fk5[1] - ry * v_fk5[0];

    let norm = (x * x + y * y + z * z).sqrt();
    [x / norm, y / norm, z / norm]
}

#[inline]
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn norm3(a: [f64; 3]) -> f64 {
    dot3(a, a).sqrt()
}

/// Breakdown of all relativistic & geometric barycentric timing terms for one event.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BaryCorrection {
    /// Earth SSB Roemer delay `(r_earth . n) / c` (seconds).
    pub roemer_earth_s: f64,
    /// Spacecraft geocentric Roemer delay `(r_sc . n) / c` (seconds).
    pub roemer_sc_s: f64,
    /// Relativistic Einstein clock term `TDB - TT` (seconds).
    pub einstein_s: f64,
    /// Spacecraft-in-Earth-frame Lorentz term `(r_sc . v_earth) / c^2` (seconds).
    pub lorentz_sc_s: f64,
    /// Solar gravitational Shapiro term `+(2 GM_sun / c^3) ln(1 + n . r_hat)` (seconds).
    pub shapiro_s: f64,
    /// Sum of geometric + relativistic corrections added to `TT` (`t_tdb_s - t_tt_s`).
    pub total_delay_s: f64,
    /// Resulting barycentric `TDB` timestamp in seconds since AstroSat MJDREF (`55197.0`).
    pub t_tdb_s: f64,
}

/// Solar System Barycentric correction calculator for a fixed target `(RA, DEC)` and orbit.
#[derive(Debug, Clone)]
pub struct BaryCalculator {
    source_icrs: [f64; 3],
    orbit: OrbitTable,
    ephem: BaryEphemeris,
}

impl BaryCalculator {
    pub fn new(
        ra_deg: f64,
        dec_deg: f64,
        frame: CoordFrame,
        orbit: OrbitTable,
        ephem: BaryEphemeris,
    ) -> Self {
        // JPL DE200 is natively in the FK5/J2000 dynamical frame, so skip FK5->ICRS bias rotation
        // when `frame == CoordFrame::Fk5` and `ephem` is DE200 (matching as1bary/axBary.c:339-354).
        let effective_frame = match (&ephem, frame) {
            (BaryEphemeris::JplDe(de), CoordFrame::Fk5) if de.numde == 200 => CoordFrame::Icrs,
            _ => frame,
        };
        Self {
            source_icrs: source_unit_vector_icrs(ra_deg, dec_deg, effective_frame),
            orbit,
            ephem,
        }
    }

    pub fn unit_vector_icrs(&self) -> [f64; 3] {
        self.source_icrs
    }

    pub fn ephemeris(&self) -> &BaryEphemeris {
        &self.ephem
    }

    /// Compute full barycentric correction for an AstroSat UTC timestamp (seconds since MJD 55197).
    pub fn correct_time(&self, t_utc_sec: f64) -> Result<BaryCorrection> {
        let t_tt_s = astrosat_utc_to_tt_sec(t_utc_sec);
        let (jd_tt_hi, jd_tt_lo) = astrosat_utc_to_jd_tt(t_utc_sec);
        let einstein_s = tdb_minus_tt(jd_tt_hi, jd_tt_lo);

        let jd_tdb_lo = jd_tt_lo + einstein_s / SECS_PER_DAY;
        let ss = self.ephem.earth_and_sun_state(jd_tt_hi, jd_tdb_lo)?;
        let sc = self.orbit.interpolate(t_utc_sec)?;

        let n = self.source_icrs;
        let roemer_earth_s = dot3(ss.earth_pos_m, n) / C_LIGHT;
        let roemer_sc_s = dot3(sc.pos_m, n) / C_LIGHT;
        let lorentz_sc_s = dot3(sc.pos_m, ss.earth_vel_ms) / (C_LIGHT * C_LIGHT);

        // Vector from Sun to spacecraft: r_sc_sun = (r_earth_ssb + r_sc) - r_sun_ssb.
        let r_sun_sc = [
            ss.earth_pos_m[0] + sc.pos_m[0] - ss.sun_pos_m[0],
            ss.earth_pos_m[1] + sc.pos_m[1] - ss.sun_pos_m[1],
            ss.earth_pos_m[2] + sc.pos_m[2] - ss.sun_pos_m[2],
        ];
        let r_mag = norm3(r_sun_sc);
        let cos_theta = if r_mag > 0.0 {
            (dot3(n, r_sun_sc) / r_mag).clamp(-1.0 + 1e-12, 1.0)
        } else {
            0.0
        };
        // Shapiro correction added to arrival time: +(2 GM_sun / c^3) ln(1 + n . r_hat).
        // At conjunction (n . r_hat -> -1), ln(...) is negative, subtracting the gravitational delay.
        let shapiro_s = 2.0 * GM_SUN_OVER_C3 * (1.0 + cos_theta).ln();

        let total_delay_s =
            einstein_s + roemer_earth_s + roemer_sc_s + lorentz_sc_s + shapiro_s;
        let t_tdb_s = t_tt_s + total_delay_s;

        Ok(BaryCorrection {
            roemer_earth_s,
            roemer_sc_s,
            einstein_s,
            lorentz_sc_s,
            shapiro_s,
            total_delay_s,
            t_tdb_s,
        })
    }

    /// Parallel batch barycentric correction from AstroSat UTC seconds to TDB seconds.
    pub fn correct_times(&self, times_utc: &[f64]) -> Result<Vec<f64>> {
        if times_utc.len() >= 8_192 {
            times_utc
                .par_iter()
                .map(|&t| self.correct_time(t).map(|c| c.t_tdb_s))
                .collect()
        } else {
            times_utc
                .iter()
                .map(|&t| self.correct_time(t).map(|c| c.t_tdb_s))
                .collect()
        }
    }
}

/// Options for [`barycenter_fits`].
#[derive(Debug, Clone)]
pub struct BaryFitsOptions {
    pub frame: CoordFrame,
    pub ephemeris: EphemerisKind,
    pub write_checksums: bool,
}

impl Default for BaryFitsOptions {
    fn default() -> Self {
        Self {
            frame: CoordFrame::Icrs,
            ephemeris: EphemerisKind::Analytical,
            write_checksums: true,
        }
    }
}

/// Summary returned by [`barycenter_fits`].
#[derive(Debug, Clone)]
pub struct BaryFitsSummary {
    pub hdus_updated: usize,
    pub events_corrected: usize,
    pub gti_rows_corrected: usize,
    pub tstart_tdb: Option<f64>,
    pub tstop_tdb: Option<f64>,
    pub plephem: String,
}

/// Apply Solar System Barycentric correction to any AstroSat CZTI, LAXPC Format B, or LAXPC Format A FITS file.
///
/// Automatically resolves HDU-0 vs HDU-1 timing keywords (no `preplaxbary` workaround needed!),
/// updates `TIME` columns in event/rate extensions and `START`/`STOP` columns in all GTI extensions,
/// updates `TIMESYS = 'TDB'`, `TIMEREF = 'SOLARSYSTEM'`, `TSTART`/`TSTOP`, `PLEPHEM`, `RA_OBJ`, `DEC_OBJ`,
/// and recomputes `DATASUM`/`CHECKSUM`.
pub fn barycenter_fits(
    input_path: impl AsRef<Path>,
    output_path: impl AsRef<Path>,
    orbit_path: impl AsRef<Path>,
    ra_deg: f64,
    dec_deg: f64,
    options: &BaryFitsOptions,
) -> Result<BaryFitsSummary> {
    let fits = FitsFile::open(input_path.as_ref())?;
    let orbit = OrbitTable::from_fits(orbit_path.as_ref())?;
    let ephem = BaryEphemeris::from_kind(&options.ephemeris)?;
    let plephem = ephem.plephem_name();
    let calc = BaryCalculator::new(ra_deg, dec_deg, options.frame, orbit, ephem);

    // Resolve MJDREFI/MJDREFF from HDU 1 if present, else HDU 0 (handles LAXPC Format A seamlessly).
    let mut mjdrefi = ASTROSAT_MJDREFI;
    let mut mjdreff = 0.0f64;
    for candidate_idx in [1usize, 0usize] {
        if let Ok(h) = fits.hdu(candidate_idx) {
            if let Some(mi) = h.header().get_i64("MJDREFI") {
                mjdrefi = mi;
                mjdreff = h.header().get_f64("MJDREFF").unwrap_or(0.0);
                break;
            }
        }
    }

    // Guard against re-barycentering a file whose principal HDU is already TDB/TCB.
    let primary_timesys = fits
        .hdu(0)
        .ok()
        .and_then(|h| h.header().get_str("TIMESYS").map(|s| s.trim().to_ascii_uppercase()));
    let principal_timesys = fits
        .hdu(1)
        .ok()
        .and_then(|h| h.header().get_str("TIMESYS").map(|s| s.trim().to_ascii_uppercase()))
        .or_else(|| primary_timesys.clone());
    if matches!(principal_timesys.as_deref(), Some("TDB" | "TCB")) {
        return Err(CoreError::InvalidInput(format!(
            "input FITS is already barycenter-corrected (TIMESYS={})",
            principal_timesys.unwrap_or_default()
        )));
    }

    // Track whether Primary HDU already had its own TSTART/TSTOP (axBary.c:580-609 preserves
    // Primary's own barycentered TSTART/TSTOP rather than letting auxiliary HDUs inflate them).
    let primary_had_tstart = fits
        .hdu(0)
        .ok()
        .and_then(|h| h.header().get_f64("TSTART"))
        .is_some();
    let primary_had_tstop = fits
        .hdu(0)
        .ok()
        .and_then(|h| h.header().get_f64("TSTOP"))
        .is_some();

    let mut out_hdus: Vec<Hdu> = Vec::with_capacity(fits.len());
    let mut hdus_updated = 0usize;
    let mut events_corrected = 0usize;
    let mut gti_rows_corrected = 0usize;
    let mut global_tstart: Option<f64> = None;
    let mut global_tstop: Option<f64> = None;
    let mut gti_fallback_tstart: Option<f64> = None;
    let mut gti_fallback_tstop: Option<f64> = None;

    for hdu_ref in fits.iter_hdus() {
        let mut hdu = Hdu::from_ref(&hdu_ref);
        let mut modified = false;
        let mut updated_bounds = false;

        // Skip any individual extension already in TDB/TCB (matching axBary.c:443-449).
        let hdu_timesys = hdu
            .header
            .get_str("TIMESYS")
            .map(|s| s.trim().to_ascii_uppercase());
        if matches!(hdu_timesys.as_deref(), Some("TDB" | "TCB")) {
            out_hdus.push(hdu);
            continue;
        }

        let timezero = hdu.header.get_f64("TIMEZERO").unwrap_or(0.0);

        if *hdu_ref.kind() == HduKind::BinTable {
            let bt = hdu_ref.as_bintable()?;
            let nrows = bt.nrows();
            let row_len = bt.row_len();

            // 1. Update TIME column if present and scalar f64.
            if let Some(time_spec) = bt.column_spec("TIME") {
                if time_spec.repeat == 1 && time_spec.col_type == astrosat_fits::ColumnType::F64 {
                    let mut times = bt.read_col_f64("TIME")?;
                    if timezero != 0.0 {
                        for t in &mut times {
                            *t += timezero;
                        }
                    }
                    let corrected = calc.correct_times(&times)?;
                    let off = time_spec.offset;
                    for (r, &t_tdb) in corrected.iter().enumerate() {
                        let start = r * row_len + off;
                        hdu.data[start..start + 8].copy_from_slice(&t_tdb.to_be_bytes());
                    }
                    if let (Some(&first), Some(&last)) = (corrected.first(), corrected.last()) {
                        let lo = corrected.iter().copied().fold(first, f64::min);
                        let hi = corrected.iter().copied().fold(last, f64::max);
                        global_tstart = Some(global_tstart.map_or(lo, |g| g.min(lo)));
                        global_tstop = Some(global_tstop.map_or(hi, |g| g.max(hi)));
                        hdu.header.set_f64("TSTART", lo, "barycentric start time (TDB)");
                        hdu.header.set_f64("TSTOP", hi, "barycentric stop time (TDB)");
                        updated_bounds = true;
                    }
                    events_corrected += nrows;
                    modified = true;
                }
            }

            // 2. Update START & STOP columns in GTI extensions.
            if let (Some(start_spec), Some(stop_spec)) =
                (bt.column_spec("START"), bt.column_spec("STOP"))
            {
                if start_spec.repeat == 1
                    && stop_spec.repeat == 1
                    && start_spec.col_type == astrosat_fits::ColumnType::F64
                    && stop_spec.col_type == astrosat_fits::ColumnType::F64
                {
                    let mut starts = bt.read_col_f64("START")?;
                    let mut stops = bt.read_col_f64("STOP")?;
                    if timezero != 0.0 {
                        for t in &mut starts {
                            *t += timezero;
                        }
                        for t in &mut stops {
                            *t += timezero;
                        }
                    }
                    let starts_tdb = calc.correct_times(&starts)?;
                    let stops_tdb = calc.correct_times(&stops)?;

                    for r in 0..nrows {
                        let s_off = r * row_len + start_spec.offset;
                        let e_off = r * row_len + stop_spec.offset;
                        hdu.data[s_off..s_off + 8].copy_from_slice(&starts_tdb[r].to_be_bytes());
                        hdu.data[e_off..e_off + 8].copy_from_slice(&stops_tdb[r].to_be_bytes());
                    }
                    if let (Some(&s0), Some(&e0)) = (starts_tdb.first(), stops_tdb.first()) {
                        let lo = starts_tdb.iter().copied().fold(s0, f64::min);
                        let hi = stops_tdb.iter().copied().fold(e0, f64::max);
                        gti_fallback_tstart = Some(gti_fallback_tstart.map_or(lo, |g| g.min(lo)));
                        gti_fallback_tstop = Some(gti_fallback_tstop.map_or(hi, |g| g.max(hi)));
                        hdu.header.set_f64("TSTART", lo, "barycentric start time (TDB)");
                        hdu.header.set_f64("TSTOP", hi, "barycentric stop time (TDB)");
                        updated_bounds = true;
                    }
                    gti_rows_corrected += nrows;
                    modified = true;
                }
            }
        }

        if modified || hdu_ref.index() == 0 {
            if let Some(ts) = hdu.header.get_f64("TSTART") {
                if !updated_bounds {
                    // Match axBary.c timezero_heuristic_function for header TSTART/TSTOP.
                    let tz_hk = if timezero.abs() < 0.1 * ts.abs() {
                        timezero
                    } else {
                        0.0
                    };
                    if let Ok(c) = calc.correct_time(ts + tz_hk) {
                        hdu.header
                            .set_f64("TSTART", c.t_tdb_s, "barycentric start time (TDB)");
                    }
                }
            }
            if let Some(te) = hdu.header.get_f64("TSTOP") {
                if !updated_bounds {
                    let tz_hk = if timezero.abs() < 0.1 * te.abs() {
                        timezero
                    } else {
                        0.0
                    };
                    if let Ok(c) = calc.correct_time(te + tz_hk) {
                        hdu.header
                            .set_f64("TSTOP", c.t_tdb_s, "barycentric stop time (TDB)");
                    }
                }
            }

            hdu.header.set_str("TIMESYS", "TDB", "barycentric dynamical time");
            hdu.header.set_str("TIMEREF", "SOLARSYSTEM", "barycentric arrival times");
            hdu.header.set_str("TREFPOS", "BARYCENTER", "time reference position");
            hdu.header.set_str("TREFDIR", "RA_OBJ,DEC_OBJ", "keywords of reference direction");
            hdu.header.set_f64("TIMEZERO", 0.0, "time zero offset");
            hdu.header.set_str("PLEPHEM", &plephem, "solar system ephemeris used");
            hdu.header.set_i64("MJDREFI", mjdrefi, "integer part of MJD reference");
            hdu.header.set_f64("MJDREFF", mjdreff, "fractional part of MJD reference");
            hdu.header.set_f64("RA_OBJ", ra_deg, "[deg] target Right Ascension");
            hdu.header.set_f64("DEC_OBJ", dec_deg, "[deg] target Declination");
            let frame_str = match options.frame {
                CoordFrame::Icrs => "ICRS",
                CoordFrame::Fk5 => "FK5",
            };
            hdu.header.set_str("RADESYS", frame_str, "celestial reference system");
            hdu.header.set_str("RADECSYS", frame_str, "celestial reference system");
            if let Some(ts) = hdu.header.get_f64("TSTART") {
                hdu.header.set_str(
                    "DATE-OBS",
                    crate::time::astrosat_utc_to_iso8601(ts),
                    "start date/time (TDB)",
                );
            }
            if let Some(te) = hdu.header.get_f64("TSTOP") {
                hdu.header.set_str(
                    "DATE-END",
                    crate::time::astrosat_utc_to_iso8601(te),
                    "end date/time (TDB)",
                );
            }
            hdu.header.add_history(format!(
                "Barycentric correction applied by astrosat-rs ({plephem}, RA={ra_deg:.6}, DEC={dec_deg:.6})"
            ));
            if modified {
                hdus_updated += 1;
            }
        }

        out_hdus.push(hdu);
    }

    // Fallback to GTI bounds if no TIME column was present (GTI-only FITS).
    let global_tstart = global_tstart.or(gti_fallback_tstart);
    let global_tstop = global_tstop.or(gti_fallback_tstop);

    // Populate Primary HDU TSTART/TSTOP and DATE-OBS/DATE-END from global_tstart/global_tstop
    // ONLY when Primary did not already carry its own TSTART/TSTOP (matching axBary.c:580-609).
    if let Some(primary) = out_hdus.first_mut() {
        if !primary_had_tstart {
            if let Some(lo) = global_tstart {
                primary.header.set_f64("TSTART", lo, "barycentric start time (TDB)");
                primary.header.set_str(
                    "DATE-OBS",
                    crate::time::astrosat_utc_to_iso8601(lo),
                    "start date/time (TDB)",
                );
            }
        }
        if !primary_had_tstop {
            if let Some(hi) = global_tstop {
                primary.header.set_f64("TSTOP", hi, "barycentric stop time (TDB)");
                primary.header.set_str(
                    "DATE-END",
                    crate::time::astrosat_utc_to_iso8601(hi),
                    "end date/time (TDB)",
                );
            }
        }
    }

    write_fits(output_path, &mut out_hdus, options.write_checksums)?;

    Ok(BaryFitsSummary {
        hdus_updated,
        events_corrected,
        gti_rows_corrected,
        tstart_tdb: global_tstart,
        tstop_tdb: global_tstop,
        plephem,
    })
}

// ============================================================================
// Built-in High-Precision Analytical Solar System Barycentric Ephemeris
// ============================================================================

/// Evaluate Earth SSB position/velocity and Sun SSB position analytically at `JD (TDB)`.
///
/// Combines:
/// 1. Standish (1992) / Simon (1994) heliocentric Earth-Moon Barycenter + major planetary perturbations.
/// 2. ELP2000/Meeus geocentric lunar reflex (`1 / 82.30056`).
/// 3. Standish (1992) Table 2a giant-planet (Jupiter, Saturn, Uranus, Neptune) + terrestrial
///    Solar System Barycenter offset of the Sun.
fn analytical_earth_and_sun(jd_hi: f64, jd_lo: f64) -> SolarSystemState {
    let dt_days = (jd_hi - JD_J2000) + jd_lo;
    let (earth_pos, sun_pos) = analytical_positions_icrs_m(dt_days);

    // Central difference over h = 60 s gives Earth velocity accurate to < 1e-5 m/s.
    let h_sec = 60.0;
    let h_days = h_sec / SECS_PER_DAY;
    let (e_plus, _) = analytical_positions_icrs_m(dt_days + h_days);
    let (e_minus, _) = analytical_positions_icrs_m(dt_days - h_days);
    let earth_vel = [
        (e_plus[0] - e_minus[0]) / (2.0 * h_sec),
        (e_plus[1] - e_minus[1]) / (2.0 * h_sec),
        (e_plus[2] - e_minus[2]) / (2.0 * h_sec),
    ];

    SolarSystemState {
        earth_pos_m: earth_pos,
        earth_vel_ms: earth_vel,
        sun_pos_m: sun_pos,
    }
}

fn analytical_positions_icrs_m(dt_days: f64) -> ([f64; 3], [f64; 3]) {
    let t = dt_days / 36525.0; // Julian centuries from J2000.0

    // 1. Heliocentric EMB in J2000 ecliptic AU (Keplerian + top planetary perturbations).
    let emb_helio = emb_heliocentric_ecliptic_au(t);

    // 2. Geocentric Moon in J2000 ecliptic AU (ELP2000/Meeus leading series).
    let moon_geo = moon_geocentric_ecliptic_au(t);

    // 3. All 7 other planets heliocentric positions (Standish 1992 Table 2a, 1800-2050).
    let mercury = kepler_planet_ecliptic_au(
        t,
        0.387_099_27,
        0.000_000_37,
        0.205_635_93,
        0.000_019_06,
        7.004_979_02,
        -0.005_947_49,
        252.250_323_50,
        149_472.674_111_75,
        77.457_796_28,
        0.160_476_89,
        48.330_765_93,
        -0.125_340_81,
        0.0,
        0.0,
        0.0,
        0.0,
    );
    let venus = kepler_planet_ecliptic_au(
        t,
        0.723_335_66,
        0.000_003_90,
        0.006_776_72,
        -0.000_041_07,
        3.394_676_05,
        -0.000_788_90,
        181.979_099_50,
        58_517.815_387_29,
        131.602_467_18,
        0.002_683_29,
        76.679_842_55,
        -0.277_694_18,
        0.0,
        0.0,
        0.0,
        0.0,
    );
    let mars = kepler_planet_ecliptic_au(
        t,
        1.523_710_34,
        0.000_018_47,
        0.093_394_10,
        0.000_078_82,
        1.849_691_42,
        -0.008_131_31,
        -4.553_432_05,
        19_140.302_684_99,
        -23.943_629_59,
        0.444_410_88,
        49.559_538_91,
        -0.292_573_43,
        0.0,
        0.0,
        0.0,
        0.0,
    );
    let jupiter = kepler_planet_ecliptic_au(
        t,
        5.202_887_00,
        -0.000_116_07,
        0.048_386_24,
        -0.000_132_53,
        1.304_396_95,
        -0.001_837_14,
        34.396_440_51,
        3_034.746_127_75,
        14.728_479_83,
        0.212_526_68,
        100.473_909_09,
        0.204_691_06,
        -0.000_124_52,
        0.060_640_60,
        -0.356_354_38,
        38.351_250_00,
    );
    let saturn = kepler_planet_ecliptic_au(
        t,
        9.536_675_94,
        -0.001_250_60,
        0.053_861_79,
        -0.000_509_91,
        2.485_991_87,
        0.001_936_09,
        49.954_244_23,
        1_222.493_622_01,
        92.598_878_31,
        -0.418_972_16,
        113.662_424_48,
        -0.288_677_94,
        0.000_258_99,
        -0.134_344_69,
        0.873_201_47,
        38.351_250_00,
    );
    let uranus = kepler_planet_ecliptic_au(
        t,
        19.189_164_64,
        -0.001_961_76,
        0.047_257_44,
        -0.000_043_97,
        0.772_637_83,
        -0.002_429_39,
        313.238_104_51,
        428.482_027_85,
        170.954_276_30,
        0.408_052_81,
        74.016_925_03,
        0.042_405_89,
        0.000_583_31,
        -0.977_318_48,
        0.176_892_45,
        7.670_250_00,
    );
    let neptune = kepler_planet_ecliptic_au(
        t,
        30.069_922_76,
        0.000_262_91,
        0.008_590_48,
        0.000_051_05,
        1.770_043_47,
        0.000_353_72,
        -55.120_029_69,
        218.459_453_25,
        44.964_762_27,
        -0.322_414_64,
        131.784_225_74,
        -0.005_086_64,
        -0.000_413_48,
        0.683_463_18,
        -0.101_625_47,
        7.670_250_00,
    );

    // Reciprocal solar mass ratios (DE405, including Mercury).
    const M_MER: f64 = 1.0 / 6_023_600.0;
    const M_VEN: f64 = 1.0 / 408_523.71;
    const M_EMB: f64 = 1.0 / 328_900.56;
    const M_MAR: f64 = 1.0 / 3_098_708.0;
    const M_JUP: f64 = 1.0 / 1_047.348_6;
    const M_SAT: f64 = 1.0 / 3_497.898;
    const M_URA: f64 = 1.0 / 22_902.98;
    const M_NEP: f64 = 1.0 / 19_412.24;
    const M_TOT: f64 = 1.0 + M_MER + M_VEN + M_EMB + M_MAR + M_JUP + M_SAT + M_URA + M_NEP;

    let mut sun_ssb_ecl = [0.0f64; 3];
    for k in 0..3 {
        sun_ssb_ecl[k] = -(M_MER * mercury[k]
            + M_VEN * venus[k]
            + M_EMB * emb_helio[k]
            + M_MAR * mars[k]
            + M_JUP * jupiter[k]
            + M_SAT * saturn[k]
            + M_URA * uranus[k]
            + M_NEP * neptune[k])
            / M_TOT;
    }

    let moon_fac = 1.0 / (1.0 + EMRAT_DE405);
    let earth_ssb_ecl = [
        emb_helio[0] - moon_fac * moon_geo[0] + sun_ssb_ecl[0],
        emb_helio[1] - moon_fac * moon_geo[1] + sun_ssb_ecl[1],
        emb_helio[2] - moon_fac * moon_geo[2] + sun_ssb_ecl[2],
    ];

    (
        ecliptic_to_icrs_m(earth_ssb_ecl),
        ecliptic_to_icrs_m(sun_ssb_ecl),
    )
}

fn ecliptic_to_icrs_m(ecl_au: [f64; 3]) -> [f64; 3] {
    let ce = EPS_J2000_RAD.cos();
    let se = EPS_J2000_RAD.sin();
    [
        ecl_au[0] * AU_M,
        (ecl_au[1] * ce - ecl_au[2] * se) * AU_M,
        (ecl_au[1] * se + ecl_au[2] * ce) * AU_M,
    ]
}

#[allow(clippy::too_many_arguments)]
fn kepler_planet_ecliptic_au(
    t: f64,
    a0: f64,
    adot: f64,
    e0: f64,
    edot: f64,
    i0_deg: f64,
    idot_deg: f64,
    l0_deg: f64,
    ldot_deg: f64,
    wbar0_deg: f64,
    wbardot_deg: f64,
    om0_deg: f64,
    omdot_deg: f64,
    b_deg: f64,
    c_deg: f64,
    s_deg: f64,
    f_deg: f64,
) -> [f64; 3] {
    let a = a0 + adot * t;
    let e = e0 + edot * t;
    let inc = (i0_deg + idot_deg * t) * DEG2RAD;
    let l_deg = l0_deg + ldot_deg * t;
    let wbar_deg = wbar0_deg + wbardot_deg * t;
    let om_deg = om0_deg + omdot_deg * t;
    let arg_peri = (wbar_deg - om_deg) * DEG2RAD;
    let om = om_deg * DEG2RAD;

    let mut m_deg = l_deg - wbar_deg + b_deg * t * t;
    if f_deg != 0.0 {
        let ft = (f_deg * t) * DEG2RAD;
        m_deg += c_deg * ft.cos() + s_deg * ft.sin();
    }
    let m = (m_deg.rem_euclid(360.0)) * DEG2RAD;

    // Solve Kepler's equation E - e sin(E) = M via Newton-Raphson.
    let mut ea = m + e * m.sin();
    for _ in 0..8 {
        let d = (ea - e * ea.sin() - m) / (1.0 - e * ea.cos());
        ea -= d;
        if d.abs() < 1e-13 {
            break;
        }
    }

    let x_orb = a * (ea.cos() - e);
    let y_orb = a * (1.0 - e * e).sqrt() * ea.sin();

    let cw = arg_peri.cos();
    let sw = arg_peri.sin();
    let co = om.cos();
    let so = om.sin();
    let ci = inc.cos();
    let si = inc.sin();

    let x = (cw * co - sw * so * ci) * x_orb + (-sw * co - cw * so * ci) * y_orb;
    let y = (cw * so + sw * co * ci) * x_orb + (-sw * so + cw * co * ci) * y_orb;
    let z = (sw * si) * x_orb + (cw * si) * y_orb;
    [x, y, z]
}

/// Heliocentric Earth-Moon Barycenter in J2000 ecliptic AU with leading planetary perturbations.
fn emb_heliocentric_ecliptic_au(t: f64) -> [f64; 3] {
    // Base J2000 EMB Keplerian orbit (Standish Table 2a).
    let base = kepler_planet_ecliptic_au(
        t,
        1.000_002_61,
        0.000_005_62,
        0.016_711_23,
        -0.000_043_92,
        -0.000_015_31,
        -0.012_946_68,
        100.464_571_66,
        35_999.372_449_81,
        102.937_681_93,
        0.323_273_64,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
    );

    // Apply 34-term first-order planetary longitude (arcsec) & radius (1e-6 AU) perturbations
    // (Venus, Mars, Jupiter, Saturn) derived from exact frequency-domain Hill equations on
    // the eccentric Standish Table 2a orbits.
    let l_e = (100.464_571_66 + 35_999.372_449_81 * t) * DEG2RAD;
    let l_v = (181.979_099_50 + 58_517.815_387_29 * t) * DEG2RAD;
    let l_m = (-4.553_432_05 + 19_140.302_684_99 * t) * DEG2RAD;
    let l_j = (34.396_440_51 + 3_034.746_127_75 * t) * DEG2RAD;
    let l_s = (49.954_244_23 + 1_222.493_622_01 * t) * DEG2RAD;

    // Each tuple: (angle, dL_cos_arcsec, dL_sin_arcsec, dR_cos_uAU, dR_sin_uAU).
    let pert_terms: [(f64, f64, f64, f64, f64); 34] = [
        // Venus (13 terms)
        (l_v - 3.0 * l_e, -0.0649, 0.0392, 0.1139, 0.1954),
        (l_v - 2.0 * l_e, -0.5085, 0.2535, 0.2576, 0.4737),
        (l_v - l_e, 0.0200, 4.8905, -5.5177, 0.0414),
        (2.0 * l_v - 4.0 * l_e, 0.1985, -0.0869, -0.1401, -0.3495),
        (2.0 * l_v - 3.0 * l_e, 1.9711, 0.1817, -0.1751, 1.4642),
        (2.0 * l_v - 2.0 * l_e, 0.0079, -5.5932, 15.9381, 0.0311),
        (3.0 * l_v - 5.0 * l_e, 0.3581, -1.1861, -0.4454, -0.0862),
        (3.0 * l_v - 4.0 * l_e, 1.3472, -0.1873, 0.4063, 2.9462),
        (3.0 * l_v - 3.0 * l_e, 0.0021, -0.6900, 2.5124, 0.0077),
        (4.0 * l_v - 6.0 * l_e, -0.0427, -0.1117, 0.1985, -0.0030),
        (4.0 * l_v - 5.0 * l_e, -0.1386, 0.0160, -0.0497, -0.4463),
        (4.0 * l_v - 4.0 * l_e, 0.0010, -0.2089, 0.8670, 0.0046),
        (5.0 * l_v - 5.0 * l_e, 0.0004, -0.0844, 0.3781, 0.0012),
        // Mars (7 terms)
        (l_m - l_e, -0.0167, 0.2687, 0.3392, 0.0226),
        (2.0 * l_m - 2.0 * l_e, -0.0436, 2.0542, 4.7540, 0.1010),
        (2.0 * l_m - l_e, 1.4770, 1.6504, -0.3504, 0.3137),
        (3.0 * l_m - 3.0 * l_e, 0.0027, -0.1220, -0.3676, -0.0079),
        (3.0 * l_m - 2.0 * l_e, 0.2233, 0.3655, 0.4289, -0.2611),
        (4.0 * l_m - 3.0 * l_e, 0.2941, 0.4216, 0.9280, -0.6472),
        (4.0 * l_m - 2.0 * l_e, 0.5935, 0.3043, -0.1274, 0.2488),
        // Jupiter (9 terms)
        (l_j - l_e, 0.2882, 7.0653, 15.9526, -0.6639),
        (l_j, 1.5217, -1.7123, 0.2144, 0.3426),
        (l_j + l_e, -0.0578, 0.3084, -0.7880, -0.1540),
        (2.0 * l_j - 2.0 * l_e, -0.0132, -2.7124, -9.2278, 0.0433),
        (2.0 * l_j - l_e, 1.9818, 1.2433, 2.6207, -4.2334),
        (2.0 * l_j, -0.0120, -0.3905, 0.1847, -0.0289),
        (3.0 * l_j - 3.0 * l_e, -0.0025, -0.1646, -0.6497, 0.0097),
        (3.0 * l_j - 2.0 * l_e, 0.1536, -0.6474, -2.1571, -0.5098),
        (3.0 * l_j - l_e, 0.1876, 0.1862, 0.3550, -0.3752),
        // Saturn (5 terms)
        (l_s - l_e, -0.0380, 0.4712, 1.1115, 0.0901),
        (l_s, 0.1346, 0.2127, -0.0167, 0.0136),
        (2.0 * l_s - 2.0 * l_e, -0.0007, -0.1047, -0.3640, 0.0025),
        (2.0 * l_s - l_e, 0.1818, 0.0273, 0.0634, -0.4205),
        (2.0 * l_s, 0.0777, -0.1185, 0.0264, 0.0145),
    ];

    let mut d_lon_arcsec = 0.0f64;
    let mut d_r_uau = 0.0f64;
    for (ang, dlc, dls, drc, drs) in pert_terms {
        let (s, c) = ang.sin_cos();
        d_lon_arcsec += dlc * c + dls * s;
        d_r_uau += drc * c + drs * s;
    }

    let d_lon = d_lon_arcsec * ARCSEC2RAD;
    let d_r_au = d_r_uau * 1.0e-6;
    let r = (base[0] * base[0] + base[1] * base[1]).sqrt();
    let lon = base[1].atan2(base[0]) + d_lon;
    let r_new = r + d_r_au;

    [r_new * lon.cos(), r_new * lon.sin(), base[2]]
}

/// Geocentric Moon position in J2000 ecliptic AU (Meeus Ch. 47 / ELP2000 leading terms).
fn moon_geocentric_ecliptic_au(t: f64) -> [f64; 3] {
    let t2 = t * t;
    let lp = (218.316_447_7 + 481_267.881_234_21 * t - 0.001_578_6 * t2) * DEG2RAD;
    let d = (297.850_192_1 + 445_267.111_403_4 * t - 0.001_881_9 * t2) * DEG2RAD;
    let m = (357.529_109_2 + 35_999.050_290_9 * t - 0.000_153_6 * t2) * DEG2RAD;
    let mp = (134.963_396_4 + 477_198.867_505_5 * t + 0.008_741_4 * t2) * DEG2RAD;
    let f = (93.272_095_0 + 483_202.017_523_3 * t - 0.003_653_9 * t2) * DEG2RAD;

    // Longitude perturbation (degrees).
    let d_lon_deg = 6.288_774 * mp.sin()
        + 1.274_027 * (2.0 * d - mp).sin()
        + 0.658_314 * (2.0 * d).sin()
        + 0.213_618 * (2.0 * mp).sin()
        - 0.185_116 * m.sin()
        - 0.114_332 * (2.0 * f).sin()
        + 0.058_793 * (2.0 * d - 2.0 * mp).sin()
        + 0.057_066 * (2.0 * d - m - mp).sin()
        + 0.053_322 * (2.0 * d + mp).sin()
        + 0.045_758 * (2.0 * d - m).sin()
        - 0.040_923 * (m - mp).sin()
        - 0.034_720 * d.sin()
        - 0.030_383 * (m + mp).sin()
        + 0.015_327 * (2.0 * d - 2.0 * f).sin()
        - 0.012_528 * (mp + 2.0 * f).sin()
        + 0.010_980 * (mp - 2.0 * f).sin()
        + 0.010_675 * (4.0 * d - mp).sin()
        + 0.010_034 * (3.0 * mp).sin();

    // Latitude (degrees).
    let lat_deg = 5.128_122 * f.sin()
        + 0.280_602 * (mp + f).sin()
        + 0.277_693 * (mp - f).sin()
        + 0.173_237 * (2.0 * d - f).sin()
        + 0.055_413 * (2.0 * d - mp + f).sin()
        + 0.046_271 * (2.0 * d - mp - f).sin()
        + 0.032_573 * (2.0 * d + f).sin()
        + 0.017_198 * (2.0 * mp + f).sin()
        + 0.009_266 * (2.0 * d + mp - f).sin()
        + 0.008_822 * (2.0 * mp - f).sin();

    // Distance (km).
    let dist_km = 385_000.56
        - 20_905.355 * mp.cos()
        - 3_699.111 * (2.0 * d - mp).cos()
        - 2_955.968 * (2.0 * d).cos()
        - 569.925 * (2.0 * mp).cos()
        + 48.888 * m.cos()
        - 3.149 * (2.0 * f).cos()
        + 246.158 * (2.0 * d - 2.0 * mp).cos()
        - 152.138 * (2.0 * d - m - mp).cos()
        - 170.733 * (2.0 * d + mp).cos()
        - 204.586 * (2.0 * d - m).cos()
        - 129.620 * (m - mp).cos()
        + 108.743 * d.cos()
        + 104.755 * (m + mp).cos();

    // Remove general precession in longitude from mean equinox of date back to J2000.
    let precess_rad = (1.396_971 * t) * DEG2RAD;
    let lon = lp + d_lon_deg * DEG2RAD - precess_rad;
    let lat = lat_deg * DEG2RAD;
    let r_au = (dist_km * 1000.0) / AU_M;

    let cl = lat.cos();
    [r_au * cl * lon.cos(), r_au * cl * lon.sin(), r_au * lat.sin()]
}

// ============================================================================
// JPL Binary DE200 / DE405 / DE421 / DE430 Chebyshev Ephemeris Reader
// ============================================================================

/// Reader for official JPL Solar System ephemerides supporting BOTH:
/// - Raw JPL binary files (`lnxp1600p2200.405`, `linux_p1550p2650.430`, etc.)
/// - NASA HEASARC FITS files (`JPLEPH.405`, `JPLEPH.200`, matching `as1bary/dpleph.c:530-605`).
#[derive(Debug, Clone)]
pub struct JplDeEphemeris {
    pub numde: i32,
    pub jd_start: f64,
    pub jd_end: f64,
    pub span_days: f64,
    pub au_km: f64,
    pub emrat: f64,
    ipt: [[usize; 3]; 13],
    ncoeff: usize,
    big_endian: bool,
    /// Number of header records preceding data record 0 (`2` for raw JPL binary, `0` for FITS `DE3`).
    data_record_offset: usize,
    data: Vec<u8>,
}

impl JplDeEphemeris {
    pub fn open(path: impl AsRef<Path>) -> Result<Self> {
        let data = std::fs::read(path.as_ref())?;
        if data.starts_with(b"SIMPLE  =") {
            return Self::from_heasarc_fits_bytes(data);
        }
        Self::from_raw_binary_bytes(data)
    }

    /// Parse a NASA HEASARC FITS JPL ephemeris (`JPLEPH.405`, `JPLEPH.200`) with `DE1`/`DE2`/`DE3` HDUs.
    fn from_heasarc_fits_bytes(bytes: Vec<u8>) -> Result<Self> {
        let fits = FitsFile::from_bytes(bytes)?;

        // HDU DE1: Cname / Cvalue -> DENUM, AU, EMRAT.
        let de1_hdu = fits
            .hdu_by_name("DE1")
            .or_else(|_| fits.hdu(1))
            .map_err(|e| CoreError::InvalidEphemeris(format!("HEASARC FITS missing DE1: {e}")))?;
        let de1 = de1_hdu.as_bintable()?;
        let cnames = de1.read_col_strings("Cname")?;
        let cvals = de1.read_col_f64("Cvalue")?;

        let mut numde: Option<i32> = None;
        let mut au_km: Option<f64> = None;
        let mut emrat: Option<f64> = None;
        for (name, &val) in cnames.iter().zip(cvals.iter()) {
            let k = name.trim().to_ascii_uppercase();
            if k == "DENUM" {
                numde = Some(val.round() as i32);
            } else if k == "AU" {
                au_km = Some(val);
            } else if k == "EMRAT" {
                emrat = Some(val);
            }
        }
        let numde = numde.ok_or_else(|| {
            CoreError::InvalidEphemeris("HEASARC FITS DE1 missing DENUM constant".into())
        })?;
        let au_km = au_km.ok_or_else(|| {
            CoreError::InvalidEphemeris("HEASARC FITS DE1 missing AU constant".into())
        })?;
        let emrat = emrat.ok_or_else(|| {
            CoreError::InvalidEphemeris("HEASARC FITS DE1 missing EMRAT constant".into())
        })?;

        // HDU DE2: Pointer / NumCoeff / NumSubIntv -> ipt[0..13].
        let de2_hdu = fits
            .hdu_by_name("DE2")
            .or_else(|_| fits.hdu(2))
            .map_err(|e| CoreError::InvalidEphemeris(format!("HEASARC FITS missing DE2: {e}")))?;
        let de2 = de2_hdu.as_bintable()?;
        let ptrs = de2.read_col_i64("Pointer")?;
        let ncoeffs = de2.read_col_i64("NumCoeff")?;
        let nsubs = de2.read_col_i64("NumSubIntv")?;
        if ptrs.len() < 12 || ncoeffs.len() < 12 || nsubs.len() < 12 {
            return Err(CoreError::InvalidEphemeris(format!(
                "HEASARC FITS DE2 has only {} rows (need >= 12)",
                ptrs.len()
            )));
        }
        let mut ipt = [[0usize; 3]; 13];
        for i in 0..ptrs.len().min(13) {
            ipt[i][0] = ptrs[i].max(0) as usize;
            ipt[i][1] = ncoeffs[i].max(0) as usize;
            ipt[i][2] = nsubs[i].max(0) as usize;
        }

        // HDU DE3: TSTART, TSTOP, TIMEDEL, NAXIS1, ChebCoeffs (big-endian f64 rows).
        let de3_hdu = fits
            .hdu_by_name("DE3")
            .or_else(|_| fits.hdu(3))
            .map_err(|e| CoreError::InvalidEphemeris(format!("HEASARC FITS missing DE3: {e}")))?;
        let de3_hdr = de3_hdu.header();
        let de3_bt = de3_hdu.as_bintable()?;

        let row_len = de3_bt.row_len();
        if row_len < 16 || !row_len.is_multiple_of(8) {
            return Err(CoreError::InvalidEphemeris(format!(
                "HEASARC FITS DE3 invalid NAXIS1={row_len}"
            )));
        }
        let ncoeff = row_len / 8;
        let jd_start = de3_hdr
            .get_f64("TSTART")
            .ok_or_else(|| CoreError::InvalidEphemeris("HEASARC FITS DE3 missing TSTART".into()))?;
        let jd_end = de3_hdr
            .get_f64("TSTOP")
            .ok_or_else(|| CoreError::InvalidEphemeris("HEASARC FITS DE3 missing TSTOP".into()))?;
        let span_days = de3_hdr
            .get_f64("TIMEDEL")
            .ok_or_else(|| CoreError::InvalidEphemeris("HEASARC FITS DE3 missing TIMEDEL".into()))?;

        for &req_body in &[2usize, 9usize, 10usize] {
            let [start_1, ncf, nsub] = ipt[req_body];
            if start_1 == 0 || ncf == 0 || nsub == 0 {
                return Err(CoreError::InvalidEphemeris(format!(
                    "HEASARC JPL DE{numde} missing required body {req_body} ipt entry"
                )));
            }
            let end = (start_1 - 1)
                .checked_add(
                    3usize
                        .checked_mul(ncf)
                        .and_then(|v| v.checked_mul(nsub))
                        .ok_or_else(|| {
                            CoreError::InvalidEphemeris("HEASARC JPL DE ipt overflow".into())
                        })?,
                )
                .ok_or_else(|| CoreError::InvalidEphemeris("HEASARC JPL DE ipt overflow".into()))?;
            if end > ncoeff {
                return Err(CoreError::InvalidEphemeris(format!(
                    "HEASARC JPL DE{numde} body {req_body} exceeds ncoeff={ncoeff}"
                )));
            }
        }

        let raw_table = de3_bt.raw_data().to_vec();
        if raw_table.len() < row_len
            || !span_days.is_finite()
            || span_days <= 0.0
            || !jd_start.is_finite()
            || !jd_end.is_finite()
            || jd_end <= jd_start
            || !emrat.is_finite()
            || emrat <= 0.0
        {
            return Err(CoreError::InvalidEphemeris(format!(
                "invalid HEASARC JPL DE{numde} geometry (ncoeff={ncoeff}, span={span_days})"
            )));
        }

        Ok(Self {
            numde,
            jd_start,
            jd_end,
            span_days,
            au_km,
            emrat,
            ipt,
            ncoeff,
            big_endian: true,
            data_record_offset: 0,
            data: raw_table,
        })
    }

    fn from_raw_binary_bytes(data: Vec<u8>) -> Result<Self> {
        if data.len() < 4096 {
            return Err(CoreError::InvalidEphemeris(
                "JPL DE file is too small (< 4096 bytes)".into(),
            ));
        }

        // Detect endianness from `numde` at byte offset 2840.
        let numde_le = i32::from_le_bytes(data[2840..2844].try_into().unwrap());
        let numde_be = i32::from_be_bytes(data[2840..2844].try_into().unwrap());
        let (big_endian, numde) = if (100..=999).contains(&numde_le) {
            (false, numde_le)
        } else if (100..=999).contains(&numde_be) {
            (true, numde_be)
        } else {
            return Err(CoreError::InvalidEphemeris(format!(
                "unrecognized JPL DE number (le={numde_le}, be={numde_be})"
            )));
        };

        let read_f64 = |off: usize| -> f64 {
            let b: [u8; 8] = data[off..off + 8].try_into().unwrap();
            if big_endian {
                f64::from_be_bytes(b)
            } else {
                f64::from_le_bytes(b)
            }
        };
        let read_i32 = |off: usize| -> i32 {
            let b: [u8; 4] = data[off..off + 4].try_into().unwrap();
            if big_endian {
                i32::from_be_bytes(b)
            } else {
                i32::from_le_bytes(b)
            }
        };

        let jd_start = read_f64(2652);
        let jd_end = read_f64(2660);
        let span_days = read_f64(2668);
        let au_km = read_f64(2680);
        let emrat = read_f64(2688);

        let mut ipt = [[0usize; 3]; 13];
        for (i, row) in ipt.iter_mut().enumerate().take(12) {
            for (j, slot) in row.iter_mut().enumerate() {
                let off = 2696 + (i * 3 + j) * 4;
                *slot = read_i32(off).max(0) as usize;
            }
        }
        for (j, slot) in ipt[12].iter_mut().enumerate() {
            let off = 2844 + j * 4;
            *slot = read_i32(off).max(0) as usize;
        }

        // Derive record size (`ncoeff` f64s) from maximum index in `ipt`.
        let mut ncoeff = 2usize;
        for (idx, entry) in ipt.iter().enumerate() {
            if entry[0] > 0 && entry[1] > 0 && entry[2] > 0 {
                let ndim: usize = if idx == 11 { 2 } else { 3 };
                let end = (entry[0] - 1)
                    .checked_add(
                        ndim.checked_mul(entry[1])
                            .and_then(|v| v.checked_mul(entry[2]))
                            .ok_or_else(|| {
                                CoreError::InvalidEphemeris("JPL DE ipt table overflow".into())
                            })?,
                    )
                    .ok_or_else(|| {
                        CoreError::InvalidEphemeris("JPL DE ipt table overflow".into())
                    })?;
                if end > ncoeff {
                    ncoeff = end;
                }
            }
        }

        // Verify EMB (2), Moon (9), and Sun (10) pointer entries are valid.
        for &req_body in &[2usize, 9usize, 10usize] {
            let [start_1, ncf, nsub] = ipt[req_body];
            if start_1 == 0 || ncf == 0 || nsub == 0 {
                return Err(CoreError::InvalidEphemeris(format!(
                    "JPL DE{numde} missing required body {req_body} ipt entry"
                )));
            }
        }

        let rec_bytes = ncoeff
            .checked_mul(8)
            .ok_or_else(|| CoreError::InvalidEphemeris("JPL DE record size overflow".into()))?;
        let min_file_bytes = rec_bytes
            .checked_mul(3)
            .ok_or_else(|| CoreError::InvalidEphemeris("JPL DE file size overflow".into()))?;
        if data.len() < min_file_bytes
            || !span_days.is_finite()
            || span_days <= 0.0
            || !jd_start.is_finite()
            || !jd_end.is_finite()
            || jd_end <= jd_start
            || !emrat.is_finite()
            || emrat <= 0.0
        {
            return Err(CoreError::InvalidEphemeris(format!(
                "invalid JPL DE{numde} geometry (ncoeff={ncoeff}, span={span_days})"
            )));
        }

        Ok(Self {
            numde,
            jd_start,
            jd_end,
            span_days,
            au_km,
            emrat,
            ipt,
            ncoeff,
            big_endian,
            data_record_offset: 2,
            data,
        })
    }

    pub fn earth_and_sun_state(&self, jd_hi: f64, jd_lo: f64) -> Result<SolarSystemState> {
        let jd = jd_hi + jd_lo;
        if !jd.is_finite() || jd < self.jd_start || jd > self.jd_end {
            return Err(CoreError::TimeOutOfRange(jd));
        }
        let rec_bytes = self.ncoeff * 8;
        let max_data_rec = (self.data.len() / rec_bytes).saturating_sub(self.data_record_offset + 1);
        let raw_rec_idx = ((jd - self.jd_start) / self.span_days).floor() as usize;
        // Clamp exact `jd == jd_end` boundary to the last record (matching dpleph.c:182).
        let rec_idx = if jd == self.jd_end {
            raw_rec_idx.min(max_data_rec)
        } else {
            raw_rec_idx
        };
        let byte_start = rec_idx
            .checked_add(self.data_record_offset)
            .and_then(|r| r.checked_mul(rec_bytes))
            .ok_or(CoreError::TimeOutOfRange(jd))?;
        let byte_end = byte_start
            .checked_add(rec_bytes)
            .ok_or(CoreError::TimeOutOfRange(jd))?;
        if byte_end > self.data.len() {
            return Err(CoreError::TimeOutOfRange(jd));
        }

        let rec_slice = &self.data[byte_start..byte_end];
        let mut coeffs = Vec::with_capacity(self.ncoeff);
        for &b in rec_slice.as_chunks::<8>().0 {
            coeffs.push(if self.big_endian {
                f64::from_be_bytes(b)
            } else {
                f64::from_le_bytes(b)
            });
        }

        let t_rec_start = coeffs[0];
        let t_rec_end = coeffs[1];
        let span = t_rec_end - t_rec_start;
        if !span.is_finite() || span <= 0.0 {
            return Err(CoreError::InvalidEphemeris(format!(
                "invalid JPL DE record span ({t_rec_start}..{t_rec_end})"
            )));
        }
        // Preserve full precision by computing offset from split JD.
        let dt_from_rec = (jd_hi - t_rec_start) + jd_lo;
        let frac = (dt_from_rec / span).clamp(0.0, 1.0);

        let (emb_pos_km, emb_vel_kmd) = self.eval_body(&coeffs, 2, frac, span)?;
        let (moon_pos_km, moon_vel_kmd) = self.eval_body(&coeffs, 9, frac, span)?;
        let (sun_pos_km, _) = self.eval_body(&coeffs, 10, frac, span)?;

        let moon_fac = 1.0 / (1.0 + self.emrat);
        let km_to_m = 1000.0;
        let kmd_to_ms = 1000.0 / SECS_PER_DAY;

        let mut earth_pos_m = [0.0f64; 3];
        let mut earth_vel_ms = [0.0f64; 3];
        let mut sun_pos_m = [0.0f64; 3];
        for k in 0..3 {
            earth_pos_m[k] = (emb_pos_km[k] - moon_fac * moon_pos_km[k]) * km_to_m;
            earth_vel_ms[k] = (emb_vel_kmd[k] - moon_fac * moon_vel_kmd[k]) * kmd_to_ms;
            sun_pos_m[k] = sun_pos_km[k] * km_to_m;
        }

        Ok(SolarSystemState {
            earth_pos_m,
            earth_vel_ms,
            sun_pos_m,
        })
    }

    fn eval_body(
        &self,
        coeffs: &[f64],
        body_idx: usize,
        frac: f64,
        span_days: f64,
    ) -> Result<([f64; 3], [f64; 3])> {
        let [start_1, ncf, nsub] = self.ipt[body_idx];
        if start_1 == 0 || ncf == 0 || nsub == 0 || !span_days.is_finite() || span_days <= 0.0 {
            return Err(CoreError::InvalidEphemeris(format!(
                "invalid JPL DE body {body_idx} descriptor"
            )));
        }
        let sub_f = frac * (nsub as f64);
        let sub_idx = (sub_f.floor() as usize).min(nsub - 1);
        let tc = 2.0 * (sub_f - (sub_idx as f64)) - 1.0; // normalized [-1, 1]
        let vfac = 2.0 * (nsub as f64) / span_days;

        let base = (start_1 - 1)
            .checked_add(
                sub_idx
                    .checked_mul(3)
                    .and_then(|v| v.checked_mul(ncf))
                    .ok_or_else(|| CoreError::InvalidEphemeris("JPL DE offset overflow".into()))?,
            )
            .ok_or_else(|| CoreError::InvalidEphemeris("JPL DE offset overflow".into()))?;
        let end = base
            .checked_add(
                3usize
                    .checked_mul(ncf)
                    .ok_or_else(|| CoreError::InvalidEphemeris("JPL DE offset overflow".into()))?,
            )
            .ok_or_else(|| CoreError::InvalidEphemeris("JPL DE offset overflow".into()))?;
        if end > coeffs.len() {
            return Err(CoreError::InvalidEphemeris(format!(
                "JPL DE body {body_idx} slice [{base}..{end}] exceeds record length {}",
                coeffs.len()
            )));
        }

        let mut pc = vec![0.0f64; ncf];
        let mut vc = vec![0.0f64; ncf];
        pc[0] = 1.0;
        if ncf > 1 {
            pc[1] = tc;
            vc[1] = 1.0;
        }
        let twot = 2.0 * tc;
        for i in 2..ncf {
            pc[i] = twot * pc[i - 1] - pc[i - 2];
            vc[i] = twot * vc[i - 1] + 2.0 * pc[i - 1] - vc[i - 2];
        }

        let mut pos = [0.0f64; 3];
        let mut vel = [0.0f64; 3];
        for k in 0..3 {
            let off = base + k * ncf;
            let mut p = 0.0;
            let mut v = 0.0;
            for i in (0..ncf).rev() {
                p += pc[i] * coeffs[off + i];
                v += vc[i] * coeffs[off + i];
            }
            pos[k] = p;
            vel[k] = v * vfac;
        }
        Ok((pos, vel))
    }
}
