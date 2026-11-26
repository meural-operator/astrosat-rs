//! Stage 2e: `cztevtencor` — In-flight pixel gain correction, secular gain drift,
//! deterministic sub-bin PI randomization, and `[LLD, ULD]` energy bounds cuts.
//!
//! Strictly idempotent via the `GAINCOR = 'YES'` guard.

use serde::{Deserialize, Serialize};

use crate::types::{QuadEvents, GAIN_DRIFT_EPOCH_SEC, PI_CHANNELS, PIXELS_PER_QUAD};

/// Configuration for `cztevtencor`.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EvtEnCorConfig {
    /// Apply CALDB per-pixel Lower Level Discriminator cut (`PI >= LLD[pix]`).
    pub apply_lld: bool,
    /// Apply CALDB per-pixel Upper Level Discriminator cut (`PI <= ULD[pix]`).
    pub apply_uld: bool,
    /// Apply secular time-dependent gain drift factor `tgaincorpar` (default: `true`).
    pub time_gcor: bool,
    /// Deterministic Wichmann-Hill PRNG seed `(ix, iy, iz)` (default: `(5, 11, 17)` matching `cztevtencor.cpp`).
    pub rng_seed: (i32, i32, i32),
}

impl Default for EvtEnCorConfig {
    fn default() -> Self {
        Self {
            apply_lld: true,
            apply_uld: true,
            time_gcor: true,
            rng_seed: (5, 11, 17),
        }
    }
}

/// CALDB per-quadrant gain correction parameters (`GAINCORPAR`, `GAINCORFLAG`) and `[LLD, ULD]`.
#[derive(Debug, Clone)]
pub struct QuadGainCorCaldb {
    pub gaincor_par: Vec<f32>,
    pub gaincor_flag: Vec<i32>,
    pub lld: Vec<i16>,
    pub uld: Vec<i16>,
}

impl Default for QuadGainCorCaldb {
    fn default() -> Self {
        Self {
            gaincor_par: vec![1.0; PIXELS_PER_QUAD],
            gaincor_flag: vec![1; PIXELS_PER_QUAD],
            lld: vec![0; PIXELS_PER_QUAD],
            uld: vec![(PI_CHANNELS - 1) as i16; PIXELS_PER_QUAD],
        }
    }
}

/// Wichmann-Hill PRNG matching `cztevtencor::urand()` & `cztbindata::urand()`
/// (`Rand_ix=5, Rand_iy=11, Rand_iz=17`, moduli `30269, 30307, 30319`).
#[derive(Debug, Clone)]
pub struct WichmannHillRng {
    ix: i32,
    iy: i32,
    iz: i32,
}

impl WichmannHillRng {
    pub fn new(ix: i32, iy: i32, iz: i32) -> Self {
        Self {
            ix: ix.max(1) % 30269,
            iy: iy.max(1) % 30307,
            iz: iz.max(1) % 30319,
        }
    }

    /// Generate a deterministic uniform variate in `[0.0, 1.0)`.
    #[inline]
    pub fn next_f64(&mut self) -> f64 {
        self.ix = (171 * self.ix) % 30269;
        self.iy = (172 * self.iy) % 30307;
        self.iz = (170 * self.iz) % 30319;
        let temp = (self.ix as f64) / 30269.0
            + (self.iy as f64) / 30307.0
            + (self.iz as f64) / 30319.0;
        temp - temp.floor()
    }
}

/// Compute the CZTI Pipeline v3.0.1 secular time-dependent gain correction factor:
///
/// `tgaincorpar = 0.994 + (0.00014925 / 30.8955) * ((tstart - 181094400.0) / 86400.0)`
#[inline]
pub fn secular_gain_drift_factor(tstart: f64, enabled: bool) -> f32 {
    if !enabled {
        return 1.0;
    }
    let ndays = ((tstart - GAIN_DRIFT_EPOCH_SEC) / 86400.0) as f32;
    0.994 + (0.000_149_25 / 30.8955) * ndays
}

/// Result of applying `evtencor` to a quadrant.
#[derive(Debug, Clone)]
pub struct QuadEvtEnCorResult {
    pub events: QuadEvents,
    pub tgaincorpar: f32,
    pub skipped_already_corrected: bool,
    pub removed_by_lld_uld: usize,
}

/// Apply in-flight `GAINCOR` + secular drift + deterministic sub-bin PI randomization +
/// `[LLD, ULD]` cuts to a quadrant's events (`cztevtencor.cpp`).
///
/// If `already_gain_corrected` is `true` (`GAINCOR == 'YES'` in FITS header), returns the events
/// unchanged (idempotent guard).
pub fn evtencor_quad(
    events: &QuadEvents,
    caldb: &QuadGainCorCaldb,
    head_tstart: f64,
    already_gain_corrected: bool,
    cfg: &EvtEnCorConfig,
) -> QuadEvtEnCorResult {
    if already_gain_corrected {
        return QuadEvtEnCorResult {
            events: events.clone(),
            tgaincorpar: 1.0,
            skipped_already_corrected: true,
            removed_by_lld_uld: 0,
        };
    }

    let tgaincorpar = secular_gain_drift_factor(head_tstart, cfg.time_gcor);
    let mut out = events.clone();
    let n = out.len();
    if n == 0 || out.pi.len() != n || out.energy.len() != n {
        return QuadEvtEnCorResult {
            events: out,
            tgaincorpar,
            skipped_already_corrected: false,
            removed_by_lld_uld: 0,
        };
    }

    let mut rng = WichmannHillRng::new(cfg.rng_seed.0, cfg.rng_seed.1, cfg.rng_seed.2);
    let mut keep = vec![true; n];
    let mut removed = 0usize;

    for j in 0..n {
        let dpix = (out.detid[j] as usize) * 256 + (out.pixid[j] as usize);
        if dpix >= PIXELS_PER_QUAD {
            keep[j] = false;
            removed += 1;
            continue;
        }
        let new_pi: i16 = if caldb.gaincor_flag[dpix] == 1 {
            let u = rng.next_f64();
            let e_cor = (((out.pi[j] as f64) + u) * 0.5 + 5.0)
                * (caldb.gaincor_par[dpix] as f64)
                * (tgaincorpar as f64);
            let e_cor_f32 = e_cor as f32;
            out.energy[j] = e_cor_f32;
            let pi_val = (((e_cor_f32 - 5.0) / 0.5) as i32).clamp(0, (PI_CHANNELS - 1) as i32);
            pi_val as i16
        } else {
            out.pi[j]
        };
        out.pi[j] = new_pi;

        let fail_lld = cfg.apply_lld && new_pi < caldb.lld[dpix];
        let fail_uld = cfg.apply_uld && new_pi > caldb.uld[dpix];
        if fail_lld || fail_uld {
            keep[j] = false;
            removed += 1;
        }
    }

    let filtered = if cfg.apply_lld || cfg.apply_uld {
        out.filter_by_mask(&keep)
    } else {
        out
    };

    QuadEvtEnCorResult {
        events: filtered,
        tgaincorpar,
        skipped_already_corrected: false,
        removed_by_lld_uld: removed,
    }
}
