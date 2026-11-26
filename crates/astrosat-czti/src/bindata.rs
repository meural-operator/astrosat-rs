//! Stage 3c: `cztbindata` & `cztrspgen` — `radec2CZTI`, coded-mask ray-tracing,
//! `BKGSCAL`-renormalized mask weights, OGIP `.pha`/`.lc`/`.rsp`, and mask-weighted `.evt`.

use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use crate::error::{CztiError, Result};
use crate::evtencor::{secular_gain_drift_factor, QuadGainCorCaldb, WichmannHillRng};
use crate::types::{
    detid_pixid_to_detxy, BadpixTable, GtiType, QuadEvents, QuadLivetime, MASK_HEIGHT_MM,
    MAX_FOV_DEG, NUM_DET_PER_QUAD, PI_CHANNELS, PIXELS_PER_QUAD, QUAD_DIM,
};

/// Oversampled compressed mask dimensions (`cztbindata.h`).
pub const COMP_COLS: usize = 327;
pub const TOTAL_MASK_ROWS: usize = 8175;
pub const MODULE_MASK_CELLS: usize = 1950;
/// Oversampled mask grid resolution in mm (20 microns).
pub const MASK_CELL_MM: f64 = 0.02;

/// Per-module additional LLD floor (PI channels) applied in `cztbindata` & `cztrspgen`
/// (`cztbindata.h:107` & `cztbindata.cpp:1930-1931`).
pub const ADDITIONAL_LLD: [[i16; 16]; 4] = [
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 58, 58, 58, 0, 62, 60, 60, 0, 54],
    [0, 0, 0, 58, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 50],
];

/// Clamp per-pixel CALDB `LLD` values against `ADDITIONAL_LLD[quad_id][det_id]`
/// (`cztbindata.cpp:1930-1931`).
pub fn apply_additional_lld(quad_id: usize, lld: &mut [i16]) {
    if quad_id >= 4 {
        return;
    }
    let floors = &ADDITIONAL_LLD[quad_id];
    for (i, slot) in lld.iter_mut().take(PIXELS_PER_QUAD).enumerate() {
        let detnum = i / 256;
        if floors[detnum] > *slot {
            *slot = floors[detnum];
        }
    }
}

/// Measured per-module boresight X offsets (`deg`) from CZTI Pipeline v3.0.1 (`cztbindata.h`).
pub const DETWISE_TX_SHIFT: [[f64; 16]; 4] = [
    [
        0.030, 0.010, 0.030, 0.020, 0.025, 0.020, 0.050, 0.020, -0.030, -0.010, 0.010, 0.000,
        0.005, -0.005, -0.035, -0.010,
    ],
    [
        -0.120, -0.095, -0.130, -0.140, -0.115, -0.140, -0.125, -0.145, -0.125, -0.136, -0.145,
        -0.155, -0.135, -0.150, -0.155, -0.170,
    ],
    [
        -0.065, -0.065, -0.100, -0.050, -0.120, -0.080, -0.085, -0.089, -0.055, -0.130, -0.140,
        -0.110, -0.125, -0.105, -0.045, -0.065,
    ],
    [
        0.070, -0.045, 0.060, 0.020, 0.000, 0.055, 0.040, -0.030, 0.000, 0.020, 0.075, 0.035,
        -0.030, -0.015, 0.010, 0.050,
    ],
];

/// Measured per-module boresight Y offsets (`deg`) from CZTI Pipeline v3.0.1 (`cztbindata.h`).
pub const DETWISE_TY_SHIFT: [[f64; 16]; 4] = [
    [
        -0.020, -0.010, -0.045, -0.050, -0.020, -0.040, -0.035, -0.055, -0.020, -0.055, -0.020,
        -0.030, -0.037, -0.030, -0.040, -0.045,
    ],
    [
        -0.080, -0.040, -0.070, -0.090, -0.030, -0.060, -0.035, 0.055, -0.030, -0.030, -0.020,
        -0.065, -0.020, -0.020, 0.070, -0.020,
    ],
    [
        0.230, 0.170, 0.120, 0.180, 0.180, 0.200, 0.170, 0.182, 0.175, 0.195, 0.200, 0.195,
        0.200, 0.155, 0.190, 0.170,
    ],
    [
        0.230, 0.200, 0.220, 0.165, 0.140, 0.180, 0.190, 0.170, 0.195, 0.165, 0.130, 0.150,
        0.135, 0.140, 0.120, 0.110,
    ],
];

/// Output product selection for `cztbindata`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum BinDataOutputType {
    Spec,
    Lc,
    Both,
}

impl BinDataOutputType {
    pub fn parse(s: &str) -> Option<Self> {
        match s.trim().to_ascii_lowercase().as_str() {
            "spec" => Some(Self::Spec),
            "lc" => Some(Self::Lc),
            "both" => Some(Self::Both),
            _ => None,
        }
    }
}

/// Configuration for `cztbindata`.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BinDataConfig {
    /// Source RA (degrees). Legacy sentinel: if `source_coords` is `None` and `(0.0, 0.0)`,
    /// falls back to `RA_PNT`/`DEC_PNT` from MKF header.
    pub ra_src_deg: f64,
    /// Source DEC (degrees).
    pub dec_src_deg: f64,
    /// Explicit optional source `(RA_deg, DEC_deg)` override (supports targeting `(0.0, 0.0)`
    /// unambiguously; takes precedence over `ra_src_deg`/`dec_src_deg` when `Some`).
    #[serde(default)]
    pub source_coords: Option<(f64, f64)>,
    /// Apply coded-mask weighting (`true` for background-subtracted coded-mask products).
    pub mask_weight: bool,
    /// Which products to generate (`Spec`, `Lc`, or `Both`).
    pub output_type: BinDataOutputType,
    /// Light curve time bin size in seconds (default: `1.0`).
    pub time_binsize: f64,
    /// Light curve lower energy bound in keV (default: `20.0`).
    pub emin_kev: f32,
    /// Light curve upper energy bound in keV (default: `100.0`).
    pub emax_kev: f32,
    /// Generate `.rsp` response matrix alongside `.pha` spectrum (default: `true`).
    pub gen_resp: bool,
    /// Generate mask-weighted event file with `OPEN_FRACTION` and `WEIGHT` columns.
    pub generate_eventfile: bool,
    /// Maximum bad-pixel flag accepted (default: `0`).
    pub badpix_threshold: u8,
    /// Apply LLD cut (always enabled when `mask_weight == true`).
    pub apply_lld: bool,
    /// Apply ULD cut (always enabled when `mask_weight == true`).
    pub apply_uld: bool,
    /// Apply default CZTI spectral channel grouping flag.
    pub do_grouping: bool,
    /// User-supplied boresight offset adjustments in degrees.
    pub delta_tx: f64,
    pub delta_ty: f64,
    /// Minimum active pixels in a 16x16 module required to renormalize mask weights (default: `100`).
    pub min_active_pix_per_module: usize,
    /// Minimum `FRACEXP` to keep a lightcurve time bin (default: `0.05`).
    pub min_fracexp: f32,
    /// GTI mode (`Quad` or `Common`).
    pub gti_type: GtiType,
    /// Quadrants to process (default: `[0, 1, 2, 3]`).
    pub quads: Vec<usize>,
}

impl Default for BinDataConfig {
    fn default() -> Self {
        Self {
            ra_src_deg: 0.0,
            dec_src_deg: 0.0,
            source_coords: None,
            mask_weight: true,
            output_type: BinDataOutputType::Both,
            time_binsize: 1.0,
            emin_kev: 20.0,
            emax_kev: 100.0,
            gen_resp: true,
            generate_eventfile: false,
            badpix_threshold: 0,
            apply_lld: true,
            apply_uld: true,
            do_grouping: true,
            delta_tx: 0.0,
            delta_ty: 0.0,
            min_active_pix_per_module: 100,
            min_fracexp: 0.05,
            gti_type: GtiType::Quad,
            quads: vec![0, 1, 2, 3],
        }
    }
}

impl BinDataConfig {
    /// Resolve target `(RA_deg, DEC_deg)` from explicit `source_coords`, legacy `(ra_src_deg, dec_src_deg)`,
    /// or MKF pointing `(ra_pnt, dec_pnt)`.
    pub fn resolve_source_radec(&self, ra_pnt: f64, dec_pnt: f64) -> (f64, f64) {
        if let Some((ra, dec)) = self.source_coords {
            if ra.is_finite() && dec.is_finite() {
                return (ra, dec);
            }
        }
        if self.ra_src_deg.is_finite()
            && self.dec_src_deg.is_finite()
            && (self.ra_src_deg != 0.0 || self.dec_src_deg != 0.0)
        {
            (self.ra_src_deg, self.dec_src_deg)
        } else {
            (ra_pnt, dec_pnt)
        }
    }
}

// ---------------------------------------------------------------------------
// 1. Exact `radec2CZTI.f` camera coordinate transformations
// ---------------------------------------------------------------------------

#[inline]
fn dotp_norm(v1: [f64; 3], v2: [f64; 3]) -> f64 {
    let m1 = (v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2]).sqrt();
    let m2 = (v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2]).sqrt();
    if m1 <= 0.0 || m2 <= 0.0 {
        0.0
    } else {
        (v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]) / (m1 * m2)
    }
}

/// Convert AstroSat `(roll_RA, roll_Dec, roll_Rot)` (degrees) to `(roll, pitch, yaw)` `(RA, Dec)`
/// in degrees (`radec2CZTI.f::roll_rot2rpy`).
pub fn roll_rot_to_rpy(
    roll_ra_deg: f64,
    roll_dec_deg: f64,
    roll_rot_deg: f64,
) -> ([f64; 2], [f64; 2], [f64; 2]) {
    let ra = roll_ra_deg.to_radians();
    let dec = roll_dec_deg.to_radians();
    let rot = roll_rot_deg.to_radians();

    let roll = [roll_ra_deg, roll_dec_deg];
    let alr = ra.cos() * dec.cos();
    let amr = ra.sin() * dec.cos();
    let anr = dec.sin();

    let amp = -rot.cos() * dec.sin() * ra.sin() - rot.sin() * ra.cos();
    let alp = -rot.cos() * dec.sin() * ra.cos() + rot.sin() * ra.sin();
    let mut anp = (1.0 - alp * alp - amp * amp).max(0.0).sqrt();
    let temp = alp * alr + amp * amr;
    anp = -anp.copysign(temp);
    if dec < 0.0 {
        anp = -anp;
    }

    let pitch_dec = anp.clamp(-1.0, 1.0).asin().to_degrees();
    let pitch_ra = if alp.abs() < 1e-15 && amp.abs() < 1e-15 {
        0.0
    } else {
        let mut r = amp.atan2(alp).to_degrees();
        if r < 0.0 {
            r += 360.0;
        }
        r
    };
    let pitch = [pitch_ra, pitch_dec];

    // Yaw = Roll x Pitch
    let aly = amr * anp - anr * amp;
    let amy = anr * alp - alr * anp;
    let any = alr * amp - amr * alp;

    let yaw_dec = any.clamp(-1.0, 1.0).asin().to_degrees();
    let yaw_ra = if aly.abs() < 1e-15 && amy.abs() < 1e-15 {
        0.0
    } else {
        let mut r = amy.atan2(aly).to_degrees();
        if r < 0.0 {
            r += 360.0;
        }
        r
    };
    let yaw = [yaw_ra, yaw_dec];

    (roll, pitch, yaw)
}

/// Convert attitude quaternion `q = [q1, q2, q3, q4]` to `(roll, pitch, yaw)` `(RA, Dec)` in
/// degrees (`radec2CZTI.f::quat2rpy`).
pub fn quat_to_rpy(q: [f64; 4]) -> ([f64; 2], [f64; 2], [f64; 2]) {
    let q11 = q[0] * q[0];
    let q22 = q[1] * q[1];
    let q33 = q[2] * q[2];
    let q44 = q[3] * q[3];
    let q12 = q[0] * q[1];
    let q13 = q[0] * q[2];
    let q23 = q[1] * q[2];
    let q34 = q[2] * q[3];
    let q24 = q[1] * q[3];
    let q14 = q[0] * q[3];

    let r11 = q44 + q11 - q22 - q33;
    let r12 = 2.0 * (q12 + q34);
    let r13 = 2.0 * (q13 - q24);
    let r21 = 2.0 * (q12 - q34);
    let r22 = q44 - q11 + q22 - q33;
    let r23 = 2.0 * (q23 + q14);
    let r31 = 2.0 * (q13 + q24);
    let r32 = 2.0 * (q23 - q14);
    let r33 = q44 - q11 - q22 + q33;

    let to_ra_dec = |a: f64, b: f64, c: f64| -> [f64; 2] {
        let dec = c.clamp(-1.0, 1.0).asin().to_degrees();
        let ra = if a.abs() < 1e-15 && b.abs() < 1e-15 {
            0.0
        } else {
            let mut v = b.atan2(a).to_degrees();
            if v < 0.0 {
                v += 360.0;
            }
            v
        };
        [ra, dec]
    };

    let yaw = to_ra_dec(r11, r12, r13);
    let roll = to_ra_dec(r21, r22, r23);
    let pitch = to_ra_dec(r31, r32, r33);
    (roll, pitch, yaw)
}

/// Compute CZTI camera coordinates `(tx_deg, ty_deg)` for a source `(ra_src_deg, dec_src_deg)`
/// given AstroSat pointing `(roll_ra_deg, roll_dec_deg, roll_rot_deg)` (`radec2CZTI.f::fcompute_tx_ty`).
pub fn radec_to_czti_tx_ty(
    roll_ra_deg: f64,
    roll_dec_deg: f64,
    roll_rot_deg: f64,
    ra_src_deg: f64,
    dec_src_deg: f64,
) -> (f64, f64) {
    const CZTI_PCS_X: [f64; 3] = [0.0, 0.0, -1.0];
    const CZTI_PCS_Y: [f64; 3] = [0.0, 1.0, 0.0];
    const CZTI_PCS_Z: [f64; 3] = [1.0, 0.0, 0.0];

    let (radec_r, radec_p, radec_y) = roll_rot_to_rpy(roll_ra_deg, roll_dec_deg, roll_rot_deg);

    let to_vec = |rd: [f64; 2]| -> [f64; 3] {
        let ra = rd[0].to_radians();
        let dec = rd[1].to_radians();
        [dec.cos() * ra.cos(), dec.cos() * ra.sin(), dec.sin()]
    };

    let roll = to_vec(radec_r);
    let pitch = to_vec(radec_p);
    let yaw = to_vec(radec_y);
    let src = to_vec([ra_src_deg, dec_src_deg]);

    let src_mid = [
        dotp_norm(roll, src),
        dotp_norm(pitch, src),
        dotp_norm(yaw, src),
    ];
    let src_fin = [
        dotp_norm(CZTI_PCS_X, src_mid),
        dotp_norm(CZTI_PCS_Y, src_mid),
        dotp_norm(CZTI_PCS_Z, src_mid),
    ];

    let mut tx = src_fin[0].atan2(src_fin[2]).to_degrees();
    let mut ty = src_fin[1].atan2(src_fin[2]).to_degrees();
    if tx.abs() < 1e-10 {
        tx = 0.0;
    }
    if ty.abs() < 1e-10 {
        ty = 0.0;
    }
    (tx, ty)
}

/// Validate that `(tx_deg, ty_deg)` lies within the CZTI coded-mask FOV (`|tx|, |ty| <= 5.8 deg`).
pub fn validate_fov(tx_deg: f64, ty_deg: f64) -> Result<()> {
    if !tx_deg.is_finite()
        || !ty_deg.is_finite()
        || tx_deg.abs() > MAX_FOV_DEG
        || ty_deg.abs() > MAX_FOV_DEG
    {
        return Err(CztiError::OutOfFov {
            tx: tx_deg,
            ty: ty_deg,
        });
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// 2. Coded-mask ray-tracing & BKGSCAL-renormalized mask weights
// ---------------------------------------------------------------------------

/// Bit lookup into the CALDB compressed 20-micron oversampled mask (`cztshadowgen.cpp::checkbit`).
#[inline]
pub fn check_compressed_mask_bit(mask_pattern: &[i32], i: usize, j: usize) -> bool {
    let jj = j / 25;
    let bit_no = j % 25;
    let idx = i * COMP_COLS + jj;
    if idx >= mask_pattern.len() {
        return false;
    }
    (mask_pattern[idx] & (1i32 << bit_no)) != 0
}

/// Ray-trace a quadrant's 64x64 pixel open fractions (`open_frac[dety * 64 + detx]`) and raw mask
/// weights (`2 * open_frac - 1`) through the compressed 20 um oversampled mask at height `478 mm`,
/// with `2.46 mm` interior / `2.28 mm` edge pixel pitch and `2.5 mm` module gap
/// (`cztshadowgen.cpp::calculateTrans` + `getShadow`).
pub fn raytrace_quad_shadow(
    tx_deg: f64,
    ty_deg: f64,
    mask_pattern: &[i32],
    tx_shift: &[f64; 16],
    ty_shift: &[f64; 16],
) -> (Vec<f32>, Vec<f32>) {
    let mut open_frac = vec![0.0f32; PIXELS_PER_QUAD];
    let mut raw_weights = vec![0.0f32; PIXELS_PER_QUAD];

    // Precompute module-local lower_left positions along both axes (mm)
    let mut pos_mm = [0.0f64; QUAD_DIM];
    let mut cells = [0usize; QUAD_DIM];
    let mut cur = 0.0f64;
    for blk in 0..4 {
        for p in 0..16 {
            let idx = blk * 16 + p;
            pos_mm[idx] = cur;
            let is_edge = p == 0 || p == 15;
            cells[idx] = if is_edge { 114 } else { 123 };
            cur += if is_edge { 2.28 } else { 2.46 };
        }
        cur += 2.5;
    }

    // Parallelize across the 64 rows `j` (x-axis in cztshadowgen convention)
    let trans_rows: Vec<(usize, Vec<f32>)> = (0..QUAD_DIM)
        .into_par_iter()
        .map(|j| {
            let yy = (j / 16) * 16;
            let yapitch = ((j / 16) as f64) * 2.5;
            let starty = ((MODULE_MASK_CELLS * (yy / 16)) as f64) + (yapitch / MASK_CELL_MM);
            let endy = ((MODULE_MASK_CELLS * (yy / 16 + 1)) as f64) + (yapitch / MASK_CELL_MM);
            let row_counter = cells[j];
            let lower_left_x = pos_mm[j];

            let mut col_vals = vec![0.0f32; QUAD_DIM];
            for i in 0..QUAD_DIM {
                let xx = (i / 16) * 16;
                let xpitch = ((i / 16) as f64) * 2.5;
                let startx = ((MODULE_MASK_CELLS * (xx / 16)) as f64) + (xpitch / MASK_CELL_MM);
                let endx = ((MODULE_MASK_CELLS * (xx / 16 + 1)) as f64) + (xpitch / MASK_CELL_MM);
                let col_counter = cells[i];
                let lower_left_y = pos_mm[i];

                let mod_idx = (j / 16) * 4 + (i / 16);
                let theta_x = (ty_deg + ty_shift[mod_idx]).to_radians();
                let theta_y = (tx_deg + tx_shift[mod_idx]).to_radians();

                let x_mask = (lower_left_y + MASK_HEIGHT_MM * theta_x.tan()) / MASK_CELL_MM;
                let y_mask = (lower_left_x + MASK_HEIGHT_MM * theta_y.tan()) / MASK_CELL_MM;

                let jj0 = y_mask as isize;
                let ii0 = x_mask as isize;
                let sy = starty.ceil() as isize;
                let ey = endy.ceil() as isize;
                let sx = startx.ceil() as isize;
                let ex = endx.ceil() as isize;

                let jj_lo = jj0.max(sy).max(0);
                let jj_hi = (jj0 + row_counter as isize).min(ey);
                let ii_lo = ii0.max(sx).max(0);
                let ii_hi = (ii0 + col_counter as isize).min(ex);

                let mut total_open = 0u32;
                if jj_lo < jj_hi && ii_lo < ii_hi {
                    for jj in jj_lo..jj_hi {
                        let jju = jj as usize;
                        for ii in ii_lo..ii_hi {
                            if check_compressed_mask_bit(mask_pattern, ii as usize, jju) {
                                total_open += 1;
                            }
                        }
                    }
                }
                col_vals[i] = (total_open as f32) / ((col_counter * row_counter) as f32);
            }
            (j, col_vals)
        })
        .collect();

    // `transValues[i * 64 + j]` -> `open_frac[i * 64 + j]` (indexed by dety * 64 + detx)
    for (j, col_vals) in trans_rows {
        for (i, &frac) in col_vals.iter().enumerate() {
            let seq = i * QUAD_DIM + j;
            open_frac[seq] = frac;
            raw_weights[seq] = 2.0 * frac - 1.0;
        }
    }

    (open_frac, raw_weights)
}

/// Ray-trace a quadrant shadow from a 64x64 binary mask pattern (synthetic/fallback helper when
/// the 20-micron oversampled CALDB mask is not supplied).
pub fn raytrace_quad_shadow_simple(
    tx_deg: f64,
    ty_deg: f64,
    mask_64x64: &[Vec<f32>],
) -> (Vec<f32>, Vec<f32>) {
    let mut open_frac = vec![0.0f32; PIXELS_PER_QUAD];
    let mut raw_weights = vec![0.0f32; PIXELS_PER_QUAD];
    let shift_x = (MASK_HEIGHT_MM * tx_deg.to_radians().tan()) / 2.46;
    let shift_y = (MASK_HEIGHT_MM * ty_deg.to_radians().tan()) / 2.46;

    for dety in 0..QUAD_DIM {
        for detx in 0..QUAD_DIM {
            let mx = (detx as f64 + shift_x).round() as isize;
            let my = (dety as f64 + shift_y).round() as isize;
            let f = if (0..QUAD_DIM as isize).contains(&mx) && (0..QUAD_DIM as isize).contains(&my)
            {
                mask_64x64[my as usize][mx as usize].clamp(0.0, 1.0)
            } else {
                0.0
            };
            let idx = dety * QUAD_DIM + detx;
            open_frac[idx] = f;
            raw_weights[idx] = 2.0 * f - 1.0;
        }
    }
    (open_frac, raw_weights)
}

/// Renormalized energy-dependent mask weights for one quadrant (`4096 * 512` floats indexed by
/// `(dety * 64 + detx) * 512 + pi`).
#[derive(Debug, Clone)]
pub struct RenormalizedMaskWeights {
    pub open_frac: Vec<f32>,
    pub raw_weights: Vec<f32>,
    /// Flattened `4096 * 512` renormalized weights.
    pub weights_512: Vec<f32>,
}

impl RenormalizedMaskWeights {
    #[inline]
    pub fn weight(&self, detx: usize, dety: usize, pi: usize) -> f32 {
        let idx = (dety * QUAD_DIM + detx) * PI_CHANNELS + pi;
        self.weights_512[idx]
    }
}

/// Compute `BKGSCAL`-renormalized module-wise mask weights (`cztshadowgen.cpp::calculate_renormalized_weights`).
///
/// - `bkgscal`: flattened `4096 * 512` array indexed as `(dety * 64 + detx) * 512 + pi`
/// - `lld_detmap` / `uld_detmap`: `64 x 64` PI channel bounds indexed as `[detx][dety]`
/// - `pix_exp_frac`: 4096 pixel exposure fractions indexed as `detx * 64 + dety`
#[allow(clippy::too_many_arguments)]
pub fn calculate_renormalized_weights(
    quad_id: usize,
    open_frac: &[f32],
    raw_weights_in: &[f32],
    badpix: &BadpixTable,
    badpix_threshold: u8,
    bkgscal: &[f32],
    lld_detmap: &[[i16; QUAD_DIM]; QUAD_DIM],
    uld_detmap: &[[i16; QUAD_DIM]; QUAD_DIM],
    pix_exp_frac: &[f64],
    min_npix: usize,
) -> RenormalizedMaskWeights {
    let mut raw_weights = raw_weights_in.to_vec();
    let mut quad_bp = [[0u8; QUAD_DIM]; QUAD_DIM];
    for detx in 0..QUAD_DIM {
        for dety in 0..QUAD_DIM {
            let flag = badpix.get_detxy(quad_id, detx as u8, dety as u8);
            quad_bp[detx][dety] = flag;
            if flag > badpix_threshold {
                raw_weights[detx + dety * QUAD_DIM] = 0.0;
            }
        }
    }

    let mut new_weights = vec![0.0f32; PIXELS_PER_QUAD * PI_CHANNELS];

    // Renormalize independently per 16x16 detector module (4x4 grid)
    for mi in 0..4 {
        for mj in 0..4 {
            let mut det_renorm_d = [0.0f64; PI_CHANNELS];
            let mut det_renorm_a = [0.0f64; PI_CHANNELS];
            let mut n_act_pix = [0usize; PI_CHANNELS];
            let mut f2a_sum = 0.0f64;
            let mut fa_sum = 0.0f64;

            for detx in (mi * 16)..((mi + 1) * 16) {
                for dety in (mj * 16)..((mj + 1) * 16) {
                    if quad_bp[detx][dety] > badpix_threshold {
                        continue;
                    }
                    let index = detx + dety * QUAD_DIM;
                    let pimin = (lld_detmap[detx][dety].max(0) as usize).min(PI_CHANNELS - 1);
                    let pimax = (uld_detmap[detx][dety].max(0) as usize).min(PI_CHANNELS - 1);
                    let w = raw_weights[index] as f64;
                    let f = open_frac[index] as f64;
                    let t_exp = pix_exp_frac[detx * QUAD_DIM + dety];

                    for pich in pimin..=pimax {
                        let b = bkgscal[index * PI_CHANNELS + pich] as f64;
                        det_renorm_d[pich] += w * b * t_exp;
                        det_renorm_a[pich] += b * t_exp;
                        f2a_sum += f * f;
                        fa_sum += f;
                        n_act_pix[pich] += 1;
                    }
                }
            }

            let mut denom = [1.0f64; PI_CHANNELS];
            let ratio = if fa_sum != 0.0 {
                2.0 * f2a_sum / fa_sum
            } else {
                1.0
            };
            for pich in 0..PI_CHANNELS {
                if det_renorm_a[pich] != 0.0 {
                    det_renorm_d[pich] /= det_renorm_a[pich];
                } else {
                    det_renorm_d[pich] = 0.0;
                }
                let d = ratio - (1.0 + det_renorm_d[pich]);
                denom[pich] = if d.abs() < 1e-12 { 1.0 } else { d };
            }

            for detx in (mi * 16)..((mi + 1) * 16) {
                for dety in (mj * 16)..((mj + 1) * 16) {
                    if quad_bp[detx][dety] > badpix_threshold {
                        continue;
                    }
                    let index = detx + dety * QUAD_DIM;
                    let pimin = (lld_detmap[detx][dety].max(0) as usize).min(PI_CHANNELS - 1);
                    let pimax = (uld_detmap[detx][dety].max(0) as usize).min(PI_CHANNELS - 1);
                    let w = raw_weights[index] as f64;

                    for pich in pimin..=pimax {
                        if n_act_pix[pich] > min_npix {
                            new_weights[index * PI_CHANNELS + pich] =
                                ((w - det_renorm_d[pich]) / denom[pich]) as f32;
                        }
                    }
                }
            }
        }
    }

    RenormalizedMaskWeights {
        open_frac: open_frac.to_vec(),
        raw_weights,
        weights_512: new_weights,
    }
}

// ---------------------------------------------------------------------------
// 3. Tantalum coded-mask absorption coefficient (`Exposure.cpp::absco_Ta`)
// ---------------------------------------------------------------------------

/// Linear attenuation coefficient of Tantalum (`cm^-1`) at energy `ekev` (keV),
/// exact piecewise polynomial from `Exposure.cpp::absco_Ta`.
#[inline]
pub fn absco_ta(ekev: f32) -> f32 {
    let rho_ta = 16.65f64;
    let e = ekev as f64;
    if (1.0..9.88).contains(&e) {
        let lge = e.log10();
        let lgabs = 4.09765 + lge * (-1.44692 + lge * (-0.700248));
        (rho_ta * 10.0f64.powf(lgabs)) as f32
    } else if (9.88..11.14).contains(&e) {
        (rho_ta * (754.71 + e * (-51.6789))) as f32
    } else if (11.14..11.68).contains(&e) {
        (rho_ta * (801.833 + e * (-50.0))) as f32
    } else if (11.68..67.42).contains(&e) {
        let lge = e.log10();
        let lgabs = 5.13391 + lge * (-2.5274 + lge * (-0.0274466));
        (rho_ta * 10.0f64.powf(lgabs)) as f32
    } else if (67.42..=1000.0).contains(&e) {
        let lge = e.log10();
        let lgabs = 2.74145 + lge * (2.52457 + lge * (-2.81744 + lge * 0.513496));
        (rho_ta * 10.0f64.powf(lgabs)) as f32
    } else {
        0.0
    }
}

/// Transmission of the 0.5 mm (`0.05 cm`) Tantalum coded mask at energy `ekev`.
#[inline]
pub fn ta_mask_transmission(ekev: f32) -> f64 {
    (-0.05 * (absco_ta(ekev) as f64)).exp()
}

// ---------------------------------------------------------------------------
// 4. Spectrum, Lightcurve, Weighted Event List, and Response Matrix Products
// ---------------------------------------------------------------------------

/// OGIP-compliant CZTI Type-I SPECTRUM product (512 PI channels).
#[derive(Debug, Clone)]
pub struct CztiSpectrum {
    pub quad_id: Option<u8>,
    pub channel: Vec<i32>,
    pub counts: Vec<f32>,
    pub stat_err: Vec<f32>,
    pub sys_err: Vec<f32>,
    pub quality: Vec<i16>,
    pub grouping: Vec<i16>,
    pub exposure: f64,
    pub thetax_deg: f64,
    pub thetay_deg: f64,
    pub mask_weighted: bool,
    pub badpix_threshold: u8,
    pub gti_type: GtiType,
    pub respfile: String,
    pub tstart: f64,
    pub tstop: f64,
}

/// OGIP-compliant CZTI LIGHTCURVE product (`FRACEXP > min_fracexp`, livetime-corrected rate).
#[derive(Debug, Clone)]
pub struct CztiLightCurve {
    pub quad_id: Option<u8>,
    pub time: Vec<f64>,
    pub rate: Vec<f32>,
    pub error: Vec<f32>,
    pub fracexp: Vec<f32>,
    pub timedel: f64,
    pub tstart: f64,
    pub tstop: f64,
    pub emin_kev: f32,
    pub emax_kev: f32,
    pub thetax_deg: f64,
    pub thetay_deg: f64,
    pub mask_weighted: bool,
    pub badpix_threshold: u8,
}

/// OGIP-compliant CZTI RESPONSE (`MATRIX` + `EBOUNDS`, OGIP v1.3.0) product.
#[derive(Debug, Clone)]
pub struct CztiResponse {
    pub quad_id: Option<u8>,
    pub n_ebins: usize,
    pub n_pibins: usize,
    pub energ_lo: Vec<f32>,
    pub energ_hi: Vec<f32>,
    /// Flattened `n_ebins * n_pibins` full RMF*ARF matrix in `cm^2`.
    pub matrix: Vec<f32>,
    pub ebounds_channel: Vec<i16>,
    pub ebounds_emin: Vec<f32>,
    pub ebounds_emax: Vec<f32>,
}

/// In-memory CALDB pixel response bundle needed for fast `rayon`-parallelized `.rsp` generation.
#[derive(Debug, Clone)]
pub struct QuadRspCaldb {
    pub n_ebins: usize,
    pub n_pibins: usize,
    pub emin: f32,
    pub emax: f32,
    /// `pixresp`: `n_groups * n_ebins * n_pibins` flattened redistribution table.
    pub pixresp_matrix: Vec<f32>,
    /// `resp_id[pix_row][temp_bin]` (`4096 x 5`)
    pub resp_id: Vec<[usize; 5]>,
    /// `eff_area`: `4096 * n_ebins`
    pub eff_area: Vec<f32>,
    /// `eff_area_cor1`: `16 * n_pibins` (module-wise PI correction)
    pub eff_area_cor1: Vec<f32>,
    /// `eff_area_cor2`: `n_ebins` (incident-energy correction)
    pub eff_area_cor2: Vec<f32>,
    /// Per-pixel LLD/ULD PI bounds (`4096`).
    pub lld: Vec<i16>,
    pub uld: Vec<i16>,
}

/// Bin a quadrant's events into spectrum counts/errors, lightcurve bins, and optional mask-weighted
/// event list (`cztbindata.cpp`).
#[allow(clippy::too_many_arguments)]
pub fn bindata_quad(
    events: &QuadEvents,
    badpix: &BadpixTable,
    livetime: &QuadLivetime,
    weights: Option<&RenormalizedMaskWeights>,
    gaincor: &QuadGainCorCaldb,
    already_gain_corrected: bool,
    syserr: &[f32; PI_CHANNELS],
    head_tstart: f64,
    head_tstop: f64,
    avg_tx_deg: f64,
    avg_ty_deg: f64,
    respfile: String,
    cfg: &BinDataConfig,
) -> (CztiSpectrum, CztiLightCurve, Option<QuadEvents>) {
    let qid = events.quad_id as usize;
    let tgaincorpar = secular_gain_drift_factor(head_tstart, true);
    let mut rng = WichmannHillRng::new(5, 11, 17);

    let pistart = (((cfg.emin_kev - 5.0) * 2.0) as i32).clamp(0, (PI_CHANNELS - 1) as i32) as usize;
    let piend = (((cfg.emax_kev - 5.0) * 2.0) as i32).clamp(0, (PI_CHANNELS - 1) as i32) as usize;

    let time_binsize = if cfg.time_binsize.is_finite() && cfg.time_binsize > 0.0 {
        cfg.time_binsize
    } else {
        1.0
    };
    let lc_start_t = if head_tstart.is_finite() {
        head_tstart.trunc()
    } else {
        0.0
    };
    let lc_stop_t = if head_tstop.is_finite() && head_tstop >= lc_start_t {
        head_tstop.trunc() + 1.0
    } else {
        lc_start_t + time_binsize
    };
    const MAX_LC_BINS: usize = 10_000_000;
    let raw_bins = (((lc_stop_t - lc_start_t).max(time_binsize)) / time_binsize).ceil();
    let n_lcbins = if raw_bins.is_finite() && raw_bins >= 1.0 {
        (raw_bins as usize).clamp(1, MAX_LC_BINS)
    } else {
        1
    };

    let mut sum_wc = [0.0f64; PI_CHANNELS];
    let mut sum_w2c = [0.0f64; PI_CHANNELS];
    let mut lc_counts = vec![0.0f64; n_lcbins];
    let mut lc_var = vec![0.0f64; n_lcbins];

    let n = events.len();
    let mut evt_open_frac = if cfg.generate_eventfile && cfg.mask_weight {
        Some(vec![0.0f32; n])
    } else {
        None
    };
    let mut evt_weight = if cfg.generate_eventfile && cfg.mask_weight {
        Some(vec![0.0f32; n])
    } else {
        None
    };
    let mut keep_evt = vec![false; n];

    for i in 0..n {
        let detx = events.detx[i] as usize;
        let dety = events.dety[i] as usize;
        if detx >= QUAD_DIM || dety >= QUAD_DIM {
            continue;
        }
        let flag = badpix.get_detxy(qid, detx as u8, dety as u8);
        let detpix_idx = (events.detid[i] as usize) * 256 + (events.pixid[i] as usize);
        if detpix_idx >= PIXELS_PER_QUAD {
            continue;
        }

        if flag > cfg.badpix_threshold {
            if let (Some(of), Some(wt), Some(mw)) = (&mut evt_open_frac, &mut evt_weight, weights) {
                of[i] = mw.open_frac[detx + dety * QUAD_DIM];
                wt[i] = 0.0;
            }
            continue;
        }

        let pi_val: i16 = if already_gain_corrected {
            events.pi[i]
        } else if gaincor.gaincor_flag[detpix_idx] == 1 {
            let u = rng.next_f64();
            let new_e = (((events.pi[i] as f64) + u) * 0.5 + 5.0)
                * (gaincor.gaincor_par[detpix_idx] as f64)
                * (tgaincorpar as f64);
            (((new_e - 5.0) / 0.5) as i32).clamp(0, (PI_CHANNELS - 1) as i32) as i16
        } else {
            events.pi[i]
        };

        let pi_u = (pi_val.max(0) as usize).min(PI_CHANNELS - 1);
        let lld_ok = pi_val >= gaincor.lld[detpix_idx];
        let uld_ok = pi_val <= gaincor.uld[detpix_idx];

        let (pass, w) = if cfg.mask_weight {
            if let Some(mw) = weights {
                let pix_w = mw
                    .weights_512
                    .get((dety * QUAD_DIM + detx) * PI_CHANNELS + pi_u)
                    .copied()
                    .unwrap_or(0.0) as f64;
                if let (Some(of), Some(wt)) = (&mut evt_open_frac, &mut evt_weight) {
                    of[i] = mw.open_frac.get(detx + dety * QUAD_DIM).copied().unwrap_or(0.0);
                    wt[i] = pix_w as f32;
                }
                (lld_ok && uld_ok, pix_w)
            } else {
                (false, 0.0)
            }
        } else {
            let ok = (lld_ok || !cfg.apply_lld) && (uld_ok || !cfg.apply_uld);
            (ok, 1.0)
        };

        if !pass {
            continue;
        }
        keep_evt[i] = true;

        sum_wc[pi_u] += w;
        sum_w2c[pi_u] += w * w;

        if pi_u >= pistart && pi_u <= piend && events.time[i].is_finite() && events.time[i] >= lc_start_t {
            let tbin = ((events.time[i] - lc_start_t) / time_binsize) as usize;
            if tbin < n_lcbins {
                lc_counts[tbin] += w;
                lc_var[tbin] += w * w;
            }
        }
    }

    // Build SPECTRUM
    let channel: Vec<i32> = (0..PI_CHANNELS as i32).collect();
    let counts: Vec<f32> = sum_wc.iter().map(|&v| v as f32).collect();
    let stat_err: Vec<f32> = sum_w2c.iter().map(|&v| v.sqrt() as f32).collect();
    let (quality, grouping) =
        compute_spectrum_quality_and_grouping(&counts, &stat_err, cfg.do_grouping && cfg.mask_weight);

    let spec = CztiSpectrum {
        quad_id: Some(events.quad_id),
        channel,
        counts,
        stat_err,
        sys_err: syserr.to_vec(),
        quality,
        grouping,
        exposure: events.exposure,
        thetax_deg: avg_tx_deg,
        thetay_deg: avg_ty_deg,
        mask_weighted: cfg.mask_weight,
        badpix_threshold: cfg.badpix_threshold,
        gti_type: cfg.gti_type,
        respfile,
        tstart: head_tstart,
        tstop: head_tstop,
    };

    // Compute lightcurve FRACEXP from livetime table
    let mut fracexp = vec![1.0f32; n_lcbins];
    if cfg.gti_type == GtiType::Quad && livetime.time.len() >= 2 {
        let frac_binsize = livetime.time[1] - livetime.time[0];
        let rem = if frac_binsize > 0.0 {
            (time_binsize % frac_binsize).abs()
        } else {
            f64::INFINITY
        };
        if frac_binsize > 0.0
            && time_binsize + 1e-4 >= frac_binsize
            && (rem < 1e-4 || (frac_binsize - rem).abs() < 1e-4)
        {
            let bin_ratio = (time_binsize / frac_binsize).round().max(1.0) as isize;
            let fr_startbin =
                ((lc_start_t - livetime.time[0] - 0.5 * frac_binsize) / frac_binsize).round()
                    as isize
                    + 1;
            if fr_startbin >= 0 {
                let nlv = livetime.fracexp.len() as isize;
                for (k, fe) in fracexp.iter_mut().enumerate() {
                    let mut acc = 0.0f64;
                    let mut cnt = 0usize;
                    let l0 = fr_startbin + (k as isize) * bin_ratio;
                    let l1 = fr_startbin + ((k as isize) + 1) * bin_ratio;
                    for l in l0..l1 {
                        if l >= 0 && l < nlv {
                            acc += livetime.fracexp[l as usize];
                        }
                        cnt += 1;
                    }
                    *fe = if cnt > 0 && acc > 0.0 {
                        (acc / time_binsize).clamp(0.0, 1.0) as f32
                    } else {
                        0.0
                    };
                }
            }
        }
    }

    // Apply livetime correction, divide by binsize, and keep bins with FRACEXP > min_fracexp (> 0)
    let min_fe = cfg.min_fracexp.max(0.0);
    let mut lc_time = Vec::new();
    let mut lc_rate = Vec::new();
    let mut lc_err = Vec::new();
    let mut lc_fe = Vec::new();
    for k in 0..n_lcbins {
        let fe = fracexp[k];
        if fe.is_finite() && fe > min_fe && fe > 0.0 {
            let scale = 1.0 / ((fe as f64) * time_binsize);
            lc_time.push(lc_start_t + (k as f64) * time_binsize + 0.5 * time_binsize);
            lc_rate.push((lc_counts[k] * scale) as f32);
            lc_err.push((lc_var[k].sqrt() * scale) as f32);
            lc_fe.push(fe);
        }
    }

    let lc = CztiLightCurve {
        quad_id: Some(events.quad_id),
        time: lc_time,
        rate: lc_rate,
        error: lc_err,
        fracexp: lc_fe,
        timedel: time_binsize,
        tstart: head_tstart,
        tstop: head_tstop,
        emin_kev: cfg.emin_kev,
        emax_kev: cfg.emax_kev,
        thetax_deg: avg_tx_deg,
        thetay_deg: avg_ty_deg,
        mask_weighted: cfg.mask_weight,
        badpix_threshold: cfg.badpix_threshold,
    };

    let weighted_evt = if cfg.generate_eventfile && cfg.mask_weight {
        let mut wevt = events.clone();
        wevt.open_fraction = evt_open_frac;
        wevt.weight = evt_weight;
        Some(wevt.filter_by_mask(&keep_evt))
    } else {
        None
    };

    (spec, lc, weighted_evt)
}

/// Compute legacy CZTI `.pha` `QUALITY` and SNR-based logarithmic `GROUPING` flags
/// matching `level2validation.cpp:4867-5004`.
pub fn compute_spectrum_quality_and_grouping(
    counts: &[f32],
    stat_err: &[f32],
    do_grouping: bool,
) -> (Vec<i16>, Vec<i16>) {
    let nchan = counts.len();
    let mut grouping = vec![1i16; nchan];
    let mut quality = vec![0i16; nchan];

    let mut tot_cnt = 0.0f32;
    let mut tot_err_sq = 0.0f32;
    let mut fchan: i32 = -1;
    let mut lchan: i32 = -1;

    for i in 0..nchan {
        let flux = counts[i];
        let err = stat_err.get(i).copied().unwrap_or(0.0);
        tot_cnt += flux;
        tot_err_sq += err * err;

        if flux == 0.0 && err == 0.0 {
            if fchan != -1 && lchan == -1 {
                lchan = (i as i32) - 1;
            }
            quality[i] = 5;
        } else {
            if fchan == -1 {
                fchan = i as i32;
            }
            quality[i] = 0;
        }
    }

    if tot_cnt == 0.0 {
        quality.fill(0);
        return (quality, grouping);
    }

    if fchan < 0 {
        fchan = 0;
    }
    if lchan == -1 && nchan > 0 {
        lchan = (nchan as i32) - 1;
    }

    if fchan < 34 {
        let end = 34.min(nchan);
        for slot in quality.iter_mut().take(end).skip(fchan as usize) {
            *slot = 5;
        }
        fchan = 34;
    }

    if lchan > 390 {
        let end = ((lchan as usize) + 1).min(nchan);
        for slot in quality.iter_mut().take(end).skip(391.min(end)) {
            *slot = 5;
        }
        lchan = 390;
    }

    let tot_err = tot_err_sq.sqrt();
    if do_grouping && tot_cnt != 0.0 && tot_err > 0.0 && lchan > fchan {
        let snr = tot_cnt / tot_err;
        let nomgrp_bins = (snr * snr).clamp(10.0, 100.0) as i32;

        let min_e = (fchan as f32) * 0.5 + 5.0;
        let max_e = (lchan as f32) * 0.5 + 5.0;
        let lmin_e = min_e.log10();
        let lgebin = (max_e.log10() - lmin_e) / (nomgrp_bins as f32);

        let mut change_chan: Vec<i32> = Vec::new();
        let mut num_chan: Vec<i32> = Vec::new();
        let mut prev_dc: i32 = -1;
        let mut prev_chan: i32 = -1;

        for i in 0..nomgrp_bins {
            let e1 = 10.0f32.powf(lmin_e + (i as f32) * lgebin);
            let e2 = 10.0f32.powf(lmin_e + ((i + 1) as f32) * lgebin);
            let mut cur_chan = ((e1 - 5.0) / 0.5) as i32;
            cur_chan = cur_chan.clamp(fchan, lchan);

            let de = e2 - e1;
            let dc = ((de / 0.5).round() as i32).max(1);

            if dc != prev_dc {
                if prev_chan != -1 {
                    let rem = (cur_chan - prev_chan) % prev_dc;
                    cur_chan += prev_dc - rem;
                }
                change_chan.push(cur_chan);
                num_chan.push(dc);
                prev_dc = dc;
                prev_chan = cur_chan;
            }
        }

        change_chan.push(lchan + 1);
        num_chan.push(1);

        for i in 0..change_chan.len().saturating_sub(1) {
            let step = num_chan[i].max(1);
            let mut ngrp = (change_chan[i + 1] - change_chan[i]) / step;
            if ngrp == 0 {
                ngrp = 1;
            }
            let start = change_chan[i].max(0) as usize;
            let end = (change_chan[i + 1].max(0) as usize).min(nchan);
            for slot in grouping.iter_mut().take(end).skip(start) {
                *slot = -1;
            }
            for j in 0..ngrp.max(0) {
                let idx = change_chan[i] + j * step;
                if idx >= 0 && (idx as usize) < nchan {
                    grouping[idx as usize] = 1;
                }
            }
        }
    }

    (quality, grouping)
}

/// Generate the OGIP RMF*ARF `.rsp` response matrix for one quadrant in parallel using `rayon`
/// (`cztrspgen.cpp`, zero inner-loop disk I/O!).
#[allow(clippy::too_many_arguments)]
pub fn generate_quad_response(
    quad_id: u8,
    caldb: &QuadRspCaldb,
    badpix: &BadpixTable,
    badpix_threshold: u8,
    avg_open_frac: &[f32],
    renorm_weights: Option<&RenormalizedMaskWeights>,
    pix_exp_frac: &[f64],
    ebounds_emin: Vec<f32>,
    ebounds_emax: Vec<f32>,
) -> CztiResponse {
    let n_ebins = caldb.n_ebins;
    let n_pibins = caldb.n_pibins;
    if n_ebins == 0 || n_pibins == 0 {
        return CztiResponse {
            quad_id: Some(quad_id),
            n_ebins: 0,
            n_pibins: 0,
            energ_lo: Vec::new(),
            energ_hi: Vec::new(),
            matrix: Vec::new(),
            ebounds_channel: Vec::new(),
            ebounds_emin,
            ebounds_emax,
        };
    }
    let ebinsize = (caldb.emax - caldb.emin) / (n_ebins as f32);

    // Precompute Tantalum transmission per incident energy bin
    let tr_ta: Vec<f64> = (0..n_ebins)
        .map(|eb| {
            let e_mid = caldb.emin + ((eb as f32) + 0.5) * ebinsize;
            ta_mask_transmission(e_mid)
        })
        .collect();

    // Parallelize across the 16 detector modules
    let det_matrices: Vec<Vec<f64>> = (0..NUM_DET_PER_QUAD)
        .into_par_iter()
        .map(|det| {
            let mut mat = vec![0.0f64; n_ebins * n_pibins];
            for pix in 0..256usize {
                let row = det * 256 + pix;
                if row >= caldb.resp_id.len()
                    || row >= caldb.lld.len()
                    || row >= caldb.uld.len()
                    || badpix.get(quad_id as usize, det as u8, pix as u8) > badpix_threshold
                {
                    continue;
                }
                let (detx, dety) = detid_pixid_to_detxy(quad_id, det as u8, pix as u8);
                let detx_u = detx as usize;
                let dety_u = dety as usize;
                let openfrac = avg_open_frac
                    .get(detx_u + dety_u * QUAD_DIM)
                    .copied()
                    .unwrap_or(0.0) as f64;
                let pexp = pix_exp_frac
                    .get(detx_u * QUAD_DIM + dety_u)
                    .copied()
                    .unwrap_or(0.0);
                if pexp <= 0.0 || !pexp.is_finite() {
                    continue;
                }

                // Temperature bin 2 (10 C nominal) matching `cztrspgen.cpp` line 601
                let resp_grp = caldb.resp_id[row][2];
                let Some(grp_base) = resp_grp
                    .checked_mul(n_ebins)
                    .and_then(|v| v.checked_mul(n_pibins))
                else {
                    continue;
                };
                if grp_base.saturating_add(n_ebins * n_pibins) > caldb.pixresp_matrix.len() {
                    continue;
                }
                let ea_base = row * n_ebins;
                if ea_base + n_ebins > caldb.eff_area.len() {
                    continue;
                }
                let lld = (caldb.lld[row].max(0) as usize).min(n_pibins - 1);
                let uld = (caldb.uld[row].max(0) as usize).min(n_pibins - 1);
                if lld > uld {
                    continue;
                }

                for eb in 0..n_ebins {
                    let ea = caldb.eff_area[ea_base + eb] as f64;
                    if ea == 0.0 {
                        continue;
                    }
                    let mask_eff = openfrac + tr_ta[eb] * (1.0 - openfrac);
                    let factor = ea * mask_eff * pexp;
                    let src_row = grp_base + eb * n_pibins;
                    let dst_row = eb * n_pibins;

                    if let Some(mw) = renorm_weights {
                        let w_base = (detx_u + dety_u * QUAD_DIM) * PI_CHANNELS;
                        if w_base + uld < mw.weights_512.len() {
                            for pi in lld..=uld {
                                let r = caldb.pixresp_matrix[src_row + pi] as f64;
                                if r != 0.0 {
                                    let w = mw.weights_512[w_base + pi] as f64;
                                    mat[dst_row + pi] += r * factor * w;
                                }
                            }
                        }
                    } else {
                        for pi in lld..=uld {
                            let r = caldb.pixresp_matrix[src_row + pi] as f64;
                            if r != 0.0 {
                                mat[dst_row + pi] += r * factor;
                            }
                        }
                    }
                }
            }
            mat
        })
        .collect();

    // Combine detector matrices with AREACOR1 (det x pi) and AREACOR2 (ebin)
    let mut total_matrix = vec![0.0f32; n_ebins * n_pibins];
    for eb in 0..n_ebins {
        let cor2 = *caldb.eff_area_cor2.get(eb).unwrap_or(&1.0) as f64;
        let row_off = eb * n_pibins;
        for pi in 0..n_pibins {
            let mut sum = 0.0f64;
            for (det, dmat) in det_matrices.iter().enumerate().take(NUM_DET_PER_QUAD) {
                let cor1 = *caldb.eff_area_cor1.get(det * n_pibins + pi).unwrap_or(&1.0) as f64;
                sum += dmat[row_off + pi] * cor1;
            }
            total_matrix[row_off + pi] = (sum * cor2) as f32;
        }
    }

    let energ_lo: Vec<f32> = (0..n_ebins)
        .map(|i| caldb.emin + (i as f32) * ebinsize)
        .collect();
    let energ_hi: Vec<f32> = (0..n_ebins)
        .map(|i| caldb.emin + ((i + 1) as f32) * ebinsize)
        .collect();
    let ebounds_channel: Vec<i16> = (0..n_pibins as i16).collect();

    CztiResponse {
        quad_id: Some(quad_id),
        n_ebins,
        n_pibins,
        energ_lo,
        energ_hi,
        matrix: total_matrix,
        ebounds_channel,
        ebounds_emin,
        ebounds_emax,
    }
}

/// Compute the 5x5 `(tx, ty)` attitude-histogram weighted average open fraction (`cztrspgen.cpp:441-505`)
/// for a quadrant, ray-tracing each populated `(tx, ty)` bin in parallel.
pub fn compute_5x5_avg_open_frac(
    gtisel_tx_ty: &[(f64, f64)],
    mask_ov: Option<&[i32]>,
    mask_64x64: Option<&[Vec<f32>]>,
    tx_shift: &[f64; 16],
    ty_shift: &[f64; 16],
    fallback_open_frac: &[f32],
) -> Vec<f32> {
    if gtisel_tx_ty.len() < 2 {
        return fallback_open_frac.to_vec();
    }
    let mut mintx = f64::INFINITY;
    let mut maxtx = f64::NEG_INFINITY;
    let mut minty = f64::INFINITY;
    let mut maxty = f64::NEG_INFINITY;
    for &(tx, ty) in gtisel_tx_ty {
        if tx.is_finite() && ty.is_finite() {
            mintx = mintx.min(tx);
            maxtx = maxtx.max(tx);
            minty = minty.min(ty);
            maxty = maxty.max(ty);
        }
    }
    if !mintx.is_finite() || !minty.is_finite() {
        return fallback_open_frac.to_vec();
    }

    const NBINS: usize = 5;
    let txbin = (maxtx - mintx + 0.005) / (NBINS as f64);
    let tybin = (maxty - minty + 0.005) / (NBINS as f64);
    if txbin <= 0.0 || tybin <= 0.0 {
        return fallback_open_frac.to_vec();
    }

    let mut sum_tx = [[0.0f64; NBINS]; NBINS];
    let mut sum_ty = [[0.0f64; NBINS]; NBINS];
    let mut cnt = [[0usize; NBINS]; NBINS];

    for &(tx, ty) in gtisel_tx_ty {
        if !tx.is_finite() || !ty.is_finite() {
            continue;
        }
        let ix = (((tx - mintx) / txbin) as isize).clamp(0, (NBINS - 1) as isize) as usize;
        let iy = (((ty - minty) / tybin) as isize).clamp(0, (NBINS - 1) as isize) as usize;
        sum_tx[ix][iy] += tx;
        sum_ty[ix][iy] += ty;
        cnt[ix][iy] += 1;
    }

    let mut populated = Vec::new();
    for ix in 0..NBINS {
        for iy in 0..NBINS {
            let n = cnt[ix][iy];
            if n > 0 {
                populated.push((sum_tx[ix][iy] / (n as f64), sum_ty[ix][iy] / (n as f64), n));
            }
        }
    }
    if populated.is_empty() {
        return fallback_open_frac.to_vec();
    }

    // Ray-trace each populated attitude bin (typically 1..25 bins)
    let bin_shadows: Vec<(Vec<f32>, usize)> = populated
        .into_par_iter()
        .map(|(mtx, mty, weight)| {
            let (of, _) = if let Some(ov) = mask_ov {
                raytrace_quad_shadow(mtx, mty, ov, tx_shift, ty_shift)
            } else if let Some(m64) = mask_64x64 {
                raytrace_quad_shadow_simple(mtx, mty, m64)
            } else {
                (fallback_open_frac.to_vec(), Vec::new())
            };
            (of, weight)
        })
        .collect();

    let mut avg = vec![0.0f64; PIXELS_PER_QUAD];
    let mut ntot = 0.0f64;
    for (of, w) in bin_shadows {
        let wf = w as f64;
        for (dst, &src) in avg.iter_mut().zip(of.iter()) {
            *dst += (src as f64) * wf;
        }
        ntot += wf;
    }
    if ntot > 0.0 {
        avg.into_iter().map(|v| (v / ntot) as f32).collect()
    } else {
        fallback_open_frac.to_vec()
    }
}

/// Co-add multiple CZTI spectra (`cztbindata` `COMMON` mode & `cztaddspec.cpp`).
///
/// - `counts`: summed across spectra
/// - `stat_err`: quadrature sum `sqrt(sum(e_i^2))`
/// - `sys_err`: exposure-weighted average
/// - `exposure`: summed if `add_exposure == true` (multi-observation `cztaddspec`), or averaged if
///   `add_exposure == false` (simultaneous quadrants in `COMMON` GTI mode).
pub fn combine_spectra(
    specs: &[CztiSpectrum],
    add_exposure: bool,
    respfile: String,
) -> CztiSpectrum {
    let Some(first) = specs.first() else {
        return CztiSpectrum {
            quad_id: None,
            channel: (0..PI_CHANNELS as i32).collect(),
            counts: vec![0.0; PI_CHANNELS],
            stat_err: vec![0.0; PI_CHANNELS],
            sys_err: vec![0.0; PI_CHANNELS],
            quality: vec![0; PI_CHANNELS],
            grouping: vec![1; PI_CHANNELS],
            exposure: 0.0,
            thetax_deg: 0.0,
            thetay_deg: 0.0,
            mask_weighted: true,
            badpix_threshold: 0,
            gti_type: GtiType::Common,
            respfile,
            tstart: 0.0,
            tstop: 0.0,
        };
    };

    let n = first.channel.len();
    let mut counts = vec![0.0f64; n];
    let mut var = vec![0.0f64; n];
    let mut sys_num = vec![0.0f64; n];
    let mut exp_sum = 0.0f64;
    let mut tstart = first.tstart;
    let mut tstop = first.tstop;

    for s in specs {
        let w = s.exposure.max(0.0);
        exp_sum += w;
        if s.tstart > 0.0 && (tstart == 0.0 || s.tstart < tstart) {
            tstart = s.tstart;
        }
        if s.tstop > tstop {
            tstop = s.tstop;
        }
        for i in 0..n.min(s.counts.len()) {
            let c = s.counts[i] as f64;
            counts[i] += c;
            let e = *s.stat_err.get(i).unwrap_or(&0.0) as f64;
            var[i] += e * e;
            let sy = *s.sys_err.get(i).unwrap_or(&0.0) as f64;
            // Weight systematic error by per-channel counts (`cztaddspec.cpp:216,238`)
            sys_num[i] += sy * c;
        }
    }

    let exposure = if add_exposure {
        exp_sum
    } else if !specs.is_empty() {
        exp_sum / (specs.len() as f64)
    } else {
        0.0
    };
    let sys_err: Vec<f32> = sys_num
        .into_iter()
        .zip(counts.iter())
        .map(|(num, &tot)| {
            if tot.abs() > 1e-12 {
                (num / tot) as f32
            } else {
                0.0
            }
        })
        .collect();

    let counts_f32: Vec<f32> = counts.into_iter().map(|v| v as f32).collect();
    let stat_err_f32: Vec<f32> = var.into_iter().map(|v| v.sqrt() as f32).collect();
    let do_group =
        add_exposure || (first.mask_weighted && specs.iter().any(|s| s.grouping.contains(&-1)));
    let (quality, grouping) =
        compute_spectrum_quality_and_grouping(&counts_f32, &stat_err_f32, do_group);

    CztiSpectrum {
        quad_id: None,
        channel: first.channel.clone(),
        counts: counts_f32,
        stat_err: stat_err_f32,
        sys_err,
        quality,
        grouping,
        exposure,
        thetax_deg: first.thetax_deg,
        thetay_deg: first.thetay_deg,
        mask_weighted: first.mask_weighted,
        badpix_threshold: first.badpix_threshold,
        gti_type: if add_exposure {
            first.gti_type
        } else {
            GtiType::Common
        },
        respfile,
        tstart,
        tstop,
    }
}

/// Combine multiple CZTI response matrices (`cztbindata` `COMMON` mode & `cztaddspec.cpp`).
///
/// - When `add_exposure == false` (simultaneous quadrants in `COMMON` GTI): responses are summed
///   (total collecting area across combined quadrants).
/// - When `add_exposure == true` (multi-observation `cztaddspec`): exposure-weighted average
///   `sum(exp_k * R_k) / sum(exp_k)`.
pub fn combine_responses(
    rsps: &[CztiResponse],
    exposures: &[f64],
    add_exposure: bool,
) -> CztiResponse {
    let Some(first) = rsps.first() else {
        return CztiResponse {
            quad_id: None,
            n_ebins: 0,
            n_pibins: 0,
            energ_lo: Vec::new(),
            energ_hi: Vec::new(),
            matrix: Vec::new(),
            ebounds_channel: Vec::new(),
            ebounds_emin: Vec::new(),
            ebounds_emax: Vec::new(),
        };
    };

    let len = first.n_ebins * first.n_pibins;
    let mut acc = vec![0.0f64; len];
    let mut tot_exp = 0.0f64;

    for (k, r) in rsps.iter().enumerate() {
        if r.matrix.len() != len {
            continue;
        }
        let w = if add_exposure {
            let e = exposures.get(k).copied().unwrap_or(1.0).max(0.0);
            tot_exp += e;
            e
        } else {
            1.0
        };
        for (dst, &src) in acc.iter_mut().zip(r.matrix.iter()) {
            *dst += (src as f64) * w;
        }
    }

    let inv = if add_exposure && tot_exp > 0.0 {
        1.0 / tot_exp
    } else {
        1.0
    };

    CztiResponse {
        quad_id: None,
        n_ebins: first.n_ebins,
        n_pibins: first.n_pibins,
        energ_lo: first.energ_lo.clone(),
        energ_hi: first.energ_hi.clone(),
        matrix: acc.into_iter().map(|v| (v * inv) as f32).collect(),
        ebounds_channel: first.ebounds_channel.clone(),
        ebounds_emin: first.ebounds_emin.clone(),
        ebounds_emax: first.ebounds_emax.clone(),
    }
}

/// Combine per-quadrant lightcurves into a single summed all-quadrant lightcurve
/// (`cztbindata` `COMMON` GTI mode).
pub fn combine_lightcurves(lcs: &[CztiLightCurve]) -> CztiLightCurve {
    let Some(first) = lcs.first() else {
        return CztiLightCurve {
            quad_id: None,
            time: Vec::new(),
            rate: Vec::new(),
            error: Vec::new(),
            fracexp: Vec::new(),
            timedel: 1.0,
            tstart: 0.0,
            tstop: 0.0,
            emin_kev: 20.0,
            emax_kev: 100.0,
            thetax_deg: 0.0,
            thetay_deg: 0.0,
            mask_weighted: true,
            badpix_threshold: 0,
        };
    };

    // Align bins by rounded bin key `(time / timedel).round() as i64`
    let dt = if first.timedel.is_finite() && first.timedel > 0.0 {
        first.timedel
    } else {
        1.0
    };
    let mut bins: std::collections::BTreeMap<i64, (f64, f64, f64, f64, usize)> =
        std::collections::BTreeMap::new();

    for lc in lcs {
        for i in 0..lc.time.len() {
            let t = lc.time[i];
            let key = (t / dt).round() as i64;
            let entry = bins.entry(key).or_insert((t, 0.0, 0.0, 0.0, 0));
            entry.1 += lc.rate[i] as f64;
            let e = lc.error[i] as f64;
            entry.2 += e * e;
            entry.3 += lc.fracexp[i] as f64;
            entry.4 += 1;
        }
    }

    let mut time = Vec::with_capacity(bins.len());
    let mut rate = Vec::with_capacity(bins.len());
    let mut error = Vec::with_capacity(bins.len());
    let mut fracexp = Vec::with_capacity(bins.len());
    for (_, (t, r, v, fe_sum, cnt)) in bins {
        time.push(t);
        rate.push(r as f32);
        error.push(v.sqrt() as f32);
        fracexp.push((fe_sum / (cnt.max(1) as f64)) as f32);
    }

    CztiLightCurve {
        quad_id: None,
        time,
        rate,
        error,
        fracexp,
        timedel: first.timedel,
        tstart: first.tstart,
        tstop: first.tstop,
        emin_kev: first.emin_kev,
        emax_kev: first.emax_kev,
        thetax_deg: first.thetax_deg,
        thetay_deg: first.thetay_deg,
        mask_weighted: first.mask_weighted,
        badpix_threshold: first.badpix_threshold,
    }
}
