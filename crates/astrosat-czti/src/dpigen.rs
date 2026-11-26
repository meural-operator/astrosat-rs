//! Stage 3a: `cztdpigen` — 64x64 per-quadrant and 128x128 full-plane Detector Plane Histogram (DPH)
//! and effective-area-corrected Detector Plane Image (DPI) generation.

use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use crate::types::{
    detxy_to_detid_pixid, detxy_to_locxy, BadpixTable, QuadEvents, FULL_DIM, NUM_QUADS,
    PIXELS_PER_QUAD, QUAD_DIM,
};

/// Configuration for `cztdpigen`.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DpiGenConfig {
    /// Maximum allowed bad-pixel flag value to include (default: `0`, good pixels only).
    pub badpix_threshold: u8,
    /// Optional energy filter `(emin_kev, emax_kev)`. `None` includes all energies.
    pub energy_range: Option<(f32, f32)>,
    /// Optional time filter `(tstart, tstop)`. `None` includes all times.
    pub time_range: Option<(f64, f64)>,
    /// Quadrants to process (default: `[0, 1, 2, 3]`).
    pub quads: Vec<usize>,
}

impl Default for DpiGenConfig {
    fn default() -> Self {
        Self {
            badpix_threshold: 0,
            energy_range: None,
            time_range: None,
            quads: vec![0, 1, 2, 3],
        }
    }
}

/// Generated DPH (integer counts) and DPI (effective-area corrected float image) products.
#[derive(Debug, Clone)]
pub struct DphDpiProduct {
    /// Per-quadrant 64x64 DPH indexed as `quad_dph[q][dety][detx]`.
    pub quad_dph: [Vec<Vec<i32>>; NUM_QUADS],
    /// Per-quadrant 64x64 DPI indexed as `quad_dpi[q][dety][detx]`.
    pub quad_dpi: [Vec<Vec<f32>>; NUM_QUADS],
    /// Full-plane 128x128 DPH indexed as `full_dph[locy][locx]`.
    pub full_dph: Vec<Vec<i32>>,
    /// Full-plane 128x128 DPI indexed as `full_dpi[locy][locx]`.
    pub full_dpi: Vec<Vec<f32>>,
}

/// Effective area table for one quadrant: 4096 pixels (`detid*256 + pixid`) x `n_ebins`.
#[derive(Debug, Clone)]
pub struct QuadEffArea {
    pub e_min: f32,
    pub e_max: f32,
    pub n_ebins: usize,
    /// Flattened `4096 * n_ebins` array (`area[pix_idx * n_ebins + ebin]`).
    pub area: Vec<f32>,
}

impl Default for QuadEffArea {
    /// Uniform unit effective area fallback when no CALDB `EFF_AREA` file is supplied.
    fn default() -> Self {
        Self {
            e_min: 10.0,
            e_max: 250.0,
            n_ebins: 1,
            area: vec![1.0; PIXELS_PER_QUAD],
        }
    }
}

impl QuadEffArea {
    /// Compute the 64x64 normalized effective area map (`[dety][detx]`) averaged over `[emin, emax]`,
    /// divided by the maximum pixel effective area in the quadrant (`EffArea::calculate_normalized_effarea`).
    pub fn normalized_map(&self, quad_id: u8, emin: f32, emax: f32) -> Vec<Vec<f32>> {
        let mut map = vec![vec![1.0f32; QUAD_DIM]; QUAD_DIM];
        if self.n_ebins == 0 || self.area.len() < PIXELS_PER_QUAD * self.n_ebins {
            return map;
        }
        let step = if self.n_ebins > 1 {
            (self.e_max - self.e_min) / (self.n_ebins as f32)
        } else {
            (self.e_max - self.e_min).max(1.0)
        };
        // Inclusive energy bin range matching `caldbHandler.cpp::EffArea::add_effarea` (`floor(eMin)-1 ..= floor(eMax)-1`)
        let b0 = (((emin - self.e_min) / step).floor().max(0.0) as usize).min(self.n_ebins - 1);
        let b1_incl =
            ((((emax - self.e_min) / step).floor().max(0.0) as usize).min(self.n_ebins - 1)).max(b0);
        let nb = (b1_incl - b0 + 1) as f32;

        let mut max_val = 0.0f32;
        for (dety, row_mut) in map.iter_mut().enumerate().take(QUAD_DIM) {
            for (detx, cell) in row_mut.iter_mut().enumerate().take(QUAD_DIM) {
                let (detid, pixid) = detxy_to_detid_pixid(quad_id, detx as u8, dety as u8);
                let row = (detid as usize) * 256 + (pixid as usize);
                let slice = &self.area[row * self.n_ebins + b0..=row * self.n_ebins + b1_incl];
                let avg = slice.iter().copied().sum::<f32>() / nb;
                *cell = avg;
                if avg > max_val {
                    max_val = avg;
                }
            }
        }
        if max_val > 0.0 {
            for row in &mut map {
                for v in row.iter_mut() {
                    *v /= max_val;
                }
            }
        }
        map
    }
}

/// Generate a single quadrant's 64x64 DPH (`[dety][detx]`) and DPI (`[dety][detx]`).
pub fn generate_quad_dph_dpi(
    events: &QuadEvents,
    badpix: &BadpixTable,
    eff_area: &QuadEffArea,
    cfg: &DpiGenConfig,
) -> (Vec<Vec<i32>>, Vec<Vec<f32>>) {
    let qid = events.quad_id as usize;
    let mut dph = vec![vec![0i32; QUAD_DIM]; QUAD_DIM];
    let has_e = events.energy.len() == events.len();

    for i in 0..events.len() {
        if let Some((t0, t1)) = cfg.time_range {
            let t = events.time[i];
            if t < t0 || t > t1 {
                continue;
            }
        }
        if let Some((e0, e1)) = cfg.energy_range {
            if has_e {
                let e = events.energy[i];
                if e < e0 || e > e1 {
                    continue;
                }
            }
        }
        let detx = events.detx[i] as usize;
        let dety = events.dety[i] as usize;
        if detx >= QUAD_DIM || dety >= QUAD_DIM {
            continue;
        }
        let flag = badpix.get(qid, events.detid[i], events.pixid[i]);
        if flag <= cfg.badpix_threshold {
            dph[dety][detx] += 1;
        }
    }

    let (emin, emax) = cfg.energy_range.unwrap_or((eff_area.e_min, eff_area.e_max));
    let norm_ea = eff_area.normalized_map(events.quad_id, emin, emax);
    let mut dpi = vec![vec![0.0f32; QUAD_DIM]; QUAD_DIM];
    for dety in 0..QUAD_DIM {
        for detx in 0..QUAD_DIM {
            let ea = norm_ea[dety][detx];
            if ea > 0.0 && dph[dety][detx] > 0 {
                dpi[dety][detx] = (dph[dety][detx] as f32) / ea;
            }
        }
    }
    (dph, dpi)
}

/// Generate 64x64 per-quadrant and 128x128 full-plane DPH and DPI in parallel across quadrants.
pub fn generate_dph_dpi(
    quad_events: &[QuadEvents; NUM_QUADS],
    badpix: &BadpixTable,
    eff_areas: &[QuadEffArea; NUM_QUADS],
    cfg: &DpiGenConfig,
) -> DphDpiProduct {
    #[allow(clippy::type_complexity)]
    let per_quad: Vec<(Vec<Vec<i32>>, Vec<Vec<f32>>)> = (0..NUM_QUADS)
        .into_par_iter()
        .map(|q| {
            if cfg.quads.contains(&q) {
                generate_quad_dph_dpi(&quad_events[q], badpix, &eff_areas[q], cfg)
            } else {
                (
                    vec![vec![0i32; QUAD_DIM]; QUAD_DIM],
                    vec![vec![0.0f32; QUAD_DIM]; QUAD_DIM],
                )
            }
        })
        .collect();

    let mut full_dph = vec![vec![0i32; FULL_DIM]; FULL_DIM];
    let mut full_dpi = vec![vec![0.0f32; FULL_DIM]; FULL_DIM];

    for (q, (qdph, qdpi)) in per_quad.iter().enumerate() {
        for dety in 0..QUAD_DIM {
            for detx in 0..QUAD_DIM {
                let (locx, locy) = detxy_to_locxy(q as u8, detx as u8, dety as u8);
                full_dph[locy][locx] = qdph[dety][detx];
                full_dpi[locy][locx] = qdpi[dety][detx];
            }
        }
    }

    let mut iter = per_quad.into_iter();
    let (d0, i0) = iter.next().unwrap();
    let (d1, i1) = iter.next().unwrap();
    let (d2, i2) = iter.next().unwrap();
    let (d3, i3) = iter.next().unwrap();

    DphDpiProduct {
        quad_dph: [d0, d1, d2, d3],
        quad_dpi: [i0, i1, i2, i3],
        full_dph,
        full_dpi,
    }
}
