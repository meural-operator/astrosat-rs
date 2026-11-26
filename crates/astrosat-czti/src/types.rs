//! Core AstroSat CZTI domain types, detector geometry helpers, and columnar tables.

use serde::{Deserialize, Serialize};

/// Number of independent CZTI quadrants (Q0, Q1, Q2, Q3).
pub const NUM_QUADS: usize = 4;
/// Number of CZT detector modules (4x4) per quadrant.
pub const NUM_DET_PER_QUAD: usize = 16;
/// Number of pixels (16x16) per CZT detector module.
pub const NUM_PIX_PER_DET: usize = 256;
/// Total pixels per quadrant (64x64 = 4096).
pub const PIXELS_PER_QUAD: usize = 4096;
/// Number of pixels along one quadrant edge.
pub const QUAD_DIM: usize = 64;
/// Number of pixels along full CZTI plane edge (2x2 quadrants).
pub const FULL_DIM: usize = 128;
/// Number of Pulse-Invariant (PI) energy channels in CZTI Level-2 products.
pub const PI_CHANNELS: usize = 512;
/// Number of Veto spectrum channels.
pub const VETOSPEC_SIZE: usize = 256;
/// Number of SSM CZT spectrum channels.
pub const CZTSPEC_SIZE: usize = 512;
/// Number of temperature calibration bins in CALDB GAIN / RESP_PAR (0, 5, 10, 15, 20 C).
pub const NUM_TEMP_BINS: usize = 5;
/// CZTI coded-mask height above detector plane in mm.
pub const MASK_HEIGHT_MM: f64 = 478.0;
/// CZTI detector module active side length (mm) for WCS scale (`XLDET` in pipeline v3.0.1).
pub const XLDET_MM: f64 = 39.0;
/// CZTI maximum field of view half-angle in degrees for mask-weighted spectroscopy/LC.
pub const MAX_FOV_DEG: f64 = 5.8;
/// Reference epoch for CZTI secular gain drift correction (AstroSat seconds).
pub const GAIN_DRIFT_EPOCH_SEC: f64 = 181_094_400.0;
/// CZTI onboard clock frequency for fractional ticks (50 kHz => 20 us per tick).
pub const CZT_TICK_HZ: f64 = 50_000.0;

/// Pixel quality flag definitions matching CZTI CALDB / `badpixTemplate`.
pub mod pix_flag {
    pub const GOOD: u8 = 0;
    pub const SPEC_BAD: u8 = 1;
    pub const FLICKERING: u8 = 2;
    pub const NOISY: u8 = 3;
    pub const DEAD: u8 = 4;
}

/// GTI selection mode for `datasel`, `pixclean`, `dpigen`, `image`, and `bindata`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum GtiType {
    /// Separate GTI per quadrant (`Q0_GTI..Q3_GTI`).
    Quad,
    /// Common GTI shared across all quadrants (`GTI`).
    Common,
}

impl GtiType {
    pub fn as_fits_str(self) -> &'static str {
        match self {
            GtiType::Quad => "QUAD",
            GtiType::Common => "COMMON",
        }
    }

    pub fn gti_extname(self, qid: usize) -> String {
        match self {
            GtiType::Common => "GTI".to_string(),
            GtiType::Quad => format!("Q{qid}_GTI"),
        }
    }

    pub fn parse(s: &str) -> Option<Self> {
        match s.trim().to_ascii_uppercase().as_str() {
            "QUAD" => Some(GtiType::Quad),
            "COMMON" => Some(GtiType::Common),
            _ => None,
        }
    }
}

/// Stage 2 event selection & noise cleaning algorithm.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum EventSelectionAlgo {
    /// Classic pipeline (`datasel` -> `pixclean` -> `evtencor` -> `evtclean`).
    Default,
    /// Generalized v3.0 noise cleaning (`noisypixclean` -> `datasel` -> `superbunchclean` ->
    /// `heavybunchclean` -> `flickpixclean` -> `eventsep` -> `evtencor` -> `evtclean`).
    Generalized,
}

/// Map `(quad_id, det_id, pix_id)` to quadrant-local `(detx, dety)` in `0..64`.
///
/// Exact CZTI Pipeline v3.0.1 convention (`generate_detx_dety`):
/// - `dx = det_id % 4`, `dy = det_id / 4`
/// - `px = pix_id % 16`, `py = pix_id / 16`
/// - `detx = dx * 16 + px`, `dety = dy * 16 + py`
/// - For Q0 and Q3: `dety = 63 - dety`
/// - For Q1 and Q2: `detx = 63 - detx`
#[inline]
pub fn detid_pixid_to_detxy(quad_id: u8, det_id: u8, pix_id: u8) -> (u8, u8) {
    let dx = det_id % 4;
    let dy = det_id / 4;
    let px = pix_id % 16;
    let py = pix_id / 16;
    let mut detx = dx * 16 + px;
    let mut dety = dy * 16 + py;
    if quad_id == 0 || quad_id == 3 {
        dety = 63 - dety;
    } else {
        detx = 63 - detx;
    }
    (detx, dety)
}

/// Map quadrant-local `(detx, dety)` in `0..64` to full-plane `(locx, locy)` in `0..128`.
///
/// Matches `cztshadowgen.cpp` and `Badpix` layout:
/// - Q0: `(detx, 64 + dety)`
/// - Q1: `(64 + detx, 64 + dety)`
/// - Q2: `(64 + detx, dety)`
/// - Q3: `(detx, dety)`
#[inline]
pub fn detxy_to_locxy(quad_id: u8, detx: u8, dety: u8) -> (usize, usize) {
    let x = detx as usize;
    let y = dety as usize;
    match quad_id {
        0 => (x, 64 + y),
        1 => (64 + x, 64 + y),
        2 => (64 + x, y),
        3 => (x, y),
        _ => (x, y),
    }
}

/// Inverse mapping from `(quad_id, detx, dety)` back to `(det_id, pix_id)`.
#[inline]
pub fn detxy_to_detid_pixid(quad_id: u8, detx: u8, dety: u8) -> (u8, u8) {
    let (ux, uy) = if quad_id == 0 || quad_id == 3 {
        (detx, 63 - dety)
    } else {
        (63 - detx, dety)
    };
    let dx = ux / 16;
    let dy = uy / 16;
    let px = ux % 16;
    let py = uy % 16;
    (dy * 4 + dx, py * 16 + px)
}

/// Columnar (SoA) event table for a single CZTI quadrant (`Q0..Q3`).
///
/// Before `pha2energy`, `pi` and `energy` may be empty; after `pha2energy` they have the same
/// length as `time`. Optional `open_fraction` and `weight` are populated by `bindata` when
/// writing mask-weighted event files.
#[derive(Debug, Clone, Default)]
pub struct QuadEvents {
    pub quad_id: u8,
    pub time: Vec<f64>,
    pub cztseccnt: Vec<f64>,
    pub cztntick: Vec<u16>,
    pub pha: Vec<u16>,
    pub detid: Vec<u8>,
    pub pixid: Vec<u8>,
    pub detx: Vec<u8>,
    pub dety: Vec<u8>,
    pub veto: Vec<u16>,
    pub alpha: Vec<u8>,
    pub pi: Vec<i16>,
    pub energy: Vec<f32>,
    pub open_fraction: Option<Vec<f32>>,
    pub weight: Option<Vec<f32>>,
    pub tstart: f64,
    pub tstop: f64,
    pub exposure: f64,
}

impl QuadEvents {
    pub fn new(quad_id: u8) -> Self {
        Self {
            quad_id,
            ..Default::default()
        }
    }

    pub fn with_capacity(quad_id: u8, cap: usize) -> Self {
        Self {
            quad_id,
            time: Vec::with_capacity(cap),
            cztseccnt: Vec::with_capacity(cap),
            cztntick: Vec::with_capacity(cap),
            pha: Vec::with_capacity(cap),
            detid: Vec::with_capacity(cap),
            pixid: Vec::with_capacity(cap),
            detx: Vec::with_capacity(cap),
            dety: Vec::with_capacity(cap),
            veto: Vec::with_capacity(cap),
            alpha: Vec::with_capacity(cap),
            pi: Vec::new(),
            energy: Vec::new(),
            open_fraction: None,
            weight: None,
            tstart: 0.0,
            tstop: 0.0,
            exposure: 0.0,
        }
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.time.len()
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.time.is_empty()
    }

    #[inline]
    pub fn has_energy(&self) -> bool {
        self.pi.len() == self.time.len() && self.energy.len() == self.time.len() && !self.time.is_empty()
    }

    /// Filter events by a boolean keep mask, returning a compact `QuadEvents`.
    pub fn filter_by_mask(&self, keep: &[bool]) -> Self {
        debug_assert_eq!(keep.len(), self.len());
        let n_keep = keep.iter().filter(|&&b| b).count();
        let has_e = self.pi.len() == self.len() && self.energy.len() == self.len();
        let mut out = Self {
            quad_id: self.quad_id,
            time: Vec::with_capacity(n_keep),
            cztseccnt: Vec::with_capacity(n_keep),
            cztntick: Vec::with_capacity(n_keep),
            pha: Vec::with_capacity(n_keep),
            detid: Vec::with_capacity(n_keep),
            pixid: Vec::with_capacity(n_keep),
            detx: Vec::with_capacity(n_keep),
            dety: Vec::with_capacity(n_keep),
            veto: Vec::with_capacity(n_keep),
            alpha: Vec::with_capacity(n_keep),
            pi: if has_e { Vec::with_capacity(n_keep) } else { Vec::new() },
            energy: if has_e { Vec::with_capacity(n_keep) } else { Vec::new() },
            open_fraction: self.open_fraction.as_ref().map(|_| Vec::with_capacity(n_keep)),
            weight: self.weight.as_ref().map(|_| Vec::with_capacity(n_keep)),
            tstart: self.tstart,
            tstop: self.tstop,
            exposure: self.exposure,
        };
        for (i, &k) in keep.iter().enumerate() {
            if !k {
                continue;
            }
            out.time.push(self.time[i]);
            out.cztseccnt.push(self.cztseccnt[i]);
            out.cztntick.push(self.cztntick[i]);
            out.pha.push(self.pha[i]);
            out.detid.push(self.detid[i]);
            out.pixid.push(self.pixid[i]);
            out.detx.push(self.detx[i]);
            out.dety.push(self.dety[i]);
            out.veto.push(self.veto[i]);
            out.alpha.push(self.alpha[i]);
            if has_e {
                out.pi.push(self.pi[i]);
                out.energy.push(self.energy[i]);
            }
            if let (Some(src), Some(dst)) = (&self.open_fraction, &mut out.open_fraction) {
                dst.push(src[i]);
            }
            if let (Some(src), Some(dst)) = (&self.weight, &mut out.weight) {
                dst.push(src[i]);
            }
        }
        out
    }
}

/// Cosmic-ray bunch record table for a single quadrant (`Q0..Q3` of `*_bunch.fits`).
#[derive(Debug, Clone, Default)]
pub struct BunchTable {
    pub quad_id: u8,
    pub time: Vec<f64>,
    pub time_dfs: Vec<u8>,
    pub time_dsl: Vec<u8>,
    /// Number of events in the bunch minus 1 (`NumEvent` column in FITS).
    pub num_event: Vec<u8>,
    pub detid1: Vec<u8>,
    pub detid2: Vec<u8>,
    pub detid3: Vec<u8>,
    pub detid4: Vec<u8>,
    /// 0-based event row index of the first event of the bunch (`RowNum` column).
    pub row_num: Vec<i32>,
    pub detid_fevt: Vec<u8>,
    pub pixid_fevt: Vec<u8>,
    pub detid_sevt: Vec<u8>,
    pub pixid_sevt: Vec<u8>,
    pub detid_levt: Vec<u8>,
    pub pixid_levt: Vec<u8>,
}

impl BunchTable {
    pub fn new(quad_id: u8) -> Self {
        Self {
            quad_id,
            ..Default::default()
        }
    }

    pub fn len(&self) -> usize {
        self.time.len()
    }

    pub fn is_empty(&self) -> bool {
        self.time.is_empty()
    }
}

/// Veto spectrum table (`VETOSPECTRUM` HDU).
#[derive(Debug, Clone, Default)]
pub struct VetoSpectrumTable {
    pub time: Vec<f64>,
    pub cztseccnt: Vec<f64>,
    /// Flattened `nrows * 256` veto spectra (`256U`).
    pub veto_spec: Vec<[u16; VETOSPEC_SIZE]>,
    pub quad_id: Vec<u8>,
}

/// SSM data table (`SSM Data` HDU).
#[derive(Debug, Clone, Default)]
pub struct SsmDataTable {
    pub time: Vec<f64>,
    pub cztseccnt: Vec<f64>,
    pub veto_spec: Vec<[u16; VETOSPEC_SIZE]>,
    pub czt_spec: Vec<[u16; CZTSPEC_SIZE]>,
    pub czt_spec_veto: Vec<[u16; CZTSPEC_SIZE]>,
    pub czt_spec_alpha: Vec<[u16; CZTSPEC_SIZE]>,
    pub quad_id: Vec<u8>,
}

/// Detector module temperature table (`TEMP` HDU, 16 modules per row).
#[derive(Debug, Clone, Default)]
pub struct TempTable {
    pub time: Vec<f64>,
    pub cztseccnt: Vec<f64>,
    pub temperature: Vec<[f32; NUM_DET_PER_QUAD]>,
    pub quad_flag: Vec<u8>,
}

/// Frame header housekeeping table (`*.hdr` FITS file).
#[derive(Debug, Clone, Default)]
pub struct FrameHeaderTable {
    pub time: Vec<f64>,
    pub cztseccnt: Vec<f64>,
    pub data_id: Vec<u8>,
    pub mode_id: Vec<u8>,
    pub word_count: Vec<i32>,
    pub frame_no: Vec<u16>,
    pub status_flag: Vec<u16>,
    pub write_pkt_no: Vec<u16>,
    pub read_pkt_no: Vec<u16>,
    pub command: Vec<i32>,
    pub laxpc_time: Vec<f64>,
    pub dcnt: Vec<u16>,
    pub error_count: Vec<u8>,
    pub error_flag: Vec<u8>,
    pub boot_page_no: Vec<u8>,
    pub czt_no: Vec<u8>,
    pub czt_status: Vec<u8>,
    pub evt_read_mode: Vec<u8>,
    pub cmd_status: Vec<u8>,
    pub buffer_no: Vec<u8>,
    pub base_add: Vec<u8>,
    pub veto_spec_range: Vec<u8>,
    pub channel_no: Vec<u8>,
    pub adc_output: Vec<u16>,
    pub cmd_1sec: Vec<u16>,
    pub alpha_count: Vec<u16>,
    pub veto_count: Vec<u16>,
    pub czt_count_lt_uld: Vec<u16>,
    pub czt_count_gt_uld: Vec<u16>,
    pub czt_data_read: Vec<u16>,
}

/// Binned livetime series for a single quadrant (`Q0..Q3` of `*_livetime.fits`).
#[derive(Debug, Clone, Default)]
pub struct QuadLivetime {
    pub quad_id: u8,
    pub binsize: f64,
    pub time: Vec<f64>,
    pub fracexp: Vec<f64>,
}

/// 4-quadrant bad pixel quality flags (`[quad][det_id * 256 + pix_id]`) and 128x128 full-plane map.
#[derive(Debug, Clone)]
pub struct BadpixTable {
    /// Indexed as `flags[quad][det_id * 256 + pix_id]` (4096 pixels per quadrant).
    pub flags: [Vec<u8>; NUM_QUADS],
}

impl Default for BadpixTable {
    fn default() -> Self {
        Self {
            flags: [
                vec![pix_flag::GOOD; PIXELS_PER_QUAD],
                vec![pix_flag::GOOD; PIXELS_PER_QUAD],
                vec![pix_flag::GOOD; PIXELS_PER_QUAD],
                vec![pix_flag::GOOD; PIXELS_PER_QUAD],
            ],
        }
    }
}

impl BadpixTable {
    /// Get the flag for `(quad_id, det_id, pix_id)`.
    #[inline]
    pub fn get(&self, quad_id: usize, det_id: u8, pix_id: u8) -> u8 {
        self.flags[quad_id][(det_id as usize) * 256 + (pix_id as usize)]
    }

    /// Get the flag indexed by quadrant-local `(detx, dety)`.
    #[inline]
    pub fn get_detxy(&self, quad_id: usize, detx: u8, dety: u8) -> u8 {
        let (det_id, pix_id) = detxy_to_detid_pixid(quad_id as u8, detx, dety);
        self.get(quad_id, det_id, pix_id)
    }

    /// Build the 128x128 full-plane badpixel map (`map[locy][locx]`).
    pub fn to_full_map(&self) -> Vec<Vec<u8>> {
        let mut map = vec![vec![0u8; FULL_DIM]; FULL_DIM];
        for q in 0..NUM_QUADS {
            for det_id in 0..16u8 {
                for pix_id in 0..=255u8 {
                    let (detx, dety) = detid_pixid_to_detxy(q as u8, det_id, pix_id);
                    let (locx, locy) = detxy_to_locxy(q as u8, detx, dety);
                    map[locy][locx] = self.get(q, det_id, pix_id);
                }
            }
        }
        map
    }

    /// Merge multiple badpixel tables by taking the worst (maximum) quality flag per pixel.
    pub fn union_all(tables: &[BadpixTable]) -> Self {
        let mut out = Self::default();
        for t in tables {
            for q in 0..NUM_QUADS {
                for idx in 0..PIXELS_PER_QUAD {
                    if t.flags[q][idx] > out.flags[q][idx] {
                        out.flags[q][idx] = t.flags[q][idx];
                    }
                }
            }
        }
        out
    }
}

/// Per-quadrant aspect time series and average pointing (`RA`, `DEC`, `TWIST` in degrees).
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct QuadAspect {
    pub quad_id: u8,
    pub teldef_name: String,
    pub time: Vec<f64>,
    pub nx: Vec<f32>,
    pub ny: Vec<f32>,
    pub nz: Vec<f32>,
    pub nxt: Vec<f32>,
    pub nyt: Vec<f32>,
    pub nzt: Vec<f32>,
    /// Average RA in degrees.
    pub ra_deg: f64,
    /// Average DEC in degrees.
    pub dec_deg: f64,
    /// Average TWIST in degrees.
    pub twist_deg: f64,
}
