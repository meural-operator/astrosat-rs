//! FITS and CALDB I/O helpers for all AstroSat CZTI Level-1, CALDB, and Level-2 products.

use std::fs::File;
use std::io::Read;
use std::path::{Path, PathBuf};

use astrosat_core::{CaldbIndex, MkfTable};
use astrosat_fits::{
    apply_ogip_meta, build_gti_hdu, build_lightcurve_hdu, build_response_hdus, write_fits,
    BinTableBuilder, BinTableRef, FitsFile, Hdu, HduKind, ImageBuilder, OgipMeta,
    ResponseMatrixData,
};
use flate2::read::GzDecoder;

use crate::bindata::{
    CztiLightCurve, CztiResponse, CztiSpectrum, QuadRspCaldb, COMP_COLS, TOTAL_MASK_ROWS,
};
use crate::dpigen::QuadEffArea;
use crate::error::{CztiError, Result};
use crate::evtencor::QuadGainCorCaldb;
use crate::gaas::{RpySampleRad, TeldefMatrix};
use crate::gtigen::CztiGtiOutput;
use crate::image::SkyImage;
use crate::pha2energy::{EboundsTable, QuadGainOffset};
use crate::types::{
    detid_pixid_to_detxy, BadpixTable, BunchTable, FrameHeaderTable, GtiType, QuadAspect,
    QuadEvents, QuadLivetime, SsmDataTable, TempTable, VetoSpectrumTable, CZTSPEC_SIZE, FULL_DIM,
    NUM_DET_PER_QUAD, NUM_QUADS, PI_CHANNELS, PIXELS_PER_QUAD, QUAD_DIM, VETOSPEC_SIZE,
};

/// Metadata copied from Level-1 FITS primary/event headers into Level-2 products.
#[derive(Debug, Clone)]
pub struct CztiObsMeta {
    pub obs_id: Option<String>,
    pub object: Option<String>,
    pub ra_obj: Option<f64>,
    pub dec_obj: Option<f64>,
    pub ra_pnt: Option<f64>,
    pub dec_pnt: Option<f64>,
    pub tstart: f64,
    pub tstop: f64,
    pub date_obs: Option<String>,
    pub date_end: Option<String>,
    pub caldb_ver: String,
}

impl Default for CztiObsMeta {
    fn default() -> Self {
        Self {
            obs_id: None,
            object: None,
            ra_obj: None,
            dec_obj: None,
            ra_pnt: None,
            dec_pnt: None,
            tstart: 0.0,
            tstop: 0.0,
            date_obs: None,
            date_end: None,
            caldb_ver: "astrosat-rs".to_string(),
        }
    }
}

impl CztiObsMeta {
    pub fn from_fits(path: impl AsRef<Path>) -> Result<Self> {
        let fits = FitsFile::open(path)?;
        let h0 = fits.hdu(0)?.header();
        let h1 = fits.hdu(1).ok().map(|h| h.header());
        let get_s = |k: &str| {
            h0.get_str(k)
                .or_else(|| h1.and_then(|h| h.get_str(k)))
                .map(|s| s.trim().to_string())
                .filter(|s| !s.is_empty())
        };
        let get_f = |k: &str| h0.get_f64(k).or_else(|| h1.and_then(|h| h.get_f64(k)));

        Ok(Self {
            obs_id: get_s("OBS_ID"),
            object: get_s("OBJECT"),
            ra_obj: get_f("RA_OBJ"),
            dec_obj: get_f("DEC_OBJ"),
            ra_pnt: get_f("RA_PNT"),
            dec_pnt: get_f("DEC_PNT"),
            tstart: get_f("TSTART").unwrap_or(0.0),
            tstop: get_f("TSTOP").unwrap_or(0.0),
            date_obs: get_s("DATE-OBS"),
            date_end: get_s("DATE-END"),
            caldb_ver: get_s("CALDBVER").unwrap_or_else(|| "astrosat-rs".to_string()),
        })
    }

    pub fn to_ogip(&self, detnam: Option<String>, exposure: f64) -> OgipMeta {
        OgipMeta {
            telescop: "ASTROSAT".to_string(),
            instrume: "CZTI".to_string(),
            detnam,
            obs_id: self.obs_id.clone(),
            object: self.object.clone(),
            ra_obj: self.ra_obj.or(self.ra_pnt),
            dec_obj: self.dec_obj.or(self.dec_pnt),
            tstart: self.tstart,
            tstop: self.tstop,
            exposure,
            ontime: Some(exposure),
            livetime: Some(exposure),
            mjdrefi: 55197,
            mjdreff: 0.0,
            timesys: "UTC".to_string(),
            timeref: "LOCAL".to_string(),
            timeunit: "s".to_string(),
            date_obs: self.date_obs.clone(),
            date_end: self.date_end.clone(),
            creator: "astrosat-czti 0.1.0".to_string(),
            origin: "CZTI-POC".to_string(),
        }
    }
}

// ---------------------------------------------------------------------------
// CALDB Wrapper
// ---------------------------------------------------------------------------

/// High-level CZTI CALDB resolver backed by `astrosat_core::CaldbIndex`.
#[derive(Debug, Clone)]
pub struct CztiCaldb {
    index: CaldbIndex,
    caldb_ver: String,
}

impl CztiCaldb {
    /// Open CZTI CALDB from a root directory (supports `data/as1/czti/caldb.indx`,
    /// `data/astrosat/czti/caldb.indx`, or direct `caldb.indx`).
    pub fn open(caldb_root: impl AsRef<Path>) -> Result<Self> {
        let root = caldb_root.as_ref();
        let candidates = [
            root.join("data/as1/czti/caldb.indx"),
            root.join("data/astrosat/czti/caldb.indx"),
            root.join("caldb.indx"),
        ];
        let found = candidates.iter().find(|p| p.exists()).cloned();
        let index = if let Some(ref p) = found {
            CaldbIndex::open(root, Some(p.to_string_lossy().as_ref()))?
        } else {
            CaldbIndex::open(root, None)?
        };

        let caldb_ver = found
            .and_then(|p| FitsFile::open(p).ok())
            .and_then(|f| {
                f.hdu(0)
                    .ok()
                    .and_then(|h| h.header().get_str("CALDBVER").map(|s| s.trim().to_string()))
            })
            .unwrap_or_else(|| "20221209".to_string());

        Ok(Self { index, caldb_ver })
    }

    pub fn version(&self) -> &str {
        &self.caldb_ver
    }

    /// Resolve the best matching CZTI CALDB file path for `(codename, detnam, tstart)`.
    pub fn resolve(&self, codename: &str, detnam: Option<&str>, tstart: f64) -> Result<PathBuf> {
        let entry = self
            .index
            .query_best("ASTROSAT", "CZTI", detnam, codename, Some(tstart), None)
            .or_else(|| self.index.query_best("ASTROSAT", "CZTI", None, codename, Some(tstart), None))
            .ok_or_else(|| CztiError::CaldbLookup {
                codename: codename.to_string(),
                detnam: detnam.unwrap_or("-").to_string(),
                reason: "no matching CIF row".to_string(),
            })?;
        if entry.full_path.exists() {
            return Ok(entry.full_path.clone());
        }
        // Fallback: search under data/as1/czti/bcf or cpf
        let alt1 = self
            .index
            .caldb_root()
            .join("data/as1/czti/bcf")
            .join(&entry.cal_file);
        if alt1.exists() {
            return Ok(alt1);
        }
        let alt2 = self
            .index
            .caldb_root()
            .join("data/as1/czti/cpf")
            .join(&entry.cal_file);
        if alt2.exists() {
            return Ok(alt2);
        }
        Err(CztiError::MissingFile(entry.full_path.clone()))
    }

    pub fn read_ebounds(&self, tstart: f64) -> Result<EboundsTable> {
        let path = self.resolve("EBOUNDS", None, tstart)?;
        read_ebounds_fits(path)
    }

    pub fn read_gain_offset(&self, quad_id: usize, tstart: f64) -> Result<QuadGainOffset> {
        let detnam = format!("QUADRANT{quad_id}");
        let path = self.resolve("GAIN", Some(&detnam), tstart)?;
        read_gain_offset_fits(path, quad_id)
    }

    pub fn read_teldef(&self, quad_id: usize, tstart: f64) -> Result<(String, TeldefMatrix)> {
        let detnam = format!("QUADRANT{quad_id}");
        let path = self.resolve("TELDEF", Some(&detnam), tstart)?;
        let name = path
            .file_name()
            .map(|s| s.to_string_lossy().to_string())
            .unwrap_or_default();
        let mat = read_teldef_fits(&path)?;
        Ok((name, mat))
    }

    pub fn read_badpix(&self, tstart: f64) -> Result<BadpixTable> {
        let path = self.resolve("BADPIX", Some("QUADRANT0"), tstart)?;
        read_badpix_fits(path)
    }

    pub fn read_gaincor_lld_uld(&self, quad_id: usize, tstart: f64) -> Result<QuadGainCorCaldb> {
        let detnam = format!("QUADRANT{quad_id}");
        let gc_path = self.resolve("GAINCOR", Some(&detnam), tstart)?;
        let lld_path = self.resolve("LLD", Some(&detnam), tstart)?;
        let uld_path = self.resolve("ULD", Some(&detnam), tstart)?;
        read_gaincor_lld_uld_fits(gc_path, lld_path, uld_path, quad_id)
    }

    pub fn read_eff_area(&self, quad_id: usize, tstart: f64) -> Result<QuadEffArea> {
        let detnam = format!("QUADRANT{quad_id}");
        let path = self.resolve("EFF_AREA", Some(&detnam), tstart)?;
        read_eff_area_fits(path, quad_id)
    }

    pub fn read_bkgscal(&self, quad_id: usize, tstart: f64) -> Result<Vec<f32>> {
        let detnam = format!("QUADRANT{quad_id}");
        let path = self.resolve("BKGSCAL", Some(&detnam), tstart)?;
        read_bkgscal_fits(path, quad_id)
    }

    pub fn read_syserr(&self, quad_id: usize, tstart: f64) -> Result<[f32; PI_CHANNELS]> {
        let detnam = format!("QUADRANT{quad_id}");
        let path = self.resolve("SYSERR", Some(&detnam), tstart)?;
        read_syserr_fits(path, quad_id)
    }

    pub fn read_mask_pattern(&self, quad_id: usize, tstart: f64) -> Result<Vec<Vec<f32>>> {
        let detnam = format!("QUADRANT{quad_id}");
        let path = self.resolve("MASK_PATTERN", Some(&detnam), tstart)?;
        read_mask_pattern_fits(path, quad_id)
    }

    pub fn read_mask_oversampled(&self, quad_id: usize, tstart: f64) -> Result<Vec<i32>> {
        let path = self.resolve("MASK_OVERSAMPLED", None, tstart)?;
        read_mask_oversampled_fits(path, quad_id)
    }

    pub fn read_shadowlib(&self, quad_id: usize, tstart: f64) -> Result<(usize, Vec<Vec<f32>>)> {
        let detnam = format!("QUADRANT{quad_id}");
        let path = self.resolve("SHADOWLIB", Some(&detnam), tstart)?;
        read_shadowlib_fits(path, quad_id)
    }

    /// Decompress/open CALDB `SHADOWLIB` once and return `(grid_dim, [Q0..Q3 shadows])`.
    pub fn read_shadowlib_all(&self, tstart: f64) -> Result<(usize, [Vec<Vec<f32>>; NUM_QUADS])> {
        let path = self.resolve("SHADOWLIB", Some("QUADRANT0"), tstart)?;
        read_shadowlib_all_fits(path)
    }

    /// Load the shared CALDB `PIXRESP` redistribution table once across all quadrants.
    pub fn read_pixresp(&self, tstart: f64) -> Result<SharedPixResp> {
        let path = self.resolve("PIXRESP", None, tstart)?;
        read_pixresp_fits(path)
    }

    /// Build a quadrant's `QuadRspCaldb` reusing an already-loaded `SharedPixResp`.
    pub fn read_rsp_caldb_with_pixresp(
        &self,
        quad_id: usize,
        tstart: f64,
        pixresp: &SharedPixResp,
    ) -> Result<QuadRspCaldb> {
        let detnam = format!("QUADRANT{quad_id}");
        let resppar_path = self.resolve("RESP_PAR", Some(&detnam), tstart)?;
        let effarea_path = self.resolve("EFF_AREA", Some(&detnam), tstart)?;
        let effareacor_path = self.resolve("EFF_AREACOR", Some(&detnam), tstart)?;
        let lld_path = self.resolve("LLD", Some(&detnam), tstart)?;
        let uld_path = self.resolve("ULD", Some(&detnam), tstart)?;
        read_rsp_caldb_with_pixresp_fits(
            pixresp,
            resppar_path,
            effarea_path,
            effareacor_path,
            lld_path,
            uld_path,
            quad_id,
        )
    }

    pub fn read_rsp_caldb(&self, quad_id: usize, tstart: f64) -> Result<QuadRspCaldb> {
        let pixresp = self.read_pixresp(tstart)?;
        self.read_rsp_caldb_with_pixresp(quad_id, tstart, &pixresp)
    }
}

/// Open a FITS file transparently handling gzip-compressed CALDB files (magic `0x1F 0x8B`),
/// using zero-copy `mmap` for uncompressed FITS files.
pub fn open_fits_maybe_gz(path: impl AsRef<Path>) -> Result<FitsFile> {
    let p = path.as_ref();
    let mut magic = [0u8; 2];
    let n = {
        let mut f = File::open(p)?;
        f.read(&mut magic)?
    };
    if n == 2 && magic == [0x1F, 0x8B] {
        let f = File::open(p)?;
        let mut decoder = GzDecoder::new(f);
        let mut bytes = Vec::new();
        decoder.read_to_end(&mut bytes)?;
        Ok(FitsFile::from_bytes(bytes)?)
    } else {
        Ok(FitsFile::open(p)?)
    }
}

/// Shared CALDB `PIXRESP` redistribution matrix loaded once for all 4 quadrants.
#[derive(Debug, Clone)]
pub struct SharedPixResp {
    pub n_ebins: usize,
    pub n_pibins: usize,
    pub emin: f32,
    pub emax: f32,
    pub matrix: Vec<f32>,
}

/// Read an integer vector/scalar column (supporting `nI`, `nJ`, scalar `I`/`J`/`B`), returning
/// the `elem_idx`-th element per row (clamped to `[0, max_val]`).
fn read_int_col_elem_i16(
    bt: &BinTableRef<'_>,
    col_name: &str,
    elem_idx: usize,
    max_val: i32,
) -> Result<Vec<i16>> {
    if let Ok((rep, vals)) = bt.read_vec_i16(col_name) {
        let k = elem_idx.min(rep.saturating_sub(1));
        let nrows = bt.nrows();
        let mut out = Vec::with_capacity(nrows);
        for r in 0..nrows {
            let v = vals.get(r * rep + k).copied().unwrap_or(0) as i32;
            out.push(v.clamp(0, max_val) as i16);
        }
        return Ok(out);
    }
    if let Ok((rep, vals)) = bt.read_vec_i32(col_name) {
        let k = elem_idx.min(rep.saturating_sub(1));
        let nrows = bt.nrows();
        let mut out = Vec::with_capacity(nrows);
        for r in 0..nrows {
            let v = vals.get(r * rep + k).copied().unwrap_or(0);
            out.push(v.clamp(0, max_val) as i16);
        }
        return Ok(out);
    }
    let vals = bt.read_col_i32(col_name)?;
    Ok(vals
        .into_iter()
        .map(|v| v.clamp(0, max_val) as i16)
        .collect())
}

// ---------------------------------------------------------------------------
// Direct CALDB FITS Readers
// ---------------------------------------------------------------------------

pub fn read_ebounds_fits(path: impl AsRef<Path>) -> Result<EboundsTable> {
    let fits = FitsFile::open(path)?;
    let hdu = fits
        .hdu_by_name("EBOUNDS")
        .or_else(|_| fits.hdu(1))?;
    let bt = hdu.as_bintable()?;
    let ch = bt.read_col_i16("CHANNEL")?;
    let e_min = bt.read_col_f32("E_MIN")?;
    let e_max = bt.read_col_f32("E_MAX")?;
    Ok(EboundsTable {
        channel: ch,
        e_min,
        e_max,
    })
}

pub fn read_gain_offset_fits(path: impl AsRef<Path>, quad_id: usize) -> Result<QuadGainOffset> {
    let fits = FitsFile::open(path)?;
    let extname = format!("Q{quad_id}");
    let hdu = fits
        .hdu_by_name(&extname)
        .or_else(|_| fits.hdu(quad_id + 1))?;
    let bt = hdu.as_bintable()?;
    let detid = bt.read_col_u8("DETID")?;
    let pixid = bt.read_col_u8("PIXID")?;
    let col_names = ["GAIN_00", "GAIN_05", "GAIN_10", "GAIN_15", "GAIN_20"];

    let mut out = QuadGainOffset::default();
    for (tbin, &cname) in col_names.iter().enumerate() {
        let (rep, vals) = bt.read_vec_f32(cname)?;
        if rep >= 2 {
            for (row, (&d, &p)) in detid.iter().zip(pixid.iter()).enumerate() {
                let idx = (d as usize) * 256 + (p as usize);
                if idx < PIXELS_PER_QUAD && row * rep + 1 < vals.len() {
                    out.gain[tbin][idx] = vals[row * rep];
                    out.offset[tbin][idx] = vals[row * rep + 1];
                }
            }
        }
    }
    Ok(out)
}

pub fn read_teldef_fits(path: impl AsRef<Path>) -> Result<TeldefMatrix> {
    let fits = FitsFile::open(path)?;
    let hdr = fits.hdu(0)?.header();
    let get = |k: &str, default: f64| hdr.get_f64(k).unwrap_or(default);
    Ok(TeldefMatrix {
        m: [
            [
                get("ALIGNM11", 1.0),
                get("ALIGNM12", 0.0),
                get("ALIGNM13", 0.0),
            ],
            [
                get("ALIGNM21", 0.0),
                get("ALIGNM22", 1.0),
                get("ALIGNM23", 0.0),
            ],
            [
                get("ALIGNM31", 0.0),
                get("ALIGNM32", 0.0),
                get("ALIGNM33", 1.0),
            ],
        ],
    })
}

pub fn read_badpix_fits(path: impl AsRef<Path>) -> Result<BadpixTable> {
    let fits = FitsFile::open(path)?;
    let mut table = BadpixTable::default();
    for (q, q_flags) in table.flags.iter_mut().enumerate() {
        let extname = format!("Q{q}");
        let hdu = fits.hdu_by_name(&extname).or_else(|_| fits.hdu(q + 1))?;
        let bt = hdu.as_bintable()?;
        let flags = bt.read_col_u8("PIX_FLAG")?;
        if let (Ok(detid), Ok(pixid)) = (bt.read_col_u8("DETID"), bt.read_col_u8("PIXID")) {
            for ((&d, &p), &f) in detid.iter().zip(pixid.iter()).zip(flags.iter()) {
                let idx = (d as usize) * 256 + (p as usize);
                if idx < PIXELS_PER_QUAD {
                    q_flags[idx] = f;
                }
            }
        } else {
            for (idx, &f) in flags.iter().take(PIXELS_PER_QUAD).enumerate() {
                q_flags[idx] = f;
            }
        }
    }
    Ok(table)
}

pub fn read_gaincor_lld_uld_fits(
    gc_path: impl AsRef<Path>,
    lld_path: impl AsRef<Path>,
    uld_path: impl AsRef<Path>,
    quad_id: usize,
) -> Result<QuadGainCorCaldb> {
    let mut out = QuadGainCorCaldb::default();
    let extname = format!("Q{quad_id}");

    {
        let fits = FitsFile::open(gc_path)?;
        let hdu = fits
            .hdu_by_name(&extname)
            .or_else(|_| fits.hdu(quad_id + 1))?;
        let bt = hdu.as_bintable()?;
        let par = bt.read_col_f32("GAINCORPAR")?;
        let flg = bt.read_col_i32("GAINCORFLAG")?;
        let n_copy = PIXELS_PER_QUAD.min(par.len()).min(flg.len());
        out.gaincor_par[..n_copy].copy_from_slice(&par[..n_copy]);
        out.gaincor_flag[..n_copy].copy_from_slice(&flg[..n_copy]);
    }
    {
        let fits = open_fits_maybe_gz(lld_path)?;
        let hdu = fits
            .hdu_by_name(&extname)
            .or_else(|_| fits.hdu(quad_id + 1))?;
        let bt = hdu.as_bintable()?;
        let vals = read_int_col_elem_i16(&bt, "LLD", 0, (PI_CHANNELS - 1) as i32)?;
        for (i, &v) in vals.iter().take(PIXELS_PER_QUAD).enumerate() {
            out.lld[i] = v;
        }
    }
    {
        let fits = open_fits_maybe_gz(uld_path)?;
        let hdu = fits
            .hdu_by_name(&extname)
            .or_else(|_| fits.hdu(quad_id + 1))?;
        let bt = hdu.as_bintable()?;
        let vals = read_int_col_elem_i16(&bt, "ULD", 0, (PI_CHANNELS - 1) as i32)?;
        for (i, &v) in vals.iter().take(PIXELS_PER_QUAD).enumerate() {
            out.uld[i] = v;
        }
    }

    Ok(out)
}

pub fn read_eff_area_fits(path: impl AsRef<Path>, quad_id: usize) -> Result<QuadEffArea> {
    let fits = open_fits_maybe_gz(path)?;
    let extname = format!("Q{quad_id}");
    let hdu = fits
        .hdu_by_name(&extname)
        .or_else(|_| fits.hdu(quad_id + 1))?;
    let (e_min, e_max) = hdu
        .header()
        .get_str("CBD10001")
        .and_then(parse_cbd_energy_range)
        .unwrap_or((1.0, 301.0));
    let bt = hdu.as_bintable()?;
    let (n_ebins, vals) = bt.read_vec_f32("AREA")?;
    let mut area = vec![1.0f32; PIXELS_PER_QUAD * n_ebins.max(1)];
    if let (Ok(detid), Ok(pixid)) = (bt.read_col_u8("DETID"), bt.read_col_u8("PIXID")) {
        for (row, (&d, &p)) in detid.iter().zip(pixid.iter()).enumerate() {
            let dst = (d as usize) * 256 + (p as usize);
            if dst < PIXELS_PER_QUAD && (row + 1) * n_ebins <= vals.len() {
                area[dst * n_ebins..(dst + 1) * n_ebins]
                    .copy_from_slice(&vals[row * n_ebins..(row + 1) * n_ebins]);
            }
        }
    } else {
        let copy_len = area.len().min(vals.len());
        area[..copy_len].copy_from_slice(&vals[..copy_len]);
    }
    Ok(QuadEffArea {
        e_min,
        e_max,
        n_ebins: n_ebins.max(1),
        area,
    })
}

pub fn read_bkgscal_fits(path: impl AsRef<Path>, quad_id: usize) -> Result<Vec<f32>> {
    let fits = FitsFile::open(path)?;
    let extname = format!("Q{quad_id}");
    let hdu = fits
        .hdu_by_name(&extname)
        .or_else(|_| fits.hdu(quad_id + 1))?;
    let bt = hdu.as_bintable()?;
    let (rep, vals) = bt.read_vec_f32("BKGSCAL")?;
    let mut out = vec![1.0f32; PIXELS_PER_QUAD * PI_CHANNELS];
    let ch = rep.min(PI_CHANNELS);
    for row in 0..PIXELS_PER_QUAD.min(bt.nrows()) {
        let src = row * rep;
        let dst = row * PI_CHANNELS;
        out[dst..dst + ch].copy_from_slice(&vals[src..src + ch]);
    }
    Ok(out)
}

pub fn read_syserr_fits(path: impl AsRef<Path>, quad_id: usize) -> Result<[f32; PI_CHANNELS]> {
    let fits = FitsFile::open(path)?;
    let extname = format!("Q{quad_id}");
    let hdu = fits
        .hdu_by_name(&extname)
        .or_else(|_| fits.hdu(quad_id + 1))?;
    let bt = hdu.as_bintable()?;
    let vals = bt.read_col_f32("SYSERR")?;
    let mut out = [0.0f32; PI_CHANNELS];
    for (i, &v) in vals.iter().take(PI_CHANNELS).enumerate() {
        out[i] = v;
    }
    Ok(out)
}

pub fn read_mask_pattern_fits(path: impl AsRef<Path>, quad_id: usize) -> Result<Vec<Vec<f32>>> {
    let fits = FitsFile::open(path)?;
    let extname = format!("Q{quad_id}");
    let hdu = fits
        .hdu_by_name(&extname)
        .or_else(|_| fits.hdu(quad_id + 1))?;
    let img = hdu.as_image()?;
    let flat = img.read_f32()?;
    let mut out = vec![vec![0.0f32; QUAD_DIM]; QUAD_DIM];
    for y in 0..QUAD_DIM {
        for x in 0..QUAD_DIM {
            if y * QUAD_DIM + x < flat.len() {
                out[y][x] = flat[y * QUAD_DIM + x];
            }
        }
    }
    Ok(out)
}

pub fn read_mask_oversampled_fits(path: impl AsRef<Path>, quad_id: usize) -> Result<Vec<i32>> {
    let fits = FitsFile::open(path)?;
    let hdu = fits.hdu(quad_id + 1)?;
    let img = hdu.as_image()?;
    let vals = img.read_f64()?;
    let mut out = vec![0i32; TOTAL_MASK_ROWS * COMP_COLS];
    for (dst, &v) in out.iter_mut().zip(vals.iter()) {
        *dst = v as i32;
    }
    Ok(out)
}

pub fn read_shadowlib_fits(path: impl AsRef<Path>, quad_id: usize) -> Result<(usize, Vec<Vec<f32>>)> {
    let fits = open_fits_maybe_gz(path)?;
    let hdu = fits
        .hdu_by_name("ShadowLibrary")
        .or_else(|_| fits.hdu(1))?;
    let bt = hdu.as_bintable()?;
    let nrows = bt.nrows();
    let grid_dim = (nrows as f64).sqrt().round() as usize;
    let col_name = format!("Q{quad_id}");
    let (rep, flat) = bt.read_vec_f32(&col_name)?;
    let mut shadows = Vec::with_capacity(nrows);
    for r in 0..nrows {
        shadows.push(flat[r * rep..(r + 1) * rep].to_vec());
    }
    Ok((grid_dim, shadows))
}

pub fn read_shadowlib_all_fits(
    path: impl AsRef<Path>,
) -> Result<(usize, [Vec<Vec<f32>>; NUM_QUADS])> {
    let fits = open_fits_maybe_gz(path)?;
    let hdu = fits
        .hdu_by_name("ShadowLibrary")
        .or_else(|_| fits.hdu(1))?;
    let bt = hdu.as_bintable()?;
    let nrows = bt.nrows();
    let grid_dim = (nrows as f64).sqrt().round() as usize;
    let mut per_quad = Vec::with_capacity(NUM_QUADS);
    for q in 0..NUM_QUADS {
        let (rep, flat) = bt.read_vec_f32(&format!("Q{q}"))?;
        let mut shadows = Vec::with_capacity(nrows);
        for r in 0..nrows {
            shadows.push(flat[r * rep..(r + 1) * rep].to_vec());
        }
        per_quad.push(shadows);
    }
    let mut it = per_quad.into_iter();
    Ok((
        grid_dim,
        [
            it.next().unwrap(),
            it.next().unwrap(),
            it.next().unwrap(),
            it.next().unwrap(),
        ],
    ))
}

pub fn read_pixresp_fits(path: impl AsRef<Path>) -> Result<SharedPixResp> {
    let pix_fits = open_fits_maybe_gz(path)?;
    let pix_hdu = pix_fits
        .hdu_by_name("PIXRESP")
        .or_else(|_| pix_fits.hdu(1))?;
    let raw_ebins = pix_hdu.header().get_i64("N_EBIN").unwrap_or(300);
    let raw_pibins = pix_hdu.header().get_i64("N_PIBIN").unwrap_or(512);
    if raw_ebins <= 0 || raw_pibins <= 0 {
        return Err(CztiError::InvalidConfig(format!(
            "Invalid PIXRESP dimensions: N_EBIN={raw_ebins}, N_PIBIN={raw_pibins}"
        )));
    }
    let n_ebins = raw_ebins as usize;
    let n_pibins = raw_pibins as usize;
    let (emin, emax) = pix_hdu
        .header()
        .get_str("CBD20001")
        .and_then(parse_cbd_energy_range)
        .unwrap_or((1.0, 301.0));
    let pix_bt = pix_hdu.as_bintable()?;
    let (_rep, matrix) = pix_bt.read_vec_f32("MATRIX")?;
    Ok(SharedPixResp {
        n_ebins,
        n_pibins,
        emin,
        emax,
        matrix,
    })
}

#[allow(clippy::too_many_arguments)]
pub fn read_rsp_caldb_with_pixresp_fits(
    pixresp: &SharedPixResp,
    resppar_path: impl AsRef<Path>,
    effarea_path: impl AsRef<Path>,
    effareacor_path: impl AsRef<Path>,
    lld_path: impl AsRef<Path>,
    uld_path: impl AsRef<Path>,
    quad_id: usize,
) -> Result<QuadRspCaldb> {
    let extname = format!("Q{quad_id}");
    let n_ebins = pixresp.n_ebins;
    let n_pibins = pixresp.n_pibins;

    // 2. RESP_PAR (`RESP_ID` is `5I` in CALDB)
    let rp_fits = open_fits_maybe_gz(resppar_path)?;
    let rp_hdu = rp_fits
        .hdu_by_name(&extname)
        .or_else(|_| rp_fits.hdu(quad_id + 1))?;
    let rp_bt = rp_hdu.as_bintable()?;
    let mut resp_id = vec![[0usize; 5]; PIXELS_PER_QUAD];
    if let Ok((rp_rep, rp_vals)) = rp_bt.read_vec_i16("RESP_ID") {
        for (r, slot) in resp_id.iter_mut().enumerate().take(rp_bt.nrows()) {
            for (t, cell) in slot.iter_mut().enumerate().take(rp_rep) {
                *cell = rp_vals[r * rp_rep + t].max(0) as usize;
            }
        }
    } else {
        let (rp_rep, rp_vals) = rp_bt.read_vec_i32("RESP_ID")?;
        for (r, slot) in resp_id.iter_mut().enumerate().take(rp_bt.nrows()) {
            for (t, cell) in slot.iter_mut().enumerate().take(rp_rep) {
                *cell = rp_vals[r * rp_rep + t].max(0) as usize;
            }
        }
    }

    // 3. EFF_AREA
    let ea = read_eff_area_fits(effarea_path, quad_id)?;

    // 4. EFF_AREACOR
    let mut eff_area_cor1 = vec![1.0f32; NUM_DET_PER_QUAD * n_pibins];
    let mut eff_area_cor2 = vec![1.0f32; n_ebins];
    if let Ok(ec_fits) = open_fits_maybe_gz(effareacor_path) {
        if let Ok(ec_hdu) = ec_fits
            .hdu_by_name(&extname)
            .or_else(|_| ec_fits.hdu(quad_id + 1))
        {
            if let Ok(ec_bt) = ec_hdu.as_bintable() {
                if let Ok((r1, v1)) = ec_bt.read_vec_f32("AREACOR1") {
                    for det in 0..NUM_DET_PER_QUAD.min(ec_bt.nrows()) {
                        let c = r1.min(n_pibins);
                        eff_area_cor1[det * n_pibins..det * n_pibins + c]
                            .copy_from_slice(&v1[det * r1..det * r1 + c]);
                    }
                }
                if let Ok((r2, v2)) = ec_bt.read_vec_f32("AREACOR2") {
                    let c = r2.min(n_ebins);
                    eff_area_cor2[..c].copy_from_slice(&v2[..c]);
                }
            }
        }
    }

    // 5. LLD & ULD (`5I` in CALDB)
    let mut lld = vec![0i16; PIXELS_PER_QUAD];
    let mut uld = vec![(n_pibins - 1) as i16; PIXELS_PER_QUAD];
    {
        let fits = open_fits_maybe_gz(lld_path)?;
        let hdu = fits
            .hdu_by_name(&extname)
            .or_else(|_| fits.hdu(quad_id + 1))?;
        let bt = hdu.as_bintable()?;
        let vals = read_int_col_elem_i16(&bt, "LLD", 0, (n_pibins - 1) as i32)?;
        for (i, &v) in vals.iter().take(PIXELS_PER_QUAD).enumerate() {
            lld[i] = v;
        }
    }
    {
        let fits = open_fits_maybe_gz(uld_path)?;
        let hdu = fits
            .hdu_by_name(&extname)
            .or_else(|_| fits.hdu(quad_id + 1))?;
        let bt = hdu.as_bintable()?;
        let vals = read_int_col_elem_i16(&bt, "ULD", 0, (n_pibins - 1) as i32)?;
        for (i, &v) in vals.iter().take(PIXELS_PER_QUAD).enumerate() {
            uld[i] = v;
        }
    }

    Ok(QuadRspCaldb {
        n_ebins,
        n_pibins,
        emin: pixresp.emin,
        emax: pixresp.emax,
        pixresp_matrix: pixresp.matrix.clone(),
        resp_id,
        eff_area: ea.area,
        eff_area_cor1,
        eff_area_cor2,
        lld,
        uld,
    })
}

#[allow(clippy::too_many_arguments)]
pub fn read_rsp_caldb_fits(
    pixresp_path: impl AsRef<Path>,
    resppar_path: impl AsRef<Path>,
    effarea_path: impl AsRef<Path>,
    effareacor_path: impl AsRef<Path>,
    lld_path: impl AsRef<Path>,
    uld_path: impl AsRef<Path>,
    quad_id: usize,
) -> Result<QuadRspCaldb> {
    let pixresp = read_pixresp_fits(pixresp_path)?;
    read_rsp_caldb_with_pixresp_fits(
        &pixresp,
        resppar_path,
        effarea_path,
        effareacor_path,
        lld_path,
        uld_path,
        quad_id,
    )
}

fn parse_cbd_energy_range(s: &str) -> Option<(f32, f32)> {
    // Format like `ENERG("1.0-301.0")keV` or `ENE("1.0-301.0")keV` or `ENER(10.0-250.0)keV`
    let start = s.find('(')? + 1;
    let end = s[start..].find(')')? + start;
    let inner = s[start..end].trim().trim_matches('"').trim_matches('\'');
    let (a, b) = inner.split_once('-')?;
    Some((a.trim().parse().ok()?, b.trim().parse().ok()?))
}

// ---------------------------------------------------------------------------
// Level-1 & MKF FITS Readers
// ---------------------------------------------------------------------------

/// Read raw 2048-byte packets for quadrant `quad_id` (`CZT_QUAD{quad_id+1}` `DataArray` column).
pub fn read_l1_quad_packets(l1_path: impl AsRef<Path>, quad_id: usize) -> Result<Vec<u8>> {
    let fits = FitsFile::open(l1_path)?;
    let extname = format!("CZT_QUAD{}", quad_id + 1);
    let hdu = fits.hdu_by_name(&extname)?;
    let bt = hdu.as_bintable()?;
    let (_rep, data) = bt.read_vec_u8("DataArray")?;
    Ok(data)
}

/// Read common `GTI`, quadrant `Q0_GTI..Q3_GTI`, and optional `BTI` from a CZTI FITS file.
pub fn read_czti_gtis(path: impl AsRef<Path>) -> Result<(CztiGtiOutput, Vec<(f64, f64)>)> {
    let fits = FitsFile::open(path)?;
    let read_ext = |name: &str| -> Vec<(f64, f64)> {
        let Ok(hdu) = fits.hdu_by_name(name) else {
            return Vec::new();
        };
        let Ok(bt) = hdu.as_bintable() else {
            return Vec::new();
        };
        let Ok(starts) = bt.read_col_f64("START") else {
            return Vec::new();
        };
        let Ok(stops) = bt.read_col_f64("STOP") else {
            return Vec::new();
        };
        starts.into_iter().zip(stops).collect()
    };

    let common = read_ext("GTI");
    let quad = [
        read_ext("Q0_GTI"),
        read_ext("Q1_GTI"),
        read_ext("Q2_GTI"),
        read_ext("Q3_GTI"),
    ];
    let bti = read_ext("BTI");
    Ok((CztiGtiOutput { common, quad }, bti))
}

/// Extract Roll/Pitch/Yaw `(RA, DEC)` in radians from an AstroSat MKF table for `cztgaas`.
pub fn extract_mkf_rpy_rad(mkf: &MkfTable) -> Result<Vec<RpySampleRad>> {
    let times = mkf.times();
    let get = |names: &[&str]| -> Result<&[f64]> {
        mkf.column_any(names).ok_or_else(|| {
            CztiError::InvalidConfig(format!("MKF missing required attitude column {:?}", names))
        })
    };
    let r_ra = get(&["ROLL_RA", "Roll_RA"])?;
    let r_dec = get(&["ROLL_DEC", "Roll_DEC"])?;
    let p_ra = get(&["PITCH_RA", "Pitch_RA"])?;
    let p_dec = get(&["PITCH_DEC", "Pitch_DEC"])?;
    let y_ra = get(&["YAW_RA", "Yaw_RA"])?;
    let y_dec = get(&["YAW_DEC", "Yaw_DEC"])?;

    let mut out = Vec::with_capacity(times.len());
    for i in 0..times.len() {
        out.push(RpySampleRad {
            time: times[i],
            roll_ra: r_ra[i].to_radians(),
            roll_dec: r_dec[i].to_radians(),
            pitch_ra: p_ra[i].to_radians(),
            pitch_dec: p_dec[i].to_radians(),
            yaw_ra: y_ra[i].to_radians(),
            yaw_dec: y_dec[i].to_radians(),
        });
    }
    Ok(out)
}

// ---------------------------------------------------------------------------
// Level-2 FITS Writers & Readers
// ---------------------------------------------------------------------------

/// Options for writing a CZTI Level-2 event file (`*_level2.fits`, `*_bc.evt`, `*_clean.evt`).
#[derive(Debug, Clone, Default)]
pub struct EventFileWriteOpts {
    pub gti_type: Option<GtiType>,
    pub ncompat: Option<i32>,
    pub gaincor_applied: bool,
    pub lld_cut: bool,
    pub uld_cut: bool,
    pub pixel_exposure: Option<[Vec<f64>; NUM_QUADS]>,
}

/// Write a complete multi-HDU CZTI Level-2 event FITS file (`Primary`, `Q0..Q3`, `VETOSPECTRUM`,
/// `SSM Data`, `TEMP`, `GTI`, `Q0_GTI..Q3_GTI`, and optional `EXPOSURE`).
#[allow(clippy::too_many_arguments)]
pub fn write_l2_event_file(
    path: impl AsRef<Path>,
    meta: &CztiObsMeta,
    quads: &[QuadEvents; NUM_QUADS],
    veto_spec: &VetoSpectrumTable,
    ssm_data: &SsmDataTable,
    temp: &TempTable,
    gtis: &CztiGtiOutput,
    opts: &EventFileWriteOpts,
) -> Result<()> {
    let mut hdus: Vec<Hdu> = Vec::with_capacity(16);

    let mut primary = Hdu::empty_primary();
    let total_exp = quads.iter().map(|q| q.exposure).fold(0.0f64, f64::max);
    apply_ogip_meta(&mut primary.header, &meta.to_ogip(None, total_exp));
    primary
        .header
        .set_str("MISSION", "ASTROSAT", "Name of mission/satellite");
    primary
        .header
        .set_str("CALDBVER", &meta.caldb_ver, "CALDB version");
    if let Some(gt) = opts.gti_type {
        primary
            .header
            .set_str("GTITYPE", gt.as_fits_str(), "GTI selection mode");
    }
    if let Some(nc) = opts.ncompat {
        primary
            .header
            .set_i64("NCOMPAT", nc as i64, "Generalized noise clean compatibility");
    }
    if opts.gaincor_applied {
        primary
            .header
            .set_str("GAINCOR", "YES", "In-flight gain correction applied");
    }
    if opts.lld_cut {
        primary.header.set_str("LLDCUT", "YES", "LLD cut applied");
    }
    if opts.uld_cut {
        primary.header.set_str("ULDCUT", "YES", "ULD cut applied");
    }
    hdus.push(primary);

    // Q0..Q3 event tables
    for (q, ev) in quads.iter().enumerate() {
        let mut b = BinTableBuilder::new(format!("Q{q}"))
            .add_f64("Time", Some("s"), ev.time.clone())
            .add_f64("CZTSECCNT", Some("s"), ev.cztseccnt.clone())
            .add_u16("CZTNTICK", Some("micro_sec"), ev.cztntick.clone())
            .add_u16("PHA", Some("counts"), ev.pha.clone())
            .add_u8("DetID", None, ev.detid.clone())
            .add_u8("pixID", None, ev.pixid.clone())
            .add_u8("DETX", None, ev.detx.clone())
            .add_u8("DETY", None, ev.dety.clone())
            .add_u16("veto", Some("counts"), ev.veto.clone())
            .add_u8("alpha", Some("count"), ev.alpha.clone());

        if ev.pi.len() == ev.len() && ev.energy.len() == ev.len() {
            b = b
                .add_i16("PI", None, ev.pi.clone())
                .add_f32("ENERGY", Some("keV"), ev.energy.clone());
        }
        if let (Some(of), Some(wt)) = (&ev.open_fraction, &ev.weight) {
            if of.len() == ev.len() && wt.len() == ev.len() {
                b = b
                    .add_f32("OPEN_FRACTION", None, of.clone())
                    .add_f32("WEIGHT", None, wt.clone());
            }
        }

        let mut hdu = b.header_i64("QUADID", q as i64, "Quadrant Number").build()?;
        let mut qmeta = meta.to_ogip(Some(format!("QUADRANT{q}")), ev.exposure);
        if ev.tstart > 0.0 {
            qmeta.tstart = ev.tstart;
        }
        if ev.tstop > 0.0 {
            qmeta.tstop = ev.tstop;
        }
        apply_ogip_meta(&mut hdu.header, &qmeta);
        if let Some(gt) = opts.gti_type {
            hdu.header.set_str("GTITYPE", gt.as_fits_str(), "GTI type");
        }
        if let Some(nc) = opts.ncompat {
            hdu.header.set_i64("NCOMPAT", nc as i64, "Noise clean compat");
        }
        if opts.gaincor_applied {
            hdu.header.set_str("GAINCOR", "YES", "Gain correction applied");
        }
        hdus.push(hdu);
    }

    // VETOSPECTRUM
    let flat_vs: Vec<u16> = veto_spec
        .veto_spec
        .iter()
        .flat_map(|r| r.iter().copied())
        .collect();
    let mut vs_hdu = BinTableBuilder::new("VETOSPECTRUM")
        .add_f64("Time", Some("s"), veto_spec.time.clone())
        .add_f64("CZTSECCNT", Some("s"), veto_spec.cztseccnt.clone())
        .add_vec_u16("VetoSpec", None, VETOSPEC_SIZE, flat_vs)
        .add_u8("QuadID", None, veto_spec.quad_id.clone())
        .build()?;
    apply_ogip_meta(&mut vs_hdu.header, &meta.to_ogip(None, total_exp));
    hdus.push(vs_hdu);

    // SSM Data
    let flat_ssm_vs: Vec<u16> = ssm_data
        .veto_spec
        .iter()
        .flat_map(|r| r.iter().copied())
        .collect();
    let flat_ssm_cs: Vec<u16> = ssm_data
        .czt_spec
        .iter()
        .flat_map(|r| r.iter().copied())
        .collect();
    let flat_ssm_csv: Vec<u16> = ssm_data
        .czt_spec_veto
        .iter()
        .flat_map(|r| r.iter().copied())
        .collect();
    let flat_ssm_csa: Vec<u16> = ssm_data
        .czt_spec_alpha
        .iter()
        .flat_map(|r| r.iter().copied())
        .collect();
    let mut ssm_hdu = BinTableBuilder::new("SSM Data")
        .add_f64("Time", Some("s"), ssm_data.time.clone())
        .add_f64("CZTSECCNT", Some("s"), ssm_data.cztseccnt.clone())
        .add_vec_u16("VetoSpec", None, VETOSPEC_SIZE, flat_ssm_vs)
        .add_vec_u16("CZTSpec", None, CZTSPEC_SIZE, flat_ssm_cs)
        .add_vec_u16("CZTSpecVeto", None, CZTSPEC_SIZE, flat_ssm_csv)
        .add_vec_u16("CZTSpecAlpha", None, CZTSPEC_SIZE, flat_ssm_csa)
        .add_u8("QuadID", None, ssm_data.quad_id.clone())
        .build()?;
    apply_ogip_meta(&mut ssm_hdu.header, &meta.to_ogip(None, total_exp));
    hdus.push(ssm_hdu);

    // TEMP
    let flat_temp: Vec<f32> = temp
        .temperature
        .iter()
        .flat_map(|r| r.iter().copied())
        .collect();
    let mut temp_hdu = BinTableBuilder::new("TEMP")
        .add_f64("TIME", Some("s"), temp.time.clone())
        .add_f64("CZTSECCNT", Some("s"), temp.cztseccnt.clone())
        .add_vec_f32("TEMPERATURE", Some("C"), NUM_DET_PER_QUAD, flat_temp)
        .add_u8("QUAD_FLAG", None, temp.quad_flag.clone())
        .build()?;
    apply_ogip_meta(&mut temp_hdu.header, &meta.to_ogip(None, total_exp));
    hdus.push(temp_hdu);

    // GTI & Q0_GTI..Q3_GTI
    let push_gti = |hdus: &mut Vec<Hdu>, name: &str, ivs: &[(f64, f64)]| -> Result<()> {
        let starts: Vec<f64> = ivs.iter().map(|&(a, _)| a).collect();
        let stops: Vec<f64> = ivs.iter().map(|&(_, b)| b).collect();
        let exp: f64 = ivs.iter().map(|&(a, b)| (b - a).max(0.0)).sum();
        let hdu = build_gti_hdu(name, &starts, &stops, &meta.to_ogip(None, exp))?;
        hdus.push(hdu);
        Ok(())
    };
    push_gti(&mut hdus, "GTI", &gtis.common)?;
    for (q, q_ivs) in gtis.quad.iter().enumerate() {
        push_gti(&mut hdus, &format!("Q{q}_GTI"), q_ivs)?;
    }

    // Optional EXPOSURE extension (4096 rows x 4 columns)
    if let Some(pexp) = &opts.pixel_exposure {
        let mut exp_hdu = BinTableBuilder::new("EXPOSURE")
            .add_f64("EXPOSURE_Q0", None, pexp[0].to_vec())
            .add_f64("EXPOSURE_Q1", None, pexp[1].to_vec())
            .add_f64("EXPOSURE_Q2", None, pexp[2].to_vec())
            .add_f64("EXPOSURE_Q3", None, pexp[3].to_vec())
            .build()?;
        if let Some(nc) = opts.ncompat {
            exp_hdu.header.set_i64("NCOMPAT", nc as i64, "Noise clean compat");
        }
        hdus.push(exp_hdu);
    }

    write_fits(path, &mut hdus, true)?;
    Ok(())
}

/// Read `Q0..Q3` event tables, optional `EXPOSURE`, `GTI`s, and `GAINCOR` flag from a Level-2 FITS file.
#[allow(clippy::type_complexity)]
pub fn read_l2_event_file(
    path: impl AsRef<Path>,
) -> Result<(
    CztiObsMeta,
    [QuadEvents; NUM_QUADS],
    CztiGtiOutput,
    Option<[Vec<f64>; NUM_QUADS]>,
    bool,
)> {
    let meta = CztiObsMeta::from_fits(&path)?;
    let fits = FitsFile::open(&path)?;
    let gaincor_applied = fits
        .hdu(0)?
        .header()
        .get_str("GAINCOR")
        .map(|s| s.trim().eq_ignore_ascii_case("YES"))
        .unwrap_or(false);

    let mut q_vec = Vec::with_capacity(NUM_QUADS);
    for q in 0..NUM_QUADS {
        let hdu = fits.hdu_by_name(&format!("Q{q}"))?;
        let bt = hdu.as_bintable()?;
        let time = bt.read_col_f64("TIME")?;
        let cztseccnt = bt.read_col_f64("CZTSECCNT")?;
        let cztntick = bt.read_col_u16("CZTNTICK")?;
        let pha = bt.read_col_u16("PHA")?;
        let detid = bt.read_col_u8("DETID")?;
        let pixid = bt.read_col_u8("PIXID")?;
        let detx = bt.read_col_u8("DETX")?;
        let dety = bt.read_col_u8("DETY")?;
        let veto = bt.read_col_u16("VETO")?;
        let alpha = bt.read_col_u8("ALPHA")?;
        let pi = bt.read_col_i16("PI").unwrap_or_default();
        let energy = bt.read_col_f32("ENERGY").unwrap_or_default();
        let exposure = hdu.header().get_f64("EXPOSURE").unwrap_or(0.0);
        let tstart = hdu
            .header()
            .get_f64("TSTART")
            .unwrap_or_else(|| time.first().copied().unwrap_or(meta.tstart));
        let tstop = hdu
            .header()
            .get_f64("TSTOP")
            .unwrap_or_else(|| time.last().copied().unwrap_or(meta.tstop));

        q_vec.push(QuadEvents {
            quad_id: q as u8,
            time,
            cztseccnt,
            cztntick,
            pha,
            detid,
            pixid,
            detx,
            dety,
            veto,
            alpha,
            pi,
            energy,
            open_fraction: bt.read_col_f32("OPEN_FRACTION").ok(),
            weight: bt.read_col_f32("WEIGHT").ok(),
            tstart,
            tstop,
            exposure,
        });
    }

    let (gtis, _) = read_czti_gtis(&path)?;

    let pixel_exposure = fits.hdu_by_name("EXPOSURE").ok().and_then(|h| {
        let bt = h.as_bintable().ok()?;
        let mut arr: [Vec<f64>; NUM_QUADS] = [
            vec![1.0f64; PIXELS_PER_QUAD],
            vec![1.0f64; PIXELS_PER_QUAD],
            vec![1.0f64; PIXELS_PER_QUAD],
            vec![1.0f64; PIXELS_PER_QUAD],
        ];
        for (q, q_arr) in arr.iter_mut().enumerate() {
            if let Ok(col) = bt.read_col_f64(&format!("EXPOSURE_Q{q}")) {
                for (i, &v) in col.iter().take(PIXELS_PER_QUAD).enumerate() {
                    q_arr[i] = v;
                }
            }
        }
        Some(arr)
    });

    let mut it = q_vec.into_iter();
    let quads = [
        it.next().unwrap(),
        it.next().unwrap(),
        it.next().unwrap(),
        it.next().unwrap(),
    ];
    Ok((meta, quads, gtis, pixel_exposure, gaincor_applied))
}

/// Write `*_bunch.fits` (`Q0..Q3` bunch tables).
pub fn write_bunch_file(
    path: impl AsRef<Path>,
    meta: &CztiObsMeta,
    bunches: &[BunchTable; NUM_QUADS],
) -> Result<()> {
    let mut hdus = Vec::with_capacity(5);
    let mut primary = Hdu::empty_primary();
    apply_ogip_meta(&mut primary.header, &meta.to_ogip(None, 0.0));
    hdus.push(primary);

    for (q, b) in bunches.iter().enumerate() {
        let hdu = BinTableBuilder::new(format!("Q{q}"))
            .add_f64("Time", Some("s"), b.time.clone())
            .add_u8("Time_dfs", None, b.time_dfs.clone())
            .add_u8("Time_dsl", None, b.time_dsl.clone())
            .add_u8("NumEvent", None, b.num_event.clone())
            .add_u8("DetId1", None, b.detid1.clone())
            .add_u8("DetId2", None, b.detid2.clone())
            .add_u8("DetId3", None, b.detid3.clone())
            .add_u8("DetId4", None, b.detid4.clone())
            .add_i32("RowNum", None, b.row_num.clone())
            .add_u8("DetId_fevt", None, b.detid_fevt.clone())
            .add_u8("PixId_fevt", None, b.pixid_fevt.clone())
            .add_u8("DetId_sevt", None, b.detid_sevt.clone())
            .add_u8("PixId_sevt", None, b.pixid_sevt.clone())
            .add_u8("DetId_levt", None, b.detid_levt.clone())
            .add_u8("PixId_levt", None, b.pixid_levt.clone())
            .header_i64("QUADID", q as i64, "Quadrant Number")
            .build()?;
        hdus.push(hdu);
    }

    write_fits(path, &mut hdus, true)?;
    Ok(())
}

/// Write `*_level2.hdr` frame housekeeping file.
pub fn write_hdr_file(
    path: impl AsRef<Path>,
    meta: &CztiObsMeta,
    hdr: &FrameHeaderTable,
) -> Result<()> {
    let mut primary = Hdu::empty_primary();
    apply_ogip_meta(&mut primary.header, &meta.to_ogip(None, 0.0));

    let hdu = BinTableBuilder::new("HEADER")
        .add_f64("TIME", Some("s"), hdr.time.clone())
        .add_f64("CZTSECCNT", Some("s"), hdr.cztseccnt.clone())
        .add_u8("DATAID", None, hdr.data_id.clone())
        .add_u8("MODEID", None, hdr.mode_id.clone())
        .add_i32("WORDCOUNT", None, hdr.word_count.clone())
        .add_u16("FRAMENO", None, hdr.frame_no.clone())
        .add_u16("STATUSFLAG", None, hdr.status_flag.clone())
        .add_u16("WRITEPKTNO", None, hdr.write_pkt_no.clone())
        .add_u16("READPKTNO", None, hdr.read_pkt_no.clone())
        .add_i32("COMMAND", None, hdr.command.clone())
        .add_f64("LAXPCTIME", Some("s"), hdr.laxpc_time.clone())
        .add_u16("DCNT", None, hdr.dcnt.clone())
        .add_u8("ERRORCOUNT", None, hdr.error_count.clone())
        .add_u8("ERRORFLAG", None, hdr.error_flag.clone())
        .add_u8("BOOTPAGENO", None, hdr.boot_page_no.clone())
        .add_u8("CZTNO", None, hdr.czt_no.clone())
        .add_u8("CZTSTATUS", None, hdr.czt_status.clone())
        .add_u8("EVTREADMODE", None, hdr.evt_read_mode.clone())
        .add_u8("CMDSTATUS", None, hdr.cmd_status.clone())
        .add_u8("BUFFERNO", None, hdr.buffer_no.clone())
        .add_u8("BASEADD", None, hdr.base_add.clone())
        .add_u8("VETOSPECRANGE", None, hdr.veto_spec_range.clone())
        .add_u8("CHANNELNO", None, hdr.channel_no.clone())
        .add_u16("ADCOUTPUT", None, hdr.adc_output.clone())
        .add_u16("CMD1SEC", None, hdr.cmd_1sec.clone())
        .add_u16("ALPHACOUNT", None, hdr.alpha_count.clone())
        .add_u16("VETOCOUNT", None, hdr.veto_count.clone())
        .add_u16("CZTCOUNT_LT_ULD", None, hdr.czt_count_lt_uld.clone())
        .add_u16("CZTCOUNT_GT_ULD", None, hdr.czt_count_gt_uld.clone())
        .add_u16("CZTDATAREAD", None, hdr.czt_data_read.clone())
        .build()?;

    write_fits(path, &mut [primary, hdu], true)?;
    Ok(())
}

/// Write `*_livetime.fits` (`Q0..Q3` with `TIME` and `FRACEXP`).
pub fn write_livetime_file(
    path: impl AsRef<Path>,
    meta: &CztiObsMeta,
    livetimes: &[QuadLivetime; NUM_QUADS],
) -> Result<()> {
    let mut hdus = Vec::with_capacity(5);
    let mut primary = Hdu::empty_primary();
    apply_ogip_meta(&mut primary.header, &meta.to_ogip(None, 0.0));
    primary
        .header
        .set_f64("LV_BINSIZE", livetimes[0].binsize, "Livetime bin size [s]");
    hdus.push(primary);

    for (q, lv) in livetimes.iter().enumerate() {
        let hdu = BinTableBuilder::new(format!("Q{q}"))
            .add_f64("TIME", Some("s"), lv.time.clone())
            .add_f64("FRACEXP", None, lv.fracexp.clone())
            .header_f64("LV_BINSIZE", lv.binsize, "Livetime bin size [s]")
            .build()?;
        hdus.push(hdu);
    }

    write_fits(path, &mut hdus, true)?;
    Ok(())
}

/// Write `*_level2.gti` (`GTI`, `Q0_GTI..Q3_GTI`).
pub fn write_gti_file(
    path: impl AsRef<Path>,
    meta: &CztiObsMeta,
    gtis: &CztiGtiOutput,
) -> Result<()> {
    let mut hdus = Vec::with_capacity(6);
    let mut primary = Hdu::empty_primary();
    apply_ogip_meta(&mut primary.header, &meta.to_ogip(None, 0.0));
    hdus.push(primary);

    let add = |hdus: &mut Vec<Hdu>, name: &str, ivs: &[(f64, f64)]| -> Result<()> {
        let starts: Vec<f64> = ivs.iter().map(|&(a, _)| a).collect();
        let stops: Vec<f64> = ivs.iter().map(|&(_, b)| b).collect();
        let exp: f64 = ivs.iter().map(|&(a, b)| (b - a).max(0.0)).sum();
        hdus.push(build_gti_hdu(name, &starts, &stops, &meta.to_ogip(None, exp))?);
        Ok(())
    };
    add(&mut hdus, "GTI", &gtis.common)?;
    for (q, q_ivs) in gtis.quad.iter().enumerate() {
        add(&mut hdus, &format!("Q{q}_GTI"), q_ivs)?;
    }
    write_fits(path, &mut hdus, true)?;
    Ok(())
}

/// Write per-quadrant `*.aspect_Q{qid}` FITS file (`ASPECT` HDU + `RA`/`DEC`/`TWIST` cards).
pub fn write_aspect_file(
    path: impl AsRef<Path>,
    meta: &CztiObsMeta,
    aspect: &QuadAspect,
) -> Result<()> {
    let mut primary = Hdu::empty_primary();
    apply_ogip_meta(&mut primary.header, &meta.to_ogip(Some(format!("QUADRANT{}", aspect.quad_id)), 0.0));

    let mut hdu = BinTableBuilder::new("ASPECT")
        .add_f64("TIME", Some("s"), aspect.time.clone())
        .add_f32("Nx", None, aspect.nx.clone())
        .add_f32("Ny", None, aspect.ny.clone())
        .add_f32("Nz", None, aspect.nz.clone())
        .add_f32("Nxt", None, aspect.nxt.clone())
        .add_f32("Nyt", None, aspect.nyt.clone())
        .add_f32("Nzt", None, aspect.nzt.clone())
        .header_i64("QUADID", aspect.quad_id as i64, "Quadrant Number")
        .header_str("TELDEF", &aspect.teldef_name, "Teldef file name")
        .header_f64("RA", aspect.ra_deg, "[deg] Average pointing RA")
        .header_f64("DEC", aspect.dec_deg, "[deg] Average pointing DEC")
        .header_f64("TWIST", aspect.twist_deg, "[deg] Average pointing TWIST")
        .header_str("HDUCLASS", "OGIP", "OGIP standard")
        .build()?;
    apply_ogip_meta(&mut hdu.header, &meta.to_ogip(Some(format!("QUADRANT{}", aspect.quad_id)), 0.0));

    write_fits(path, &mut [primary, hdu], true)?;
    Ok(())
}

/// Write `*_badpix.fits` (`Q0..Q3` tables + 128x128 `BADPIX` image extension).
pub fn write_badpix_file(
    path: impl AsRef<Path>,
    meta: &CztiObsMeta,
    badpix: &BadpixTable,
) -> Result<()> {
    let mut hdus = Vec::with_capacity(6);
    let mut primary = Hdu::empty_primary();
    apply_ogip_meta(&mut primary.header, &meta.to_ogip(None, 0.0));
    hdus.push(primary);

    for q in 0..NUM_QUADS {
        let mut detid = Vec::with_capacity(PIXELS_PER_QUAD);
        let mut pixid = Vec::with_capacity(PIXELS_PER_QUAD);
        let mut pixx = Vec::with_capacity(PIXELS_PER_QUAD);
        let mut pixy = Vec::with_capacity(PIXELS_PER_QUAD);
        let mut flags = Vec::with_capacity(PIXELS_PER_QUAD);
        for d in 0..16u8 {
            for p in 0..=255u8 {
                let (dx, dy) = detid_pixid_to_detxy(q as u8, d, p);
                detid.push(d);
                pixid.push(p);
                pixx.push(dx);
                pixy.push(dy);
                flags.push(badpix.get(q, d, p));
            }
        }
        let hdu = BinTableBuilder::new(format!("Q{q}"))
            .add_u8("DetID", None, detid)
            .add_u8("PixID", None, pixid)
            .add_u8("PixX", None, pixx)
            .add_u8("PixY", None, pixy)
            .add_u8("PIX_FLAG", None, flags)
            .build()?;
        hdus.push(hdu);
    }

    let full_map = badpix.to_full_map();
    let flat_map: Vec<u8> = full_map.into_iter().flatten().collect();
    let img_hdu = ImageBuilder::from_u8(Some("BADPIX"), vec![FULL_DIM, FULL_DIM], &flat_map)?
        .build_extension();
    hdus.push(img_hdu);

    write_fits(path, &mut hdus, true)?;
    Ok(())
}

/// Write `*.dph` (`Q0..Q3` 64x64 + `FULL_DPH` 128x128 I32 IMAGE HDUs).
pub fn write_dph_file(
    path: impl AsRef<Path>,
    meta: &CztiObsMeta,
    gti_type: GtiType,
    quad_dph: &[Vec<Vec<i32>>; NUM_QUADS],
    full_dph: &[Vec<i32>],
) -> Result<()> {
    let mut hdus = Vec::with_capacity(6);
    let mut primary = Hdu::empty_primary();
    apply_ogip_meta(&mut primary.header, &meta.to_ogip(None, 0.0));
    primary
        .header
        .set_str("GTITYPE", gti_type.as_fits_str(), "GTI selection mode");
    hdus.push(primary);

    for (q, qimg) in quad_dph.iter().enumerate() {
        let flat: Vec<i32> = qimg.iter().flat_map(|r| r.iter().copied()).collect();
        let mut hdu = ImageBuilder::from_i32(Some(&format!("Q{q}")), vec![QUAD_DIM, QUAD_DIM], &flat)?
            .build_extension();
        apply_ogip_meta(&mut hdu.header, &meta.to_ogip(Some(format!("QUADRANT{q}")), 0.0));
        hdus.push(hdu);
    }
    let flat_full: Vec<i32> = full_dph.iter().flat_map(|r| r.iter().copied()).collect();
    let mut full_hdu =
        ImageBuilder::from_i32(Some("FULL_DPH"), vec![FULL_DIM, FULL_DIM], &flat_full)?
            .build_extension();
    apply_ogip_meta(&mut full_hdu.header, &meta.to_ogip(None, 0.0));
    hdus.push(full_hdu);

    write_fits(path, &mut hdus, true)?;
    Ok(())
}

/// Write `*.dpi` (`Q0..Q3` 64x64 + `FULL_DPI` 128x128 F32 IMAGE HDUs).
pub fn write_dpi_file(
    path: impl AsRef<Path>,
    meta: &CztiObsMeta,
    gti_type: GtiType,
    quad_dpi: &[Vec<Vec<f32>>; NUM_QUADS],
    full_dpi: &[Vec<f32>],
) -> Result<()> {
    let mut hdus = Vec::with_capacity(6);
    let mut primary = Hdu::empty_primary();
    apply_ogip_meta(&mut primary.header, &meta.to_ogip(None, 0.0));
    primary
        .header
        .set_str("GTITYPE", gti_type.as_fits_str(), "GTI selection mode");
    hdus.push(primary);

    for (q, qimg) in quad_dpi.iter().enumerate() {
        let flat: Vec<f32> = qimg.iter().flat_map(|r| r.iter().copied()).collect();
        let mut hdu = ImageBuilder::from_f32(Some(&format!("Q{q}")), vec![QUAD_DIM, QUAD_DIM], &flat)?
            .build_extension();
        apply_ogip_meta(&mut hdu.header, &meta.to_ogip(Some(format!("QUADRANT{q}")), 0.0));
        hdus.push(hdu);
    }
    let flat_full: Vec<f32> = full_dpi.iter().flat_map(|r| r.iter().copied()).collect();
    let mut full_hdu =
        ImageBuilder::from_f32(Some("FULL_DPI"), vec![FULL_DIM, FULL_DIM], &flat_full)?
            .build_extension();
    apply_ogip_meta(&mut full_hdu.header, &meta.to_ogip(None, 0.0));
    hdus.push(full_hdu);

    write_fits(path, &mut hdus, true)?;
    Ok(())
}

/// Write a reconstructed sky image (`*.img_Q{q}` or `*.img`) with full CZTI ICRS + local camera WCS.
pub fn write_sky_image_file(
    path: impl AsRef<Path>,
    extname: &str,
    meta: &CztiObsMeta,
    img: &SkyImage,
) -> Result<()> {
    let mut primary = Hdu::empty_primary();
    apply_ogip_meta(&mut primary.header, &meta.to_ogip(None, 0.0));

    let w = &img.wcs;
    let mut ext = ImageBuilder::from_f32(Some(extname), vec![img.width, img.height], &img.data)?
        .header_str("CTYPE1", &w.ctype1, "ICRS RA TAN projection")
        .header_str("CTYPE2", &w.ctype2, "ICRS DEC TAN projection")
        .header_f64("CRPIX1", w.crpix1, "Reference pixel X")
        .header_f64("CRPIX2", w.crpix2, "Reference pixel Y")
        .header_f64("CRVAL1", w.crval1, "[deg] Reference RA")
        .header_f64("CRVAL2", w.crval2, "[deg] Reference DEC")
        .header_f64("CDELT1", w.cdelt1, "[deg/pix] Pixel scale X")
        .header_f64("CDELT2", w.cdelt2, "[deg/pix] Pixel scale Y")
        .header_f64("CROTA1", w.crota1, "[deg] Rotation angle X")
        .header_f64("CROTA2", w.crota2, "[deg] Rotation angle Y")
        .header_f64("CD1_1", w.cd1_1, "WCS CD matrix 1,1")
        .header_f64("CD1_2", w.cd1_2, "WCS CD matrix 1,2")
        .header_f64("CD2_1", w.cd2_1, "WCS CD matrix 2,1")
        .header_f64("CD2_2", w.cd2_2, "WCS CD matrix 2,2")
        .header_str("RADESYS", "ICRS", "Coordinate system")
        .header_f64("EQUINOX", 2000.0, "Equinox")
        .header_str("CTYPE1B", "TX---TAN", "Local camera tangent X")
        .header_str("CTYPE2B", "TY---TAN", "Local camera tangent Y")
        .header_f64("CRPIX1B", w.crpix1b, "Camera ref pixel X")
        .header_f64("CRPIX2B", w.crpix2b, "Camera ref pixel Y")
        .header_f64("CRVAL1B", w.crval1b, "[deg] Camera ref X")
        .header_f64("CRVAL2B", w.crval2b, "[deg] Camera ref Y")
        .header_f64("CDELT1B", w.cdelt1b, "[deg/pix] Camera scale X")
        .header_f64("CDELT2B", w.cdelt2b, "[deg/pix] Camera scale Y")
        .header_f64("CD1_1B", w.cd1_1b, "Camera CD 1,1")
        .header_f64("CD1_2B", w.cd1_2b, "Camera CD 1,2")
        .header_f64("CD2_1B", w.cd2_1b, "Camera CD 2,1")
        .header_f64("CD2_2B", w.cd2_2b, "Camera CD 2,2")
        .build_extension();
    apply_ogip_meta(&mut ext.header, &meta.to_ogip(None, 0.0));

    write_fits(path, &mut [primary, ext], true)?;
    Ok(())
}

/// Write an OGIP 1.2.0 CZTI `.pha` SPECTRUM file (matching `spectrumTemplate`).
pub fn write_spectrum_file(
    path: impl AsRef<Path>,
    meta: &CztiObsMeta,
    spec: &CztiSpectrum,
) -> Result<()> {
    let detnam = spec.quad_id.map(|q| format!("QUADRANT{q}"));
    let ogip = meta.to_ogip(detnam, spec.exposure);

    let mut primary = Hdu::empty_primary();
    apply_ogip_meta(&mut primary.header, &ogip);
    primary.header.set_str("CONTENT", "SPECTRUM", "File content");

    let mut hdu = BinTableBuilder::new("SPECTRUM")
        .add_i32("CHANNEL", None, spec.channel.clone())
        .with_last_col_range(0, (PI_CHANNELS - 1) as i64)
        .add_f32("COUNTS", Some("count"), spec.counts.clone())
        .add_f32("STAT_ERR", Some("count"), spec.stat_err.clone())
        .add_f32("SYS_ERR", None, spec.sys_err.clone())
        .add_i16("QUALITY", None, spec.quality.clone())
        .add_i16("GROUPING", None, spec.grouping.clone())
        .header_str("HDUCLASS", "OGIP", "OGIP Standard")
        .header_str("HDUCLAS1", "SPECTRUM", "PHA dataset")
        .header_str("HDUVERS1", "1.2.0", "OGIP SPECTRUM version")
        .header_str("HDUVERS", "1.2.0", "OGIP SPECTRUM version")
        .header_str("HDUCLAS2", "TOTAL", "total spectrum")
        .header_str("HDUCLAS3", "COUNT", "units in counts")
        .header_str("CHANTYPE", "PI", "Pulse-Invariant channel")
        .header_i64("DETCHANS", PI_CHANNELS as i64, "Number of PI channels")
        .header_f64("AREASCAL", 1.0, "Area scaling factor")
        .header_f64("BACKSCAL", 1.0, "Background scaling factor")
        .header_f64("CORRSCAL", 1.0, "Correction scaling factor")
        .header_str("BACKFILE", "none", "Background file")
        .header_str("RESPFILE", &spec.respfile, "Response matrix file")
        .header_str("ANCRFILE", "none", "Ancillary response file")
        .header_bool("POISSERR", false, "Error in STAT_ERR column")
        .header_f64("THETAX", spec.thetax_deg, "[deg] Camera coordinate X")
        .header_f64("THETAY", spec.thetay_deg, "[deg] Camera coordinate Y")
        .header_i64("MASKWT", if spec.mask_weighted { 1 } else { 0 }, "Mask weighted")
        .header_i64("PIXTHRES", spec.badpix_threshold as i64, "Badpix threshold")
        .header_i64("QUADID", spec.quad_id.map(|q| q as i64).unwrap_or(-1), "Quadrant ID")
        .header_str("GTITYPE", spec.gti_type.as_fits_str(), "GTI type")
        .build()?;
    apply_ogip_meta(&mut hdu.header, &ogip);

    write_fits(path, &mut [primary, hdu], true)?;
    Ok(())
}

/// Write an OGIP CZTI `.lc` LIGHTCURVE file (matching `lcTemplate`).
pub fn write_lightcurve_file(
    path: impl AsRef<Path>,
    meta: &CztiObsMeta,
    lc: &CztiLightCurve,
) -> Result<()> {
    let detnam = lc.quad_id.map(|q| format!("QUADRANT{q}"));
    let mut ogip = meta.to_ogip(detnam, (lc.tstop - lc.tstart).max(0.0));
    ogip.tstart = lc.tstart;
    ogip.tstop = lc.tstop;

    let mut primary = Hdu::empty_primary();
    apply_ogip_meta(&mut primary.header, &ogip);
    primary
        .header
        .set_str("CONTENT", "LIGHT CURVE", "File content");

    let mut lc_hdu = build_lightcurve_hdu(
        "LIGHTCURVE",
        &lc.time,
        &lc.rate,
        &lc.error,
        Some(&lc.fracexp),
        lc.timedel,
        &ogip,
    )?;
    lc_hdu
        .header
        .set_f64("THETAX", lc.thetax_deg, "[deg] Camera coordinate X");
    lc_hdu
        .header
        .set_f64("THETAY", lc.thetay_deg, "[deg] Camera coordinate Y");
    lc_hdu.header.set_i64(
        "MASKWT",
        if lc.mask_weighted { 1 } else { 0 },
        "Maskweight applied (1) or not (0)",
    );
    lc_hdu
        .header
        .set_i64("PIXTHRES", lc.badpix_threshold as i64, "Badpix threshold");
    lc_hdu.header.set_str(
        "ERANGE",
        format!("{:.1}-{:.1}", lc.emin_kev, lc.emax_kev),
        "Energy range [keV]",
    );

    write_fits(path, &mut [primary, lc_hdu], true)?;
    Ok(())
}

/// Write an OGIP 1.3.0 CZTI `.rsp` RESPONSE file (`EBOUNDS` + `MATRIX` HDUs, compressed at `lo_thresh = 1e-10`).
pub fn write_response_file(
    path: impl AsRef<Path>,
    meta: &CztiObsMeta,
    rsp: &CztiResponse,
) -> Result<()> {
    let detnam = rsp.quad_id.map(|q| format!("QUADRANT{q}"));
    let ogip = meta.to_ogip(detnam, 0.0);
    let lo_thresh = 1.0e-10f64;

    let mut n_grp = Vec::with_capacity(rsp.n_ebins);
    let mut f_chan = Vec::with_capacity(rsp.n_ebins);
    let mut n_chan = Vec::with_capacity(rsp.n_ebins);
    let mut matrix_rows = Vec::with_capacity(rsp.n_ebins);

    for eb in 0..rsp.n_ebins {
        let row = &rsp.matrix[eb * rsp.n_pibins..(eb + 1) * rsp.n_pibins];
        // Run-length sparse compression above lo_thresh
        let mut fc = Vec::new();
        let mut nc = Vec::new();
        let mut vals = Vec::new();
        let mut ch = 0usize;
        while ch < rsp.n_pibins {
            if (row[ch].abs() as f64) > lo_thresh {
                let start = ch;
                while ch < rsp.n_pibins && (row[ch].abs() as f64) > lo_thresh {
                    vals.push(row[ch]);
                    ch += 1;
                }
                fc.push(start as i16);
                nc.push((ch - start) as i16);
            } else {
                ch += 1;
            }
        }
        if fc.is_empty() {
            // Keep at least one 0-value channel entry for XSPEC compatibility
            fc.push(0);
            nc.push(1);
            vals.push(0.0);
        }
        n_grp.push(fc.len() as i16);
        f_chan.push(fc);
        n_chan.push(nc);
        matrix_rows.push(vals);
    }

    let rmf = ResponseMatrixData {
        chantype: "PI".to_string(),
        detchans: rsp.n_pibins as i64,
        tlmin4: 0,
        energ_lo: rsp.energ_lo.clone(),
        energ_hi: rsp.energ_hi.clone(),
        n_grp,
        f_chan,
        n_chan,
        matrix: matrix_rows,
        channels: rsp.ebounds_channel.clone(),
        e_min: rsp.ebounds_emin.clone(),
        e_max: rsp.ebounds_emax.clone(),
        lo_thresh,
    };

    let mut primary = Hdu::empty_primary();
    apply_ogip_meta(&mut primary.header, &ogip);
    let (matrix_hdu, ebounds_hdu) = build_response_hdus(&rmf, &ogip)?;

    write_fits(path, &mut [primary, ebounds_hdu, matrix_hdu], true)?;
    Ok(())
}

/// Read `TEMP` extension from a `modeSS` Level-2 FITS file.
pub fn read_temp_extension(path: impl AsRef<Path>) -> Result<TempTable> {
    let fits = FitsFile::open(path)?;
    let hdu = fits.hdu_by_name("TEMP")?;
    let bt = hdu.as_bintable()?;
    let time = bt.read_col_f64("TIME")?;
    let cztseccnt = bt.read_col_f64("CZTSECCNT").unwrap_or_else(|_| vec![0.0; time.len()]);
    let quad_flag = bt.read_col_u8("QUAD_FLAG")?;
    let (rep, flat) = bt.read_vec_f32("TEMPERATURE")?;
    let mut temperature = Vec::with_capacity(time.len());
    for r in 0..time.len() {
        let mut arr = [10.0f32; NUM_DET_PER_QUAD];
        for d in 0..NUM_DET_PER_QUAD.min(rep) {
            arr[d] = flat[r * rep + d];
        }
        temperature.push(arr);
    }
    Ok(TempTable {
        time,
        cztseccnt,
        temperature,
        quad_flag,
    })
}

/// Read auxiliary `VETOSPECTRUM`, `SSM Data`, and `TEMP` HDUs from an existing CZTI Level-2 FITS
/// file so standalone stage commands preserve all telemetry extensions.
pub fn read_aux_tables(path: impl AsRef<Path>) -> (VetoSpectrumTable, SsmDataTable, TempTable) {
    let Ok(fits) = FitsFile::open(&path) else {
        return (
            VetoSpectrumTable::default(),
            SsmDataTable::default(),
            TempTable::default(),
        );
    };

    let mut veto = VetoSpectrumTable::default();
    if let Ok(hdu) = fits.hdu_by_name("VETOSPECTRUM") {
        if let Ok(bt) = hdu.as_bintable() {
            if let (Ok(time), Ok(czt), Ok((rep, flat)), Ok(qid)) = (
                bt.read_col_f64("Time").or_else(|_| bt.read_col_f64("TIME")),
                bt.read_col_f64("CZTSECCNT"),
                bt.read_vec_u16("VetoSpec"),
                bt.read_col_u8("QuadID"),
            ) {
                let nrows = time.len().min(czt.len()).min(qid.len());
                veto.time = time[..nrows].to_vec();
                veto.cztseccnt = czt[..nrows].to_vec();
                veto.quad_id = qid[..nrows].to_vec();
                veto.veto_spec.reserve(nrows);
                for r in 0..nrows {
                    let mut row = [0u16; VETOSPEC_SIZE];
                    let n = rep.min(VETOSPEC_SIZE);
                    if (r + 1) * rep <= flat.len() {
                        row[..n].copy_from_slice(&flat[r * rep..r * rep + n]);
                    }
                    veto.veto_spec.push(row);
                }
            }
        }
    }

    let mut ssm = SsmDataTable::default();
    if let Ok(hdu) = fits.hdu_by_name("SSM Data") {
        if let Ok(bt) = hdu.as_bintable() {
            if let (
                Ok(time),
                Ok(czt),
                Ok((rv, fv)),
                Ok((rc, fc)),
                Ok((rcv, fcv)),
                Ok((rca, fca)),
                Ok(qid),
            ) = (
                bt.read_col_f64("Time").or_else(|_| bt.read_col_f64("TIME")),
                bt.read_col_f64("CZTSECCNT"),
                bt.read_vec_u16("VetoSpec"),
                bt.read_vec_u16("CZTSpec"),
                bt.read_vec_u16("CZTSpecVeto"),
                bt.read_vec_u16("CZTSpecAlpha"),
                bt.read_col_u8("QuadID"),
            ) {
                let nrows = time.len().min(czt.len()).min(qid.len());
                ssm.time = time[..nrows].to_vec();
                ssm.cztseccnt = czt[..nrows].to_vec();
                ssm.quad_id = qid[..nrows].to_vec();
                for r in 0..nrows {
                    let mut vs = [0u16; VETOSPEC_SIZE];
                    if (r + 1) * rv <= fv.len() {
                        let n = rv.min(VETOSPEC_SIZE);
                        vs[..n].copy_from_slice(&fv[r * rv..r * rv + n]);
                    }
                    let mut cs = [0u16; CZTSPEC_SIZE];
                    if (r + 1) * rc <= fc.len() {
                        let n = rc.min(CZTSPEC_SIZE);
                        cs[..n].copy_from_slice(&fc[r * rc..r * rc + n]);
                    }
                    let mut csv = [0u16; CZTSPEC_SIZE];
                    if (r + 1) * rcv <= fcv.len() {
                        let n = rcv.min(CZTSPEC_SIZE);
                        csv[..n].copy_from_slice(&fcv[r * rcv..r * rcv + n]);
                    }
                    let mut csa = [0u16; CZTSPEC_SIZE];
                    if (r + 1) * rca <= fca.len() {
                        let n = rca.min(CZTSPEC_SIZE);
                        csa[..n].copy_from_slice(&fca[r * rca..r * rca + n]);
                    }
                    ssm.veto_spec.push(vs);
                    ssm.czt_spec.push(cs);
                    ssm.czt_spec_veto.push(csv);
                    ssm.czt_spec_alpha.push(csa);
                }
            }
        }
    }

    let temp = read_temp_extension(path).unwrap_or_default();
    (veto, ssm, temp)
}

/// Read an OGIP Type-I CZTI `.pha` spectrum file (`SPECTRUM` HDU).
pub fn read_spectrum_file(path: impl AsRef<Path>) -> Result<(CztiObsMeta, CztiSpectrum)> {
    let meta = CztiObsMeta::from_fits(&path)?;
    let fits = FitsFile::open(&path)?;
    let hdu = fits.hdu_by_name("SPECTRUM").or_else(|_| fits.hdu(1))?;
    let hdr = hdu.header();
    let bt = hdu.as_bintable()?;
    let channel = bt.read_col_i32("CHANNEL")?;
    let n = channel.len();
    let counts = bt
        .read_col_f32("COUNTS")
        .or_else(|_| {
            bt.read_col_f32("RATE").map(|r| {
                let exp = hdr.get_f64("EXPOSURE").unwrap_or(1.0) as f32;
                r.into_iter().map(|v| v * exp).collect()
            })
        })?;
    let stat_err = bt
        .read_col_f32("STAT_ERR")
        .unwrap_or_else(|_| vec![0.0; n]);
    let sys_err = bt.read_col_f32("SYS_ERR").unwrap_or_else(|_| vec![0.0; n]);
    let quality = bt.read_col_i16("QUALITY").unwrap_or_else(|_| vec![0; n]);
    let grouping = bt.read_col_i16("GROUPING").unwrap_or_else(|_| vec![1; n]);
    let exposure = hdr.get_f64("EXPOSURE").unwrap_or(0.0);
    let thetax_deg = hdr.get_f64("THETAX").unwrap_or(0.0);
    let thetay_deg = hdr.get_f64("THETAY").unwrap_or(0.0);
    let mask_weighted = hdr
        .get_str("MASKWT")
        .map(|s| s.trim().eq_ignore_ascii_case("YES"))
        .unwrap_or(true);
    let badpix_threshold = hdr.get_i64("BADPIXTH").unwrap_or(0).clamp(0, 255) as u8;
    let gti_type = hdr
        .get_str("GTITYPE")
        .and_then(GtiType::parse)
        .unwrap_or(GtiType::Quad);
    let respfile = hdr
        .get_str("RESPFILE")
        .unwrap_or("NONE")
        .trim()
        .to_string();
    let tstart = hdr.get_f64("TSTART").unwrap_or(meta.tstart);
    let tstop = hdr.get_f64("TSTOP").unwrap_or(meta.tstop);
    let quad_id = hdr.get_i64("QUADID").and_then(|v| u8::try_from(v).ok());

    Ok((
        meta,
        CztiSpectrum {
            quad_id,
            channel,
            counts,
            stat_err,
            sys_err,
            quality,
            grouping,
            exposure,
            thetax_deg,
            thetay_deg,
            mask_weighted,
            badpix_threshold,
            gti_type,
            respfile,
            tstart,
            tstop,
        },
    ))
}

/// Read an OGIP CZTI `.rsp` response file (`MATRIX` + `EBOUNDS` HDUs).
pub fn read_response_file(path: impl AsRef<Path>) -> Result<(CztiObsMeta, CztiResponse)> {
    let meta = CztiObsMeta::from_fits(&path)?;
    let fits = FitsFile::open(&path)?;
    let mhdu = fits
        .hdu_by_name("MATRIX")
        .or_else(|_| fits.hdu_by_name("SPECRESP MATRIX"))
        .or_else(|_| fits.hdu(1))?;
    let mbt = mhdu.as_bintable()?;
    let energ_lo = mbt.read_col_f32("ENERG_LO")?;
    let energ_hi = mbt.read_col_f32("ENERG_HI")?;
    let n_ebins = energ_lo.len();
    let raw_detchans = mhdu.header().get_i64("DETCHANS").unwrap_or(PI_CHANNELS as i64);
    if raw_detchans <= 0 {
        return Err(CztiError::InvalidConfig(format!(
            "Invalid DETCHANS={raw_detchans} in response file"
        )));
    }
    let n_pibins = raw_detchans as usize;
    const MAX_RSP_MATRIX_ELEMS: usize = 64_000_000;
    let total_elems = n_ebins.checked_mul(n_pibins).ok_or_else(|| {
        CztiError::InvalidConfig(format!(
            "Response matrix size overflow: {n_ebins} x {n_pibins}"
        ))
    })?;
    if total_elems > MAX_RSP_MATRIX_ELEMS {
        return Err(CztiError::InvalidConfig(format!(
            "Response matrix too large ({n_ebins} x {n_pibins} = {total_elems} > {MAX_RSP_MATRIX_ELEMS})"
        )));
    }

    let mut matrix = vec![0.0f32; total_elems];
    if let Ok((rep, flat)) = mbt.read_vec_f32("MATRIX") {
        let cols = rep.min(n_pibins);
        for r in 0..n_ebins {
            if (r + 1) * rep <= flat.len() {
                matrix[r * n_pibins..r * n_pibins + cols]
                    .copy_from_slice(&flat[r * rep..r * rep + cols]);
            }
        }
    } else if let Ok(ngrp) = mbt.read_col_i16("N_GRP") {
        for r in 0..n_ebins {
            let ng = ngrp.get(r).copied().unwrap_or(0).max(0) as usize;
            let fcs = mbt.read_var_i16("F_CHAN", r).unwrap_or_default();
            let ncs = mbt.read_var_i16("N_CHAN", r).unwrap_or_default();
            let mvs = mbt.read_var_f32("MATRIX", r).unwrap_or_default();
            let mut off = 0usize;
            for g in 0..ng.min(fcs.len()).min(ncs.len()) {
                let f = fcs[g].max(0) as usize;
                let nc = ncs[g].max(0) as usize;
                for c in 0..nc {
                    if f + c < n_pibins && off < mvs.len() {
                        matrix[r * n_pibins + f + c] = mvs[off];
                    }
                    off += 1;
                }
            }
        }
    }

    let (ebounds_channel, ebounds_emin, ebounds_emax) = if let Ok(ehdu) = fits.hdu_by_name("EBOUNDS") {
        let ebt = ehdu.as_bintable()?;
        (
            ebt.read_col_i16("CHANNEL")?,
            ebt.read_col_f32("E_MIN")?,
            ebt.read_col_f32("E_MAX")?,
        )
    } else {
        let eb = EboundsTable::default();
        (eb.channel, eb.e_min, eb.e_max)
    };
    let quad_id = mhdu.header().get_i64("QUADID").and_then(|v| u8::try_from(v).ok());

    Ok((
        meta,
        CztiResponse {
            quad_id,
            n_ebins,
            n_pibins,
            energ_lo,
            energ_hi,
            matrix,
            ebounds_channel,
            ebounds_emin,
            ebounds_emax,
        },
    ))
}

/// Read MKF header `RA_PNT`, `DEC_PNT`, plus attitude columns (`TIME`, `Roll_RA`, `Roll_DEC`, `ROLL_ROT`)
/// needed for `cztbindata`.
#[allow(clippy::type_complexity)]
pub fn read_mkf_bindata_attitude(
    mkf_path: impl AsRef<Path>,
) -> Result<(f64, f64, Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>)> {
    let fits = FitsFile::open(&mkf_path)?;
    let hdu = fits
        .hdu_by_name("MKF")
        .or_else(|_| fits.hdu_by_name("MKF_DATA"))
        .or_else(|_| {
            fits.iter_hdus()
                .find(|h| *h.kind() == HduKind::BinTable)
                .ok_or_else(|| astrosat_fits::FitsError::HduNotFound("MKF".into()))
        })?;
    let ra_pnt = hdu
        .header()
        .get_f64("RA_PNT")
        .or_else(|| fits.hdu(0).ok().and_then(|h| h.header().get_f64("RA_PNT")))
        .unwrap_or(0.0);
    let dec_pnt = hdu
        .header()
        .get_f64("DEC_PNT")
        .or_else(|| fits.hdu(0).ok().and_then(|h| h.header().get_f64("DEC_PNT")))
        .unwrap_or(0.0);

    let mkf = MkfTable::from_fits(mkf_path)?;
    let times = mkf.times().to_vec();
    let roll_ra = mkf
        .column_any(&["ROLL_RA", "Roll_RA"])
        .ok_or_else(|| CztiError::InvalidConfig("MKF missing Roll_RA".into()))?
        .to_vec();
    let roll_dec = mkf
        .column_any(&["ROLL_DEC", "Roll_DEC"])
        .ok_or_else(|| CztiError::InvalidConfig("MKF missing Roll_DEC".into()))?
        .to_vec();
    let roll_rot = mkf
        .column_any(&["ROLL_ROT", "Roll_ROT"])
        .ok_or_else(|| CztiError::InvalidConfig("MKF missing ROLL_ROT".into()))?
        .to_vec();

    Ok((ra_pnt, dec_pnt, times, roll_ra, roll_dec, roll_rot))
}
