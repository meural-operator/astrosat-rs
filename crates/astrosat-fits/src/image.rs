//! FITS `IMAGE` HDU reading and building (`BITPIX` 8, 16, 32, 64, -32, -64, `BZERO`/`BSCALE`).

use crate::card::Header;
use crate::error::{FitsError, Result};
use crate::io::Hdu;
use crate::ogip::{apply_wcs, WcsParams};

/// FITS `BITPIX` pixel format.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Bitpix {
    U8 = 8,
    I16 = 16,
    I32 = 32,
    I64 = 64,
    F32 = -32,
    F64 = -64,
}

impl Bitpix {
    pub fn from_i64(v: i64) -> Result<Self> {
        match v {
            8 => Ok(Self::U8),
            16 => Ok(Self::I16),
            32 => Ok(Self::I32),
            64 => Ok(Self::I64),
            -32 => Ok(Self::F32),
            -64 => Ok(Self::F64),
            other => Err(FitsError::UnsupportedBitpix(other)),
        }
    }

    pub fn byte_size(self) -> usize {
        ((self as i64).unsigned_abs() / 8) as usize
    }

    pub fn as_i64(self) -> i64 {
        self as i64
    }
}

/// Zero-copy reference to a FITS `IMAGE` HDU payload.
#[derive(Debug, Clone)]
pub struct ImageRef<'a> {
    bitpix: Bitpix,
    shape: Vec<usize>,
    bzero: f64,
    bscale: f64,
    data: &'a [u8],
}

impl<'a> ImageRef<'a> {
    pub(crate) fn from_header_and_data(header: &Header, data: &'a [u8]) -> Result<Self> {
        let bp_raw = header
            .get_i64("BITPIX")
            .ok_or_else(|| FitsError::InvalidHeader("missing BITPIX".into()))?;
        let bitpix = Bitpix::from_i64(bp_raw)?;
        let naxis_i64 = header.get_i64("NAXIS").unwrap_or(0);
        if !(0..=999).contains(&naxis_i64) {
            return Err(FitsError::InvalidHeader(format!(
                "invalid IMAGE NAXIS: {naxis_i64}"
            )));
        }
        let naxis = naxis_i64 as usize;
        let mut shape = Vec::with_capacity(naxis);
        for i in 1..=naxis {
            let dim = header
                .get_i64(&format!("NAXIS{i}"))
                .ok_or_else(|| FitsError::InvalidHeader(format!("missing NAXIS{i}")))?;
            let dim_u = usize::try_from(dim)
                .map_err(|_| FitsError::InvalidHeader(format!("invalid NAXIS{i}: {dim}")))?;
            shape.push(dim_u);
        }
        let num_pixels = if shape.is_empty() {
            0usize
        } else {
            shape.iter().try_fold(1usize, |acc, &d| acc.checked_mul(d)).ok_or_else(|| {
                FitsError::OutOfBounds(format!("image pixel count overflow for shape {shape:?}"))
            })?
        };
        let need_bytes = num_pixels.checked_mul(bitpix.byte_size()).ok_or_else(|| {
            FitsError::OutOfBounds(format!("image byte size overflow for shape {shape:?}"))
        })?;
        if data.len() < need_bytes {
            return Err(FitsError::OutOfBounds(format!(
                "image data truncated: need {need_bytes} bytes, got {}",
                data.len()
            )));
        }
        let bzero = header.get_f64("BZERO").unwrap_or(0.0);
        let bscale = header.get_f64("BSCALE").unwrap_or(1.0);
        Ok(Self {
            bitpix,
            shape,
            bzero,
            bscale,
            data,
        })
    }

    pub fn bitpix(&self) -> Bitpix {
        self.bitpix
    }

    /// Axis dimensions in FITS order `[NAXIS1, NAXIS2, ...]` where `NAXIS1` is fastest-varying.
    pub fn shape(&self) -> &[usize] {
        &self.shape
    }

    pub fn num_pixels(&self) -> usize {
        if self.shape.is_empty() {
            0
        } else {
            self.shape
                .iter()
                .try_fold(1usize, |acc, &d| acc.checked_mul(d))
                .unwrap_or(0)
        }
    }

    fn checked_data_bytes(&self) -> Result<usize> {
        let n = if self.shape.is_empty() {
            0usize
        } else {
            self.shape
                .iter()
                .try_fold(1usize, |acc, &d| acc.checked_mul(d))
                .ok_or_else(|| FitsError::OutOfBounds("image pixel count overflow".into()))?
        };
        let need = n
            .checked_mul(self.bitpix.byte_size())
            .ok_or_else(|| FitsError::OutOfBounds("image byte size overflow".into()))?;
        if self.data.len() < need {
            return Err(FitsError::OutOfBounds(format!(
                "image data truncated: need {need} bytes, got {}",
                self.data.len()
            )));
        }
        Ok(need)
    }

    pub fn bzero(&self) -> f64 {
        self.bzero
    }

    pub fn bscale(&self) -> f64 {
        self.bscale
    }

    pub fn raw_data(&self) -> &'a [u8] {
        self.data
    }

    /// Decode all pixels as `f64`, applying `BSCALE * raw + BZERO`.
    pub fn read_f64(&self) -> Result<Vec<f64>> {
        let need = self.checked_data_bytes()?;
        let raw = &self.data[..need];
        let bscale = self.bscale;
        let bzero = self.bzero;
        let scaled = bscale != 1.0 || bzero != 0.0;

        let out = match self.bitpix {
            Bitpix::U8 => raw
                .iter()
                .map(|&b| {
                    let v = b as f64;
                    if scaled { bscale * v + bzero } else { v }
                })
                .collect(),
            Bitpix::I16 => raw
                .as_chunks::<2>()
                .0
                .iter()
                .map(|&c| {
                    let v = i16::from_be_bytes(c) as f64;
                    if scaled { bscale * v + bzero } else { v }
                })
                .collect(),
            Bitpix::I32 => raw
                .as_chunks::<4>()
                .0
                .iter()
                .map(|&c| {
                    let v = i32::from_be_bytes(c) as f64;
                    if scaled { bscale * v + bzero } else { v }
                })
                .collect(),
            Bitpix::I64 => raw
                .as_chunks::<8>()
                .0
                .iter()
                .map(|&c| {
                    let v = i64::from_be_bytes(c) as f64;
                    if scaled { bscale * v + bzero } else { v }
                })
                .collect(),
            Bitpix::F32 => raw
                .as_chunks::<4>()
                .0
                .iter()
                .map(|&c| {
                    let v = f32::from_be_bytes(c) as f64;
                    if scaled { bscale * v + bzero } else { v }
                })
                .collect(),
            Bitpix::F64 => raw
                .as_chunks::<8>()
                .0
                .iter()
                .map(|&c| {
                    let v = f64::from_be_bytes(c);
                    if scaled { bscale * v + bzero } else { v }
                })
                .collect(),
        };
        Ok(out)
    }

    pub fn read_f32(&self) -> Result<Vec<f32>> {
        if self.bitpix == Bitpix::F32 && self.bscale == 1.0 && self.bzero == 0.0 {
            let n = self.num_pixels();
            let need = n * 4;
            if self.data.len() < need {
                return Err(FitsError::OutOfBounds("image data truncated".into()));
            }
            return Ok(self.data[..need]
                .as_chunks::<4>()
                .0
                .iter()
                .map(|&c| f32::from_be_bytes(c))
                .collect());
        }
        Ok(self.read_f64()?.into_iter().map(|v| v as f32).collect())
    }

    pub fn read_i32(&self) -> Result<Vec<i32>> {
        if self.bitpix == Bitpix::I32 && self.bscale == 1.0 && self.bzero == 0.0 {
            let n = self.num_pixels();
            let need = n * 4;
            if self.data.len() < need {
                return Err(FitsError::OutOfBounds("image data truncated".into()));
            }
            return Ok(self.data[..need]
                .as_chunks::<4>()
                .0
                .iter()
                .map(|&c| i32::from_be_bytes(c))
                .collect());
        }
        Ok(self
            .read_f64()?
            .into_iter()
            .map(|v| v.round() as i32)
            .collect())
    }

    pub fn read_i16(&self) -> Result<Vec<i16>> {
        if self.bitpix == Bitpix::I16 && self.bscale == 1.0 && self.bzero == 0.0 {
            let n = self.num_pixels();
            let need = n * 2;
            if self.data.len() < need {
                return Err(FitsError::OutOfBounds("image data truncated".into()));
            }
            return Ok(self.data[..need]
                .as_chunks::<2>()
                .0
                .iter()
                .map(|&c| i16::from_be_bytes(c))
                .collect());
        }
        Ok(self
            .read_f64()?
            .into_iter()
            .map(|v| v.round() as i16)
            .collect())
    }

    pub fn read_u8(&self) -> Result<Vec<u8>> {
        if self.bitpix == Bitpix::U8 && self.bscale == 1.0 && self.bzero == 0.0 {
            let n = self.num_pixels();
            if self.data.len() < n {
                return Err(FitsError::OutOfBounds("image data truncated".into()));
            }
            return Ok(self.data[..n].to_vec());
        }
        Ok(self
            .read_f64()?
            .into_iter()
            .map(|v| v.clamp(0.0, 255.0).round() as u8)
            .collect())
    }
}

/// Builder for Primary or Extension `IMAGE` HDUs.
pub struct ImageBuilder {
    extname: Option<String>,
    bitpix: Bitpix,
    shape: Vec<usize>,
    data: Vec<u8>,
    extra_cards: Header,
    bzero: Option<f64>,
    bscale: Option<f64>,
}

impl ImageBuilder {
    pub fn from_f32(extname: Option<&str>, shape: Vec<usize>, pixels: &[f32]) -> Result<Self> {
        check_shape(&shape, pixels.len())?;
        let mut data = Vec::with_capacity(pixels.len() * 4);
        for &v in pixels {
            data.extend_from_slice(&v.to_be_bytes());
        }
        Ok(Self {
            extname: extname.map(|s| s.to_string()),
            bitpix: Bitpix::F32,
            shape,
            data,
            extra_cards: Header::new(),
            bzero: None,
            bscale: None,
        })
    }

    pub fn from_f64(extname: Option<&str>, shape: Vec<usize>, pixels: &[f64]) -> Result<Self> {
        check_shape(&shape, pixels.len())?;
        let mut data = Vec::with_capacity(pixels.len() * 8);
        for &v in pixels {
            data.extend_from_slice(&v.to_be_bytes());
        }
        Ok(Self {
            extname: extname.map(|s| s.to_string()),
            bitpix: Bitpix::F64,
            shape,
            data,
            extra_cards: Header::new(),
            bzero: None,
            bscale: None,
        })
    }

    pub fn from_i32(extname: Option<&str>, shape: Vec<usize>, pixels: &[i32]) -> Result<Self> {
        check_shape(&shape, pixels.len())?;
        let mut data = Vec::with_capacity(pixels.len() * 4);
        for &v in pixels {
            data.extend_from_slice(&v.to_be_bytes());
        }
        Ok(Self {
            extname: extname.map(|s| s.to_string()),
            bitpix: Bitpix::I32,
            shape,
            data,
            extra_cards: Header::new(),
            bzero: None,
            bscale: None,
        })
    }

    pub fn from_i16(extname: Option<&str>, shape: Vec<usize>, pixels: &[i16]) -> Result<Self> {
        check_shape(&shape, pixels.len())?;
        let mut data = Vec::with_capacity(pixels.len() * 2);
        for &v in pixels {
            data.extend_from_slice(&v.to_be_bytes());
        }
        Ok(Self {
            extname: extname.map(|s| s.to_string()),
            bitpix: Bitpix::I16,
            shape,
            data,
            extra_cards: Header::new(),
            bzero: None,
            bscale: None,
        })
    }

    pub fn from_u8(extname: Option<&str>, shape: Vec<usize>, pixels: &[u8]) -> Result<Self> {
        check_shape(&shape, pixels.len())?;
        Ok(Self {
            extname: extname.map(|s| s.to_string()),
            bitpix: Bitpix::U8,
            shape,
            data: pixels.to_vec(),
            extra_cards: Header::new(),
            bzero: None,
            bscale: None,
        })
    }

    pub fn with_bzero_bscale(mut self, bzero: f64, bscale: f64) -> Self {
        self.bzero = Some(bzero);
        self.bscale = Some(bscale);
        self
    }

    pub fn with_wcs(mut self, wcs: &WcsParams) -> Self {
        apply_wcs(&mut self.extra_cards, wcs);
        self
    }

    pub fn header_str(mut self, key: &str, val: impl AsRef<str>, comment: &str) -> Self {
        self.extra_cards.set_str(key, val, comment);
        self
    }

    pub fn header_i64(mut self, key: &str, val: i64, comment: &str) -> Self {
        self.extra_cards.set_i64(key, val, comment);
        self
    }

    pub fn header_f64(mut self, key: &str, val: f64, comment: &str) -> Self {
        self.extra_cards.set_f64(key, val, comment);
        self
    }

    pub fn header_bool(mut self, key: &str, val: bool, comment: &str) -> Self {
        self.extra_cards.set_bool(key, val, comment);
        self
    }

    pub fn build_primary(self) -> Hdu {
        let mut hdr = Header::new();
        hdr.set_bool("SIMPLE", true, "file does conform to FITS standard");
        hdr.set_i64("BITPIX", self.bitpix.as_i64(), "number of bits per data pixel");
        hdr.set_i64("NAXIS", self.shape.len() as i64, "number of data axes");
        for (i, &dim) in self.shape.iter().enumerate() {
            hdr.set_i64(&format!("NAXIS{}", i + 1), dim as i64, "length of data axis");
        }
        hdr.set_bool("EXTEND", true, "FITS dataset may contain extensions");
        if let Some(bz) = self.bzero {
            hdr.set_f64("BZERO", bz, "zero point in scaling equation");
        }
        if let Some(bs) = self.bscale {
            hdr.set_f64("BSCALE", bs, "linear factor in scaling equation");
        }
        for c in self.extra_cards.cards() {
            hdr.set(c.clone());
        }
        Hdu::new(hdr, self.data)
    }

    pub fn build_extension(self) -> Hdu {
        let mut hdr = Header::new();
        hdr.set_str("XTENSION", "IMAGE", "IMAGE extension");
        hdr.set_i64("BITPIX", self.bitpix.as_i64(), "number of bits per data pixel");
        hdr.set_i64("NAXIS", self.shape.len() as i64, "number of data axes");
        for (i, &dim) in self.shape.iter().enumerate() {
            hdr.set_i64(&format!("NAXIS{}", i + 1), dim as i64, "length of data axis");
        }
        hdr.set_i64("PCOUNT", 0, "required keyword; must = 0");
        hdr.set_i64("GCOUNT", 1, "required keyword; must = 1");
        if let Some(name) = self.extname {
            hdr.set_str("EXTNAME", name, "extension name");
        }
        if let Some(bz) = self.bzero {
            hdr.set_f64("BZERO", bz, "zero point in scaling equation");
        }
        if let Some(bs) = self.bscale {
            hdr.set_f64("BSCALE", bs, "linear factor in scaling equation");
        }
        for c in self.extra_cards.cards() {
            hdr.set(c.clone());
        }
        Hdu::new(hdr, self.data)
    }
}

fn check_shape(shape: &[usize], len: usize) -> Result<()> {
    let expected: usize = if shape.is_empty() {
        0
    } else {
        shape
            .iter()
            .try_fold(1usize, |acc, &d| acc.checked_mul(d))
            .ok_or_else(|| {
                FitsError::FormatError(format!("image shape {shape:?} pixel count overflow"))
            })?
    };
    if expected != len {
        return Err(FitsError::FormatError(format!(
            "image pixel count {len} does not match shape {shape:?} (expected {expected})"
        )));
    }
    Ok(())
}
