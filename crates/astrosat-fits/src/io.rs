//! Memory-mapped FITS file reader and buffered streaming writer.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use memmap2::Mmap;

use crate::bintable::BinTableRef;
use crate::card::Header;
use crate::checksum::{compute_hdu_checksums, verify_hdu_checksum, FITS_BLOCK_SIZE};
use crate::error::{FitsError, Result};
use crate::image::ImageRef;

/// Classification of a FITS HDU.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum HduKind {
    PrimaryImage,
    Image,
    BinTable,
    AsciiTable,
    Other(String),
}

enum BackingStore {
    Mmap(Mmap),
    Owned(Vec<u8>),
}

impl BackingStore {
    #[inline]
    fn as_bytes(&self) -> &[u8] {
        match self {
            BackingStore::Mmap(m) => m.as_ref(),
            BackingStore::Owned(v) => v.as_slice(),
        }
    }
}

#[derive(Debug, Clone)]
struct HduIndex {
    kind: HduKind,
    header: Header,
    hdr_start: usize,
    hdr_len: usize,
    data_start: usize,
    data_len: usize,
    data_padded_len: usize,
}

/// Memory-mapped (zero-copy) FITS file reader.
pub struct FitsFile {
    store: BackingStore,
    hdus: Vec<HduIndex>,
}

impl FitsFile {
    /// Open and index a FITS file from disk using `memmap2`.
    pub fn open(path: impl AsRef<Path>) -> Result<Self> {
        let file = File::open(path.as_ref())?;
        let meta = file.metadata()?;
        if meta.len() == 0 {
            return Err(FitsError::InvalidBlockAlignment(
                "FITS file is empty (0 bytes)".into(),
            ));
        }
        // SAFETY: We treat the file as read-only; if mmap fails (e.g. special filesystem), fall back to reading bytes.
        let store = match unsafe { Mmap::map(&file) } {
            Ok(mmap) => BackingStore::Mmap(mmap),
            Err(_) => {
                let bytes = std::fs::read(path.as_ref())?;
                BackingStore::Owned(bytes)
            }
        };
        Self::from_store(store)
    }

    /// Parse a FITS file from an in-memory byte vector.
    pub fn from_bytes(bytes: Vec<u8>) -> Result<Self> {
        Self::from_store(BackingStore::Owned(bytes))
    }

    fn from_store(store: BackingStore) -> Result<Self> {
        let buf = store.as_bytes();
        let mut hdus = Vec::new();
        let mut offset = 0usize;

        while offset < buf.len() {
            // Stop gracefully if only trailing NUL/space padding remains after a valid HDU.
            if !hdus.is_empty()
                && buf[offset..].iter().all(|&b| b == 0 || b == b' ')
            {
                break;
            }
            if buf.len() - offset < FITS_BLOCK_SIZE {
                return Err(FitsError::InvalidBlockAlignment(format!(
                    "trailing {} bytes at offset {offset} (< 2880)",
                    buf.len() - offset
                )));
            }

            let (header, hdr_len) = Header::parse(&buf[offset..])?;
            let is_primary = hdus.is_empty();
            let kind = if is_primary {
                if header.get_bool("SIMPLE") != Some(true) {
                    return Err(FitsError::InvalidHeader(
                        "Primary HDU missing SIMPLE = T".into(),
                    ));
                }
                HduKind::PrimaryImage
            } else {
                match header.get_str("XTENSION").map(|s| s.trim()) {
                    Some("BINTABLE") => HduKind::BinTable,
                    Some("IMAGE") => HduKind::Image,
                    Some("TABLE") => HduKind::AsciiTable,
                    Some(other) => HduKind::Other(other.to_string()),
                    None => {
                        return Err(FitsError::InvalidHeader(format!(
                            "extension HDU {} missing XTENSION keyword",
                            hdus.len()
                        )));
                    }
                }
            };

            let data_len = compute_data_len(&header)?;
            let data_padded_len = if data_len == 0 {
                0
            } else {
                data_len
                    .checked_next_multiple_of(FITS_BLOCK_SIZE)
                    .ok_or_else(|| {
                        FitsError::OutOfBounds(format!(
                            "HDU {} padded data size overflow",
                            hdus.len()
                        ))
                    })?
            };

            let data_start = offset.checked_add(hdr_len).ok_or_else(|| {
                FitsError::OutOfBounds(format!("HDU {} offset overflow", hdus.len()))
            })?;
            let data_end = data_start.checked_add(data_len).ok_or_else(|| {
                FitsError::OutOfBounds(format!("HDU {} data length overflow", hdus.len()))
            })?;
            if data_end > buf.len() {
                return Err(FitsError::OutOfBounds(format!(
                    "HDU {} data extends past EOF (need {data_end}, file size {})",
                    hdus.len(),
                    buf.len()
                )));
            }

            let actual_padded = data_start
                .saturating_add(data_padded_len)
                .min(buf.len())
                - data_start;

            hdus.push(HduIndex {
                kind,
                header,
                hdr_start: offset,
                hdr_len,
                data_start,
                data_len,
                data_padded_len: actual_padded,
            });

            offset = data_start + actual_padded;
        }

        if hdus.is_empty() {
            return Err(FitsError::InvalidHeader("no HDUs found in FITS file".into()));
        }

        Ok(Self { store, hdus })
    }

    pub fn len(&self) -> usize {
        self.hdus.len()
    }

    pub fn is_empty(&self) -> bool {
        self.hdus.is_empty()
    }

    pub fn hdu(&self, idx: usize) -> Result<HduRef<'_>> {
        let meta = self
            .hdus
            .get(idx)
            .ok_or_else(|| FitsError::HduNotFound(format!("index {idx}")))?;
        Ok(HduRef {
            index: idx,
            meta,
            buf: self.store.as_bytes(),
        })
    }

    pub fn hdu_by_name(&self, name: &str) -> Result<HduRef<'_>> {
        let needle = name.trim();
        for (idx, h) in self.hdus.iter().enumerate() {
            if let Some(ext) = h.header.extname() {
                if ext.trim().eq_ignore_ascii_case(needle) {
                    return self.hdu(idx);
                }
            }
        }
        Err(FitsError::HduNotFound(name.to_string()))
    }

    pub fn hdu_by_name_ver(&self, name: &str, extver: i64) -> Result<HduRef<'_>> {
        let needle = name.trim();
        for (idx, h) in self.hdus.iter().enumerate() {
            if let Some(ext) = h.header.extname() {
                if ext.trim().eq_ignore_ascii_case(needle)
                    && h.header.extver().unwrap_or(1) == extver
                {
                    return self.hdu(idx);
                }
            }
        }
        Err(FitsError::HduNotFound(format!("{name}[{extver}]")))
    }

    pub fn iter_hdus(&self) -> impl Iterator<Item = HduRef<'_>> {
        let buf = self.store.as_bytes();
        self.hdus.iter().enumerate().map(move |(index, meta)| HduRef {
            index,
            meta,
            buf,
        })
    }

    pub fn hdu_names(&self) -> Vec<String> {
        self.hdus
            .iter()
            .enumerate()
            .map(|(i, h)| {
                h.header
                    .extname()
                    .map(|s| s.to_string())
                    .unwrap_or_else(|| {
                        if i == 0 {
                            "PRIMARY".to_string()
                        } else {
                            format!("HDU{i}")
                        }
                    })
            })
            .collect()
    }
}

fn compute_data_len(header: &Header) -> Result<usize> {
    let bitpix = header
        .get_i64("BITPIX")
        .ok_or_else(|| FitsError::InvalidHeader("missing BITPIX".into()))?;
    let elem_bytes = crate::image::Bitpix::from_i64(bitpix)?.byte_size();
    let naxis = header
        .get_i64("NAXIS")
        .ok_or_else(|| FitsError::InvalidHeader("missing NAXIS".into()))?;
    if !(0..=999).contains(&naxis) {
        return Err(FitsError::InvalidHeader(format!("invalid NAXIS: {naxis}")));
    }
    if naxis == 0 {
        return Ok(0);
    }

    let mut prod: usize = 1;
    for i in 1..=naxis {
        let dim = header
            .get_i64(&format!("NAXIS{i}"))
            .ok_or_else(|| FitsError::InvalidHeader(format!("missing NAXIS{i}")))?;
        let dim_u = usize::try_from(dim)
            .map_err(|_| FitsError::InvalidHeader(format!("invalid NAXIS{i}: {dim}")))?;
        prod = prod
            .checked_mul(dim_u)
            .ok_or_else(|| FitsError::OutOfBounds("HDU NAXIS product overflow".into()))?;
    }

    let pcount_i64 = header.get_i64("PCOUNT").unwrap_or(0);
    let pcount = usize::try_from(pcount_i64)
        .map_err(|_| FitsError::InvalidHeader(format!("invalid PCOUNT: {pcount_i64}")))?;
    let gcount_i64 = header.get_i64("GCOUNT").unwrap_or(1);
    let gcount = usize::try_from(gcount_i64)
        .map_err(|_| FitsError::InvalidHeader(format!("invalid GCOUNT: {gcount_i64}")))?;

    let bytes = prod
        .checked_add(pcount)
        .and_then(|v| v.checked_mul(gcount))
        .and_then(|v| v.checked_mul(elem_bytes))
        .ok_or_else(|| FitsError::OutOfBounds("HDU data byte size overflow".into()))?;
    Ok(bytes)
}

/// Zero-copy reference to a single HDU inside an open [`FitsFile`].
pub struct HduRef<'a> {
    index: usize,
    meta: &'a HduIndex,
    buf: &'a [u8],
}

impl<'a> HduRef<'a> {
    pub fn index(&self) -> usize {
        self.index
    }

    pub fn kind(&self) -> &HduKind {
        &self.meta.kind
    }

    pub fn header(&self) -> &'a Header {
        &self.meta.header
    }

    pub fn extname(&self) -> Option<&'a str> {
        self.meta.header.extname()
    }

    pub fn extver(&self) -> Option<i64> {
        self.meta.header.extver()
    }

    pub fn raw_header(&self) -> &'a [u8] {
        &self.buf[self.meta.hdr_start..self.meta.hdr_start + self.meta.hdr_len]
    }

    /// Exact logical data payload bytes (excluding trailing FITS block padding).
    pub fn raw_data(&self) -> &'a [u8] {
        &self.buf[self.meta.data_start..self.meta.data_start + self.meta.data_len]
    }

    /// Block-padded data payload bytes.
    pub fn raw_data_padded(&self) -> &'a [u8] {
        &self.buf[self.meta.data_start..self.meta.data_start + self.meta.data_padded_len]
    }

    pub fn as_bintable(&self) -> Result<BinTableRef<'a>> {
        if self.meta.kind != HduKind::BinTable {
            return Err(FitsError::TypeMismatch(format!(
                "HDU {} is {:?}, not BINTABLE",
                self.index, self.meta.kind
            )));
        }
        BinTableRef::from_header_and_data(&self.meta.header, self.raw_data())
    }

    pub fn as_image(&self) -> Result<ImageRef<'a>> {
        if !matches!(self.meta.kind, HduKind::PrimaryImage | HduKind::Image) {
            return Err(FitsError::TypeMismatch(format!(
                "HDU {} is {:?}, not IMAGE",
                self.index, self.meta.kind
            )));
        }
        ImageRef::from_header_and_data(&self.meta.header, self.raw_data())
    }

    /// Verify `DATASUM` and `CHECKSUM` on this HDU.
    pub fn verify_checksum(&self) -> (bool, bool) {
        verify_hdu_checksum(&self.meta.header, self.raw_header(), self.raw_data_padded())
    }
}

/// Owned FITS HDU ready for modification or serialization.
#[derive(Debug, Clone)]
pub struct Hdu {
    pub header: Header,
    pub data: Vec<u8>,
}

impl Hdu {
    pub fn new(header: Header, data: Vec<u8>) -> Self {
        Self { header, data }
    }

    /// Create a standard empty Primary HDU (`SIMPLE=T, BITPIX=8, NAXIS=0, EXTEND=T`).
    pub fn empty_primary() -> Self {
        let mut hdr = Header::new();
        hdr.set_bool("SIMPLE", true, "file does conform to FITS standard");
        hdr.set_i64("BITPIX", 8, "number of bits per data pixel");
        hdr.set_i64("NAXIS", 0, "number of data axes");
        hdr.set_bool("EXTEND", true, "FITS dataset may contain extensions");
        Self {
            header: hdr,
            data: Vec::new(),
        }
    }

    /// Clone an existing [`HduRef`] into an owned mutable [`Hdu`].
    pub fn from_ref(hdu: &HduRef<'_>) -> Self {
        Self {
            header: hdu.header().clone(),
            data: hdu.raw_data().to_vec(),
        }
    }

    /// Compute and update `DATASUM` and `CHECKSUM` keywords on this HDU.
    pub fn update_checksums(&mut self) -> (u32, String) {
        compute_hdu_checksums(&mut self.header, &self.data)
    }
}

/// Streaming / buffered FITS 4.0 writer.
pub struct FitsWriter<W: Write> {
    writer: W,
    hdus_written: usize,
}

impl FitsWriter<BufWriter<File>> {
    pub fn create(path: impl AsRef<Path>) -> Result<Self> {
        let file = File::create(path.as_ref())?;
        Ok(Self::new(BufWriter::with_capacity(256 * 1024, file)))
    }
}

impl<W: Write> FitsWriter<W> {
    pub fn new(writer: W) -> Self {
        Self {
            writer,
            hdus_written: 0,
        }
    }

    pub fn write_hdu(&mut self, hdu: &Hdu) -> Result<()> {
        let hdr_bytes = hdu.header.to_fits_bytes();
        self.writer.write_all(&hdr_bytes)?;

        if !hdu.data.is_empty() {
            self.writer.write_all(&hdu.data)?;
            let rem = hdu.data.len() % FITS_BLOCK_SIZE;
            if rem != 0 {
                let pad_len = FITS_BLOCK_SIZE - rem;
                let pad_byte = if hdu
                    .header
                    .get_str("XTENSION")
                    .map(|s| s.trim().eq_ignore_ascii_case("TABLE"))
                    .unwrap_or(false)
                {
                    b' '
                } else {
                    0u8
                };
                let pad = vec![pad_byte; pad_len];
                self.writer.write_all(&pad)?;
            }
        }

        self.hdus_written += 1;
        Ok(())
    }

    pub fn write_hdu_with_checksum(&mut self, hdu: &mut Hdu) -> Result<()> {
        hdu.update_checksums();
        self.write_hdu(hdu)
    }

    pub fn finish(mut self) -> Result<W> {
        self.writer.flush()?;
        Ok(self.writer)
    }
}

/// Convenience helper to write a sequence of HDUs to disk.
pub fn write_fits(
    path: impl AsRef<Path>,
    hdus: &mut [Hdu],
    with_checksums: bool,
) -> Result<()> {
    let mut writer = FitsWriter::create(path)?;
    for hdu in hdus.iter_mut() {
        if with_checksums {
            writer.write_hdu_with_checksum(hdu)?;
        } else {
            writer.write_hdu(hdu)?;
        }
    }
    writer.finish()?;
    Ok(())
}
