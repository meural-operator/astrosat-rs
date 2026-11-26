//! FITS `BINTABLE` HDU zero-copy reader and column-oriented builder.

use crate::card::{Card, Header};
use crate::error::{FitsError, Result};
use crate::io::Hdu;

/// FITS BINTABLE primitive element type (from `TFORMn`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ColumnType {
    Logical,    // L
    Bit,        // X
    Byte,       // B
    I16,        // I
    I32,        // J
    I64,        // K
    F32,        // E
    F64,        // D
    Ascii,      // A
    ComplexF32, // C
    ComplexF64, // M
}

impl ColumnType {
    pub fn from_char(c: char) -> Result<Self> {
        match c.to_ascii_uppercase() {
            'L' => Ok(Self::Logical),
            'X' => Ok(Self::Bit),
            'B' => Ok(Self::Byte),
            'I' => Ok(Self::I16),
            'J' => Ok(Self::I32),
            'K' => Ok(Self::I64),
            'E' => Ok(Self::F32),
            'D' => Ok(Self::F64),
            'A' => Ok(Self::Ascii),
            'C' => Ok(Self::ComplexF32),
            'M' => Ok(Self::ComplexF64),
            _ => Err(FitsError::InvalidTform(c.to_string())),
        }
    }

    pub fn code_char(self) -> char {
        match self {
            Self::Logical => 'L',
            Self::Bit => 'X',
            Self::Byte => 'B',
            Self::I16 => 'I',
            Self::I32 => 'J',
            Self::I64 => 'K',
            Self::F32 => 'E',
            Self::F64 => 'D',
            Self::Ascii => 'A',
            Self::ComplexF32 => 'C',
            Self::ComplexF64 => 'M',
        }
    }

    pub fn elem_size(self) -> usize {
        match self {
            Self::Logical | Self::Byte | Self::Ascii | Self::Bit => 1,
            Self::I16 => 2,
            Self::I32 | Self::F32 => 4,
            Self::I64 | Self::F64 | Self::ComplexF32 => 8,
            Self::ComplexF64 => 16,
        }
    }
}

/// Variable-length descriptor kind for `P` (32-bit) or `Q` (64-bit) heap arrays.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VarDescriptor {
    /// Fixed-length column in main table.
    Fixed,
    /// 32-bit (count, offset) descriptor (`PT`).
    P32,
    /// 64-bit (count, offset) descriptor (`QT`).
    Q64,
}

/// Parsed schema for one column of a FITS `BINTABLE`.
#[derive(Debug, Clone)]
pub struct ColumnSpec {
    pub index: usize,
    pub name: String,
    pub col_type: ColumnType,
    pub repeat: usize,
    pub var_desc: VarDescriptor,
    pub unit: Option<String>,
    pub tzero: Option<f64>,
    pub tscal: Option<f64>,
    pub tdisp: Option<String>,
    pub tform: String,
    /// Byte offset from the start of a row in the main table.
    pub offset: usize,
    /// Number of bytes occupied per row in the main table.
    pub byte_width: usize,
}

/// Parse a FITS BINTABLE `TFORMn` string such as `"1D"`, `"D"`, `"1024I"`, `"2242B"`, `"16A"`, `"PE(128)"`, `"1QJ"`.
pub fn parse_tform(tform: &str) -> Result<(usize, VarDescriptor, ColumnType, usize)> {
    let s = tform.trim();
    if s.is_empty() {
        return Err(FitsError::InvalidTform(tform.to_string()));
    }
    let bytes = s.as_bytes();
    let mut i = 0;
    while i < bytes.len() && bytes[i].is_ascii_digit() {
        i += 1;
    }
    let repeat: usize = if i == 0 {
        1
    } else {
        s[..i]
            .parse::<usize>()
            .map_err(|_| FitsError::InvalidTform(tform.to_string()))?
    };

    if i >= bytes.len() {
        return Err(FitsError::InvalidTform(tform.to_string()));
    }

    let c = (bytes[i] as char).to_ascii_uppercase();
    if c == 'P' || c == 'Q' {
        let var_desc = if c == 'P' {
            VarDescriptor::P32
        } else {
            VarDescriptor::Q64
        };
        if i + 1 >= bytes.len() {
            return Err(FitsError::InvalidTform(tform.to_string()));
        }
        let col_type = ColumnType::from_char(bytes[i + 1] as char)?;
        let width = match var_desc {
            VarDescriptor::P32 => 8,
            VarDescriptor::Q64 => 16,
            VarDescriptor::Fixed => unreachable!(),
        };
        return Ok((repeat, var_desc, col_type, width));
    }

    let col_type = ColumnType::from_char(c)?;
    let byte_width = if col_type == ColumnType::Bit {
        repeat.div_ceil(8)
    } else {
        repeat
            .checked_mul(col_type.elem_size())
            .ok_or_else(|| FitsError::InvalidTform(tform.to_string()))?
    };
    Ok((repeat, VarDescriptor::Fixed, col_type, byte_width))
}

/// Zero-copy view into a FITS `BINTABLE` HDU.
#[derive(Debug, Clone)]
pub struct BinTableRef<'a> {
    nrows: usize,
    row_len: usize,
    theap: usize,
    pcount: usize,
    columns: Vec<ColumnSpec>,
    data: &'a [u8],
}

impl<'a> BinTableRef<'a> {
    pub(crate) fn from_header_and_data(header: &Header, data: &'a [u8]) -> Result<Self> {
        let row_len_i64 = header
            .get_i64("NAXIS1")
            .ok_or_else(|| FitsError::InvalidHeader("BINTABLE missing NAXIS1".into()))?;
        let row_len = usize::try_from(row_len_i64)
            .map_err(|_| FitsError::InvalidHeader(format!("invalid NAXIS1: {row_len_i64}")))?;

        let nrows_i64 = header
            .get_i64("NAXIS2")
            .ok_or_else(|| FitsError::InvalidHeader("BINTABLE missing NAXIS2".into()))?;
        let nrows = usize::try_from(nrows_i64)
            .map_err(|_| FitsError::InvalidHeader(format!("invalid NAXIS2: {nrows_i64}")))?;

        let pcount_i64 = header.get_i64("PCOUNT").unwrap_or(0);
        let pcount = usize::try_from(pcount_i64)
            .map_err(|_| FitsError::InvalidHeader(format!("invalid PCOUNT: {pcount_i64}")))?;

        let tfields_i64 = header
            .get_i64("TFIELDS")
            .ok_or_else(|| FitsError::InvalidHeader("BINTABLE missing TFIELDS".into()))?;
        let tfields = usize::try_from(tfields_i64)
            .map_err(|_| FitsError::InvalidHeader(format!("invalid TFIELDS: {tfields_i64}")))?;
        if tfields > header.cards().len() {
            return Err(FitsError::InvalidHeader(format!(
                "BINTABLE TFIELDS ({tfields}) exceeds header card count ({})",
                header.cards().len()
            )));
        }

        let main_table_bytes = row_len.checked_mul(nrows).ok_or_else(|| {
            FitsError::InvalidHeader("BINTABLE NAXIS1 * NAXIS2 overflow".into())
        })?;
        let theap = match header.get_i64("THEAP") {
            Some(v) => usize::try_from(v)
                .map_err(|_| FitsError::InvalidHeader(format!("invalid THEAP: {v}")))?,
            None => main_table_bytes,
        };

        let heap_end = if pcount > 0 {
            theap.checked_add(pcount).ok_or_else(|| {
                FitsError::OutOfBounds("BINTABLE THEAP + PCOUNT overflow".into())
            })?
        } else {
            0
        };
        let min_bytes = main_table_bytes.max(heap_end);
        if data.len() < min_bytes {
            return Err(FitsError::OutOfBounds(format!(
                "BINTABLE data truncated: need {min_bytes} bytes, got {}",
                data.len()
            )));
        }

        let mut columns = Vec::with_capacity(tfields);
        let mut offset = 0usize;
        for idx in 1..=tfields {
            let tform = header
                .get_str(&format!("TFORM{idx}"))
                .ok_or_else(|| FitsError::InvalidHeader(format!("missing TFORM{idx}")))?
                .trim()
                .to_string();
            let name = header
                .get_str(&format!("TTYPE{idx}"))
                .map(|s| s.trim().to_string())
                .unwrap_or_else(|| format!("COL{idx}"));
            let unit = header
                .get_str(&format!("TUNIT{idx}"))
                .map(|s| s.trim().to_string());
            let tzero = header.get_f64(&format!("TZERO{idx}"));
            let tscal = header.get_f64(&format!("TSCAL{idx}"));
            let tdisp = header
                .get_str(&format!("TDISP{idx}"))
                .map(|s| s.trim().to_string());

            let (repeat, var_desc, col_type, byte_width) = parse_tform(&tform)?;
            columns.push(ColumnSpec {
                index: idx - 1,
                name,
                col_type,
                repeat,
                var_desc,
                unit,
                tzero,
                tscal,
                tdisp,
                tform,
                offset,
                byte_width,
            });
            offset = offset.checked_add(byte_width).ok_or_else(|| {
                FitsError::InvalidHeader("BINTABLE column byte offset overflow".into())
            })?;
        }

        if offset != row_len && tfields > 0 {
            return Err(FitsError::InvalidHeader(format!(
                "BINTABLE sum of column widths ({offset}) != NAXIS1 ({row_len})"
            )));
        }

        Ok(Self {
            nrows,
            row_len,
            theap,
            pcount,
            columns,
            data,
        })
    }

    pub fn nrows(&self) -> usize {
        self.nrows
    }

    pub fn row_len(&self) -> usize {
        self.row_len
    }

    pub fn columns(&self) -> &[ColumnSpec] {
        &self.columns
    }

    pub fn column_index(&self, name: &str) -> Option<usize> {
        let needle = name.trim();
        self.columns
            .iter()
            .position(|c| c.name.eq_ignore_ascii_case(needle))
    }

    pub fn column_spec(&self, name: &str) -> Option<&ColumnSpec> {
        self.column_index(name).map(|i| &self.columns[i])
    }

    fn require_col(&self, name: &str) -> Result<&ColumnSpec> {
        self.column_spec(name)
            .ok_or_else(|| FitsError::ColumnNotFound(name.to_string()))
    }

    pub fn raw_data(&self) -> &'a [u8] {
        self.data
    }

    /// Zero-copy slice for row `row` (`0..nrows`).
    #[inline]
    pub fn row_bytes(&self, row: usize) -> &'a [u8] {
        let start = row * self.row_len;
        &self.data[start..start + self.row_len]
    }

    /// Zero-copy slice for a single cell `(row, col_idx)`.
    #[inline]
    pub fn cell_bytes(&self, row: usize, col_idx: usize) -> &'a [u8] {
        let spec = &self.columns[col_idx];
        let start = row * self.row_len + spec.offset;
        &self.data[start..start + spec.byte_width]
    }

    /// Zero-copy slice for a single cell by column name.
    pub fn cell_bytes_by_name(&self, row: usize, name: &str) -> Result<&'a [u8]> {
        if row >= self.nrows {
            return Err(FitsError::OutOfBounds(format!(
                "row {row} >= nrows {}",
                self.nrows
            )));
        }
        let idx = self
            .column_index(name)
            .ok_or_else(|| FitsError::ColumnNotFound(name.to_string()))?;
        Ok(self.cell_bytes(row, idx))
    }

    /// Iterate over zero-copy row slices.
    pub fn iter_rows(&self) -> impl Iterator<Item = &'a [u8]> + '_ {
        (0..self.nrows).map(move |r| self.row_bytes(r))
    }

    /// Read a scalar numeric column as `Vec<f64>`, applying `TSCAL`/`TZERO` if present.
    pub fn read_col_f64(&self, name: &str) -> Result<Vec<f64>> {
        let col = self.require_col(name)?;
        if col.var_desc != VarDescriptor::Fixed || col.repeat != 1 {
            return Err(FitsError::TypeMismatch(format!(
                "column '{}' is not a scalar fixed column (TFORM={})",
                col.name, col.tform
            )));
        }
        let tscal = col.tscal.unwrap_or(1.0);
        let tzero = col.tzero.unwrap_or(0.0);
        let scaled = tscal != 1.0 || tzero != 0.0;
        let off = col.offset;
        let rl = self.row_len;

        let mut out = Vec::with_capacity(self.nrows);
        match col.col_type {
            ColumnType::F64 => {
                for r in 0..self.nrows {
                    let s = r * rl + off;
                    let v = f64::from_be_bytes(self.data[s..s + 8].try_into().unwrap());
                    out.push(if scaled { tscal * v + tzero } else { v });
                }
            }
            ColumnType::F32 => {
                for r in 0..self.nrows {
                    let s = r * rl + off;
                    let v = f32::from_be_bytes(self.data[s..s + 4].try_into().unwrap()) as f64;
                    out.push(if scaled { tscal * v + tzero } else { v });
                }
            }
            ColumnType::I64 => {
                for r in 0..self.nrows {
                    let s = r * rl + off;
                    let v = i64::from_be_bytes(self.data[s..s + 8].try_into().unwrap()) as f64;
                    out.push(if scaled { tscal * v + tzero } else { v });
                }
            }
            ColumnType::I32 => {
                for r in 0..self.nrows {
                    let s = r * rl + off;
                    let v = i32::from_be_bytes(self.data[s..s + 4].try_into().unwrap()) as f64;
                    out.push(if scaled { tscal * v + tzero } else { v });
                }
            }
            ColumnType::I16 => {
                for r in 0..self.nrows {
                    let s = r * rl + off;
                    let v = i16::from_be_bytes(self.data[s..s + 2].try_into().unwrap()) as f64;
                    out.push(if scaled { tscal * v + tzero } else { v });
                }
            }
            ColumnType::Byte => {
                for r in 0..self.nrows {
                    let v = self.data[r * rl + off] as f64;
                    out.push(if scaled { tscal * v + tzero } else { v });
                }
            }
            other => {
                return Err(FitsError::TypeMismatch(format!(
                    "cannot read column '{}' of type {:?} as f64",
                    col.name, other
                )));
            }
        }
        Ok(out)
    }

    pub fn read_col_f32(&self, name: &str) -> Result<Vec<f32>> {
        let col = self.require_col(name)?;
        if col.var_desc == VarDescriptor::Fixed
            && col.repeat == 1
            && col.col_type == ColumnType::F32
            && col.tscal.unwrap_or(1.0) == 1.0
            && col.tzero.unwrap_or(0.0) == 0.0
        {
            let off = col.offset;
            let rl = self.row_len;
            let mut out = Vec::with_capacity(self.nrows);
            for r in 0..self.nrows {
                let s = r * rl + off;
                out.push(f32::from_be_bytes(self.data[s..s + 4].try_into().unwrap()));
            }
            return Ok(out);
        }
        Ok(self.read_col_f64(name)?.into_iter().map(|v| v as f32).collect())
    }

    /// Read a scalar integer column as `Vec<i64>` (preserving exact integer precision for `TZERO` unsigned offsets).
    pub fn read_col_i64(&self, name: &str) -> Result<Vec<i64>> {
        let col = self.require_col(name)?;
        if col.var_desc != VarDescriptor::Fixed || col.repeat != 1 {
            return Err(FitsError::TypeMismatch(format!(
                "column '{}' is not a scalar fixed column (TFORM={})",
                col.name, col.tform
            )));
        }
        let tscal = col.tscal.unwrap_or(1.0);
        let tzero = col.tzero.unwrap_or(0.0);
        let int_shift = if tscal == 1.0 && tzero.is_finite() && tzero.fract() == 0.0 {
            if tzero == 9_223_372_036_854_775_808.0 {
                // FITS standard unsigned 64-bit TZERO = 2^63 wraps to i64::MIN in two's complement.
                Some(i64::MIN)
            } else if tzero >= (i64::MIN as f64) && tzero < 9_223_372_036_854_775_808.0 {
                Some(tzero as i64)
            } else {
                None
            }
        } else {
            None
        };
        let off = col.offset;
        let rl = self.row_len;
        let mut out = Vec::with_capacity(self.nrows);

        match col.col_type {
            ColumnType::I64 => {
                for r in 0..self.nrows {
                    let s = r * rl + off;
                    let raw = i64::from_be_bytes(self.data[s..s + 8].try_into().unwrap());
                    out.push(match int_shift {
                        Some(z) => raw.wrapping_add(z),
                        None => (tscal * (raw as f64) + tzero).round() as i64,
                    });
                }
            }
            ColumnType::I32 => {
                for r in 0..self.nrows {
                    let s = r * rl + off;
                    let raw = i32::from_be_bytes(self.data[s..s + 4].try_into().unwrap()) as i64;
                    out.push(match int_shift {
                        Some(z) => raw.wrapping_add(z),
                        None => (tscal * (raw as f64) + tzero).round() as i64,
                    });
                }
            }
            ColumnType::I16 => {
                for r in 0..self.nrows {
                    let s = r * rl + off;
                    let raw = i16::from_be_bytes(self.data[s..s + 2].try_into().unwrap()) as i64;
                    out.push(match int_shift {
                        Some(z) => raw.wrapping_add(z),
                        None => (tscal * (raw as f64) + tzero).round() as i64,
                    });
                }
            }
            ColumnType::Byte => {
                for r in 0..self.nrows {
                    let raw = self.data[r * rl + off] as i64;
                    out.push(match int_shift {
                        Some(z) => raw.wrapping_add(z),
                        None => (tscal * (raw as f64) + tzero).round() as i64,
                    });
                }
            }
            ColumnType::F32 | ColumnType::F64 => {
                return Ok(self
                    .read_col_f64(name)?
                    .into_iter()
                    .map(|v| v.round() as i64)
                    .collect());
            }
            other => {
                return Err(FitsError::TypeMismatch(format!(
                    "cannot read column '{}' of type {:?} as i64",
                    col.name, other
                )));
            }
        }
        Ok(out)
    }

    pub fn read_col_i32(&self, name: &str) -> Result<Vec<i32>> {
        Ok(self.read_col_i64(name)?.into_iter().map(|v| v as i32).collect())
    }

    pub fn read_col_i16(&self, name: &str) -> Result<Vec<i16>> {
        Ok(self.read_col_i64(name)?.into_iter().map(|v| v as i16).collect())
    }

    pub fn read_col_u16(&self, name: &str) -> Result<Vec<u16>> {
        Ok(self.read_col_i64(name)?.into_iter().map(|v| v as u16).collect())
    }

    pub fn read_col_u32(&self, name: &str) -> Result<Vec<u32>> {
        Ok(self.read_col_i64(name)?.into_iter().map(|v| v as u32).collect())
    }

    pub fn read_col_u8(&self, name: &str) -> Result<Vec<u8>> {
        let col = self.require_col(name)?;
        if col.var_desc == VarDescriptor::Fixed
            && col.repeat == 1
            && col.col_type == ColumnType::Byte
            && col.tscal.unwrap_or(1.0) == 1.0
            && col.tzero.unwrap_or(0.0) == 0.0
        {
            let off = col.offset;
            let rl = self.row_len;
            let mut out = Vec::with_capacity(self.nrows);
            for r in 0..self.nrows {
                out.push(self.data[r * rl + off]);
            }
            return Ok(out);
        }
        Ok(self.read_col_i64(name)?.into_iter().map(|v| v as u8).collect())
    }

    pub fn read_col_bool(&self, name: &str) -> Result<Vec<bool>> {
        let col = self.require_col(name)?;
        if col.var_desc != VarDescriptor::Fixed || col.repeat != 1 {
            return Err(FitsError::TypeMismatch(format!(
                "column '{}' is not a scalar logical column",
                col.name
            )));
        }
        let off = col.offset;
        let rl = self.row_len;
        let mut out = Vec::with_capacity(self.nrows);
        match col.col_type {
            ColumnType::Logical => {
                for r in 0..self.nrows {
                    out.push(self.data[r * rl + off] == b'T');
                }
            }
            ColumnType::Byte => {
                for r in 0..self.nrows {
                    out.push(self.data[r * rl + off] != 0);
                }
            }
            other => {
                return Err(FitsError::TypeMismatch(format!(
                    "cannot read column '{}' of type {:?} as bool",
                    col.name, other
                )));
            }
        }
        Ok(out)
    }

    pub fn read_col_strings(&self, name: &str) -> Result<Vec<String>> {
        let col = self.require_col(name)?;
        if col.var_desc != VarDescriptor::Fixed || col.col_type != ColumnType::Ascii {
            return Err(FitsError::TypeMismatch(format!(
                "column '{}' is not a fixed ASCII string column (TFORM={})",
                col.name, col.tform
            )));
        }
        let off = col.offset;
        let w = col.byte_width;
        let rl = self.row_len;
        let mut out = Vec::with_capacity(self.nrows);
        for r in 0..self.nrows {
            let slice = &self.data[r * rl + off..r * rl + off + w];
            // Strip trailing NULs and spaces per FITS convention.
            let end = slice
                .iter()
                .position(|&b| b == 0)
                .unwrap_or(slice.len());
            let s = String::from_utf8_lossy(&slice[..end])
                .trim_end_matches(' ')
                .to_string();
            out.push(s);
        }
        Ok(out)
    }

    /// Read a fixed-width byte vector column (`nB`, e.g., `2242B`) as `(repeat, flat_data)`.
    pub fn read_vec_u8(&self, name: &str) -> Result<(usize, Vec<u8>)> {
        let col = self.require_col(name)?;
        if col.var_desc != VarDescriptor::Fixed || col.col_type != ColumnType::Byte {
            return Err(FitsError::TypeMismatch(format!(
                "column '{}' is not a byte vector column (TFORM={})",
                col.name, col.tform
            )));
        }
        let rep = col.repeat;
        let off = col.offset;
        let rl = self.row_len;
        let mut out = Vec::with_capacity(self.nrows * rep);
        for r in 0..self.nrows {
            let s = r * rl + off;
            out.extend_from_slice(&self.data[s..s + rep]);
        }
        Ok((rep, out))
    }

    /// Read a fixed-width `i16` vector column (`nI`, e.g., `1024I`) as `(repeat, flat_data)`.
    pub fn read_vec_i16(&self, name: &str) -> Result<(usize, Vec<i16>)> {
        let col = self.require_col(name)?;
        if col.var_desc != VarDescriptor::Fixed || col.col_type != ColumnType::I16 {
            return Err(FitsError::TypeMismatch(format!(
                "column '{}' is not an I16 vector column (TFORM={})",
                col.name, col.tform
            )));
        }
        let rep = col.repeat;
        let off = col.offset;
        let rl = self.row_len;
        let mut out = Vec::with_capacity(self.nrows * rep);
        for r in 0..self.nrows {
            let s = r * rl + off;
            let row_slice = &self.data[s..s + rep * 2];
            for &c in row_slice.as_chunks::<2>().0 {
                out.push(i16::from_be_bytes(c));
            }
        }
        Ok((rep, out))
    }

    /// Read a fixed-width `u16` vector column (`nI`, applying `TZERO` such as `32768`) as `(repeat, flat_data)`.
    pub fn read_vec_u16(&self, name: &str) -> Result<(usize, Vec<u16>)> {
        let col = self.require_col(name)?;
        if col.var_desc != VarDescriptor::Fixed || col.col_type != ColumnType::I16 {
            return Err(FitsError::TypeMismatch(format!(
                "column '{}' is not an I16/U16 vector column (TFORM={})",
                col.name, col.tform
            )));
        }
        let rep = col.repeat;
        let off = col.offset;
        let rl = self.row_len;
        let tzero = col.tzero.unwrap_or(0.0).round() as i32;
        let mut out = Vec::with_capacity(self.nrows * rep);
        for r in 0..self.nrows {
            let s = r * rl + off;
            let row_slice = &self.data[s..s + rep * 2];
            for &c in row_slice.as_chunks::<2>().0 {
                let raw = i16::from_be_bytes(c) as i32;
                out.push(raw.wrapping_add(tzero) as u16);
            }
        }
        Ok((rep, out))
    }

    /// Read a fixed-width `i32` vector column (`nJ`) as `(repeat, flat_data)`.
    pub fn read_vec_i32(&self, name: &str) -> Result<(usize, Vec<i32>)> {
        let col = self.require_col(name)?;
        if col.var_desc != VarDescriptor::Fixed || col.col_type != ColumnType::I32 {
            return Err(FitsError::TypeMismatch(format!(
                "column '{}' is not an I32 vector column (TFORM={})",
                col.name, col.tform
            )));
        }
        let rep = col.repeat;
        let off = col.offset;
        let rl = self.row_len;
        let mut out = Vec::with_capacity(self.nrows * rep);
        for r in 0..self.nrows {
            let s = r * rl + off;
            let row_slice = &self.data[s..s + rep * 4];
            for &c in row_slice.as_chunks::<4>().0 {
                out.push(i32::from_be_bytes(c));
            }
        }
        Ok((rep, out))
    }

    /// Read a fixed-width `f32` vector column (`nE`, e.g., `512E`, `16E`) as `(repeat, flat_data)`.
    pub fn read_vec_f32(&self, name: &str) -> Result<(usize, Vec<f32>)> {
        let col = self.require_col(name)?;
        if col.var_desc != VarDescriptor::Fixed {
            return Err(FitsError::TypeMismatch(format!(
                "column '{}' is variable-length (TFORM={}); use read_var_f32",
                col.name, col.tform
            )));
        }
        let rep = col.repeat;
        let off = col.offset;
        let rl = self.row_len;
        let mut out = Vec::with_capacity(self.nrows * rep);
        match col.col_type {
            ColumnType::F32 => {
                for r in 0..self.nrows {
                    let s = r * rl + off;
                    let row_slice = &self.data[s..s + rep * 4];
                    for &c in row_slice.as_chunks::<4>().0 {
                        out.push(f32::from_be_bytes(c));
                    }
                }
            }
            ColumnType::F64 => {
                for r in 0..self.nrows {
                    let s = r * rl + off;
                    let row_slice = &self.data[s..s + rep * 8];
                    for &c in row_slice.as_chunks::<8>().0 {
                        out.push(f64::from_be_bytes(c) as f32);
                    }
                }
            }
            other => {
                return Err(FitsError::TypeMismatch(format!(
                    "cannot read column '{}' of type {:?} as vec f32",
                    col.name, other
                )));
            }
        }
        Ok((rep, out))
    }

    /// Read a fixed-width `f64` vector column (`nD`) as `(repeat, flat_data)`.
    pub fn read_vec_f64(&self, name: &str) -> Result<(usize, Vec<f64>)> {
        let col = self.require_col(name)?;
        if col.var_desc != VarDescriptor::Fixed {
            return Err(FitsError::TypeMismatch(format!(
                "column '{}' is variable-length (TFORM={}); use read_var_f64",
                col.name, col.tform
            )));
        }
        let rep = col.repeat;
        let off = col.offset;
        let rl = self.row_len;
        let mut out = Vec::with_capacity(self.nrows * rep);
        match col.col_type {
            ColumnType::F64 => {
                for r in 0..self.nrows {
                    let s = r * rl + off;
                    let row_slice = &self.data[s..s + rep * 8];
                    for &c in row_slice.as_chunks::<8>().0 {
                        out.push(f64::from_be_bytes(c));
                    }
                }
            }
            ColumnType::F32 => {
                for r in 0..self.nrows {
                    let s = r * rl + off;
                    let row_slice = &self.data[s..s + rep * 4];
                    for &c in row_slice.as_chunks::<4>().0 {
                        out.push(f32::from_be_bytes(c) as f64);
                    }
                }
            }
            other => {
                return Err(FitsError::TypeMismatch(format!(
                    "cannot read column '{}' of type {:?} as vec f64",
                    col.name, other
                )));
            }
        }
        Ok((rep, out))
    }

    /// Resolve a variable-length (`P`/`Q`) OR fixed-length cell for row `row` into `(count, raw_bytes)`.
    fn cell_array_bytes(&self, col: &ColumnSpec, row: usize) -> Result<(usize, &'a [u8])> {
        if row >= self.nrows {
            return Err(FitsError::OutOfBounds(format!(
                "row {row} >= nrows {}",
                self.nrows
            )));
        }
        let cell = self.cell_bytes(row, col.index);
        match col.var_desc {
            VarDescriptor::Fixed => Ok((col.repeat, cell)),
            VarDescriptor::P32 => {
                let count = u32::from_be_bytes(cell[0..4].try_into().unwrap()) as usize;
                let heap_off = u32::from_be_bytes(cell[4..8].try_into().unwrap()) as usize;
                let nbytes = count
                    .checked_mul(col.col_type.elem_size())
                    .ok_or_else(|| FitsError::OutOfBounds(format!("heap count overflow for '{}'", col.name)))?;
                let heap_end = heap_off
                    .checked_add(nbytes)
                    .ok_or_else(|| FitsError::OutOfBounds(format!("heap offset overflow for '{}'", col.name)))?;
                let start = self
                    .theap
                    .checked_add(heap_off)
                    .ok_or_else(|| FitsError::OutOfBounds(format!("THEAP + offset overflow for '{}'", col.name)))?;
                let end = start
                    .checked_add(nbytes)
                    .ok_or_else(|| FitsError::OutOfBounds(format!("heap end overflow for '{}'", col.name)))?;
                if heap_end > self.pcount || end > self.data.len() {
                    return Err(FitsError::OutOfBounds(format!(
                        "heap slice [{start}..{end}] out of bounds for column '{}'",
                        col.name
                    )));
                }
                Ok((count, &self.data[start..end]))
            }
            VarDescriptor::Q64 => {
                let count_u64 = u64::from_be_bytes(cell[0..8].try_into().unwrap());
                let heap_off_u64 = u64::from_be_bytes(cell[8..16].try_into().unwrap());
                let count = usize::try_from(count_u64)
                    .map_err(|_| FitsError::OutOfBounds(format!("Q64 count overflow for '{}'", col.name)))?;
                let heap_off = usize::try_from(heap_off_u64)
                    .map_err(|_| FitsError::OutOfBounds(format!("Q64 offset overflow for '{}'", col.name)))?;
                let nbytes = count
                    .checked_mul(col.col_type.elem_size())
                    .ok_or_else(|| FitsError::OutOfBounds(format!("heap count overflow for '{}'", col.name)))?;
                let heap_end = heap_off
                    .checked_add(nbytes)
                    .ok_or_else(|| FitsError::OutOfBounds(format!("heap offset overflow for '{}'", col.name)))?;
                let start = self
                    .theap
                    .checked_add(heap_off)
                    .ok_or_else(|| FitsError::OutOfBounds(format!("THEAP + offset overflow for '{}'", col.name)))?;
                let end = start
                    .checked_add(nbytes)
                    .ok_or_else(|| FitsError::OutOfBounds(format!("heap end overflow for '{}'", col.name)))?;
                if heap_end > self.pcount || end > self.data.len() {
                    return Err(FitsError::OutOfBounds(format!(
                        "heap slice [{start}..{end}] out of bounds for column '{}'",
                        col.name
                    )));
                }
                Ok((count, &self.data[start..end]))
            }
        }
    }

    /// Read a variable-length or fixed-length `i16`/`i32` array cell at `row` (e.g. OGIP RMF `F_CHAN`, `N_CHAN`).
    pub fn read_var_i16(&self, name: &str, row: usize) -> Result<Vec<i16>> {
        let col = self.require_col(name)?;
        let (count, bytes) = self.cell_array_bytes(col, row)?;
        let mut out = Vec::with_capacity(count);
        match col.col_type {
            ColumnType::I16 => {
                for &c in bytes.as_chunks::<2>().0 {
                    out.push(i16::from_be_bytes(c));
                }
            }
            ColumnType::I32 => {
                for &c in bytes.as_chunks::<4>().0 {
                    out.push(i32::from_be_bytes(c) as i16);
                }
            }
            other => {
                return Err(FitsError::TypeMismatch(format!(
                    "cannot read var_i16 from column '{}' of type {:?}",
                    col.name, other
                )));
            }
        }
        Ok(out)
    }

    /// Read a variable-length or fixed-length `i32` array cell at `row`.
    pub fn read_var_i32(&self, name: &str, row: usize) -> Result<Vec<i32>> {
        let col = self.require_col(name)?;
        let (count, bytes) = self.cell_array_bytes(col, row)?;
        let mut out = Vec::with_capacity(count);
        match col.col_type {
            ColumnType::I32 => {
                for &c in bytes.as_chunks::<4>().0 {
                    out.push(i32::from_be_bytes(c));
                }
            }
            ColumnType::I16 => {
                for &c in bytes.as_chunks::<2>().0 {
                    out.push(i16::from_be_bytes(c) as i32);
                }
            }
            other => {
                return Err(FitsError::TypeMismatch(format!(
                    "cannot read var_i32 from column '{}' of type {:?}",
                    col.name, other
                )));
            }
        }
        Ok(out)
    }

    /// Read a variable-length or fixed-length `f32` array cell at `row` (e.g. OGIP RMF `MATRIX`).
    pub fn read_var_f32(&self, name: &str, row: usize) -> Result<Vec<f32>> {
        let col = self.require_col(name)?;
        let (count, bytes) = self.cell_array_bytes(col, row)?;
        let mut out = Vec::with_capacity(count);
        match col.col_type {
            ColumnType::F32 => {
                for &c in bytes.as_chunks::<4>().0 {
                    out.push(f32::from_be_bytes(c));
                }
            }
            ColumnType::F64 => {
                for &c in bytes.as_chunks::<8>().0 {
                    out.push(f64::from_be_bytes(c) as f32);
                }
            }
            other => {
                return Err(FitsError::TypeMismatch(format!(
                    "cannot read var_f32 from column '{}' of type {:?}",
                    col.name, other
                )));
            }
        }
        Ok(out)
    }

    /// Read a variable-length or fixed-length `f64` array cell at `row`.
    pub fn read_var_f64(&self, name: &str, row: usize) -> Result<Vec<f64>> {
        let col = self.require_col(name)?;
        let (count, bytes) = self.cell_array_bytes(col, row)?;
        let mut out = Vec::with_capacity(count);
        match col.col_type {
            ColumnType::F64 => {
                for &c in bytes.as_chunks::<8>().0 {
                    out.push(f64::from_be_bytes(c));
                }
            }
            ColumnType::F32 => {
                for &c in bytes.as_chunks::<4>().0 {
                    out.push(f32::from_be_bytes(c) as f64);
                }
            }
            other => {
                return Err(FitsError::TypeMismatch(format!(
                    "cannot read var_f64 from column '{}' of type {:?}",
                    col.name, other
                )));
            }
        }
        Ok(out)
    }
}

/// Owned column payload for `BinTableBuilder`.
#[derive(Debug, Clone)]
pub enum ColumnData {
    Logical(Vec<bool>),
    Bytes { repeat: usize, data: Vec<u8> },
    I16 { repeat: usize, data: Vec<i16> },
    I32 { repeat: usize, data: Vec<i32> },
    I64 { repeat: usize, data: Vec<i64> },
    F32 { repeat: usize, data: Vec<f32> },
    F64 { repeat: usize, data: Vec<f64> },
    Strings { width: usize, data: Vec<String> },
    VarI16(Vec<Vec<i16>>),
    VarF32(Vec<Vec<f32>>),
}

struct PendingColumn {
    name: String,
    unit: Option<String>,
    tzero: Option<f64>,
    tscal: Option<f64>,
    tdisp: Option<String>,
    tlmin: Option<i64>,
    tlmax: Option<i64>,
    data: ColumnData,
}

impl PendingColumn {
    fn nrows(&self) -> Result<usize> {
        match &self.data {
            ColumnData::Logical(v) => Ok(v.len()),
            ColumnData::Bytes { repeat, data } => check_div(data.len(), *repeat, &self.name),
            ColumnData::I16 { repeat, data } => check_div(data.len(), *repeat, &self.name),
            ColumnData::I32 { repeat, data } => check_div(data.len(), *repeat, &self.name),
            ColumnData::I64 { repeat, data } => check_div(data.len(), *repeat, &self.name),
            ColumnData::F32 { repeat, data } => check_div(data.len(), *repeat, &self.name),
            ColumnData::F64 { repeat, data } => check_div(data.len(), *repeat, &self.name),
            ColumnData::Strings { data, .. } => Ok(data.len()),
            ColumnData::VarI16(v) => Ok(v.len()),
            ColumnData::VarF32(v) => Ok(v.len()),
        }
    }

    fn tform(&self) -> String {
        match &self.data {
            ColumnData::Logical(_) => "1L".into(),
            ColumnData::Bytes { repeat, .. } => format!("{repeat}B"),
            ColumnData::I16 { repeat, .. } => format!("{repeat}I"),
            ColumnData::I32 { repeat, .. } => format!("{repeat}J"),
            ColumnData::I64 { repeat, .. } => format!("{repeat}K"),
            ColumnData::F32 { repeat, .. } => format!("{repeat}E"),
            ColumnData::F64 { repeat, .. } => format!("{repeat}D"),
            ColumnData::Strings { width, .. } => format!("{width}A"),
            ColumnData::VarI16(rows) => {
                let max_len = rows.iter().map(|r| r.len()).max().unwrap_or(0);
                format!("PI({max_len})")
            }
            ColumnData::VarF32(rows) => {
                let max_len = rows.iter().map(|r| r.len()).max().unwrap_or(0);
                format!("PE({max_len})")
            }
        }
    }

    fn byte_width(&self) -> usize {
        match &self.data {
            ColumnData::Logical(_) => 1,
            ColumnData::Bytes { repeat, .. } => *repeat,
            ColumnData::I16 { repeat, .. } => *repeat * 2,
            ColumnData::I32 { repeat, .. } => *repeat * 4,
            ColumnData::I64 { repeat, .. } => *repeat * 8,
            ColumnData::F32 { repeat, .. } => *repeat * 4,
            ColumnData::F64 { repeat, .. } => *repeat * 8,
            ColumnData::Strings { width, .. } => *width,
            ColumnData::VarI16(_) | ColumnData::VarF32(_) => 8,
        }
    }
}

fn check_div(len: usize, repeat: usize, name: &str) -> Result<usize> {
    if repeat == 0 {
        return Err(FitsError::FormatError(format!(
            "column '{name}' has repeat=0"
        )));
    }
    if !len.is_multiple_of(repeat) {
        return Err(FitsError::FormatError(format!(
            "column '{name}' length {len} not divisible by repeat {repeat}"
        )));
    }
    Ok(len / repeat)
}

/// Column-oriented builder for FITS `BINTABLE` HDUs.
pub struct BinTableBuilder {
    extname: String,
    columns: Vec<PendingColumn>,
    extra_cards: Header,
}

impl BinTableBuilder {
    pub fn new(extname: impl Into<String>) -> Self {
        Self {
            extname: extname.into(),
            columns: Vec::new(),
            extra_cards: Header::new(),
        }
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

    pub fn header_card(mut self, card: Card) -> Self {
        self.extra_cards.set(card);
        self
    }

    pub fn add_comment(mut self, text: impl Into<String>) -> Self {
        self.extra_cards.add_comment(text);
        self
    }

    pub fn add_history(mut self, text: impl Into<String>) -> Self {
        self.extra_cards.add_history(text);
        self
    }

    fn push_col(mut self, col: PendingColumn) -> Self {
        self.columns.push(col);
        self
    }

    pub fn add_f64(self, name: &str, unit: Option<&str>, data: Vec<f64>) -> Self {
        self.push_col(PendingColumn {
            name: name.to_string(),
            unit: unit.map(|s| s.to_string()),
            tzero: None,
            tscal: None,
            tdisp: None,
            tlmin: None,
            tlmax: None,
            data: ColumnData::F64 { repeat: 1, data },
        })
    }

    pub fn add_f32(self, name: &str, unit: Option<&str>, data: Vec<f32>) -> Self {
        self.push_col(PendingColumn {
            name: name.to_string(),
            unit: unit.map(|s| s.to_string()),
            tzero: None,
            tscal: None,
            tdisp: None,
            tlmin: None,
            tlmax: None,
            data: ColumnData::F32 { repeat: 1, data },
        })
    }

    pub fn add_i64(self, name: &str, unit: Option<&str>, data: Vec<i64>) -> Self {
        self.push_col(PendingColumn {
            name: name.to_string(),
            unit: unit.map(|s| s.to_string()),
            tzero: None,
            tscal: None,
            tdisp: None,
            tlmin: None,
            tlmax: None,
            data: ColumnData::I64 { repeat: 1, data },
        })
    }

    pub fn add_i32(self, name: &str, unit: Option<&str>, data: Vec<i32>) -> Self {
        self.push_col(PendingColumn {
            name: name.to_string(),
            unit: unit.map(|s| s.to_string()),
            tzero: None,
            tscal: None,
            tdisp: None,
            tlmin: None,
            tlmax: None,
            data: ColumnData::I32 { repeat: 1, data },
        })
    }

    pub fn add_u32(self, name: &str, unit: Option<&str>, data: Vec<u32>) -> Self {
        let shifted: Vec<i32> = data
            .into_iter()
            .map(|v| (v as i64 - 2_147_483_648i64) as i32)
            .collect();
        self.push_col(PendingColumn {
            name: name.to_string(),
            unit: unit.map(|s| s.to_string()),
            tzero: Some(2_147_483_648.0),
            tscal: Some(1.0),
            tdisp: None,
            tlmin: None,
            tlmax: None,
            data: ColumnData::I32 {
                repeat: 1,
                data: shifted,
            },
        })
    }

    pub fn add_i16(self, name: &str, unit: Option<&str>, data: Vec<i16>) -> Self {
        self.push_col(PendingColumn {
            name: name.to_string(),
            unit: unit.map(|s| s.to_string()),
            tzero: None,
            tscal: None,
            tdisp: None,
            tlmin: None,
            tlmax: None,
            data: ColumnData::I16 { repeat: 1, data },
        })
    }

    pub fn add_u16(self, name: &str, unit: Option<&str>, data: Vec<u16>) -> Self {
        let shifted: Vec<i16> = data
            .into_iter()
            .map(|v| (v as i32 - 32768) as i16)
            .collect();
        self.push_col(PendingColumn {
            name: name.to_string(),
            unit: unit.map(|s| s.to_string()),
            tzero: Some(32768.0),
            tscal: Some(1.0),
            tdisp: None,
            tlmin: None,
            tlmax: None,
            data: ColumnData::I16 {
                repeat: 1,
                data: shifted,
            },
        })
    }

    pub fn add_u8(self, name: &str, unit: Option<&str>, data: Vec<u8>) -> Self {
        self.push_col(PendingColumn {
            name: name.to_string(),
            unit: unit.map(|s| s.to_string()),
            tzero: None,
            tscal: None,
            tdisp: None,
            tlmin: None,
            tlmax: None,
            data: ColumnData::Bytes { repeat: 1, data },
        })
    }

    pub fn add_bool(self, name: &str, data: Vec<bool>) -> Self {
        self.push_col(PendingColumn {
            name: name.to_string(),
            unit: None,
            tzero: None,
            tscal: None,
            tdisp: None,
            tlmin: None,
            tlmax: None,
            data: ColumnData::Logical(data),
        })
    }

    pub fn add_strings(self, name: &str, width: usize, data: Vec<String>) -> Self {
        let w = width.max(1);
        self.push_col(PendingColumn {
            name: name.to_string(),
            unit: None,
            tzero: None,
            tscal: None,
            tdisp: None,
            tlmin: None,
            tlmax: None,
            data: ColumnData::Strings { width: w, data },
        })
    }

    pub fn add_vec_u8(self, name: &str, unit: Option<&str>, repeat: usize, data: Vec<u8>) -> Self {
        self.push_col(PendingColumn {
            name: name.to_string(),
            unit: unit.map(|s| s.to_string()),
            tzero: None,
            tscal: None,
            tdisp: None,
            tlmin: None,
            tlmax: None,
            data: ColumnData::Bytes { repeat, data },
        })
    }

    pub fn add_vec_i16(
        self,
        name: &str,
        unit: Option<&str>,
        repeat: usize,
        data: Vec<i16>,
    ) -> Self {
        self.push_col(PendingColumn {
            name: name.to_string(),
            unit: unit.map(|s| s.to_string()),
            tzero: None,
            tscal: None,
            tdisp: None,
            tlmin: None,
            tlmax: None,
            data: ColumnData::I16 { repeat, data },
        })
    }

    pub fn add_vec_u16(
        self,
        name: &str,
        unit: Option<&str>,
        repeat: usize,
        data: Vec<u16>,
    ) -> Self {
        let shifted: Vec<i16> = data
            .into_iter()
            .map(|v| (v as i32 - 32768) as i16)
            .collect();
        self.push_col(PendingColumn {
            name: name.to_string(),
            unit: unit.map(|s| s.to_string()),
            tzero: Some(32768.0),
            tscal: Some(1.0),
            tdisp: None,
            tlmin: None,
            tlmax: None,
            data: ColumnData::I16 {
                repeat,
                data: shifted,
            },
        })
    }

    pub fn add_vec_i32(
        self,
        name: &str,
        unit: Option<&str>,
        repeat: usize,
        data: Vec<i32>,
    ) -> Self {
        self.push_col(PendingColumn {
            name: name.to_string(),
            unit: unit.map(|s| s.to_string()),
            tzero: None,
            tscal: None,
            tdisp: None,
            tlmin: None,
            tlmax: None,
            data: ColumnData::I32 { repeat, data },
        })
    }

    pub fn add_vec_f32(
        self,
        name: &str,
        unit: Option<&str>,
        repeat: usize,
        data: Vec<f32>,
    ) -> Self {
        self.push_col(PendingColumn {
            name: name.to_string(),
            unit: unit.map(|s| s.to_string()),
            tzero: None,
            tscal: None,
            tdisp: None,
            tlmin: None,
            tlmax: None,
            data: ColumnData::F32 { repeat, data },
        })
    }

    pub fn add_vec_f64(
        self,
        name: &str,
        unit: Option<&str>,
        repeat: usize,
        data: Vec<f64>,
    ) -> Self {
        self.push_col(PendingColumn {
            name: name.to_string(),
            unit: unit.map(|s| s.to_string()),
            tzero: None,
            tscal: None,
            tdisp: None,
            tlmin: None,
            tlmax: None,
            data: ColumnData::F64 { repeat, data },
        })
    }

    pub fn add_var_i16(self, name: &str, unit: Option<&str>, rows: Vec<Vec<i16>>) -> Self {
        self.push_col(PendingColumn {
            name: name.to_string(),
            unit: unit.map(|s| s.to_string()),
            tzero: None,
            tscal: None,
            tdisp: None,
            tlmin: None,
            tlmax: None,
            data: ColumnData::VarI16(rows),
        })
    }

    pub fn add_var_f32(self, name: &str, unit: Option<&str>, rows: Vec<Vec<f32>>) -> Self {
        self.push_col(PendingColumn {
            name: name.to_string(),
            unit: unit.map(|s| s.to_string()),
            tzero: None,
            tscal: None,
            tdisp: None,
            tlmin: None,
            tlmax: None,
            data: ColumnData::VarF32(rows),
        })
    }

    /// Set `TLMINn` / `TLMAXn` on the most recently added column.
    pub fn with_last_col_range(mut self, tlmin: i64, tlmax: i64) -> Self {
        if let Some(last) = self.columns.last_mut() {
            last.tlmin = Some(tlmin);
            last.tlmax = Some(tlmax);
        }
        self
    }

    /// Build the finished `BINTABLE` `Hdu`.
    pub fn build(self) -> Result<Hdu> {
        let mut nrows_opt: Option<usize> = None;
        for col in &self.columns {
            let nr = col.nrows()?;
            if let Some(prev) = nrows_opt {
                if prev != nr {
                    return Err(FitsError::FormatError(format!(
                        "column '{}' has {nr} rows, expected {prev}",
                        col.name
                    )));
                }
            } else {
                nrows_opt = Some(nr);
            }
        }
        let nrows = nrows_opt.unwrap_or(0);
        let row_len: usize = self.columns.iter().map(|c| c.byte_width()).sum();

        let mut main_buf = vec![0u8; nrows * row_len];
        let mut heap_buf: Vec<u8> = Vec::new();

        // Fill column-by-column into row-interleaved big-endian FITS layout.
        let mut col_offset = 0usize;
        for col in &self.columns {
            let w = col.byte_width();
            match &col.data {
                ColumnData::Logical(v) => {
                    for (r, &b) in v.iter().enumerate() {
                        main_buf[r * row_len + col_offset] = if b { b'T' } else { b'F' };
                    }
                }
                ColumnData::Bytes { repeat, data } => {
                    let rep = *repeat;
                    for r in 0..nrows {
                        let dst = r * row_len + col_offset;
                        main_buf[dst..dst + rep].copy_from_slice(&data[r * rep..(r + 1) * rep]);
                    }
                }
                ColumnData::I16 { repeat, data } => {
                    let rep = *repeat;
                    for r in 0..nrows {
                        let mut dst = r * row_len + col_offset;
                        for &val in &data[r * rep..(r + 1) * rep] {
                            main_buf[dst..dst + 2].copy_from_slice(&val.to_be_bytes());
                            dst += 2;
                        }
                    }
                }
                ColumnData::I32 { repeat, data } => {
                    let rep = *repeat;
                    for r in 0..nrows {
                        let mut dst = r * row_len + col_offset;
                        for &val in &data[r * rep..(r + 1) * rep] {
                            main_buf[dst..dst + 4].copy_from_slice(&val.to_be_bytes());
                            dst += 4;
                        }
                    }
                }
                ColumnData::I64 { repeat, data } => {
                    let rep = *repeat;
                    for r in 0..nrows {
                        let mut dst = r * row_len + col_offset;
                        for &val in &data[r * rep..(r + 1) * rep] {
                            main_buf[dst..dst + 8].copy_from_slice(&val.to_be_bytes());
                            dst += 8;
                        }
                    }
                }
                ColumnData::F32 { repeat, data } => {
                    let rep = *repeat;
                    for r in 0..nrows {
                        let mut dst = r * row_len + col_offset;
                        for &val in &data[r * rep..(r + 1) * rep] {
                            main_buf[dst..dst + 4].copy_from_slice(&val.to_be_bytes());
                            dst += 4;
                        }
                    }
                }
                ColumnData::F64 { repeat, data } => {
                    let rep = *repeat;
                    for r in 0..nrows {
                        let mut dst = r * row_len + col_offset;
                        for &val in &data[r * rep..(r + 1) * rep] {
                            main_buf[dst..dst + 8].copy_from_slice(&val.to_be_bytes());
                            dst += 8;
                        }
                    }
                }
                ColumnData::Strings { width, data } => {
                    let width = *width;
                    for (r, s) in data.iter().enumerate() {
                        let dst = r * row_len + col_offset;
                        main_buf[dst..dst + width].fill(b' ');
                        let bytes = s.as_bytes();
                        let n = bytes.len().min(width);
                        main_buf[dst..dst + n].copy_from_slice(&bytes[..n]);
                    }
                }
                ColumnData::VarI16(rows) => {
                    for (r, vals) in rows.iter().enumerate() {
                        let count = vals.len() as u32;
                        let heap_off = heap_buf.len() as u32;
                        for &v in vals {
                            heap_buf.extend_from_slice(&v.to_be_bytes());
                        }
                        let dst = r * row_len + col_offset;
                        main_buf[dst..dst + 4].copy_from_slice(&count.to_be_bytes());
                        main_buf[dst + 4..dst + 8].copy_from_slice(&heap_off.to_be_bytes());
                    }
                }
                ColumnData::VarF32(rows) => {
                    for (r, vals) in rows.iter().enumerate() {
                        let count = vals.len() as u32;
                        let heap_off = heap_buf.len() as u32;
                        for &v in vals {
                            heap_buf.extend_from_slice(&v.to_be_bytes());
                        }
                        let dst = r * row_len + col_offset;
                        main_buf[dst..dst + 4].copy_from_slice(&count.to_be_bytes());
                        main_buf[dst + 4..dst + 8].copy_from_slice(&heap_off.to_be_bytes());
                    }
                }
            }
            col_offset += w;
        }

        let pcount = heap_buf.len();
        if !heap_buf.is_empty() {
            main_buf.extend_from_slice(&heap_buf);
        }

        let mut hdr = Header::new();
        hdr.set_str("XTENSION", "BINTABLE", "binary table extension");
        hdr.set_i64("BITPIX", 8, "8-bit bytes");
        hdr.set_i64("NAXIS", 2, "2-dimensional binary table");
        hdr.set_i64("NAXIS1", row_len as i64, "width of table in bytes");
        hdr.set_i64("NAXIS2", nrows as i64, "number of rows in table");
        hdr.set_i64("PCOUNT", pcount as i64, "size of special data area");
        hdr.set_i64("GCOUNT", 1, "one data group (required keyword)");
        hdr.set_i64("TFIELDS", self.columns.len() as i64, "number of fields in each row");

        for (i, col) in self.columns.iter().enumerate() {
            let n = i + 1;
            hdr.set_str(&format!("TTYPE{n}"), &col.name, "label for field");
            hdr.set_str(&format!("TFORM{n}"), col.tform(), "data format of field");
            if let Some(unit) = &col.unit {
                if !unit.is_empty() {
                    hdr.set_str(&format!("TUNIT{n}"), unit, "physical unit of field");
                }
            }
            if let Some(tz) = col.tzero {
                if tz.fract() == 0.0 && tz.abs() < (i64::MAX as f64) {
                    hdr.set_i64(&format!("TZERO{n}"), tz as i64, "offset for unsigned integers");
                } else {
                    hdr.set_f64(&format!("TZERO{n}"), tz, "zero point offset");
                }
            }
            if let Some(ts) = col.tscal {
                hdr.set_f64(&format!("TSCAL{n}"), ts, "scale factor");
            }
            if let Some(td) = &col.tdisp {
                hdr.set_str(&format!("TDISP{n}"), td, "display format");
            }
            if let Some(tlmin) = col.tlmin {
                hdr.set_i64(&format!("TLMIN{n}"), tlmin, "minimum legal value");
            }
            if let Some(tlmax) = col.tlmax {
                hdr.set_i64(&format!("TLMAX{n}"), tlmax, "maximum legal value");
            }
        }

        if !self.extname.is_empty() {
            hdr.set_str("EXTNAME", &self.extname, "name of this binary table extension");
        }

        for card in self.extra_cards.cards() {
            // Do not overwrite structural BINTABLE keywords.
            let k = card.keyword.as_str();
            if matches!(
                k,
                "XTENSION" | "BITPIX" | "NAXIS" | "NAXIS1" | "NAXIS2" | "PCOUNT" | "GCOUNT" | "TFIELDS"
            ) {
                continue;
            }
            hdr.set(card.clone());
        }

        Ok(Hdu::new(hdr, main_buf))
    }
}
