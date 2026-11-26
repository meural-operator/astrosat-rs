//! FITS 4.0 `DATASUM` and `CHECKSUM` calculation per Seaman, Pence & Rots (2012).

use crate::card::Header;
use crate::error::{FitsError, Result};

/// FITS block size in bytes.
pub const FITS_BLOCK_SIZE: usize = 2880;

/// Excluded ASCII punctuation bytes in the 16-char CHECKSUM encoding.
const fn is_excluded(b: u8) -> bool {
    (b >= 0x3A && b <= 0x40) || (b >= 0x5B && b <= 0x60)
}

/// Accumulate ones'-complement 32-bit big-endian sum over a byte slice.
///
/// If `data.len()` is not a multiple of 4, trailing bytes are treated as zero-padded
/// (consistent with FITS 2880-byte zero-padded data blocks).
pub fn accumulate_checksum(initial: u32, data: &[u8]) -> u32 {
    let mut sum: u64 = initial as u64;
    let (chunks, rem) = data.as_chunks::<4>();

    for &c in chunks {
        let word = u32::from_be_bytes(c);
        sum += word as u64;
    }

    if !rem.is_empty() {
        let mut tail = [0u8; 4];
        tail[..rem.len()].copy_from_slice(rem);
        sum += u32::from_be_bytes(tail) as u64;
    }

    while (sum >> 32) != 0 {
        sum = (sum & 0xFFFF_FFFF) + (sum >> 32);
    }

    sum as u32
}

/// Compute the 32-bit ones'-complement `DATASUM` of an HDU data payload.
///
/// Trailing zero-padding to a multiple of 2880 bytes does not change the sum.
pub fn datasum(data: &[u8]) -> u32 {
    accumulate_checksum(0, data)
}

/// Encode a 32-bit unsigned value into the 16-character FITS `CHECKSUM` ASCII representation.
///
/// Pass `!hdu_sum` to obtain the value that makes the entire HDU sum to `0xFFFF_FFFF` (`-0`).
pub fn encode_checksum(value: u32) -> String {
    let bytes = value.to_be_bytes();
    let mut unrotated = [0u8; 16];

    for (i, &b) in bytes.iter().enumerate() {
        let quot = b / 4 + 0x30;
        let rem = b % 4;
        let mut ch = [quot + rem, quot, quot, quot];

        // Eliminate excluded punctuation (0x3A..=0x40 and 0x5B..=0x60) while keeping
        // the sum ch[0]+ch[1]+ch[2]+ch[3] invariant.
        loop {
            let mut changed = false;
            for j in 0..4 {
                for k in (j + 1)..4 {
                    if is_excluded(ch[j]) || is_excluded(ch[k]) {
                        if ch[j] < 0x7A && ch[k] > 0x30 {
                            ch[j] += 1;
                            ch[k] -= 1;
                            changed = true;
                        } else if ch[k] < 0x7A && ch[j] > 0x30 {
                            ch[k] += 1;
                            ch[j] -= 1;
                            changed = true;
                        }
                    }
                }
            }
            if !changed {
                break;
            }
        }

        // Distribute into the 4 big-endian words at byte lane `i`.
        for w in 0..4 {
            unrotated[w * 4 + i] = ch[w];
        }
    }

    // Rotate right by 1 character because the CHECKSUM string begins at byte offset 11
    // (1-indexed column 12), which is 3 mod 4.
    let mut rotated = [0u8; 16];
    for idx in 0..16 {
        rotated[(idx + 1) % 16] = unrotated[idx];
    }

    String::from_utf8(rotated.to_vec()).expect("CHECKSUM encoding is valid ASCII")
}

/// Decode a 16-character FITS `CHECKSUM` string back into the 32-bit value it represents.
pub fn decode_checksum(encoded: &str) -> Result<u32> {
    let bytes = encoded.as_bytes();
    if bytes.len() != 16 {
        return Err(FitsError::FormatError(format!(
            "CHECKSUM must be 16 ASCII characters, got {}",
            bytes.len()
        )));
    }
    // Undo the right-rotation by 1 (rotate left by 1).
    let mut unrotated = [0u8; 16];
    for idx in 0..16 {
        unrotated[idx] = bytes[(idx + 1) % 16];
    }
    // Sum the 4 big-endian 32-bit words and subtract 4 * 0x30303030 = 0xC0C0C0C0 in 1's complement.
    let raw_sum = accumulate_checksum(0, &unrotated);
    let offset_comp = !0xC0C0_C0C0u32; // adding !x is subtracting x in 1's complement
    let val = add_ones_complement(raw_sum, offset_comp);
    Ok(if val == 0xFFFF_FFFF { 0 } else { val })
}

#[inline]
pub fn add_ones_complement(a: u32, b: u32) -> u32 {
    let mut s = (a as u64) + (b as u64);
    while (s >> 32) != 0 {
        s = (s & 0xFFFF_FFFF) + (s >> 32);
    }
    s as u32
}

/// Compute and insert/update `DATASUM` and `CHECKSUM` cards in `header` for the given `data`.
///
/// Returns `(datasum_u32, checksum_str)`.
pub fn compute_hdu_checksums(header: &mut Header, data: &[u8]) -> (u32, String) {
    let dsum = datasum(data);
    header.set_str(
        "DATASUM",
        dsum.to_string(),
        "data unit checksum updated by astrosat-fits",
    );
    header.set_str(
        "CHECKSUM",
        "0000000000000000",
        "HDU checksum updated by astrosat-fits",
    );

    let hdr_bytes = header.to_fits_bytes();
    let hdr_sum = accumulate_checksum(0, &hdr_bytes);
    let total = add_ones_complement(hdr_sum, dsum);
    let encoded = encode_checksum(!total);

    header.set_str(
        "CHECKSUM",
        &encoded,
        "HDU checksum updated by astrosat-fits",
    );

    (dsum, encoded)
}

/// Verify an HDU's `DATASUM` and `CHECKSUM` given its raw 2880-padded header and data bytes.
///
/// Returns `(datasum_ok, checksum_ok)`. If a keyword is absent, its corresponding flag is `true`.
pub fn verify_hdu_checksum(header: &Header, raw_header: &[u8], raw_data: &[u8]) -> (bool, bool) {
    let dsum = datasum(raw_data);

    let datasum_ok = match header.get_str("DATASUM") {
        Some(s) => s.trim().parse::<u32>().map(|v| v == dsum).unwrap_or(false),
        None => true,
    };

    let checksum_ok = match header.get_str("CHECKSUM") {
        Some(_) => {
            let hdr_sum = accumulate_checksum(0, raw_header);
            let total = add_ones_complement(hdr_sum, dsum);
            total == 0xFFFF_FFFF || total == 0
        }
        None => true,
    };

    (datasum_ok, checksum_ok)
}
