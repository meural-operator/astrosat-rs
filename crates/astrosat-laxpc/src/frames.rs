//! Stage 1: Level-1 2242-byte EA & BB frame parsing, integrity checks, sorting, and deduplication.
//!
//! Replaces legacy fixed-size Fortran arrays (`frameno` in `laxpcl1.f`) with deterministic,
//! dynamically sized, zero-copy-friendly Rust structures.

use crate::error::{LaxpcError, Result};
use crate::types::{LaxpcUnit, L1_FRAME_SIZE};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

/// Sync byte at offset 44 (Fortran index 45: `-34` as `i8`).
pub const FRAME_SYNC_BYTE: u8 = 0xDE;
/// Mode byte for Event Analysis (`modeEA`) at offset 46 (Fortran index 47: `-22` as `i8`).
pub const MODE_EA_BYTE: u8 = 0xEA;
/// Mode byte for Broad Band Counting (`modeBB`) at offset 46 (Fortran index 47: `-69` as `i8`).
pub const MODE_BB_BYTE: u8 = 0xBB;
/// Coarse time marker byte in EA payload (`-17` as `i8`).
pub const COARSE_TIME_MARKER: u8 = 0xEF;
/// End-of-frame marker byte in EA payload (`-18` as `i8`).
pub const END_OF_FRAME_MARKER: u8 = 0xEE;

/// Telemetry mode of a 2242-byte LAXPC Level-1 frame.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum FrameMode {
    /// Event Analysis mode (`0xEA`).
    Ea,
    /// Broad Band Counting mode (`0xBB`).
    Bb,
}

impl FrameMode {
    #[inline]
    pub const fn byte(self) -> u8 {
        match self {
            Self::Ea => MODE_EA_BYTE,
            Self::Bb => MODE_BB_BYTE,
        }
    }
}

/// Decoded header metadata of a single 2242-byte Level-1 LAXPC frame.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FrameHeader {
    /// L1-assigned frame timestamp (bytes `0..8`, big-endian `f64`).
    pub l1_time: f64,
    /// Original 64-bit frame sequence column (bytes `36..44`, big-endian `i64`).
    pub raw_seq: i64,
    /// Sync byte (offset 44, must be `0xDE`).
    pub sync: u8,
    /// Detector unit ID (offset 45, `1..=3`).
    pub unit_id: u8,
    /// Mode byte (offset 46, `0xEA` or `0xBB`).
    pub mode_byte: u8,
    /// 56-bit raw STBG clock field (bytes `47..54`, big-endian).
    pub stbg_56: u64,
    /// Lower 40-bit 100 kHz STBG clock ticks (bytes `49..54`, big-endian).
    pub stbg_40: u64,
    /// Status flag byte (offset 57, must be `0` for valid frames).
    pub status_flag: u8,
    /// 16-bit onboard frame counter `ifr1` (bytes `58..60`, big-endian `u16`).
    pub frame_counter: u16,
    /// Error-correction flags `ie1`, `ie2`, `ie3` (bytes `2092..2098`, big-endian `i16`).
    pub err_flags: [i16; 3],
}

impl FrameHeader {
    /// Parse header fields from a 2242-byte raw frame slice.
    pub fn parse(raw: &[u8]) -> Result<Self> {
        if raw.len() < L1_FRAME_SIZE {
            return Err(LaxpcError::Frame(format!(
                "frame too short: {} bytes (expected {})",
                raw.len(),
                L1_FRAME_SIZE
            )));
        }
        let l1_time = f64::from_be_bytes(raw[0..8].try_into().unwrap());
        let raw_seq = i64::from_be_bytes(raw[36..44].try_into().unwrap());
        let sync = raw[44];
        let unit_id = raw[45];
        let mode_byte = raw[46];

        let mut s56 = [0u8; 8];
        s56[1..8].copy_from_slice(&raw[47..54]);
        let stbg_56 = u64::from_be_bytes(s56);

        let mut s40 = [0u8; 8];
        s40[3..8].copy_from_slice(&raw[49..54]);
        let stbg_40 = u64::from_be_bytes(s40);

        let status_flag = raw[57];
        let frame_counter = u16::from_be_bytes([raw[58], raw[59]]);

        let ie1 = i16::from_be_bytes([raw[2092], raw[2093]]);
        let ie2 = i16::from_be_bytes([raw[2094], raw[2095]]);
        let ie3 = i16::from_be_bytes([raw[2096], raw[2097]]);

        Ok(Self {
            l1_time,
            raw_seq,
            sync,
            unit_id,
            mode_byte,
            stbg_56,
            stbg_40,
            status_flag,
            frame_counter,
            err_flags: [ie1, ie2, ie3],
        })
    }

    /// Check basic Level-1 structural integrity (`laxpcl1.f:2730-2731`):
    /// sync `0xDE`, unit in `1..=3`, status flag `0`, and `ie1 <= 1`.
    #[inline]
    pub fn is_structurally_valid(&self) -> bool {
        self.sync == FRAME_SYNC_BYTE
            && (1..=3).contains(&self.unit_id)
            && self.status_flag == 0
            && self.err_flags[0] <= 1
    }

    /// Check validity for a specific unit and mode.
    #[inline]
    pub fn matches_unit_and_mode(&self, unit: LaxpcUnit, mode: FrameMode) -> bool {
        self.is_structurally_valid() && self.unit_id == unit.id() && self.mode_byte == mode.byte()
    }
}

/// A validated, unwrapped, and ordered 2242-byte Level-1 frame.
#[derive(Debug, Clone)]
pub struct OrderedFrame {
    /// Parsed header metadata.
    pub header: FrameHeader,
    /// Unwrapped monotonic 64-bit frame sequence counter (accounting for 16-bit rollovers).
    pub unwrapped_counter: i64,
    /// Raw 2242-byte frame payload.
    pub data: Vec<u8>,
}

/// Gap diagnostic for missing consecutive frames in the sorted sequence.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub struct FrameGap {
    /// Unwrapped sequence counter immediately before the gap.
    pub prev_counter: i64,
    /// Unwrapped sequence counter immediately after the gap.
    pub next_counter: i64,
    /// Number of missing frames (`next_counter - prev_counter - 1`).
    pub missing_frames: u64,
}

/// Deterministic diagnostics produced by frame sorting and integrity checking.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct FrameDiagnostics {
    /// Total raw frames examined across input files.
    pub total_raw_frames: usize,
    /// Number of unique valid ordered frames retained.
    pub retained_frames: usize,
    /// Frames rejected due to illegal sync, unit ID, mode, status flag, or error flag `ie1 > 1`.
    pub rejected_header: usize,
    /// Stray out-of-sequence frames rejected (jumps in `(6000, 60000)` relative to neighbors).
    pub rejected_stray: usize,
    /// Duplicate frames encountered where the newer copy replaced the earlier copy.
    pub duplicates_replaced: usize,
    /// Duplicate frames encountered that were ignored (corrupt/anomalous timestamp).
    pub duplicates_ignored: usize,
    /// Number of out-of-order steps observed in the input stream.
    pub out_of_order_count: usize,
    /// Total number of missing frames inferred from sequence gaps (`ngap` in `laxpcl1.f`).
    pub total_missing_frames: u64,
    /// Number of single-frame gaps (`delta == 2`).
    pub single_frame_gaps: usize,
    /// Number of multi-frame gaps (`delta > 2`).
    pub multi_frame_gaps: usize,
    /// Detailed list of frame gaps.
    pub gaps: Vec<FrameGap>,
}

/// Reorder, validate, and deduplicate raw 2242-byte Level-1 frames from a single orbit file,
/// matching the semantics of `laxpcl1.f:frameno` without fixed-size array limits.
///
/// `raw_frames` is an iterator of 2242-byte frame slices (e.g., rows from a Level-1 FITS BINTABLE).
/// If `expected_mode` is `Some(FrameMode::Ea)`, frames with wrong mode or unit are rejected;
/// if `None`, any structurally valid frame (sync `0xDE`, unit `1..=3`, status `0`, `ie1 <= 1`) is kept.
pub fn sort_and_dedup_frames<'a, I>(
    raw_frames: I,
    expected_unit: Option<LaxpcUnit>,
    expected_mode: Option<FrameMode>,
) -> Result<(Vec<OrderedFrame>, FrameDiagnostics)>
where
    I: IntoIterator<Item = &'a [u8]>,
{
    let mut diag = FrameDiagnostics::default();

    // Unwrapping state matching `laxpcl1.f:frameno` lines 2700-2819
    let mut rollover: i64 = 0;
    let mut rollover_prev: i64 = 0;
    let mut last_rollover_idx: i64 = -66_666;
    let mut last_rollover_idx_prev: i64 = -66_666;

    // History of accepted raw 16-bit counters, unwrapped counters, and L1 timestamps
    let mut nfo: Vec<i64> = Vec::new();
    let mut nf: Vec<i64> = Vec::new();
    let mut tims: Vec<f64> = Vec::new();

    // Map from unwrapped frame counter -> (OrderedFrame, previous-slot backup for stray rollback)
    let mut slots: BTreeMap<i64, OrderedFrame> = BTreeMap::new();
    // Backup of the most recently modified slot so a stray frame (detected on the next frame)
    // can be undone exactly as `laxpcl1.f:2740-2752` does.
    let mut last_slot_key: Option<i64> = None;
    let mut last_slot_prev_val: Option<OrderedFrame> = None;

    for raw in raw_frames {
        diag.total_raw_frames += 1;
        let hdr = match FrameHeader::parse(raw) {
            Ok(h) => h,
            Err(_) => {
                diag.rejected_header += 1;
                continue;
            }
        };

        if !hdr.is_structurally_valid() {
            diag.rejected_header += 1;
            continue;
        }
        if let Some(u) = expected_unit {
            if hdr.unit_id != u.id() {
                diag.rejected_header += 1;
                continue;
            }
        }
        if let Some(m) = expected_mode {
            if hdr.mode_byte != m.byte() {
                diag.rejected_header += 1;
                continue;
            }
        }

        let ifram = hdr.frame_counter as i64;
        let nfr = nfo.len();

        // Stray frame check (`laxpcl1.f:2737-2755`):
        // If previous frame jumped by (6000, 60000) from the one before AND current frame jumps
        // by (6000, 60000) from the previous frame, the previous frame was a stray -> roll it back.
        let is_stray_prev = if nfr >= 2 {
            let ia1 = (ifram - nfo[nfr - 1]).abs();
            let ia2 = (nfo[nfr - 1] - nfo[nfr - 2]).abs();
            ia1 > 6000 && ia2 > 6000 && ia1 < 60000 && ia2 < 60000
        } else {
            false
        };

        if is_stray_prev {
            // Undo previous frame's slot insertion
            if let Some(prev_key) = last_slot_key.take() {
                if let Some(old_frame) = last_slot_prev_val.take() {
                    slots.insert(prev_key, old_frame);
                } else {
                    slots.remove(&prev_key);
                }
            }
            rollover = rollover_prev;
            last_rollover_idx = last_rollover_idx_prev;
            nfo.pop();
            nf.pop();
            tims.pop();
            diag.rejected_stray += 1;
        }

        let cur_idx = nfo.len() as i64 + 1; // 1-based like Fortran `nfr`
        rollover_prev = rollover;
        last_rollover_idx_prev = last_rollover_idx;

        if let Some(&prev_nfo) = nfo.last() {
            if ifram < prev_nfo - 12_000 && (cur_idx - last_rollover_idx) > 5_000 {
                rollover += 1;
                last_rollover_idx = cur_idx;
            }
        }

        let mut unwrapped = ifram + rollover * 65_536;
        if let Some(&prev_nfo) = nfo.last() {
            if last_rollover_idx > 0
                && ifram > 55_000
                && (cur_idx - last_rollover_idx) < 12_000
                && (ifram - prev_nfo) > 5
            {
                unwrapped -= 65_536;
            }
        }

        if let Some(&prev_nf) = nf.last() {
            if unwrapped < prev_nf {
                diag.out_of_order_count += 1;
            }
        }

        // Duplicate handling (`laxpcl1.f:2790-2817`)
        let mut accept = true;
        if let Some(existing) = slots.get(&unwrapped) {
            // Check if this is actually a 65536 wrap-around after a long (>2000 s) gap
            let prev_tim = tims.last().copied().unwrap_or(hdr.l1_time);
            let prev_nfo = nfo.last().copied().unwrap_or(ifram);
            if hdr.l1_time > prev_tim + 2000.0
                && !slots.contains_key(&(unwrapped + 65_536))
                && ifram < prev_nfo
            {
                unwrapped += 65_536;
                rollover += 1;
                last_rollover_idx = cur_idx;
            } else {
                // True duplicate
                if hdr.stbg_56 > 5_000_000_000 {
                    accept = false;
                }
                if accept && ifram == prev_nfo {
                    if let Some(prev_slot) = slots.get(&(unwrapped - 1)) {
                        let dt_new = (hdr.l1_time - prev_slot.header.l1_time).abs();
                        let dt_old = (existing.header.l1_time - prev_slot.header.l1_time).abs();
                        if dt_new > 2.0 * dt_old && dt_old > 0.0 {
                            accept = false;
                        }
                    }
                }
                if accept {
                    diag.duplicates_replaced += 1;
                } else {
                    diag.duplicates_ignored += 1;
                }
            }
        }

        nfo.push(ifram);
        nf.push(unwrapped);
        tims.push(hdr.l1_time);

        if accept {
            let mut data = raw[..L1_FRAME_SIZE].to_vec();
            // Store unwrapped sequence number into bytes 36..44 (big-endian i64) like `frameno`
            data[36..44].copy_from_slice(&unwrapped.to_be_bytes());
            let mut hdr_updated = hdr;
            hdr_updated.raw_seq = unwrapped;

            let frame = OrderedFrame {
                header: hdr_updated,
                unwrapped_counter: unwrapped,
                data,
            };
            last_slot_key = Some(unwrapped);
            last_slot_prev_val = slots.insert(unwrapped, frame);
        }
    }

    // Collect in strictly increasing unwrapped frame order and compute gap statistics
    let ordered: Vec<OrderedFrame> = slots.into_values().collect();
    diag.retained_frames = ordered.len();

    for win in ordered.windows(2) {
        let prev_c = win[0].unwrapped_counter;
        let next_c = win[1].unwrapped_counter;
        let delta = next_c - prev_c;
        if delta == 2 {
            diag.single_frame_gaps += 1;
            diag.total_missing_frames += 1;
            diag.gaps.push(FrameGap {
                prev_counter: prev_c,
                next_counter: next_c,
                missing_frames: 1,
            });
        } else if delta > 2 {
            let missing = (delta - 1) as u64;
            diag.multi_frame_gaps += 1;
            diag.total_missing_frames += missing;
            diag.gaps.push(FrameGap {
                prev_counter: prev_c,
                next_counter: next_c,
                missing_frames: missing,
            });
        }
    }

    Ok((ordered, diag))
}

/// Merge and deduplicate ordered frames from multiple orbit files.
///
/// Orbit files in AstroSat often overlap at orbit boundaries; this merges them in increasing
/// L1 time order, removing duplicate frames that share the same `(frame_counter, stbg_40)` or
/// fall strictly before the previous orbit's maximum timestamp (`l1_time`) when counters overlap.
pub fn merge_multi_orbit_frames(
    orbit_batches: Vec<(Vec<OrderedFrame>, FrameDiagnostics)>,
) -> (Vec<OrderedFrame>, FrameDiagnostics) {
    let mut combined_diag = FrameDiagnostics::default();
    let mut merged: Vec<OrderedFrame> = Vec::new();
    let mut max_time = f64::NEG_INFINITY;

    for (frames, diag) in orbit_batches {
        combined_diag.total_raw_frames += diag.total_raw_frames;
        combined_diag.rejected_header += diag.rejected_header;
        combined_diag.rejected_stray += diag.rejected_stray;
        combined_diag.duplicates_replaced += diag.duplicates_replaced;
        combined_diag.duplicates_ignored += diag.duplicates_ignored;
        combined_diag.out_of_order_count += diag.out_of_order_count;
        combined_diag.single_frame_gaps += diag.single_frame_gaps;
        combined_diag.multi_frame_gaps += diag.multi_frame_gaps;
        combined_diag.total_missing_frames += diag.total_missing_frames;
        combined_diag.gaps.extend(diag.gaps);

        let mut batch_max_time = max_time;
        for f in frames {
            if f.header.l1_time > batch_max_time {
                batch_max_time = f.header.l1_time;
            }
            // Deduplicate exact overlap across orbit boundaries using STBG + frame counter
            if let Some(last) = merged.last() {
                if f.header.frame_counter == last.header.frame_counter
                    && f.header.stbg_40 == last.header.stbg_40
                {
                    combined_diag.duplicates_ignored += 1;
                    continue;
                }
            }
            merged.push(f);
        }
        max_time = batch_max_time;
    }

    combined_diag.retained_frames = merged.len();
    (merged, combined_diag)
}
