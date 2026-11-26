//! Stage 2: Event Analysis (`modeEA`) 5-byte event unpacking and Xe K-escape double-event logic.
//!
//! Faithfully implements the bit-exact decoding of `laxpcl1.f:conver1`, including:
//! - 40-bit 100 kHz STBG clock reconstruction + 5th-byte rollover correction (post-Sep 2024 fix)
//! - `0xEF` coarse-time marker carry propagation, glitch rejection, and retroactive pre-marker sync
//! - 5-byte event record unpacking (`ic1`, `ic2`, `ia1`, `ia2`, `k1`, `k2`, fine time `10 us`)
//! - Single-event and Xe K-escape double-event summing + dominant anode assignment
//! - Frame-end abort/truncation (`0xEE` or `na0 > 4`)
//! - Internal frame gap detection (`>= 50 ms`) for exposure correction

use crate::frames::{OrderedFrame, COARSE_TIME_MARKER, END_OF_FRAME_MARKER, FRAME_SYNC_BYTE, MODE_EA_BYTE};
use crate::types::{EventFilter, LaxpcUnit, L1_FRAME_SIZE, STBG_CLOCK_HZ};
use serde::{Deserialize, Serialize};

/// Reference STBG-to-UTC offset constant used in `laxpcl1.f:2191` for 5th-byte rollover disambiguation.
const STBG_EPOCH_OFFSET_SEC: f64 = 1.8121823219164e8;
/// $256^4$ ticks (42,949.67296 s at 100 kHz).
const POW256_4: f64 = 4_294_967_296.0;

/// Classification of a decoded EA event.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum EventKind {
    /// Single-anode event (`ic2 == 0`).
    Single = 1,
    /// Double-anode event with Xe K-escape (`(k1 == 1 || k2 == 1) && ia1 < 8 && ia2 < 8`).
    Double = 2,
}

/// Raw unpacked fields from a 5-byte EA event word before single/double selection.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RawEaWord {
    /// Fine time offset in 10 us ticks (`0..=255`, already adjusted by `-2 mod 256`).
    pub fine_ticks: u8,
    /// Primary anode ID (`0..=15`).
    pub ia1: u8,
    /// Secondary anode ID (`0..=15`).
    pub ia2: u8,
    /// Primary 10-bit pulse height (`0..=1023`).
    pub ic1: u16,
    /// Secondary 10-bit pulse height (`0..=1023`).
    pub ic2: u16,
    /// Primary Xe K-X-ray flag (`0` or `1`).
    pub k1: u8,
    /// Secondary Xe K-X-ray flag (`0` or `1`).
    pub k2: u8,
}

impl RawEaWord {
    /// Unpack a 5-byte EA event record `[b0, b1, b2, b3, b4]` matching `laxpcl1.f:2230-2300`.
    #[inline]
    pub fn unpack(bytes: [u8; 5]) -> Self {
        let i5 = bytes[0];
        let i4 = bytes[1];
        let i3 = bytes[2] as u16;
        let i2 = bytes[3] as u16;
        let i1 = bytes[4];

        let fine_ticks = i4.wrapping_sub(2);
        let ia1 = i5 & 0x0F;
        let ia2 = (i5 >> 4) & 0x0F;
        let k1 = i1 & 0x01;
        let l1 = ((i1 >> 2) & 0x03) as u16;
        let k2 = (i1 >> 4) & 0x01;
        let l2 = ((i1 >> 6) & 0x03) as u16;
        let ic1 = (i3 << 2) | l1;
        let ic2 = (i2 << 2) | l2;

        Self {
            fine_ticks,
            ia1,
            ia2,
            ic1,
            ic2,
            k1,
            k2,
        }
    }

    /// Resolve single-event or Xe K-escape double-event `(pha, anode, kind)` matching
    /// `laxpcl1.f:2301-2314`. Returns `None` if the event is rejected (e.g. non-K double or veto double).
    #[inline]
    pub fn resolve(self) -> Option<(u16, u8, EventKind)> {
        if self.ic2 == 0 {
            Some((self.ic1, self.ia1, EventKind::Single))
        } else if (self.k1 == 1 || self.k2 == 1) && self.ia1 < 8 && self.ia2 < 8 {
            let pha = (self.ic1 + self.ic2).min(1023);
            let anode = if self.k1 == 1
                && (self.ic1 <= self.ic2 || (self.k2 == 0 && self.ic2 > 1))
            {
                self.ia2
            } else {
                self.ia1
            };
            Some((pha, anode, EventKind::Double))
        } else {
            None
        }
    }
}

/// A decoded LAXPC EA event in instrument time (before or after TCT UTC conversion).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EaEvent {
    /// Instrument timestamp (seconds since STBG switch-on, 10 us resolution).
    pub time_ins: f64,
    /// Raw 10-bit pulse height channel (`1..=999`).
    pub pha: u16,
    /// Assigned anode ID (`1..=10`).
    pub anode: u8,
    /// Detector unit ID (`1..=3`).
    pub unit_id: u8,
    /// Single or Xe K-escape double event.
    pub kind: EventKind,
}

/// Persistent cross-frame clock state for sequential `conver1` EA decoding within an orbit.
#[derive(Debug, Clone)]
pub struct EaDecoderState {
    pub tims: f64,
    pub timf_prev: f64,
    pub ifrp: i32,
    pub it5: i64,
    pub i2p: i32,
    pub i3p: i32,
    pub i4p: i32,
}

impl Default for EaDecoderState {
    fn default() -> Self {
        Self::new()
    }
}

impl EaDecoderState {
    /// Create a fresh decoder state for the start of an orbit file.
    pub fn new() -> Self {
        Self {
            tims: 0.0,
            timf_prev: 0.0,
            ifrp: -1,
            it5: 0,
            i2p: -1,
            i3p: -1,
            i4p: -1,
        }
    }

    /// Reset per-orbit coarse-time tracking (as `laxpcl1.f` does when `timf == 0`).
    pub fn reset_orbit(&mut self) {
        self.tims = 0.0;
        self.timf_prev = 0.0;
        self.ifrp = -1;
        self.i2p = -1;
        self.i3p = -1;
        self.i4p = -1;
    }
}

/// Result of decoding a single 2242-byte EA frame.
#[derive(Debug, Clone)]
pub struct DecodedEaFrame {
    /// L1 timestamp from frame header.
    pub l1_time: f64,
    /// Reconstructed 40-bit STBG header time (in seconds).
    pub header_time_ins: f64,
    /// 16-bit onboard frame counter.
    pub frame_counter: u16,
    /// Valid unpacked events from this frame (before cross-frame range/overlap filters).
    pub events: Vec<EaEvent>,
    /// Number of low/out-of-range/rejected words (`na0`) in this frame.
    pub anomalous_words: usize,
    /// Whether an `0xEE` end-of-frame marker was hit.
    pub hit_eof_marker: bool,
}

/// Decode a single 2242-byte EA frame updating `state`, bit-exact with `laxpcl1.f:conver1`.
pub fn decode_ea_frame(
    raw: &[u8],
    unit: LaxpcUnit,
    state: &mut EaDecoderState,
) -> Option<DecodedEaFrame> {
    if raw.len() < L1_FRAME_SIZE {
        return None;
    }
    if raw[44] != FRAME_SYNC_BYTE
        || raw[45] != unit.id()
        || raw[46] != MODE_EA_BYTE
        || raw[57] != 0
    {
        return None;
    }

    let l1_time = f64::from_be_bytes(raw[0..8].try_into().unwrap());
    let mut i5 = raw[49] as i64;
    let i4 = raw[50] as i64;
    let i3 = raw[51] as i64;
    let i2 = raw[52] as i64;
    let i1 = raw[53] as i64;

    if state.timf_prev == 0.0 {
        state.i2p = -1;
        state.i3p = -1;
        state.i4p = -1;
        state.tims = 0.0;
        state.ifrp = -1;
    }

    if i5 < state.it5 {
        i5 = state.it5;
    }

    let mut timf = ((((i5 * 256 + i4) * 256 + i3) * 256 + i2) * 256 + i1) as f64;
    let ifr1 = u16::from_be_bytes([raw[58], raw[59]]) as i32;

    // 5th time-byte rollover correction (`laxpcl1.f:2191-2199`, Dec 2024 fix)
    let dt = l1_time - timf / STBG_CLOCK_HZ - STBG_EPOCH_OFFSET_SEC;
    let d4 = dt * STBG_CLOCK_HZ / POW256_4;
    let n4 = if d4 > 0.0 {
        (d4 + 0.3) as i64
    } else if d4 < 0.0 {
        (d4 - 0.3) as i64
    } else {
        0
    };
    i5 += n4;
    timf += (n4 as f64) * POW256_4;
    state.timf_prev = timf;

    if state.tims == 0.0 {
        state.tims = timf - (i1 as f64);
    }
    if (timf - (i1 as f64)) > state.tims + 256.0 && ifr1 > state.ifrp + 1 {
        state.tims = timf - (i1 as f64);
        state.i2p = i2 as i32;
        state.i3p = i3 as i32;
        state.i4p = i4 as i32;
    }

    if i5 < 0 {
        return None;
    }
    state.ifrp = ifr1;
    state.it5 = i5;

    let mut ichk = 0;
    let mut na0: usize = 0;
    let mut hit_eof_marker = false;

    let mut tev_ticks: Vec<f64> = Vec::with_capacity(410);
    let mut phav: Vec<u16> = Vec::with_capacity(410);
    let mut anodev: Vec<u8> = Vec::with_capacity(410);
    let mut kindv: Vec<EventKind> = Vec::with_capacity(410);

    // Track coarse-time marker history for single-marker glitch filtering (`j2, j3, j4`)
    let mut j2_hist: Vec<i32> = Vec::with_capacity(64);
    let mut j3_hist: Vec<i32> = Vec::with_capacity(64);
    let mut j4_hist: Vec<i32> = Vec::with_capacity(64);

    // Loop over 5-byte words: Fortran `do i=61,2092,5` -> 0-based `pos = 60, 65, ..., 2090`
    let mut pos = 60usize;
    while pos <= 2091 {
        if raw[pos] == END_OF_FRAME_MARKER {
            hit_eof_marker = true;
            break;
        }
        if pos + 5 > raw.len() {
            break;
        }
        let w = [
            raw[pos],
            raw[pos + 1],
            raw[pos + 2],
            raw[pos + 3],
            raw[pos + 4],
        ];

        if w[0] == COARSE_TIME_MARKER {
            let mut m4 = w[1] as i32;
            let mut m3 = w[2] as i32;
            let mut m2 = w[3] as i32;

            // Carry propagation (`laxpcl1.f:2235-2247`)
            if m2 < state.i2p && m3 == state.i3p && m2 == 0 && state.i2p == 255 {
                m3 += 1;
            }
            if m3 > 255 {
                m3 = 0;
                m4 += 1;
            } else if m3 < state.i3p && m4 == state.i4p && m3 == 0 && state.i3p == 255 {
                m4 += 1;
            }
            if m4 > 255 {
                m4 = 0;
                state.it5 += 1;
            } else if m4 < state.i4p && m4 == 0 && state.i4p == 255 {
                state.it5 += 1;
            }

            // Single-marker glitch filter (`laxpcl1.f:2248-2259`)
            let nt = j2_hist.len(); // 0-based count before pushing (Fortran `nt - 1`)
            if nt >= 2 {
                let prev2_j2 = j2_hist[nt - 2];
                let prev2_j3 = j3_hist[nt - 2];
                let prev2_j4 = j4_hist[nt - 2];
                if m2 < state.i2p && m2 > 5 && state.i2p == prev2_j2 + 1 {
                    m2 = state.i2p + 1;
                }
                if m2 > state.i2p + 1 && m2 > 5 && state.i2p == prev2_j2 + 1 {
                    m2 = state.i2p + 1;
                }
                if m3 < state.i3p && m3 > 0 && state.i3p == prev2_j3 {
                    m3 = state.i3p;
                }
                if m3 > state.i3p && m2 > 5 && state.i3p == prev2_j3 {
                    m3 = state.i3p;
                }
                if m4 < state.i4p && m4 > 0 && state.i4p == prev2_j4 {
                    m4 = state.i4p;
                }
                if m4 > state.i4p && m3 > 5 && state.i4p == prev2_j4 {
                    m4 = state.i4p;
                }
            }
            j2_hist.push(m2);
            j3_hist.push(m3);
            j4_hist.push(m4);

            state.tims = ((((state.it5 * 256 + m4 as i64) * 256 + m3 as i64) * 256 + m2 as i64)
                * 256) as f64;

            // Retroactively align events emitted before the first 0xEF marker in this frame
            if ichk == 0
                && ((!tev_ticks.is_empty() && timf > state.tims && state.i2p == -1)
                    || m2 == state.i2p)
            {
                for t in &mut tev_ticks {
                    *t -= 256.0;
                }
            } else if ichk == 0 && m2 != state.i2p + 1 && m2 > 0 && state.i2p != -1 {
                let i2c = (state.i2p - m2 + 1) as f64;
                let i3c = (state.i3p - m3) as f64;
                let i4c = if m3 > 0 { (state.i4p - m4) as f64 } else { 0.0 };
                let shift = 256.0 * i2c + 65536.0 * i3c + 16_777_216.0 * i4c;
                for t in &mut tev_ticks {
                    *t -= shift;
                }
            }

            state.i2p = m2;
            state.i3p = m3;
            state.i4p = m4;
            ichk = 1;
        } else {
            let word = RawEaWord::unpack(w);
            let t_tick = state.tims + (word.fine_ticks as f64);
            let (pha, anode, kind) = match word.resolve() {
                Some(v) => v,
                None => (0, 0, EventKind::Single),
            };
            if !(10..=1000).contains(&pha) || anode == 0 {
                na0 += 1;
            }
            tev_ticks.push(t_tick);
            phav.push(pha);
            anodev.push(anode);
            kindv.push(kind);
        }

        pos += 5;
    }

    // Frame-end abort/truncation check (`laxpcl1.f:2330-2338`)
    let mut nad = tev_ticks.len();
    if (na0 > 4 || hit_eof_marker) && nad > 0 {
        let mut cut = nad; // 1-based `j = nad + 1` -> `j - 2 = nad - 1`
        for j in 1..nad {
            if tev_ticks[j] <= tev_ticks[j - 1] || tev_ticks[j] > tev_ticks[j - 1] + 1.0e5 {
                cut = j; // 1-based `j+1`, so `(j+1) - 2 = j - 1`
                break;
            }
        }
        nad = cut.saturating_sub(1);
    }

    let mut events = Vec::with_capacity(nad);
    for idx in 0..nad {
        events.push(EaEvent {
            time_ins: tev_ticks[idx] / STBG_CLOCK_HZ,
            pha: phav[idx],
            anode: anodev[idx],
            unit_id: unit.id(),
            kind: kindv[idx],
        });
    }

    Some(DecodedEaFrame {
        l1_time,
        header_time_ins: timf / STBG_CLOCK_HZ,
        frame_counter: ifr1 as u16,
        events,
        anomalous_words: na0,
        hit_eof_marker,
    })
}

/// Output of decoding an entire stream of ordered EA frames across one or more orbits.
#[derive(Debug, Clone)]
pub struct EaStreamResult {
    /// All valid events passing basic validity (`time_ins > 1.0`, `1 <= pha <= 999`, `1 <= anode <= 10`).
    pub events: Vec<EaEvent>,
    /// Internal frame-boundary time gaps (`>= 50 ms`) in instrument time `(start_ins, stop_ins)`.
    pub frame_gaps: Vec<(f64, f64)>,
    /// Large time gaps (`>= 10.0 s`) in instrument time `(start_ins, stop_ins)`.
    pub large_gaps: Vec<(f64, f64)>,
    /// Number of non-monotonic event timestamp steps encountered.
    pub non_monotonic_steps: usize,
    /// Total valid events decoded before user anode/channel/single filtering.
    pub total_valid_events: usize,
    /// Minimum and maximum instrument timestamps seen.
    pub t_ins_min: f64,
    pub t_ins_max: f64,
}

/// Decode a multi-orbit sequence of ordered EA frame batches (one slice per orbit),
/// applying cross-orbit overlap filtering (`tim <= t00`) and internal frame-gap detection
/// matching `laxpcl1.f:513-650`.
pub fn decode_ea_orbits(orbit_frames: &[Vec<OrderedFrame>], unit: LaxpcUnit) -> EaStreamResult {
    let mut events: Vec<EaEvent> = Vec::new();
    let mut frame_gaps: Vec<(f64, f64)> = Vec::new();
    let mut large_gaps: Vec<(f64, f64)> = Vec::new();
    let mut non_monotonic_steps: usize = 0;

    let mut t00 = 0.0_f64;
    let mut tmax = 0.0_f64;
    let mut t_ins_min = f64::INFINITY;
    let mut prev_valid_tim: Option<f64> = None;

    let mut state = EaDecoderState::new();

    for frames in orbit_frames {
        state.reset_orbit();
        for frame in frames {
            let decoded = match decode_ea_frame(&frame.data, unit, &mut state) {
                Some(d) => d,
                None => continue,
            };

            for (k, ev) in decoded.events.into_iter().enumerate() {
                let tim = ev.time_ins;
                if tim <= 1.0 {
                    continue;
                }

                // Cross-orbit overlap / gap-fill guard (`laxpcl1.f:535-545`)
                let mut filled_gap = false;
                if tim <= t00 {
                    if non_monotonic_steps == 0 {
                        for gap in &mut large_gaps {
                            if tim > gap.0 && tim < gap.1 {
                                gap.0 = tim;
                                filled_gap = true;
                                break;
                            }
                        }
                    }
                    if !filled_gap {
                        continue;
                    }
                }

                // Validity check (`laxpcl1.f:549-550`)
                if ev.pha == 0 || ev.pha >= 1000 || ev.anode == 0 || ev.anode > 10 {
                    continue;
                }

                if tim < t_ins_min {
                    t_ins_min = tim;
                }
                if tim > tmax {
                    tmax = tim;
                }

                if let Some(prev_t) = prev_valid_tim {
                    if tim <= prev_t {
                        non_monotonic_steps += 1;
                    }
                    if tim >= prev_t + 10.0 && !filled_gap {
                        let last_end = large_gaps.last().map(|g| g.1).unwrap_or(-1.0);
                        if prev_t > last_end {
                            large_gaps.push((prev_t, tim));
                        }
                    }
                    if k == 0 && tim >= prev_t + 5.0e-2 {
                        let last_end = frame_gaps.last().map(|g| g.1).unwrap_or(-1.0);
                        if prev_t > last_end {
                            frame_gaps.push((prev_t, tim));
                        }
                    }
                }

                prev_valid_tim = Some(tim);
                events.push(ev);
            }
        }
        if tmax > t00 {
            t00 = tmax;
        }
    }

    let total_valid_events = events.len();
    EaStreamResult {
        events,
        frame_gaps,
        large_gaps,
        non_monotonic_steps,
        total_valid_events,
        t_ins_min: if t_ins_min.is_finite() { t_ins_min } else { 0.0 },
        t_ins_max: tmax,
    }
}

/// Check whether an `EaEvent` passes a user `EventFilter` (anode, channel range, single/double).
#[inline]
pub fn passes_event_filter(ev: &EaEvent, filter: &EventFilter) -> bool {
    if filter.single_only && ev.kind != EventKind::Single {
        return false;
    }
    if !filter.anodes.matches(ev.anode) {
        return false;
    }
    // Note: in `laxpcl1.f:363-364,614-615`, icl/icu are 1-based shifted by +1 from 0..=1023,
    // while `ic(3)` is 1..=999; using inclusive `[pha_min, pha_max]` directly on 0..=1023 PHA.
    ev.pha >= filter.pha_min && ev.pha <= filter.pha_max
}
