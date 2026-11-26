//! Stage 2f: `cztevtclean` & `cztflagbadpix` — Alpha/Veto event selection and bad-pixel map union.

use serde::{Deserialize, Serialize};

use crate::types::{BadpixTable, QuadEvents};

/// Configuration for `cztevtclean` (default: `alpha == 0`, `veto in 0..=0`).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EvtCleanConfig {
    /// Allowed alpha values (default: `[0]`, i.e., non-alpha X-ray events).
    pub allowed_alpha: Vec<u8>,
    /// Inclusive veto range `(min_veto, max_veto)` (default: `(0, 0)`).
    pub veto_range: (u16, u16),
    /// Whether input contains double-event pairs (if true, both members of a pair must pass).
    pub is_double_event: bool,
}

impl Default for EvtCleanConfig {
    fn default() -> Self {
        Self {
            allowed_alpha: vec![0],
            veto_range: (0, 0),
            is_double_event: false,
        }
    }
}

impl EvtCleanConfig {
    /// Parse pipeline CLI strings like `alphaval="0"` and `vetorange="0-0"`.
    pub fn from_strings(alphaval: &str, vetorange: &str, is_double_event: bool) -> Self {
        let allowed_alpha: Vec<u8> = alphaval
            .split(',')
            .filter_map(|s| s.trim().parse::<u8>().ok())
            .collect();
        let veto_parts: Vec<u16> = vetorange
            .split('-')
            .filter_map(|s| s.trim().parse::<u16>().ok())
            .collect();
        let veto_range = match veto_parts.as_slice() {
            [a, b] => (*a, *b),
            [a] => (*a, *a),
            _ => (0, 0),
        };
        Self {
            allowed_alpha: if allowed_alpha.is_empty() {
                vec![0]
            } else {
                allowed_alpha
            },
            veto_range,
            is_double_event,
        }
    }
}

/// Filter events in a single quadrant by `alpha` and `veto` criteria (`cztevtclean.cpp`).
pub fn evtclean_quad(events: &QuadEvents, cfg: &EvtCleanConfig) -> QuadEvents {
    let n = events.len();
    if n == 0 {
        return events.clone();
    }

    let passes = |i: usize| -> bool {
        let a_ok = cfg.allowed_alpha.contains(&events.alpha[i]);
        let v = events.veto[i];
        let v_ok = v >= cfg.veto_range.0 && v <= cfg.veto_range.1;
        a_ok && v_ok
    };

    let mut keep = vec![false; n];
    if cfg.is_double_event {
        let mut i = 0usize;
        while i + 1 < n {
            if passes(i) && passes(i + 1) {
                keep[i] = true;
                keep[i + 1] = true;
            }
            i += 2;
        }
    } else {
        for (i, slot) in keep.iter_mut().enumerate() {
            *slot = passes(i);
        }
    }

    events.filter_by_mask(&keep)
}

/// `cztflagbadpix`: merge one or more per-quadrant `BadpixTable`s by taking the maximum (worst)
/// pixel quality flag at every `(quad, detid, pixid)`.
pub fn flag_badpix_union(tables: &[BadpixTable]) -> BadpixTable {
    BadpixTable::union_all(tables)
}
