//! Stage 2c: `cztdatasel` — Filter CZTI events against `COMMON` or `QUAD` GTI tables.

use rayon::prelude::*;

use crate::gtigen::normalize_intervals;
use crate::types::QuadEvents;

/// Filter a single quadrant's events to those falling inside `gti_intervals` (`tstart <= t <= tstop`).
///
/// Updates `tstart`, `tstop`, and `exposure` (sum of GTI interval lengths) on the returned `QuadEvents`.
pub fn datasel_quad(events: &QuadEvents, gti_intervals: &[(f64, f64)]) -> QuadEvents {
    let gti = normalize_intervals(gti_intervals);
    let exposure: f64 = gti.iter().map(|&(a, b)| (b - a).max(0.0)).sum();

    if events.is_empty() || gti.is_empty() {
        let mut empty = QuadEvents::new(events.quad_id);
        empty.exposure = exposure;
        if let Some(&(a, _)) = gti.first() {
            empty.tstart = a;
        }
        if let Some(&(_, b)) = gti.last() {
            empty.tstop = b;
        }
        return empty;
    }

    // Parallel binary-search membership check
    let keep: Vec<bool> = events
        .time
        .par_iter()
        .map(|&t| {
            let idx = gti.partition_point(|&(_, stop)| stop < t);
            idx < gti.len() && t >= gti[idx].0 && t <= gti[idx].1
        })
        .collect();

    let mut out = events.filter_by_mask(&keep);
    out.exposure = exposure;
    out.tstart = gti[0].0;
    out.tstop = gti[gti.len() - 1].1;
    out
}
