//! Pure-Rust, `rayon`-parallelized AstroSat LAXPC Level-1 to Level-2 scientific pipeline.
//!
//! Covers all three LAXPC proportional counter units (LAXPC10 / `lxp1`, LAXPC20 / `lxp2`,
//! LAXPC30 / `lxp3`) with bit-exact parity against TIFR `LaxpcSoft v3.4.6` (`laxpcl1.f` &
//! `backshiftv3.f`) and zero CFITSIO/HEASoft/Fortran build dependencies:
//!
//! - [`frames`]: 2242-byte Level-1 EA/BB frame parsing, integrity checks, sorting, deduplication, and gap diagnostics.
//! - [`ea`]: 5-byte EA event unpacking (10 us STBG clock, `0xEF` coarse markers, `0xEE` EOF) and Xe K-escape double-event handling.
//! - [`bb`]: 31 x 64-byte BBC subframe decoding (`LLD`, `JAN`, `LDU`, `IGE`), 32.768 s binning, and ULD quartile rates.
//! - [`gti`]: TCT conversion (via `astrosat-core::tct::LaxpcTct`) and standard MKF + BBC GTI generation (ELV, SAA, ULD anomalies, gap merge/drop, frame-gap exposure correction).
//! - [`products`]: Non-paralyzable deadtime correction (43 us, 1-s running mean for `tbin < 0.2 s`, 23000 cps singularity guard) and OGIP Format B + Format A FITS generation.
//! - [`background`]: 2D Earth lat/lon B-spline + 5-harmonic diurnal background model, `gaina8.dat`/`lxhv.dat` gain alignment, spectral/light-curve shifting, and automatic RMF recommendation.
//! - [`pipeline`]: High-level single-unit and parallel multi-unit orchestrator (`run_pipeline`, `run_unit_pipeline`, `discover_l1_orbits`).

pub mod background;
pub mod bb;
pub mod ea;
pub mod error;
pub mod frames;
pub mod gti;
pub mod pipeline;
pub mod products;
pub mod types;

pub use background::{
    apply_lightcurve_gain_shift, background_lightcurve_shift_fraction, compute_gain_shift,
    estimate_lxp30_density, eval_bspline_2d, eval_bspline_basis, recommend_rmf,
    secular_background_scale, shift_background_spectrum, BackgroundCatalog, BackgroundCatalogEntry,
    BackgroundModel2D, BackgroundPha, ChannelMap, GainEntry, GainShiftSolution, GainTable,
    HvJumpTable, DEFAULT_DIURNAL_PERIOD_SEC,
};
pub use bb::{decode_bb_frame, BbcRecord, BbcSeries, UldQuartileSummary};
pub use ea::{
    decode_ea_frame, decode_ea_orbits, passes_event_filter, DecodedEaFrame, EaDecoderState,
    EaEvent, EaStreamResult, EventKind, RawEaWord,
};
pub use error::{LaxpcError, Result};
pub use frames::{
    merge_multi_orbit_frames, sort_and_dedup_frames, FrameDiagnostics, FrameGap, FrameHeader,
    FrameMode, OrderedFrame, COARSE_TIME_MARKER, END_OF_FRAME_MARKER, FRAME_SYNC_BYTE,
    MODE_BB_BYTE, MODE_EA_BYTE,
};
pub use gti::{
    generate_standard_gti, interp_lon_deg, load_laxpc_tct_files, load_merged_mkf_files,
    lxp10_uld_anomaly_threshold, utc_to_ins, GtiConfig, GtiResult,
};
pub use pipeline::{
    discover_l1_orbits, read_and_sort_l1_fits_frames, run_pipeline, run_unit_pipeline,
    write_orbit_from_mkf, MultiUnitPipelineResult, OrbitInputFiles, PipelineConfig,
    UnitPipelineResult,
};
pub use products::{
    bin_events_raw, build_calibrated_events, build_lightcurve, build_spectrum,
    compute_bin_gap_fractions, compute_gti_frame_gap_loss, deadtime_correct_rate, deadtime_factor,
    smooth_all_counts_for_deadtime, write_ascii_lc, write_ascii_pha, write_format_a_event_fits,
    write_format_a_event_fits_with_channels, write_format_b_event_fits, write_gti_fits,
    write_lightcurve_fits, write_orbit_fits, write_spectrum_fits, CalibratedEvent, LaxpcLightCurve,
    LaxpcSpectrum, RawBinnedCounts,
};
pub use types::{
    utc_to_mission_day, AnodeSelection, EventFilter, LaxpcUnit, BBC_BIN_SEC,
    DEADTIME_RATE_GUARD_CPS, DEFAULT_DEADTIME_SEC, L1_FRAME_SIZE, STBG_CLOCK_HZ,
};
