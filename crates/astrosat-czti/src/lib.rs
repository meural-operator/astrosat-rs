//! Pure-Rust, `rayon`-parallelized AstroSat CZTI Level-1 -> Level-2 scientific pipeline.
//!
//! Faithfully implements all 3 stages of ISRO AstroSat CZTI Pipeline v3.0.1 with zero C/CFITSIO/HEASoft
//! build dependencies:
//! - **Stage 1**: [`science2event`], [`bunchclean`], [`pha2energy`]
//! - **Stage 2**: [`gtigen`], [`gaas`], [`datasel`], [`pixclean`], [`noiseclean`], [`evtencor`], [`evtclean`]
//! - **Stage 3**: [`dpigen`], [`image`], [`bindata`]
//! - **Orchestrator & FITS/CALDB I/O**: [`pipeline::run_pipeline`], [`fits_io`]

pub mod bindata;
pub mod bunchclean;
pub mod datasel;
pub mod dpigen;
pub mod error;
pub mod evtclean;
pub mod evtencor;
pub mod fits_io;
pub mod gaas;
pub mod gtigen;
pub mod image;
pub mod noiseclean;
pub mod pha2energy;
pub mod pipeline;
pub mod pixclean;
pub mod science2event;
pub mod stages;
pub mod types;

pub use bindata::{
    absco_ta, apply_additional_lld, bindata_quad, calculate_renormalized_weights,
    combine_lightcurves, combine_responses, combine_spectra, compute_5x5_avg_open_frac,
    compute_spectrum_quality_and_grouping, generate_quad_response, quat_to_rpy,
    radec_to_czti_tx_ty, raytrace_quad_shadow, raytrace_quad_shadow_simple, roll_rot_to_rpy,
    ta_mask_transmission, validate_fov, BinDataConfig, BinDataOutputType, CztiLightCurve,
    CztiResponse, CztiSpectrum, QuadRspCaldb, RenormalizedMaskWeights, ADDITIONAL_LLD,
};
pub use bunchclean::{
    apply_gti_to_livetime, bunchclean_quad, re_update_livetime, update_livetime, BunchCleanConfig,
    QuadBunchCleanResult,
};
pub use datasel::datasel_quad;
pub use dpigen::{generate_dph_dpi, generate_quad_dph_dpi, DphDpiProduct, DpiGenConfig, QuadEffArea};
pub use error::{CztiError, Result};
pub use evtclean::{evtclean_quad, flag_badpix_union, EvtCleanConfig};
pub use evtencor::{
    evtencor_quad, secular_gain_drift_factor, EvtEnCorConfig, QuadEvtEnCorResult, QuadGainCorCaldb,
    WichmannHillRng,
};
pub use fits_io::{CztiCaldb, CztiObsMeta};
pub use gaas::{
    calculate_ra_dec_twist_deg, compute_quad_aspect, det_to_inertial, interpolate_rpy,
    RpySampleRad, TeldefMatrix,
};
pub use gtigen::{
    default_mkf_thresholds, drop_short_intervals, generate_czti_gtis, intersect_intervals,
    intervals_from_mkf_mask, normalize_intervals, parse_mkf_threshold_text, parse_user_gti_text,
    subtract_intervals, CztiGtiOutput, GtiGenConfig, MkfThresholdRule,
};
pub use image::{
    balanced_cross_correlation, balanced_cross_correlation_from_mask, combine_sky_images,
    fft_cross_correlation, fft_cross_correlation_with_shift, CztiWcs, ImageAlgorithm, ImageConfig,
    SkyImage,
};
pub use noiseclean::{
    eventsep_quad, flickpixclean_quad, flickpixclean_quad_with_gti, heavybunchclean_quad,
    heavybunchclean_quad_with_exposure, merge_split_bunches, noisypixclean_quad, pix_threshold_fix,
    superbunchclean_quad, FlickPixCleanResult, NoiseCleanConfig, NoisyPixCleanResult,
    SuperBunchCleanResult,
};
pub use pha2energy::{
    pha2energy_quad, temp_to_bin, EboundsTable, QuadGainOffset, QuadTempSeries,
};
pub use pipeline::{run_pipeline, run_science2event, CztiPipelineConfig, CztiPipelineSummary};
pub use stages::{
    read_aspect_file, read_bunch_file, read_livetime_file, run_addspec, run_bindata,
    run_bunchclean, run_datasel, run_dpigen, run_evtclean, run_evtencor, run_gaas, run_gtigen,
    run_image, run_noiseclean, run_pha2energy, run_pixclean,
};
pub use pixclean::{find_bad_pixels, pixclean_quad, PixCleanConfig, QuadPixCleanResult};
pub use science2event::{
    assemble_frames, clean_junk_temperatures, decode_packet, decode_quad_packets,
    packet_bytes_to_words, DecodedFrame, DecodedPacket, DecodedQuadStream, FebHeader,
    PacketFrameHeader, Science2EventConfig, PACKET_SIZE_BYTES, PACKET_WORDS, SYNC1, SYNC2,
};
pub use types::{
    detid_pixid_to_detxy, detxy_to_detid_pixid, detxy_to_locxy, pix_flag, BadpixTable, BunchTable,
    EventSelectionAlgo, FrameHeaderTable, GtiType, QuadAspect, QuadEvents, QuadLivetime,
    SsmDataTable, TempTable, VetoSpectrumTable, CZTSPEC_SIZE, CZT_TICK_HZ, FULL_DIM,
    GAIN_DRIFT_EPOCH_SEC, MASK_HEIGHT_MM, MAX_FOV_DEG, NUM_DET_PER_QUAD, NUM_PIX_PER_DET,
    NUM_QUADS, NUM_TEMP_BINS, PIXELS_PER_QUAD, PI_CHANNELS, QUAD_DIM, VETOSPEC_SIZE, XLDET_MM,
};
