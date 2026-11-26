//! Complete `clap` (derive) CLI definition for the unified `astrosat` binary.

use std::path::PathBuf;

use clap::parser::ValueSource;
use clap::{ArgMatches, Args, Parser, Subcommand, ValueEnum};

/// Walk down nested subcommands to the leaf `ArgMatches` for value-source inspection.
pub fn leaf_matches(mut m: &ArgMatches) -> &ArgMatches {
    while let Some((_name, sub)) = m.subcommand() {
        m = sub;
    }
    m
}

/// True iff the flag/option `id` was explicitly supplied on the command line.
#[inline]
pub fn is_cli_arg(m: &ArgMatches, id: &str) -> bool {
    m.value_source(id) == Some(ValueSource::CommandLine)
}

/// Resolve a boolean option with symmetric CLI `--pos`/`--neg` override over config & default.
#[inline]
pub fn resolve_flag(
    m: &ArgMatches,
    pos_id: &str,
    neg_id: &str,
    cfg_val: Option<bool>,
    default: bool,
) -> bool {
    if is_cli_arg(m, pos_id) {
        true
    } else if is_cli_arg(m, neg_id) {
        false
    } else {
        cfg_val.unwrap_or(default)
    }
}

#[derive(Debug, Parser)]
#[command(
    name = "astrosat",
    version,
    about = "Pure-Rust, zero-patch AstroSat CZTI & LAXPC Level-1 -> Level-2 scientific pipeline",
    long_about = "astrosat-rs: High-throughput pure-Rust AstroSat CZTI & LAXPC pipeline, \
                  native Solar System barycentric timing, FITS 4.0/OGIP inspector & validator, \
                  and pulsar periodicity analysis with zero CFITSIO/HEASoft/glog dependencies."
)]
pub struct Cli {
    /// Number of Rayon worker threads (defaults to all available CPU cores).
    #[arg(long, global = true, value_name = "N")]
    pub threads: Option<usize>,

    /// Emit machine-readable JSON on stdout (diagnostics go to stderr).
    #[arg(long, global = true)]
    pub json: bool,

    /// Optional TOML or JSON configuration file.
    #[arg(long, global = true, value_name = "PATH")]
    pub config: Option<PathBuf>,

    /// Enable verbose debug logging.
    #[arg(short, long, global = true)]
    pub verbose: bool,

    /// Suppress non-essential human output.
    #[arg(short, long, global = true)]
    pub quiet: bool,

    #[command(subcommand)]
    pub command: TopCommand,
}

#[derive(Debug, Subcommand)]
pub enum TopCommand {
    /// Run AstroSat CZTI Level-1 -> Level-2 pipeline or individual stages.
    Czti {
        #[command(subcommand)]
        command: CztiCommand,
    },
    /// Run AstroSat LAXPC Level-1 -> Level-2 pipeline or individual stages.
    Laxpc {
        #[command(subcommand)]
        command: LaxpcCommand,
    },
    /// Apply Solar System barycentric correction (works natively on CZTI, LAXPC Format A & B).
    Bary(BaryArgs),
    /// Fast human/JSON summary of FITS HDUs, rows, timing keywords, exposure, and health.
    Inspect(InspectArgs),
    /// Built-in FITS 4.0 + CHECKSUM/DATASUM + OGIP schema validation.
    Validate(ValidateArgs),
    /// Pulsar epoch-folding and Z^2_n periodicity test on a barycentered event FITS file.
    Timing(TimingArgs),
}

// ============================================================================
// Inspect, Validate, Timing, Bary
// ============================================================================

#[derive(Debug, Args)]
pub struct InspectArgs {
    /// Input FITS file(s) to inspect.
    #[arg(required = true, value_name = "FITS")]
    pub files: Vec<PathBuf>,

    /// Scan event/GTI columns for time monotonicity and NaN/Inf values.
    #[arg(long)]
    pub full: bool,
}

#[derive(Debug, Args)]
pub struct ValidateArgs {
    /// Input FITS file(s) to validate.
    #[arg(required = true, value_name = "FITS")]
    pub files: Vec<PathBuf>,

    /// Treat warnings (e.g. missing CHECKSUM or EXTNAME) as validation failures.
    #[arg(long)]
    pub strict: bool,
}

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum FrameArg {
    Icrs,
    Fk5,
}

#[derive(Debug, Args)]
pub struct BaryArgs {
    /// Input Level-2 event or light-curve FITS file (CZTI, LAXPC Format A, or Format B).
    #[arg(value_name = "INPUT_FITS")]
    pub input: Option<PathBuf>,

    /// Input FITS file (alternative flag form).
    #[arg(short = 'i', long = "input", value_name = "FITS")]
    pub input_flag: Option<PathBuf>,

    /// Output barycentered FITS file path.
    #[arg(short = 'o', long = "output", required = true, value_name = "OUTPUT_FITS")]
    pub output: PathBuf,

    /// AstroSat orbit FITS file (`.orb` / `AstroSat.orb`).
    #[arg(long, conflicts_with = "mkf", value_name = "ORB_FITS")]
    pub orbit: Option<PathBuf>,

    /// Optional MKF file to synthesize an in-memory orbit if `--orbit` is omitted.
    #[arg(long, value_name = "MKF_FITS")]
    pub mkf: Option<PathBuf>,

    /// Target Right Ascension in degrees (defaults to `RA_OBJ`/`RA_NOM`/`RA_PNT` from input FITS).
    #[arg(long, allow_hyphen_values = true, value_name = "DEG")]
    pub ra: Option<f64>,

    /// Target Declination in degrees (defaults to `DEC_OBJ`/`DEC_NOM`/`DEC_PNT` from input FITS).
    #[arg(long, allow_hyphen_values = true, value_name = "DEG")]
    pub dec: Option<f64>,

    /// Coordinate reference frame (`icrs` or `fk5`).
    #[arg(long, value_enum, default_value_t = FrameArg::Icrs)]
    pub frame: FrameArg,

    /// Ephemeris backend: `analytical` (built-in) or path to a JPL DE200/DE405/DE421/DE430 binary or HEASARC FITS (`JPLEPH.*`) file.
    #[arg(long, default_value = "analytical", value_name = "EPHEM")]
    pub ephem: String,

    /// Recompute CHECKSUM/DATASUM headers on output FITS (overrides config `write_checksums = false`).
    #[arg(long, conflicts_with = "no_checksum")]
    pub checksum: bool,

    /// Skip recomputing CHECKSUM/DATASUM headers on output FITS.
    #[arg(long)]
    pub no_checksum: bool,
}

#[derive(Debug, Args)]
pub struct TimingArgs {
    /// Input barycentered event FITS file.
    #[arg(value_name = "BARY_EVT_FITS")]
    pub input: PathBuf,

    /// HDU name containing event timestamps (auto-detects `EVENTS`, `LAXPC_MODE_EA`, `Q0..Q3`).
    #[arg(long)]
    pub hdu: Option<String>,

    /// Optional GTI FITS file to filter event timestamps.
    #[arg(long, value_name = "GTI_FITS")]
    pub gti: Option<PathBuf>,

    /// Spin frequency in Hz.
    #[arg(long, value_name = "HZ")]
    pub freq: Option<f64>,

    /// Spin period in seconds (alternative to `--freq`).
    #[arg(long, value_name = "SEC")]
    pub period: Option<f64>,

    /// Frequency derivative df/dt in Hz/s.
    #[arg(long, default_value_t = 0.0, allow_hyphen_values = true)]
    pub fdot: f64,

    /// Second frequency derivative d2f/dt2 in Hz/s^2.
    #[arg(long, default_value_t = 0.0, allow_hyphen_values = true)]
    pub fddot: f64,

    /// Reference epoch t0 in same time units as `TIME` (defaults to midpoint of data span).
    #[arg(long)]
    pub epoch: Option<f64>,

    /// Minimum search frequency in Hz (used with `--freq-max`).
    #[arg(long)]
    pub freq_min: Option<f64>,

    /// Maximum search frequency in Hz (used with `--freq-min`).
    #[arg(long)]
    pub freq_max: Option<f64>,

    /// Half-width (+/- Hz) to search around `--freq` / `--period`.
    #[arg(long, default_value_t = 0.0)]
    pub search_span: f64,

    /// Frequency step df in Hz (defaults to 1 / (5 * T_span)).
    #[arg(long)]
    pub df: Option<f64>,

    /// Explicit number of frequency trials.
    #[arg(long)]
    pub n_trials: Option<usize>,

    /// Number of harmonics `n` for Z^2_n statistic.
    #[arg(long, default_value_t = 2)]
    pub harmonics: usize,

    /// Number of phase bins for epoch-folding chi-squared and pulse profile.
    #[arg(long, default_value_t = 32)]
    pub nbins: usize,

    /// Minimum PI/Channel value to include.
    #[arg(long)]
    pub pi_min: Option<i32>,

    /// Maximum PI/Channel value to include.
    #[arg(long)]
    pub pi_max: Option<i32>,

    /// Filter by LAXPC detector unit ID (1, 2, or 3) when Format A `LAXPC_No` is present.
    #[arg(long)]
    pub unit: Option<u8>,

    /// Filter by LAXPC anode layer when `Layer` column is present.
    #[arg(long)]
    pub layer: Option<u8>,
}

// ============================================================================
// CZTI Subcommands
// ============================================================================

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum CztiGtiTypeArg {
    Quad,
    Common,
}

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum CztiAlgoArg {
    Default,
    Generalized,
}

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum CztiImageAlgoArg {
    Balcc,
    Fft,
}

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum CztiOutputTypeArg {
    Spec,
    Lc,
    Both,
}

#[derive(Debug, Subcommand)]
pub enum CztiCommand {
    /// Run the full end-to-end AstroSat CZTI Level-1 -> Level-2 pipeline.
    Run(CztiRunArgs),
    /// Stage 1a: Decode 2048-byte Level-1 modeM0 packets into raw events, bunches, and headers.
    Decode(CztiDecodeArgs),
    /// Stage 1b: Identify and clean cosmic-ray bunches and compute livetime/pixel exposure.
    Bunchclean(CztiBunchcleanArgs),
    /// Stage 1c: Convert event PHA to calibrated Energy (keV) and 512-channel PI.
    Pha2energy(CztiPha2energyArgs),
    /// Stage 2a: Generate COMMON and per-quadrant (Q0_GTI..Q3_GTI) GTIs from M0 + MKF.
    Gtigen(CztiGtigenArgs),
    /// Stage 2b: Compute detector pointing aspect (RA, DEC, TWIST) per quadrant.
    Gaas(CztiGaasArgs),
    /// Stage 2c: Filter events against COMMON or QUAD GTIs.
    Datasel(CztiDataselArgs),
    /// Stage 2d: Bad/noisy/flickering pixel cleaning (`pixclean` or generalized `noiseclean`).
    Clean(CztiCleanArgs),
    /// Stage 2e: Apply in-flight gain correction, secular drift, and LLD/ULD cuts.
    Evtencor(CztiEvtencorArgs),
    /// Stage 2f: Filter events by alpha/veto flags and bad-pixel map.
    Evtclean(CztiEvtcleanArgs),
    /// Stage 3a: Generate 128x128 Detector Plane Histogram (DPH) and effective-area DPI.
    Dpigen(CztiDpigenArgs),
    /// Stage 3b: Reconstruct coded-mask sky images (`balcc` or `fft`) with WCS.
    Image(CztiImageArgs),
    /// Stage 3c: Coded-mask ray-traced spectra (.pha), light curves (.lc), responses (.rsp), and weighted events.
    Bindata(CztiBindataArgs),
    /// Co-add multiple CZTI OGIP `.pha` spectra and (optionally) exposure-weighted `.rsp` responses (`cztaddspec`).
    Addspec(CztiAddspecArgs),
}

#[derive(Debug, Args)]
pub struct CztiRunArgs {
    /// Level-1 CZTI observation directory (auto-discovers M0, SS, TCT, and MKF files).
    #[arg(short = 'i', long = "input", value_name = "L1_DIR")]
    pub input: Option<PathBuf>,
    /// Explicit Level-1 modeM0 FITS file (`*cztM0_level1.fits`).
    #[arg(long, value_name = "FITS")]
    pub m0: Option<PathBuf>,
    /// Explicit Level-1 modeSS FITS file (`*cztSS_level1.fits`) for temperature telemetry.
    #[arg(long, value_name = "FITS")]
    pub ss: Option<PathBuf>,
    /// Explicit Level-1 TCT file (`*.tct`).
    #[arg(long, value_name = "TCT")]
    pub tct: Option<PathBuf>,
    /// Explicit MKF housekeeping file (`*.mkf`).
    #[arg(long, value_name = "MKF")]
    pub mkf: Option<PathBuf>,
    /// CZTI CALDB directory (containing `data/as1/czti/caldb.indx` or `caldb.indx`).
    #[arg(long, value_name = "CALDB_DIR")]
    pub caldb: Option<PathBuf>,
    /// Output directory for Level-2 FITS products.
    #[arg(short = 'o', long = "output", value_name = "OUT_DIR")]
    pub output: Option<PathBuf>,
    /// Output filename prefix (auto-derived from modeM0 filename if omitted).
    #[arg(long)]
    pub prefix: Option<String>,
    /// Optional `mkfThresholds.txt` override.
    #[arg(long, value_name = "TXT")]
    pub mkf_thresholds: Option<PathBuf>,
    /// Optional user ASCII GTI file.
    #[arg(long, value_name = "TXT")]
    pub user_gti: Option<PathBuf>,
    /// Pipeline entry stage (1, 2, or 3).
    #[arg(long, default_value_t = 1)]
    pub entry_stage: u8,
    /// Pipeline exit stage (1, 2, or 3).
    #[arg(long, default_value_t = 3)]
    pub exit_stage: u8,
    /// GTI mode (`quad` or `common`).
    #[arg(long, value_enum)]
    pub gti_type: Option<CztiGtiTypeArg>,
    /// Cleaning algorithm (`default` or `generalized`).
    #[arg(long, value_enum)]
    pub algo: Option<CztiAlgoArg>,
    /// Source RA in degrees for bindata (defaults to MKF `RA_PNT` / event header).
    #[arg(long, allow_hyphen_values = true)]
    pub ra: Option<f64>,
    /// Source DEC in degrees for bindata (defaults to MKF `DEC_PNT` / event header).
    #[arg(long, allow_hyphen_values = true)]
    pub dec: Option<f64>,
    /// Light curve time bin size in seconds.
    #[arg(long)]
    pub tbin: Option<f64>,
    /// Light curve minimum energy (keV).
    #[arg(long)]
    pub emin: Option<f32>,
    /// Light curve maximum energy (keV).
    #[arg(long)]
    pub emax: Option<f32>,
    /// Bindata output products (`spec`, `lc`, or `both`).
    #[arg(long, value_enum, default_value_t = CztiOutputTypeArg::Both)]
    pub output_type: CztiOutputTypeArg,
    /// Enable coded-mask weighting in bindata (overrides config `mask_weight = false`).
    #[arg(long, conflicts_with = "no_mask_weight")]
    pub mask_weight: bool,
    /// Disable coded-mask weighting in bindata.
    #[arg(long)]
    pub no_mask_weight: bool,
    /// Generate `.rsp` response matrices (overrides config `gen_resp = false`).
    #[arg(long, conflicts_with = "no_resp")]
    pub resp: bool,
    /// Skip generating `.rsp` response matrices.
    #[arg(long)]
    pub no_resp: bool,
    /// Also write mask-weighted Level-2 event file (`*_weight.evt`).
    #[arg(long, conflicts_with = "no_gen_evt")]
    pub gen_evt: bool,
    /// Do not write mask-weighted Level-2 event file (overrides config `generate_eventfile = true`).
    #[arg(long)]
    pub no_gen_evt: bool,
    /// Sky image reconstruction algorithm (`balcc` or `fft`).
    #[arg(long, value_enum)]
    pub image_algo: Option<CztiImageAlgoArg>,
    /// Oversampling factor for FFT sky imaging.
    #[arg(long)]
    pub oversampling: Option<usize>,
    /// Comma-separated quadrants to process (e.g. `0,1,2,3`).
    #[arg(long, value_delimiter = ',')]
    pub quads: Option<Vec<usize>>,
}

#[derive(Debug, Args)]
pub struct CztiDecodeArgs {
    /// Input Level-1 modeM0 FITS file.
    #[arg(long, required = true)]
    pub m0: PathBuf,
    /// Input Level-1 TCT file.
    #[arg(long, required = true)]
    pub tct: PathBuf,
    /// Output Level-2 event FITS file.
    #[arg(short = 'o', long, required = true)]
    pub output: PathBuf,
    /// Output bunch table FITS file.
    #[arg(long)]
    pub bunch: Option<PathBuf>,
    /// Output frame header FITS file.
    #[arg(long)]
    pub hdr: Option<PathBuf>,
    /// Interpret Level-1 packets as little-endian instead of standard big-endian.
    #[arg(long)]
    pub little_endian: bool,
}

#[derive(Debug, Args)]
pub struct CztiBunchcleanArgs {
    #[arg(short = 'i', long, required = true)]
    pub input: PathBuf,
    #[arg(long, required = true)]
    pub bunch: PathBuf,
    #[arg(short = 'o', long, required = true)]
    pub output: PathBuf,
    #[arg(long)]
    pub out_bunch: Option<PathBuf>,
    #[arg(long)]
    pub out_livetime: Option<PathBuf>,
    #[arg(long)]
    pub generalized: bool,
}

#[derive(Debug, Args)]
pub struct CztiPha2energyArgs {
    #[arg(short = 'i', long, required = true)]
    pub input: PathBuf,
    #[arg(long, required = true)]
    pub caldb: PathBuf,
    /// Optional modeSS FITS file (either Level-2 with `TEMP` HDU or Level-1 `modeSS` with `--tct`).
    #[arg(long)]
    pub ss: Option<PathBuf>,
    /// Level-1 TCT file (required when `--ss` is a Level-1 `modeSS` packet FITS file).
    #[arg(long)]
    pub tct: Option<PathBuf>,
    #[arg(short = 'o', long, required = true)]
    pub output: PathBuf,
}

#[derive(Debug, Args)]
pub struct CztiGtigenArgs {
    #[arg(long, required = true)]
    pub m0: PathBuf,
    #[arg(long, required = true)]
    pub mkf: PathBuf,
    #[arg(short = 'o', long, required = true)]
    pub output: PathBuf,
    #[arg(long)]
    pub mkf_thresholds: Option<PathBuf>,
    #[arg(long)]
    pub user_gti: Option<PathBuf>,
    #[arg(long, default_value_t = 1.05)]
    pub min_exposure: f64,
}

#[derive(Debug, Args)]
pub struct CztiGaasArgs {
    #[arg(long, required = true)]
    pub mkf: PathBuf,
    #[arg(long, required = true)]
    pub caldb: PathBuf,
    #[arg(short = 'o', long, required = true)]
    pub output_dir: PathBuf,
    #[arg(long)]
    pub gti: Option<PathBuf>,
    #[arg(long, value_enum, default_value_t = CztiGtiTypeArg::Quad)]
    pub gti_type: CztiGtiTypeArg,
    #[arg(long, default_value = "czti_level2")]
    pub prefix: String,
}

#[derive(Debug, Args)]
pub struct CztiDataselArgs {
    #[arg(short = 'i', long, required = true)]
    pub input: PathBuf,
    #[arg(long, required = true)]
    pub gti: PathBuf,
    #[arg(short = 'o', long, required = true)]
    pub output: PathBuf,
    #[arg(long, value_enum, default_value_t = CztiGtiTypeArg::Quad)]
    pub gti_type: CztiGtiTypeArg,
}

#[derive(Debug, Args)]
pub struct CztiCleanArgs {
    #[arg(short = 'i', long, required = true)]
    pub input: PathBuf,
    #[arg(long, required = true)]
    pub caldb: PathBuf,
    #[arg(short = 'o', long, required = true)]
    pub output: PathBuf,
    #[arg(long)]
    pub out_badpix: Option<PathBuf>,
    #[arg(long)]
    pub out_livetime: Option<PathBuf>,
    #[arg(long)]
    pub livetime: Option<PathBuf>,
    #[arg(long)]
    pub bunch: Option<PathBuf>,
    #[arg(long, value_enum, default_value_t = CztiAlgoArg::Default)]
    pub algo: CztiAlgoArg,
    #[arg(long, value_enum, default_value_t = CztiGtiTypeArg::Quad)]
    pub gti_type: CztiGtiTypeArg,
}

#[derive(Debug, Args)]
pub struct CztiEvtencorArgs {
    #[arg(short = 'i', long, required = true)]
    pub input: PathBuf,
    #[arg(long, required = true)]
    pub caldb: PathBuf,
    #[arg(short = 'o', long, required = true)]
    pub output: PathBuf,
    #[arg(long)]
    pub no_lld: bool,
    #[arg(long)]
    pub no_uld: bool,
    #[arg(long)]
    pub no_time_gcor: bool,
    #[arg(long)]
    pub seed: Option<i32>,
}

#[derive(Debug, Args)]
pub struct CztiEvtcleanArgs {
    #[arg(short = 'i', long, required = true)]
    pub input: PathBuf,
    #[arg(short = 'o', long, required = true)]
    pub output: PathBuf,
    #[arg(long, default_value_t = 0)]
    pub alpha: u8,
    #[arg(long, default_value_t = 0)]
    pub veto_min: u16,
    #[arg(long, default_value_t = 0)]
    pub veto_max: u16,
    #[arg(long)]
    pub badpix: Option<PathBuf>,
}

#[derive(Debug, Args)]
pub struct CztiDpigenArgs {
    #[arg(short = 'i', long, required = true)]
    pub input: PathBuf,
    #[arg(long, required = true)]
    pub caldb: PathBuf,
    #[arg(long, required = true)]
    pub badpix: PathBuf,
    #[arg(long, required = true)]
    pub out_dph: PathBuf,
    #[arg(long, required = true)]
    pub out_dpi: PathBuf,
    #[arg(long, value_enum, default_value_t = CztiGtiTypeArg::Quad)]
    pub gti_type: CztiGtiTypeArg,
    #[arg(long, default_value_t = 20.0)]
    pub emin: f32,
    #[arg(long, default_value_t = 100.0)]
    pub emax: f32,
}

#[derive(Debug, Args)]
pub struct CztiImageArgs {
    #[arg(long, required = true)]
    pub dpi: PathBuf,
    /// Path to aspect file (`*.aspect_Q0`) or directory containing `*_Q0..Q3`.
    #[arg(long, required = true)]
    pub aspect: PathBuf,
    #[arg(long, required = true)]
    pub caldb: PathBuf,
    #[arg(short = 'o', long, required = true)]
    pub output_dir: PathBuf,
    #[arg(long, default_value = "czti_level2")]
    pub prefix: String,
    #[arg(long, value_enum, default_value_t = CztiImageAlgoArg::Balcc)]
    pub algo: CztiImageAlgoArg,
    #[arg(long, default_value_t = 16)]
    pub oversampling: usize,
    #[arg(long, value_delimiter = ',', default_value = "0,1,2,3")]
    pub quads: Vec<usize>,
}

#[derive(Debug, Args)]
pub struct CztiBindataArgs {
    #[arg(short = 'i', long, required = true)]
    pub input: PathBuf,
    #[arg(long, required = true)]
    pub mkf: PathBuf,
    #[arg(long, required = true)]
    pub caldb: PathBuf,
    #[arg(long, required = true)]
    pub badpix: PathBuf,
    #[arg(long, required = true)]
    pub livetime: PathBuf,
    #[arg(short = 'o', long, required = true)]
    pub output_dir: PathBuf,
    #[arg(long, default_value = "czti_level2")]
    pub prefix: String,
    /// Source RA in degrees for bindata (defaults to MKF `RA_PNT` / event header).
    #[arg(long, allow_hyphen_values = true)]
    pub ra: Option<f64>,
    /// Source DEC in degrees for bindata (defaults to MKF `DEC_PNT` / event header).
    #[arg(long, allow_hyphen_values = true)]
    pub dec: Option<f64>,
    #[arg(long, default_value_t = 1.0)]
    pub tbin: f64,
    #[arg(long, default_value_t = 20.0)]
    pub emin: f32,
    #[arg(long, default_value_t = 100.0)]
    pub emax: f32,
    #[arg(long, value_enum, default_value_t = CztiOutputTypeArg::Both)]
    pub output_type: CztiOutputTypeArg,
    #[arg(long)]
    pub no_mask_weight: bool,
    #[arg(long)]
    pub no_resp: bool,
    #[arg(long)]
    pub gen_evt: bool,
    #[arg(long, value_enum, default_value_t = CztiGtiTypeArg::Quad)]
    pub gti_type: CztiGtiTypeArg,
    #[arg(long, value_delimiter = ',', default_value = "0,1,2,3")]
    pub quads: Vec<usize>,
}

#[derive(Debug, Args)]
pub struct CztiAddspecArgs {
    /// Input CZTI `.pha` spectrum files to co-add.
    #[arg(long, required = true, num_args = 1..)]
    pub pha: Vec<PathBuf>,
    /// Optional input CZTI `.rsp` response files matching `--pha`.
    #[arg(long, num_args = 0..)]
    pub rsp: Vec<PathBuf>,
    /// Output co-added `.pha` spectrum path.
    #[arg(short = 'o', long = "out-pha", required = true)]
    pub out_pha: PathBuf,
    /// Optional output co-added `.rsp` response path.
    #[arg(long = "out-rsp")]
    pub out_rsp: Option<PathBuf>,
    /// Sum exposure across inputs (multi-orbit) instead of averaging (multi-quadrant).
    #[arg(long)]
    pub add_exposure: bool,
}

// ============================================================================
// LAXPC Subcommands
// ============================================================================

#[derive(Debug, Subcommand)]
pub enum LaxpcCommand {
    /// Run the full end-to-end AstroSat LAXPC Level-1 -> Level-2 pipeline.
    Run(LaxpcRunArgs),
    /// Stage 1+2: Sort, deduplicate, and decode Level-1 EA/BB telemetry frames.
    Decode(LaxpcDecodeArgs),
    /// Stage 4: Generate standard MKF + BBC Good Time Intervals (`GTI`).
    Gtigen(LaxpcGtigenArgs),
    /// Stage 5a: Generate a deadtime-corrected OGIP `RATE` light curve (`.lcurv`).
    Lightcurve(LaxpcLightcurveArgs),
    /// Stage 5b: Generate a deadtime-corrected OGIP `SPECTRUM` (`.spec`).
    Spectrum(LaxpcSpectrumArgs),
    /// Stage 6: Compute gain shift (`gaina8.dat` / `lxhv.dat`), shift background spectrum, and recommend RMF.
    Backshift(LaxpcBackshiftArgs),
    /// Build a 1-second cadence `AstroSat.orb` FITS file from MKF orbital columns.
    Orbit(LaxpcOrbitArgs),
}

#[derive(Debug, Args)]
pub struct LaxpcRunArgs {
    /// Level-1 LAXPC root, `laxpc/`, or orbit directory for automatic file discovery.
    #[arg(short = 'i', long = "input", value_name = "L1_DIR")]
    pub input: Option<PathBuf>,
    /// Output directory for Format B / Format A Level-2 products.
    #[arg(short = 'o', long = "output", value_name = "OUT_DIR")]
    pub output: Option<PathBuf>,
    /// Comma-separated LAXPC units (`10,20,30`, `lxp1,lxp2`, or `all`).
    #[arg(long, value_delimiter = ',')]
    pub units: Option<Vec<String>>,
    /// LAXPC CALDB directory (`gaina8.dat`, `lxhv.dat`, `lxp{10,20,30}.channel`).
    #[arg(long, value_name = "CALDB_DIR")]
    pub caldb: Option<PathBuf>,
    /// Light curve time bin width in seconds.
    #[arg(long)]
    pub tbin: Option<f64>,
    /// Non-paralyzable deadtime constant in seconds (default 43e-6).
    #[arg(long)]
    pub deadtime: Option<f64>,
    /// Spectrum channel count (`256`, `512`, or `1024`).
    #[arg(long)]
    pub channels: Option<usize>,
    /// Write combined multi-unit Format A `level2.event.fits`.
    #[arg(long, conflicts_with = "no_format_a")]
    pub format_a: bool,
    /// Skip combined Format A output (overrides config `write_format_a_combined = true`).
    #[arg(long)]
    pub no_format_a: bool,
    /// Write per-unit Format B `.evn` event FITS files (overrides config `write_format_b_events = false`).
    #[arg(long, conflicts_with = "no_format_b")]
    pub format_b: bool,
    /// Skip writing per-unit Format B `.evn` event FITS files.
    #[arg(long)]
    pub no_format_b: bool,
    /// Write `AstroSat.orb` (overrides config `write_orbit_file = false`).
    #[arg(long, conflicts_with = "no_orbit")]
    pub orbit: bool,
    /// Skip writing `AstroSat.orb`.
    #[arg(long)]
    pub no_orbit: bool,
    /// Also write ASCII `.pha` and `.lc` diagnostics.
    #[arg(long, conflicts_with = "no_ascii")]
    pub ascii: bool,
    /// Skip ASCII `.pha` and `.lc` diagnostics (overrides config `write_ascii_products = true`).
    #[arg(long)]
    pub no_ascii: bool,
    /// Retain only single-anode events (exclude Xe K-escape doubles).
    #[arg(long, conflicts_with = "no_single_only")]
    pub single_only: bool,
    /// Include both single and Xe K-escape double events (overrides config `single_only = true`).
    #[arg(long)]
    pub no_single_only: bool,
    /// Retain only Layer 1 (anodes A1/A2).
    #[arg(long, conflicts_with = "no_layer1_only")]
    pub layer1_only: bool,
    /// Use default anode selection instead of Layer 1 only (overrides config `layer1_only = true`).
    #[arg(long)]
    pub no_layer1_only: bool,
    /// Anode selection (`all`, `layer1`, or comma-separated anode IDs like `1,2,3`).
    #[arg(long)]
    pub anodes: Option<String>,
    /// Energy selection in keV as `EMIN:EMAX` (requires `lxp{10,20,30}.channel` in `--caldb`).
    #[arg(long, value_name = "EMIN:EMAX")]
    pub energy_range: Option<String>,
    /// Minimum 1024-scale raw PHA channel.
    #[arg(long)]
    pub pha_min: Option<u16>,
    /// Maximum 1024-scale raw PHA channel.
    #[arg(long)]
    pub pha_max: Option<u16>,
    /// Minimum Earth elevation angle (`ELV`) in degrees (default 3.5).
    #[arg(long, allow_hyphen_values = true)]
    pub elv_min: Option<f64>,
    /// SAA entry longitude threshold in degrees.
    #[arg(long, allow_hyphen_values = true)]
    pub saa_entry_lon: Option<f64>,
    /// SAA exit longitude threshold in degrees.
    #[arg(long, allow_hyphen_values = true)]
    pub saa_exit_lon: Option<f64>,
    /// Merge GTI gaps shorter than this duration (seconds).
    #[arg(long)]
    pub merge_gap: Option<f64>,
    /// Drop GTI intervals shorter than this duration (seconds).
    #[arg(long)]
    pub min_duration: Option<f64>,
    /// Optional user GTI file (FITS `GTI` or 2-column ASCII UTC intervals).
    #[arg(long, value_name = "GTI")]
    pub user_gti: Option<PathBuf>,
    /// Bypass automatic GTI generation and use `--user-gti` directly.
    #[arg(long, conflicts_with = "no_user_gti_only")]
    pub user_gti_only: bool,
    /// Intersect `--user-gti` with automatic GTI (overrides config `user_gti_only = true`).
    #[arg(long)]
    pub no_user_gti_only: bool,
    /// Optional 2D background model file (`backfit.lxp?0*`).
    #[arg(long, value_name = "FILE")]
    pub bkg_model: Option<PathBuf>,
    /// Optional ASCII background spectrum (`backlxp?0*.pha`).
    #[arg(long, value_name = "FILE")]
    pub bkg_pha: Option<PathBuf>,
    /// Optional explicit 30 keV peak channel (`chb`) for custom `--bkg-pha` spectra.
    #[arg(long, value_name = "CH")]
    pub bkg_peak_ch: Option<f64>,
}

#[derive(Debug, Args)]
pub struct LaxpcDecodeArgs {
    /// Level-1 modeEA FITS file(s).
    #[arg(long, required = true)]
    pub ea: Vec<PathBuf>,
    /// Optional Level-1 modeBB FITS file(s).
    #[arg(long)]
    pub bb: Vec<PathBuf>,
    /// LAXPC unit (`10`, `20`, `30`, `lxp1`, `lxp2`, `lxp3`); auto-detected if omitted.
    #[arg(long)]
    pub unit: Option<String>,
    /// Optional TCT file(s) if writing calibrated Format B `.evn` output.
    #[arg(long)]
    pub tct: Vec<PathBuf>,
    /// Optional GTI FITS file for event filtering when writing `.evn`.
    #[arg(long)]
    pub gti: Option<PathBuf>,
    /// Optional output Format B `.evn` event FITS file (requires `--tct`).
    #[arg(short = 'o', long)]
    pub output: Option<PathBuf>,
}

#[derive(Debug, Args)]
pub struct LaxpcGtigenArgs {
    #[arg(long, required = true)]
    pub ea: Vec<PathBuf>,
    #[arg(long, required = true)]
    pub mkf: Vec<PathBuf>,
    #[arg(long, required = true)]
    pub tct: Vec<PathBuf>,
    #[arg(long)]
    pub bb: Vec<PathBuf>,
    #[arg(long)]
    pub unit: Option<String>,
    #[arg(short = 'o', long, required = true)]
    pub output: PathBuf,
    #[arg(long, default_value_t = 3.5, allow_hyphen_values = true)]
    pub elv_min: f64,
    #[arg(long)]
    pub user_gti: Option<PathBuf>,
}

#[derive(Debug, Args)]
pub struct LaxpcLightcurveArgs {
    #[arg(long, required = true)]
    pub ea: Vec<PathBuf>,
    #[arg(long, required = true)]
    pub mkf: Vec<PathBuf>,
    #[arg(long, required = true)]
    pub tct: Vec<PathBuf>,
    #[arg(long)]
    pub bb: Vec<PathBuf>,
    #[arg(long)]
    pub gti: Option<PathBuf>,
    #[arg(long)]
    pub unit: Option<String>,
    #[arg(long, default_value_t = 1.0)]
    pub tbin: f64,
    #[arg(long, default_value_t = 43e-6)]
    pub deadtime: f64,
    #[arg(long)]
    pub caldb: Option<PathBuf>,
    #[arg(long)]
    pub bkg_model: Option<PathBuf>,
    #[arg(long)]
    pub single_only: bool,
    #[arg(long)]
    pub layer1_only: bool,
    /// Energy selection in keV as `EMIN:EMAX` (requires `lxp{10,20,30}.channel` in `--caldb`).
    #[arg(long, value_name = "EMIN:EMAX")]
    pub energy_range: Option<String>,
    /// Minimum 1024-scale raw PHA channel.
    #[arg(long)]
    pub pha_min: Option<u16>,
    /// Maximum 1024-scale raw PHA channel.
    #[arg(long)]
    pub pha_max: Option<u16>,
    #[arg(short = 'o', long, required = true)]
    pub output: PathBuf,
    #[arg(long)]
    pub ascii: Option<PathBuf>,
}

#[derive(Debug, Args)]
pub struct LaxpcSpectrumArgs {
    #[arg(long, required = true)]
    pub ea: Vec<PathBuf>,
    #[arg(long, required = true)]
    pub mkf: Vec<PathBuf>,
    #[arg(long, required = true)]
    pub tct: Vec<PathBuf>,
    #[arg(long)]
    pub bb: Vec<PathBuf>,
    #[arg(long)]
    pub gti: Option<PathBuf>,
    #[arg(long)]
    pub unit: Option<String>,
    #[arg(long)]
    pub channels: Option<usize>,
    #[arg(long, default_value_t = 43e-6)]
    pub deadtime: f64,
    /// LAXPC CALDB directory containing `lxp{10,20,30}.channel` (required for `--energy-range`).
    #[arg(long)]
    pub caldb: Option<PathBuf>,
    #[arg(long)]
    pub single_only: bool,
    #[arg(long)]
    pub layer1_only: bool,
    /// Energy selection in keV as `EMIN:EMAX` (requires `lxp{10,20,30}.channel` in `--caldb`).
    #[arg(long, value_name = "EMIN:EMAX")]
    pub energy_range: Option<String>,
    /// Minimum 1024-scale raw PHA channel.
    #[arg(long)]
    pub pha_min: Option<u16>,
    /// Maximum 1024-scale raw PHA channel.
    #[arg(long)]
    pub pha_max: Option<u16>,
    #[arg(short = 'o', long, required = true)]
    pub output: PathBuf,
    #[arg(long)]
    pub ascii: Option<PathBuf>,
}

#[derive(Debug, Args)]
pub struct LaxpcBackshiftArgs {
    /// LAXPC unit (`10`, `20`, `30`, `lxp1`, `lxp2`, `lxp3`).
    #[arg(long, required = true)]
    pub unit: String,
    /// LAXPC CALDB directory containing `gaina8.dat`, `lxhv.dat`, and `lxp{10,20,30}.channel`.
    #[arg(long, required = true)]
    pub caldb: PathBuf,
    /// Observation mission day (days since 2015-09-28 / AstroSat launch).
    #[arg(long)]
    pub mission_day: Option<f64>,
    /// Alternative: extract observation mission day from a FITS `TSTART` keyword.
    #[arg(long)]
    pub fits: Option<PathBuf>,
    /// Optional background `.pha` file (`backlxp?0*.pha`) to gain-shift.
    #[arg(long)]
    pub bkg_pha: Option<PathBuf>,
    /// Optional explicit 30 keV peak channel (`chb`) for custom `--bkg-pha` spectra.
    #[arg(long, value_name = "CH")]
    pub bkg_peak_ch: Option<f64>,
    /// Optional output FITS path for the shifted background spectrum.
    #[arg(short = 'o', long)]
    pub output: Option<PathBuf>,
    #[arg(long)]
    pub layer1_only: bool,
    #[arg(long)]
    pub single_only: bool,
    /// Minimum 1024-scale raw PHA channel for light-curve shift fraction (`fshbk`) calculation.
    #[arg(long, default_value_t = 0)]
    pub pha_min: u16,
    /// Maximum 1024-scale raw PHA channel for light-curve shift fraction (`fshbk`) calculation.
    #[arg(long, default_value_t = 1023)]
    pub pha_max: u16,
}

#[derive(Debug, Args)]
pub struct LaxpcOrbitArgs {
    /// Input MKF file(s) (`*.mkf`).
    #[arg(long, required = true)]
    pub mkf: Vec<PathBuf>,
    /// Output `AstroSat.orb` FITS file path.
    #[arg(short = 'o', long, required = true)]
    pub output: PathBuf,
}

#[cfg(test)]
mod tests {
    use super::Cli;
    use clap::CommandFactory;

    #[test]
    fn verify_clap_cli_invariants() {
        Cli::command().debug_assert();
    }
}
