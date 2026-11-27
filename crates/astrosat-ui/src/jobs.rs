//! Reproducible CLI/TOML preview generator (`POST /api/cli-preview`) and in-process
//! asynchronous pipeline job runner (`POST /api/jobs/run`, `GET /api/jobs`, `GET /api/jobs/{id}`).

use std::collections::HashMap;
use std::panic::{catch_unwind, AssertUnwindSafe};
use std::path::{Path, PathBuf};
use std::str::FromStr;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Arc, Mutex};
use std::time::Instant;

use anyhow::{bail, Context, Result};
use astrosat_core::{
    barycenter_fits, toml_quote, BaryFitsOptions, CoordFrame, EphemerisKind, GtiInterval, GtiSet,
};
use astrosat_czti::bindata::BinDataOutputType;
use astrosat_czti::image::ImageAlgorithm;
use astrosat_czti::pipeline::{run_pipeline as run_czti_pipeline, CztiPipelineConfig};
use astrosat_czti::stages::run_addspec;
use astrosat_czti::types::{EventSelectionAlgo, GtiType};
use astrosat_fits::{FitsFile, OgipMeta};
use astrosat_laxpc::background::{
    background_lightcurve_shift_fraction, compute_gain_shift, shift_background_spectrum,
    BackgroundCatalog, BackgroundPha, ChannelMap, GainTable, HvJumpTable,
};
use astrosat_laxpc::gti::{load_merged_mkf_files, GtiConfig};
use astrosat_laxpc::pipeline::{
    run_pipeline as run_laxpc_pipeline, write_orbit_from_mkf, PipelineConfig as LaxpcPipelineConfig,
};
use astrosat_laxpc::products::write_spectrum_fits;
use astrosat_laxpc::types::{utc_to_mission_day, AnodeSelection, EventFilter, LaxpcUnit};
use chrono::Utc;
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};

use crate::browse::classify_product;
use crate::http_util::{resolve_existing_dir, resolve_existing_file, resolve_output_path, ApiError, ApiResult};

// ---------------------------------------------------------------------------
// Request / Response Types
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Default, Deserialize, Serialize)]
#[serde(default)]
pub struct CztiJobParams {
    pub input_dir: Option<String>,
    pub m0: Option<String>,
    pub ss: Option<String>,
    pub tct: Option<String>,
    pub mkf: Option<String>,
    pub caldb_dir: Option<String>,
    pub output_dir: Option<String>,
    pub output_prefix: Option<String>,
    pub mkf_threshold_file: Option<String>,
    pub user_gti_file: Option<String>,
    pub entry_stage: Option<u8>,
    pub exit_stage: Option<u8>,
    pub gti_type: Option<String>,
    pub algo: Option<String>,
    pub ra_deg: Option<f64>,
    pub dec_deg: Option<f64>,
    pub time_binsize: Option<f64>,
    pub emin_kev: Option<f32>,
    pub emax_kev: Option<f32>,
    pub output_type: Option<String>,
    pub mask_weight: Option<bool>,
    pub gen_resp: Option<bool>,
    pub generate_eventfile: Option<bool>,
    pub image_algo: Option<String>,
    pub oversampling: Option<usize>,
    pub quads: Option<Vec<usize>>,
}

#[derive(Debug, Clone, Default, Deserialize, Serialize)]
#[serde(default)]
pub struct LaxpcJobParams {
    pub input_dir: Option<String>,
    pub output_dir: Option<String>,
    pub caldb_dir: Option<String>,
    pub units: Option<Vec<String>>,
    pub tbin: Option<f64>,
    pub deadtime_sec: Option<f64>,
    pub spectrum_channels: Option<usize>,
    pub write_format_a_combined: Option<bool>,
    pub write_format_b_events: Option<bool>,
    pub write_orbit_file: Option<bool>,
    pub write_ascii_products: Option<bool>,
    pub single_only: Option<bool>,
    pub layer1_only: Option<bool>,
    pub anodes: Option<String>,
    pub energy_range: Option<String>,
    pub pha_min: Option<u16>,
    pub pha_max: Option<u16>,
    pub min_elv_deg: Option<f64>,
    pub saa_entry_lon_deg: Option<f64>,
    pub saa_exit_lon_deg: Option<f64>,
    pub merge_gap_sec: Option<f64>,
    pub min_duration_sec: Option<f64>,
    pub user_gti_file: Option<String>,
    pub user_gti_only: Option<bool>,
    pub background_model_path: Option<String>,
    pub background_pha_path: Option<String>,
    pub background_peak_ch: Option<f64>,
}

#[derive(Debug, Clone, Default, Deserialize, Serialize)]
#[serde(default)]
pub struct BaryJobParams {
    pub input: Option<String>,
    pub output: Option<String>,
    pub orbit: Option<String>,
    pub mkf: Option<String>,
    pub ra_deg: Option<f64>,
    pub dec_deg: Option<f64>,
    pub frame: Option<String>,
    pub ephemeris: Option<String>,
    pub write_checksums: Option<bool>,
}

#[derive(Debug, Clone, Default, Deserialize, Serialize)]
#[serde(default)]
pub struct CztiAddspecJobParams {
    pub pha: Vec<String>,
    pub rsp: Vec<String>,
    pub out_pha: Option<String>,
    pub out_rsp: Option<String>,
    pub add_exposure: Option<bool>,
}

#[derive(Debug, Clone, Default, Deserialize, Serialize)]
#[serde(default)]
pub struct LaxpcBackshiftJobParams {
    pub unit: Option<String>,
    pub caldb_dir: Option<String>,
    pub mission_day: Option<f64>,
    pub fits: Option<String>,
    pub bkg_pha: Option<String>,
    pub bkg_peak_ch: Option<f64>,
    pub output: Option<String>,
    pub layer1_only: Option<bool>,
    pub single_only: Option<bool>,
    pub pha_min: Option<u16>,
    pub pha_max: Option<u16>,
}

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct PipelineJobRequest {
    pub pipeline: String,
    #[serde(default)]
    pub threads: Option<usize>,
    #[serde(default)]
    pub czti: Option<CztiJobParams>,
    #[serde(default)]
    pub laxpc: Option<LaxpcJobParams>,
    #[serde(default)]
    pub bary: Option<BaryJobParams>,
    #[serde(default)]
    pub czti_addspec: Option<CztiAddspecJobParams>,
    #[serde(default)]
    pub laxpc_backshift: Option<LaxpcBackshiftJobParams>,
}

#[derive(Debug, Clone, Serialize)]
pub struct CliPreviewResponse {
    pub pipeline: String,
    pub cli_command: String,
    pub toml_config: String,
}

#[derive(Debug, Clone, Serialize)]
pub struct JobOutputFile {
    pub label: String,
    pub path: String,
    pub product: String,
    pub size_bytes: u64,
}

#[derive(Debug, Clone, Serialize)]
pub struct JobRecord {
    pub id: String,
    pub pipeline: String,
    pub status: String,
    pub created_at: String,
    pub elapsed_ms: u64,
    pub cli_command: String,
    pub toml_config: String,
    pub summary: Option<Value>,
    pub outputs: Vec<JobOutputFile>,
    pub error: Option<String>,
}

#[derive(Debug, Clone, Serialize)]
pub struct JobListResponse {
    pub jobs: Vec<JobRecord>,
}

// ---------------------------------------------------------------------------
// Thread-Safe Job Manager
// ---------------------------------------------------------------------------

pub struct JobManager {
    counter: AtomicU64,
    jobs: Mutex<HashMap<String, (Instant, JobRecord)>>,
    order: Mutex<Vec<String>>,
}

impl Default for JobManager {
    fn default() -> Self {
        Self::new()
    }
}

impl JobManager {
    pub fn new() -> Self {
        Self {
            counter: AtomicU64::new(1),
            jobs: Mutex::new(HashMap::new()),
            order: Mutex::new(Vec::new()),
        }
    }

    pub fn list(&self) -> JobListResponse {
        let order = self.order.lock().unwrap_or_else(|e| e.into_inner());
        let map = self.jobs.lock().unwrap_or_else(|e| e.into_inner());
        let mut jobs = Vec::with_capacity(order.len());
        for id in order.iter().rev() {
            if let Some((started, rec)) = map.get(id) {
                let mut r = rec.clone();
                if r.status == "running" {
                    r.elapsed_ms = started.elapsed().as_millis() as u64;
                }
                jobs.push(r);
            }
        }
        JobListResponse { jobs }
    }

    pub fn get(&self, id: &str) -> Option<JobRecord> {
        let map = self.jobs.lock().unwrap_or_else(|e| e.into_inner());
        map.get(id).map(|(started, rec)| {
            let mut r = rec.clone();
            if r.status == "running" {
                r.elapsed_ms = started.elapsed().as_millis() as u64;
            }
            r
        })
    }

    fn insert(&self, rec: JobRecord) {
        let id = rec.id.clone();
        if let Ok(mut map) = self.jobs.lock() {
            map.insert(id.clone(), (Instant::now(), rec));
        }
        if let Ok(mut order) = self.order.lock() {
            order.push(id);
        }
    }

    fn mark_succeeded(&self, id: &str, summary: Value, outputs: Vec<JobOutputFile>) {
        if let Ok(mut map) = self.jobs.lock() {
            if let Some((started, rec)) = map.get_mut(id) {
                rec.status = "succeeded".to_string();
                rec.elapsed_ms = started.elapsed().as_millis() as u64;
                rec.summary = Some(summary);
                rec.outputs = outputs;
                rec.error = None;
            }
        }
    }

    fn mark_failed(&self, id: &str, err_msg: String) {
        if let Ok(mut map) = self.jobs.lock() {
            if let Some((started, rec)) = map.get_mut(id) {
                rec.status = "failed".to_string();
                rec.elapsed_ms = started.elapsed().as_millis() as u64;
                rec.error = Some(err_msg);
            }
        }
    }

    /// Validate request, spawn in-process background job, and return initial `running` record.
    pub fn spawn_job(
        self: &Arc<Self>,
        req: PipelineJobRequest,
        data_dir: PathBuf,
    ) -> ApiResult<JobRecord> {
        let preview = build_cli_preview(&req)?;
        let pipeline = preview.pipeline.clone();

        if let Some(t) = req.threads {
            if t == 0 || t > 512 {
                return Err(ApiError::bad_request(format!(
                    "threads must be in 1..=512, got {t}"
                )));
            }
        }

        // Pre-validate required paths synchronously so bad requests fail immediately with 4xx!
        let prepared = prepare_job(&req, &pipeline, &data_dir)?;

        let num = self.counter.fetch_add(1, Ordering::SeqCst);
        let id = format!("job-{num}");
        let initial = JobRecord {
            id: id.clone(),
            pipeline: pipeline.clone(),
            status: "running".to_string(),
            created_at: Utc::now().to_rfc3339(),
            elapsed_ms: 0,
            cli_command: preview.cli_command,
            toml_config: preview.toml_config,
            summary: None,
            outputs: Vec::new(),
            error: None,
        };
        self.insert(initial.clone());

        let mgr = Arc::clone(self);
        let job_id = id.clone();
        let threads = req.threads;
        std::thread::Builder::new()
            .name(format!("astrosat-ui-{job_id}"))
            .spawn(move || {
                let res = catch_unwind(AssertUnwindSafe(|| {
                    if let Some(n_threads) = threads {
                        match rayon::ThreadPoolBuilder::new().num_threads(n_threads).build() {
                            Ok(pool) => pool.install(|| execute_prepared_job(prepared)),
                            Err(e) => Err(anyhow::anyhow!("failed creating rayon pool: {e}")),
                        }
                    } else {
                        execute_prepared_job(prepared)
                    }
                }));
                match res {
                    Ok(Ok((summary, outputs))) => mgr.mark_succeeded(&job_id, summary, outputs),
                    Ok(Err(e)) => mgr.mark_failed(&job_id, format!("{e:#}")),
                    Err(_) => mgr.mark_failed(&job_id, "job panicked unexpectedly".to_string()),
                }
            })
            .map_err(|e| ApiError::internal(format!("failed spawning worker thread: {e}")))?;

        Ok(initial)
    }
}

// ---------------------------------------------------------------------------
// CLI & TOML Preview Generator
// ---------------------------------------------------------------------------

pub fn build_cli_preview(req: &PipelineJobRequest) -> ApiResult<CliPreviewResponse> {
    if let Some(t) = req.threads {
        if t == 0 || t > 512 {
            return Err(ApiError::bad_request(format!(
                "threads must be in 1..=512, got {t}"
            )));
        }
    }
    let pipe = req.pipeline.trim().to_ascii_lowercase();
    match pipe.as_str() {
        "czti" => {
            let p = req.czti.clone().unwrap_or_default();
            let (cli, toml) = format_czti_preview(&p, req.threads)?;
            Ok(CliPreviewResponse {
                pipeline: "czti".to_string(),
                cli_command: cli,
                toml_config: toml,
            })
        }
        "laxpc" => {
            let p = req.laxpc.clone().unwrap_or_default();
            let (cli, toml) = format_laxpc_preview(&p, req.threads)?;
            Ok(CliPreviewResponse {
                pipeline: "laxpc".to_string(),
                cli_command: cli,
                toml_config: toml,
            })
        }
        "bary" => {
            let p = req.bary.clone().unwrap_or_default();
            if non_empty(&p.orbit).is_some() && non_empty(&p.mkf).is_some() {
                return Err(ApiError::bad_request(
                    "bary orbit and mkf are mutually exclusive; provide only one",
                ));
            }
            let (cli, toml) = format_bary_preview(&p, req.threads)?;
            Ok(CliPreviewResponse {
                pipeline: "bary".to_string(),
                cli_command: cli,
                toml_config: toml,
            })
        }
        "czti-addspec" => {
            let p = req.czti_addspec.clone().unwrap_or_default();
            let (cli, toml) = format_czti_addspec_preview(&p, req.threads)?;
            Ok(CliPreviewResponse {
                pipeline: "czti-addspec".to_string(),
                cli_command: cli,
                toml_config: toml,
            })
        }
        "laxpc-backshift" => {
            let p = req.laxpc_backshift.clone().unwrap_or_default();
            if p.mission_day.is_some() && non_empty(&p.fits).is_some() {
                return Err(ApiError::bad_request(
                    "laxpc-backshift mission_day and fits are mutually exclusive; provide only one",
                ));
            }
            let (cli, toml) = format_laxpc_backshift_preview(&p, req.threads)?;
            Ok(CliPreviewResponse {
                pipeline: "laxpc-backshift".to_string(),
                cli_command: cli,
                toml_config: toml,
            })
        }
        other => Err(ApiError::bad_request(format!(
            "unsupported pipeline '{other}' (expected 'czti', 'laxpc', 'bary', 'czti-addspec', or 'laxpc-backshift')"
        ))),
    }
}

fn sh_quote(s: &str) -> String {
    if !s.is_empty()
        && s.bytes()
            .all(|b| b.is_ascii_alphanumeric() || b"-_./:=+@".contains(&b))
    {
        s.to_string()
    } else {
        format!("'{}'", s.replace('\'', "'\\''"))
    }
}

fn non_empty(opt: &Option<String>) -> Option<&str> {
    opt.as_deref().map(str::trim).filter(|s| !s.is_empty())
}

fn canonicalize_czti_output_type(raw: &str) -> ApiResult<(&'static str, BinDataOutputType)> {
    match raw.trim().to_ascii_lowercase().as_str() {
        "both" | "all" => Ok(("both", BinDataOutputType::Both)),
        "spec" | "spectrum" => Ok(("spec", BinDataOutputType::Spec)),
        "lc" | "lightcurve" => Ok(("lc", BinDataOutputType::Lc)),
        other => Err(ApiError::bad_request(format!(
            "invalid output_type '{other}' (expected 'spec', 'lc', or 'both')"
        ))),
    }
}

fn parse_and_canonicalize_laxpc_energy_range(raw: &str) -> ApiResult<(f64, f64, String)> {
    let trimmed = raw.trim();
    let (elo, ehi) = if let Some((a, b)) = trimmed.split_once(':') {
        (a.trim(), b.trim())
    } else {
        let parts: Vec<&str> = trimmed.split('-').collect();
        if parts.len() != 2 {
            return Err(ApiError::bad_request(format!(
                "invalid energy_range '{trimmed}' (expected EMIN:EMAX or EMIN-EMAX, e.g. '3:80')"
            )));
        }
        (parts[0].trim(), parts[1].trim())
    };
    let lo: f64 = elo
        .parse()
        .map_err(|_| ApiError::bad_request("invalid energy_range lower bound"))?;
    let hi: f64 = ehi
        .parse()
        .map_err(|_| ApiError::bad_request("invalid energy_range upper bound"))?;
    if !lo.is_finite() || !hi.is_finite() || lo < 0.0 || hi <= lo {
        return Err(ApiError::bad_request(
            "energy_range must have finite 0 <= emin < emax",
        ));
    }
    Ok((lo, hi, format!("{lo}:{hi}")))
}

fn format_czti_preview(p: &CztiJobParams, threads: Option<usize>) -> ApiResult<(String, String)> {
    if let Some(es) = p.entry_stage {
        if !(1..=3).contains(&es) {
            return Err(ApiError::bad_request("entry_stage must be 1..=3"));
        }
    }
    if let Some(xs) = p.exit_stage {
        if !(1..=3).contains(&xs) {
            return Err(ApiError::bad_request("exit_stage must be 1..=3"));
        }
    }
    let es_eff = p.entry_stage.unwrap_or(1);
    let xs_eff = p.exit_stage.unwrap_or(3);
    if es_eff > xs_eff {
        return Err(ApiError::bad_request(format!(
            "entry_stage ({es_eff}) must be <= exit_stage ({xs_eff})"
        )));
    }
    if let Some(emin) = p.emin_kev {
        if !emin.is_finite() || emin < 0.0 {
            return Err(ApiError::bad_request("emin_kev must be finite and >= 0"));
        }
    }
    if let Some(emax) = p.emax_kev {
        if !emax.is_finite() || emax <= 0.0 {
            return Err(ApiError::bad_request("emax_kev must be finite and > 0"));
        }
    }
    if let (Some(emin), Some(emax)) = (p.emin_kev, p.emax_kev) {
        if emin >= emax {
            return Err(ApiError::bad_request(format!(
                "emin_kev ({emin}) must be < emax_kev ({emax})"
            )));
        }
    }
    if let Some(tb) = p.time_binsize {
        if !tb.is_finite() || tb <= 0.0 {
            return Err(ApiError::bad_request("time_binsize must be finite and > 0"));
        }
    }
    if let Some(ra) = p.ra_deg {
        if !ra.is_finite() || !(-360.0..=360.0).contains(&ra) {
            return Err(ApiError::bad_request("ra_deg must be finite in [-360, 360]"));
        }
    }
    if let Some(dec) = p.dec_deg {
        if !dec.is_finite() || !(-90.0..=90.0).contains(&dec) {
            return Err(ApiError::bad_request("dec_deg must be finite in [-90, 90]"));
        }
    }
    if let Some(ov) = p.oversampling {
        if ov == 0 || ov > 64 {
            return Err(ApiError::bad_request("oversampling must be 1..=64"));
        }
    }
    if let Some(ref qs) = p.quads {
        if qs.is_empty() {
            return Err(ApiError::bad_request("quads list cannot be empty"));
        }
        if qs.iter().any(|&q| q > 3) {
            return Err(ApiError::bad_request("quads must contain values in 0..=3"));
        }
    }
    if let Some(gt) = non_empty(&p.gti_type) {
        if GtiType::parse(gt).is_none() {
            return Err(ApiError::bad_request(format!(
                "invalid gti_type '{gt}' (expected 'quad' or 'common')"
            )));
        }
    }
    if let Some(al) = non_empty(&p.algo) {
        if !al.eq_ignore_ascii_case("default") && !al.eq_ignore_ascii_case("generalized") {
            return Err(ApiError::bad_request(format!("invalid algo '{al}'")));
        }
    }
    if let Some(ia) = non_empty(&p.image_algo) {
        if !ia.eq_ignore_ascii_case("balcc") && !ia.eq_ignore_ascii_case("fft") {
            return Err(ApiError::bad_request(format!("invalid image_algo '{ia}'")));
        }
    }

    let mut args = vec!["astrosat".to_string()];
    let mut toml = String::new();
    if let Some(t) = threads {
        args.push("--threads".into());
        args.push(t.to_string());
        toml.push_str(&format!("threads = {t}\n\n"));
    }
    args.push("czti".into());
    args.push("run".into());
    toml.push_str("[czti]\n");

    if let Some(v) = non_empty(&p.input_dir) {
        args.extend(["--input".into(), sh_quote(v)]);
        toml.push_str(&format!("input_dir = {}\n", toml_quote(v)));
    }
    if let Some(v) = non_empty(&p.m0) {
        args.extend(["--m0".into(), sh_quote(v)]);
        toml.push_str(&format!("l1_m0_file = {}\n", toml_quote(v)));
    }
    if let Some(v) = non_empty(&p.ss) {
        args.extend(["--ss".into(), sh_quote(v)]);
        toml.push_str(&format!("l1_ss_file = {}\n", toml_quote(v)));
    }
    if let Some(v) = non_empty(&p.tct) {
        args.extend(["--tct".into(), sh_quote(v)]);
        toml.push_str(&format!("tct_file = {}\n", toml_quote(v)));
    }
    if let Some(v) = non_empty(&p.mkf) {
        args.extend(["--mkf".into(), sh_quote(v)]);
        toml.push_str(&format!("mkf_file = {}\n", toml_quote(v)));
    }
    if let Some(v) = non_empty(&p.caldb_dir) {
        args.extend(["--caldb".into(), sh_quote(v)]);
        toml.push_str(&format!("caldb_dir = {}\n", toml_quote(v)));
    }
    if let Some(v) = non_empty(&p.output_dir) {
        args.extend(["--output".into(), sh_quote(v)]);
        toml.push_str(&format!("output_dir = {}\n", toml_quote(v)));
    }
    if let Some(v) = non_empty(&p.output_prefix) {
        args.extend(["--prefix".into(), sh_quote(v)]);
        toml.push_str(&format!("output_prefix = {}\n", toml_quote(v)));
    }
    if let Some(v) = non_empty(&p.mkf_threshold_file) {
        args.extend(["--mkf-thresholds".into(), sh_quote(v)]);
        toml.push_str(&format!("mkf_threshold_file = {}\n", toml_quote(v)));
    }
    if let Some(v) = non_empty(&p.user_gti_file) {
        args.extend(["--user-gti".into(), sh_quote(v)]);
        toml.push_str(&format!("user_gti_file = {}\n", toml_quote(v)));
    }
    if let Some(v) = p.entry_stage {
        args.extend(["--entry-stage".into(), v.to_string()]);
        toml.push_str(&format!("entry_stage = {v}\n"));
    }
    if let Some(v) = p.exit_stage {
        args.extend(["--exit-stage".into(), v.to_string()]);
        toml.push_str(&format!("exit_stage = {v}\n"));
    }
    if let Some(v) = non_empty(&p.gti_type) {
        args.extend(["--gti-type".into(), sh_quote(v)]);
        toml.push_str(&format!("gti_type = {}\n", toml_quote(v)));
    }
    if let Some(v) = non_empty(&p.algo) {
        args.extend(["--algo".into(), sh_quote(v)]);
        toml.push_str(&format!("algo = {}\n", toml_quote(v)));
    }
    if let Some(v) = p.ra_deg.filter(|x| x.is_finite()) {
        args.extend(["--ra".into(), format!("{v}")]);
        toml.push_str(&format!("ra_deg = {v}\n"));
    }
    if let Some(v) = p.dec_deg.filter(|x| x.is_finite()) {
        args.extend(["--dec".into(), format!("{v}")]);
        toml.push_str(&format!("dec_deg = {v}\n"));
    }
    if let Some(v) = p.time_binsize.filter(|x| x.is_finite() && *x > 0.0) {
        args.extend(["--tbin".into(), format!("{v}")]);
        toml.push_str(&format!("time_binsize = {v}\n"));
    }
    if let Some(v) = p.emin_kev.filter(|x| x.is_finite()) {
        args.extend(["--emin".into(), format!("{v}")]);
        toml.push_str(&format!("emin_kev = {v}\n"));
    }
    if let Some(v) = p.emax_kev.filter(|x| x.is_finite()) {
        args.extend(["--emax".into(), format!("{v}")]);
        toml.push_str(&format!("emax_kev = {v}\n"));
    }
    if let Some(v) = non_empty(&p.output_type) {
        let (canon, _) = canonicalize_czti_output_type(v)?;
        args.extend(["--output-type".into(), sh_quote(canon)]);
        toml.push_str(&format!("output_type = {}\n", toml_quote(canon)));
    }
    if let Some(v) = p.mask_weight {
        args.push(if v { "--mask-weight" } else { "--no-mask-weight" }.into());
        toml.push_str(&format!("mask_weight = {v}\n"));
    }
    if let Some(v) = p.gen_resp {
        args.push(if v { "--resp" } else { "--no-resp" }.into());
        toml.push_str(&format!("gen_resp = {v}\n"));
    }
    if let Some(v) = p.generate_eventfile {
        args.push(if v { "--gen-evt" } else { "--no-gen-evt" }.into());
        toml.push_str(&format!("generate_eventfile = {v}\n"));
    }
    if let Some(v) = non_empty(&p.image_algo) {
        args.extend(["--image-algo".into(), sh_quote(v)]);
        toml.push_str(&format!("image_algo = {}\n", toml_quote(v)));
    }
    if let Some(v) = p.oversampling {
        args.extend(["--oversampling".into(), v.to_string()]);
        toml.push_str(&format!("oversampling = {v}\n"));
    }
    if let Some(ref qs) = p.quads {
        if !qs.is_empty() {
            let joined = qs
                .iter()
                .map(|q| q.to_string())
                .collect::<Vec<_>>()
                .join(",");
            args.extend(["--quads".into(), joined]);
            let arr = qs
                .iter()
                .map(|q| q.to_string())
                .collect::<Vec<_>>()
                .join(", ");
            toml.push_str(&format!("quads = [{arr}]\n"));
        }
    }

    Ok((args.join(" "), toml))
}

fn format_laxpc_preview(p: &LaxpcJobParams, threads: Option<usize>) -> ApiResult<(String, String)> {
    let mut args = vec!["astrosat".to_string()];
    let mut toml = String::new();
    if let Some(t) = threads {
        args.push("--threads".into());
        args.push(t.to_string());
        toml.push_str(&format!("threads = {t}\n\n"));
    }
    args.push("laxpc".into());
    args.push("run".into());
    toml.push_str("[laxpc]\n");

    if let Some(v) = non_empty(&p.input_dir) {
        args.extend(["--input".into(), sh_quote(v)]);
        toml.push_str(&format!("input_dir = {}\n", toml_quote(v)));
    }
    if let Some(v) = non_empty(&p.output_dir) {
        args.extend(["--output".into(), sh_quote(v)]);
        toml.push_str(&format!("output_dir = {}\n", toml_quote(v)));
    }
    if let Some(v) = non_empty(&p.caldb_dir) {
        args.extend(["--caldb".into(), sh_quote(v)]);
        toml.push_str(&format!("caldb_dir = {}\n", toml_quote(v)));
    }
    if let Some(ref us) = p.units {
        if us.is_empty() {
            return Err(ApiError::bad_request("units list cannot be empty"));
        }
        for u in us {
            LaxpcUnit::from_str(u.trim()).map_err(|_| {
                ApiError::bad_request(format!("invalid LAXPC unit '{u}'"))
            })?;
        }
        args.extend(["--units".into(), us.join(",")]);
        let arr = us
            .iter()
            .map(|u| toml_quote(u))
            .collect::<Vec<_>>()
            .join(", ");
        toml.push_str(&format!("units = [{arr}]\n"));
    }
    if let Some(v) = p.tbin {
        if !v.is_finite() || v <= 0.0 {
            return Err(ApiError::bad_request("tbin must be finite and > 0"));
        }
        args.extend(["--tbin".into(), format!("{v}")]);
        toml.push_str(&format!("tbin = {v}\n"));
    }
    if let Some(v) = p.deadtime_sec {
        if !v.is_finite() || v < 0.0 {
            return Err(ApiError::bad_request("deadtime_sec must be finite and >= 0"));
        }
        args.extend(["--deadtime".into(), format!("{v}")]);
        toml.push_str(&format!("deadtime_sec = {v}\n"));
    }
    if let Some(v) = p.spectrum_channels {
        if !matches!(v, 256 | 512 | 1024) {
            return Err(ApiError::bad_request(
                "spectrum_channels must be 256, 512, or 1024",
            ));
        }
        args.extend(["--channels".into(), v.to_string()]);
        toml.push_str(&format!("spectrum_channels = {v}\n"));
    }
    if let Some(v) = p.write_format_a_combined {
        args.push(if v { "--format-a" } else { "--no-format-a" }.into());
        toml.push_str(&format!("write_format_a_combined = {v}\n"));
    }
    if let Some(v) = p.write_format_b_events {
        args.push(if v { "--format-b" } else { "--no-format-b" }.into());
        toml.push_str(&format!("write_format_b_events = {v}\n"));
    }
    if let Some(v) = p.write_orbit_file {
        args.push(if v { "--orbit" } else { "--no-orbit" }.into());
        toml.push_str(&format!("write_orbit_file = {v}\n"));
    }
    if let Some(v) = p.write_ascii_products {
        args.push(if v { "--ascii" } else { "--no-ascii" }.into());
        toml.push_str(&format!("write_ascii_products = {v}\n"));
    }
    if let Some(v) = p.single_only {
        args.push(if v { "--single-only" } else { "--no-single-only" }.into());
        toml.push_str(&format!("single_only = {v}\n"));
    }
    if let Some(v) = p.layer1_only {
        args.push(if v { "--layer1-only" } else { "--no-layer1-only" }.into());
        toml.push_str(&format!("layer1_only = {v}\n"));
    }
    if let Some(v) = non_empty(&p.anodes) {
        parse_laxpc_anodes(v)?;
        args.extend(["--anodes".into(), sh_quote(v)]);
        toml.push_str(&format!("anodes = {}\n", toml_quote(v)));
    }
    if let Some(v) = non_empty(&p.energy_range) {
        let (_, _, canon) = parse_and_canonicalize_laxpc_energy_range(v)?;
        args.extend(["--energy-range".into(), sh_quote(&canon)]);
        toml.push_str(&format!("energy_range = {}\n", toml_quote(&canon)));
    }
    let pha_min_eff = p.pha_min.unwrap_or(0);
    let pha_max_eff = p.pha_max.unwrap_or(1023);
    if pha_max_eff > 1023 || pha_min_eff > pha_max_eff {
        return Err(ApiError::bad_request(format!(
            "invalid PHA range [{pha_min_eff}, {pha_max_eff}]"
        )));
    }
    if let Some(v) = p.pha_min {
        args.extend(["--pha-min".into(), v.to_string()]);
        toml.push_str(&format!("pha_min = {v}\n"));
    }
    if let Some(v) = p.pha_max {
        args.extend(["--pha-max".into(), v.to_string()]);
        toml.push_str(&format!("pha_max = {v}\n"));
    }
    if let Some(v) = p.min_elv_deg {
        if !v.is_finite() {
            return Err(ApiError::bad_request("min_elv_deg must be finite"));
        }
        args.extend(["--elv-min".into(), format!("{v}")]);
        toml.push_str(&format!("min_elv_deg = {v}\n"));
    }
    if let Some(v) = p.saa_entry_lon_deg {
        if !v.is_finite() {
            return Err(ApiError::bad_request("saa_entry_lon_deg must be finite"));
        }
        args.extend(["--saa-entry-lon".into(), format!("{v}")]);
        toml.push_str(&format!("saa_entry_lon_deg = {v}\n"));
    }
    if let Some(v) = p.saa_exit_lon_deg {
        if !v.is_finite() {
            return Err(ApiError::bad_request("saa_exit_lon_deg must be finite"));
        }
        args.extend(["--saa-exit-lon".into(), format!("{v}")]);
        toml.push_str(&format!("saa_exit_lon_deg = {v}\n"));
    }
    if let Some(v) = p.merge_gap_sec {
        if !v.is_finite() || v < 0.0 {
            return Err(ApiError::bad_request(
                "merge_gap_sec must be finite and >= 0",
            ));
        }
        args.extend(["--merge-gap".into(), format!("{v}")]);
        toml.push_str(&format!("merge_gap_sec = {v}\n"));
    }
    if let Some(v) = p.min_duration_sec {
        if !v.is_finite() || v < 0.0 {
            return Err(ApiError::bad_request(
                "min_duration_sec must be finite and >= 0",
            ));
        }
        args.extend(["--min-duration".into(), format!("{v}")]);
        toml.push_str(&format!("min_duration_sec = {v}\n"));
    }
    if let Some(v) = non_empty(&p.user_gti_file) {
        args.extend(["--user-gti".into(), sh_quote(v)]);
        toml.push_str(&format!("user_gti_file = {}\n", toml_quote(v)));
    }
    if let Some(v) = p.user_gti_only {
        args.push(if v { "--user-gti-only" } else { "--no-user-gti-only" }.into());
        toml.push_str(&format!("user_gti_only = {v}\n"));
    }
    if let Some(v) = non_empty(&p.background_model_path) {
        args.extend(["--bkg-model".into(), sh_quote(v)]);
        toml.push_str(&format!("background_model_path = {}\n", toml_quote(v)));
    }
    if let Some(v) = non_empty(&p.background_pha_path) {
        args.extend(["--bkg-pha".into(), sh_quote(v)]);
        toml.push_str(&format!("background_pha_path = {}\n", toml_quote(v)));
    }
    if let Some(v) = p.background_peak_ch {
        if !v.is_finite() || v <= 0.0 {
            return Err(ApiError::bad_request(
                "background_peak_ch must be finite and > 0",
            ));
        }
        args.extend(["--bkg-peak-ch".into(), format!("{v}")]);
        toml.push_str(&format!("background_peak_ch = {v}\n"));
    }

    Ok((args.join(" "), toml))
}

fn format_bary_preview(p: &BaryJobParams, threads: Option<usize>) -> ApiResult<(String, String)> {
    if let Some(ra) = p.ra_deg {
        if !ra.is_finite() || !(-360.0..=360.0).contains(&ra) {
            return Err(ApiError::bad_request("ra_deg must be finite in [-360, 360]"));
        }
    }
    if let Some(dec) = p.dec_deg {
        if !dec.is_finite() || !(-90.0..=90.0).contains(&dec) {
            return Err(ApiError::bad_request("dec_deg must be finite in [-90, 90]"));
        }
    }
    if let Some(f) = non_empty(&p.frame) {
        if !f.eq_ignore_ascii_case("icrs") && !f.eq_ignore_ascii_case("fk5") {
            return Err(ApiError::bad_request(format!(
                "invalid frame '{f}' (expected 'icrs' or 'fk5')"
            )));
        }
    }

    let mut args = vec!["astrosat".to_string()];
    let mut toml = String::new();
    if let Some(t) = threads {
        args.push("--threads".into());
        args.push(t.to_string());
        toml.push_str(&format!("threads = {t}\n\n"));
    }
    args.push("bary".into());
    toml.push_str("[bary]\n");

    if let Some(v) = non_empty(&p.input) {
        args.push(sh_quote(v));
    }
    if let Some(v) = non_empty(&p.orbit) {
        args.extend(["--orbit".into(), sh_quote(v)]);
        toml.push_str(&format!("orbit_file = {}\n", toml_quote(v)));
    } else if let Some(v) = non_empty(&p.mkf) {
        args.extend(["--mkf".into(), sh_quote(v)]);
        toml.push_str(&format!(
            "# mkf = {} (CLI-only --mkf; AppConfig [bary] uses orbit_file)\n",
            toml_quote(v)
        ));
    }
    if let Some(v) = non_empty(&p.output) {
        args.extend(["--output".into(), sh_quote(v)]);
    }
    if let Some(v) = p.ra_deg.filter(|x| x.is_finite()) {
        args.extend(["--ra".into(), format!("{v}")]);
        toml.push_str(&format!("ra_deg = {v}\n"));
    }
    if let Some(v) = p.dec_deg.filter(|x| x.is_finite()) {
        args.extend(["--dec".into(), format!("{v}")]);
        toml.push_str(&format!("dec_deg = {v}\n"));
    }
    if let Some(v) = non_empty(&p.frame) {
        args.extend(["--frame".into(), sh_quote(v)]);
        toml.push_str(&format!("frame = {}\n", toml_quote(v)));
    }
    if let Some(v) = non_empty(&p.ephemeris) {
        args.extend(["--ephem".into(), sh_quote(v)]);
        toml.push_str(&format!("ephemeris = {}\n", toml_quote(v)));
    }
    if let Some(v) = p.write_checksums {
        args.push(if v { "--checksum" } else { "--no-checksum" }.into());
        toml.push_str(&format!("write_checksums = {v}\n"));
    }

    Ok((args.join(" "), toml))
}

fn format_czti_addspec_preview(
    p: &CztiAddspecJobParams,
    threads: Option<usize>,
) -> ApiResult<(String, String)> {
    if p.pha.iter().any(|s| s.trim().is_empty()) {
        return Err(ApiError::bad_request(
            "czti-addspec pha list contains an empty path",
        ));
    }
    if !p.rsp.is_empty() {
        if p.rsp.len() != p.pha.len() {
            return Err(ApiError::bad_request(format!(
                "czti-addspec rsp count ({}) must be 0 or match pha count ({})",
                p.rsp.len(),
                p.pha.len()
            )));
        }
        if p.rsp.iter().any(|s| s.trim().is_empty()) {
            return Err(ApiError::bad_request(
                "czti-addspec rsp list contains an empty path",
            ));
        }
    }
    // Reject out_rsp without rsp when input .pha rows are present (allow partial UI typing when pha is empty).
    if non_empty(&p.out_rsp).is_some() && p.rsp.is_empty() && !p.pha.is_empty() {
        return Err(ApiError::bad_request(
            "out_rsp requires matching input rsp files",
        ));
    }

    let mut args = vec!["astrosat".to_string()];
    let mut toml = String::new();
    if let Some(t) = threads {
        args.push("--threads".into());
        args.push(t.to_string());
        toml.push_str(&format!("threads = {t}\n\n"));
    }
    args.extend(["czti".into(), "addspec".into()]);

    let pha_items: Vec<&str> = p.pha.iter().map(|s| s.trim()).filter(|s| !s.is_empty()).collect();
    if !pha_items.is_empty() {
        args.push("--pha".into());
        for f in pha_items {
            args.push(sh_quote(f));
        }
    }
    let rsp_items: Vec<&str> = p.rsp.iter().map(|s| s.trim()).filter(|s| !s.is_empty()).collect();
    if !rsp_items.is_empty() {
        args.push("--rsp".into());
        for f in rsp_items {
            args.push(sh_quote(f));
        }
    }
    if let Some(v) = non_empty(&p.out_pha) {
        args.extend(["-o".into(), sh_quote(v)]);
    }
    if let Some(v) = non_empty(&p.out_rsp) {
        args.extend(["--out-rsp".into(), sh_quote(v)]);
    }
    if p.add_exposure == Some(true) {
        args.push("--add-exposure".into());
    }

    toml.push_str(
        "# Standalone CLI tool: astrosat czti addspec uses CLI flags (no AppConfig TOML table).\n",
    );
    Ok((args.join(" "), toml))
}

fn format_laxpc_backshift_preview(
    p: &LaxpcBackshiftJobParams,
    threads: Option<usize>,
) -> ApiResult<(String, String)> {
    if let Some(u) = non_empty(&p.unit) {
        LaxpcUnit::from_str(u).map_err(|e| ApiError::bad_request(format!("{e}")))?;
    }
    if let Some(md) = p.mission_day {
        if !md.is_finite() || md < 0.0 {
            return Err(ApiError::bad_request(
                "mission_day must be finite and >= 0",
            ));
        }
    }
    if let Some(ch) = p.bkg_peak_ch {
        if !ch.is_finite() || ch <= 0.0 {
            return Err(ApiError::bad_request(
                "bkg_peak_ch must be finite and > 0",
            ));
        }
    }
    let pha_min_eff = p.pha_min.unwrap_or(0);
    let pha_max_eff = p.pha_max.unwrap_or(1023);
    if pha_max_eff > 1023 || pha_min_eff > pha_max_eff {
        return Err(ApiError::bad_request(format!(
            "invalid PHA range [{pha_min_eff}, {pha_max_eff}]"
        )));
    }
    if non_empty(&p.output).is_some() && non_empty(&p.bkg_pha).is_none() {
        return Err(ApiError::bad_request(
            "output requires bkg_pha to write a shifted background spectrum",
        ));
    }

    let mut args = vec!["astrosat".to_string()];
    let mut toml = String::new();
    if let Some(t) = threads {
        args.push("--threads".into());
        args.push(t.to_string());
        toml.push_str(&format!("threads = {t}\n\n"));
    }
    args.extend(["laxpc".into(), "backshift".into()]);

    if let Some(v) = non_empty(&p.unit) {
        args.extend(["--unit".into(), sh_quote(v)]);
    }
    if let Some(v) = non_empty(&p.caldb_dir) {
        args.extend(["--caldb".into(), sh_quote(v)]);
    }
    if let Some(md) = p.mission_day {
        args.extend(["--mission-day".into(), format!("{md}")]);
    } else if let Some(fp) = non_empty(&p.fits) {
        args.extend(["--fits".into(), sh_quote(fp)]);
    }
    if let Some(v) = non_empty(&p.bkg_pha) {
        args.extend(["--bkg-pha".into(), sh_quote(v)]);
    }
    if let Some(ch) = p.bkg_peak_ch {
        args.extend(["--bkg-peak-ch".into(), format!("{ch}")]);
    }
    if let Some(v) = non_empty(&p.output) {
        args.extend(["-o".into(), sh_quote(v)]);
    }
    if p.layer1_only == Some(true) {
        args.push("--layer1-only".into());
    }
    if p.single_only == Some(true) {
        args.push("--single-only".into());
    }
    if let Some(v) = p.pha_min {
        args.extend(["--pha-min".into(), v.to_string()]);
    }
    if let Some(v) = p.pha_max {
        args.extend(["--pha-max".into(), v.to_string()]);
    }

    toml.push_str(
        "# Standalone CLI tool: astrosat laxpc backshift uses CLI flags (no AppConfig TOML table).\n",
    );
    Ok((args.join(" "), toml))
}

// ---------------------------------------------------------------------------
// Job Preparation & Execution
// ---------------------------------------------------------------------------

enum PreparedJob {
    Czti(Box<CztiPipelineConfig>),
    Laxpc(Box<LaxpcPipelineConfig>),
    Bary {
        input: PathBuf,
        output: PathBuf,
        orbit: Option<PathBuf>,
        mkf: Option<PathBuf>,
        ra_deg: Option<f64>,
        dec_deg: Option<f64>,
        opts: BaryFitsOptions,
    },
    CztiAddspec {
        pha_files: Vec<PathBuf>,
        rsp_files: Vec<PathBuf>,
        out_pha: PathBuf,
        out_rsp: Option<PathBuf>,
        add_exposure: bool,
    },
    LaxpcBackshift {
        unit: LaxpcUnit,
        caldb_dir: PathBuf,
        mission_day: Option<f64>,
        fits: Option<PathBuf>,
        bkg_pha: Option<PathBuf>,
        bkg_peak_ch: Option<f64>,
        output: Option<PathBuf>,
        layer1_only: bool,
        single_only: bool,
        pha_min: u16,
        pha_max: u16,
    },
}

fn prepare_job(
    req: &PipelineJobRequest,
    pipeline: &str,
    data_dir: &Path,
) -> ApiResult<PreparedJob> {
    match pipeline {
        "czti" => {
            let p = req.czti.as_ref().ok_or_else(|| {
                ApiError::bad_request("missing 'czti' configuration object for czti job")
            })?;
            // Validate numeric, stage, energy, and enum parameters up front
            let entry_stage = p.entry_stage.unwrap_or(1);
            if !(1..=3).contains(&entry_stage) {
                return Err(ApiError::bad_request("entry_stage must be 1..=3"));
            }
            let exit_stage = p.exit_stage.unwrap_or(3);
            if !(1..=3).contains(&exit_stage) {
                return Err(ApiError::bad_request("exit_stage must be 1..=3"));
            }
            if entry_stage > exit_stage {
                return Err(ApiError::bad_request(format!(
                    "entry_stage ({entry_stage}) must be <= exit_stage ({exit_stage})"
                )));
            }
            let gti_type = if let Some(gt) = non_empty(&p.gti_type) {
                Some(GtiType::parse(gt).ok_or_else(|| {
                    ApiError::bad_request(format!(
                        "invalid gti_type '{gt}' (expected 'quad' or 'common')"
                    ))
                })?)
            } else {
                None
            };
            let evt_sel_algo = if let Some(al) = non_empty(&p.algo) {
                Some(if al.eq_ignore_ascii_case("generalized") {
                    EventSelectionAlgo::Generalized
                } else if al.eq_ignore_ascii_case("default") {
                    EventSelectionAlgo::Default
                } else {
                    return Err(ApiError::bad_request(format!("invalid algo '{al}'")));
                })
            } else {
                None
            };
            let source_coords = match (p.ra_deg, p.dec_deg) {
                (Some(ra), Some(dec)) => {
                    if !ra.is_finite() || !(-360.0..=360.0).contains(&ra) {
                        return Err(ApiError::bad_request("ra_deg must be finite in [-360, 360]"));
                    }
                    if !dec.is_finite() || !(-90.0..=90.0).contains(&dec) {
                        return Err(ApiError::bad_request("dec_deg must be finite in [-90, 90]"));
                    }
                    Some((ra, dec))
                }
                (None, None) => None,
                _ => {
                    return Err(ApiError::bad_request(
                        "ra_deg and dec_deg must be provided together",
                    ));
                }
            };
            if let Some(tb) = p.time_binsize {
                if !tb.is_finite() || tb <= 0.0 {
                    return Err(ApiError::bad_request("time_binsize must be finite and > 0"));
                }
            }
            if let Some(emin) = p.emin_kev {
                if !emin.is_finite() || emin < 0.0 {
                    return Err(ApiError::bad_request("emin_kev must be finite and >= 0"));
                }
            }
            if let Some(emax) = p.emax_kev {
                if !emax.is_finite() || emax <= 0.0 {
                    return Err(ApiError::bad_request("emax_kev must be finite and > 0"));
                }
            }
            if let (Some(emin), Some(emax)) = (p.emin_kev, p.emax_kev) {
                if emin >= emax {
                    return Err(ApiError::bad_request(format!(
                        "emin_kev ({emin}) must be < emax_kev ({emax})"
                    )));
                }
            }
            let output_type = if let Some(ot) = non_empty(&p.output_type) {
                let (_, parsed_ot) = canonicalize_czti_output_type(ot)?;
                Some(parsed_ot)
            } else {
                None
            };
            let image_algo = if let Some(ia) = non_empty(&p.image_algo) {
                Some(match ia.to_ascii_lowercase().as_str() {
                    "balcc" => ImageAlgorithm::BalCC,
                    "fft" => ImageAlgorithm::Fft,
                    other => {
                        return Err(ApiError::bad_request(format!(
                            "invalid image_algo '{other}'"
                        )))
                    }
                })
            } else {
                None
            };
            if let Some(ov) = p.oversampling {
                if ov == 0 || ov > 64 {
                    return Err(ApiError::bad_request("oversampling must be 1..=64"));
                }
            }
            if let Some(ref qs) = p.quads {
                if qs.is_empty() {
                    return Err(ApiError::bad_request("quads list cannot be empty"));
                }
                if qs.iter().any(|&q| q > 3) {
                    return Err(ApiError::bad_request("quads must contain values in 0..=3"));
                }
            }

            let discovered = if let Some(dir_s) = non_empty(&p.input_dir) {
                let d = resolve_existing_dir(Some(dir_s), data_dir)?;
                Some(discover_czti_l1(&d))
            } else {
                None
            };
            let m0 = if let Some(s) = non_empty(&p.m0) {
                resolve_existing_file(s, data_dir)?
            } else if let Some(ref f) = discovered.as_ref().and_then(|d| d.m0.clone()) {
                f.clone()
            } else {
                return Err(ApiError::bad_request(
                    "missing CZTI modeM0 L1 FITS (set input_dir or m0)",
                ));
            };
            let ss = if let Some(s) = non_empty(&p.ss) {
                Some(resolve_existing_file(s, data_dir)?)
            } else {
                discovered.as_ref().and_then(|d| d.ss.clone())
            };
            let tct = if let Some(s) = non_empty(&p.tct) {
                resolve_existing_file(s, data_dir)?
            } else if let Some(ref f) = discovered.as_ref().and_then(|d| d.tct.clone()) {
                f.clone()
            } else {
                return Err(ApiError::bad_request(
                    "missing CZTI .tct file (set input_dir or tct)",
                ));
            };
            let mkf = if let Some(s) = non_empty(&p.mkf) {
                resolve_existing_file(s, data_dir)?
            } else if let Some(ref f) = discovered.as_ref().and_then(|d| d.mkf.clone()) {
                f.clone()
            } else {
                return Err(ApiError::bad_request(
                    "missing CZTI .mkf file (set input_dir or mkf)",
                ));
            };
            let caldb_s = non_empty(&p.caldb_dir).ok_or_else(|| {
                ApiError::bad_request("missing CZTI caldb_dir")
            })?;
            let caldb = resolve_existing_dir(Some(caldb_s), data_dir)?;
            let out_s = non_empty(&p.output_dir).ok_or_else(|| {
                ApiError::bad_request("missing CZTI output_dir")
            })?;
            let output_dir = resolve_output_path(out_s, data_dir)?;

            let mut cfg = CztiPipelineConfig::new(m0, ss, tct, mkf, caldb, output_dir);
            cfg.output_prefix = non_empty(&p.output_prefix).map(str::to_string);
            cfg.mkf_threshold_file = if let Some(s) = non_empty(&p.mkf_threshold_file) {
                Some(resolve_existing_file(s, data_dir)?)
            } else {
                discovered.as_ref().and_then(|d| d.mkf_thresholds.clone())
            };
            cfg.user_gti_file = if let Some(s) = non_empty(&p.user_gti_file) {
                Some(resolve_existing_file(s, data_dir)?)
            } else {
                None
            };
            cfg.entry_stage = entry_stage;
            cfg.exit_stage = exit_stage;
            if let Some(gt) = gti_type {
                cfg.gti_type = gt;
            }
            if let Some(al) = evt_sel_algo {
                cfg.evt_sel_algo = al;
            }
            if let Some((ra, dec)) = source_coords {
                cfg.bindata.ra_src_deg = ra;
                cfg.bindata.dec_src_deg = dec;
                cfg.bindata.source_coords = Some((ra, dec));
            }
            if let Some(tb) = p.time_binsize {
                cfg.bindata.time_binsize = tb;
            }
            if let Some(emin) = p.emin_kev {
                cfg.bindata.emin_kev = emin;
            }
            if let Some(emax) = p.emax_kev {
                cfg.bindata.emax_kev = emax;
            }
            if cfg.bindata.emin_kev >= cfg.bindata.emax_kev {
                return Err(ApiError::bad_request(format!(
                    "emin_kev ({}) must be < emax_kev ({})",
                    cfg.bindata.emin_kev, cfg.bindata.emax_kev
                )));
            }
            if let Some(ot) = output_type {
                cfg.bindata.output_type = ot;
            }
            if let Some(mw) = p.mask_weight {
                cfg.bindata.mask_weight = mw;
            }
            if let Some(gr) = p.gen_resp {
                cfg.bindata.gen_resp = gr;
            }
            if let Some(ge) = p.generate_eventfile {
                cfg.bindata.generate_eventfile = ge;
            }
            if let Some(ia) = image_algo {
                cfg.image.algorithm = ia;
            }
            if let Some(ov) = p.oversampling {
                cfg.image.oversampling_factor = ov;
            }
            if let Some(ref qs) = p.quads {
                if !qs.is_empty() {
                    cfg.bindata.quads = qs.clone();
                    cfg.image.quads = qs.clone();
                }
            }
            Ok(PreparedJob::Czti(Box::new(cfg)))
        }
        "laxpc" => {
            let p = req.laxpc.as_ref().ok_or_else(|| {
                ApiError::bad_request("missing 'laxpc' configuration object for laxpc job")
            })?;
            let in_s = non_empty(&p.input_dir)
                .ok_or_else(|| ApiError::bad_request("missing LAXPC input_dir"))?;
            let input_dir = resolve_existing_dir(Some(in_s), data_dir)?;
            let out_s = non_empty(&p.output_dir)
                .ok_or_else(|| ApiError::bad_request("missing LAXPC output_dir"))?;
            let output_dir = resolve_output_path(out_s, data_dir)?;
            let caldb_dir = if let Some(cs) = non_empty(&p.caldb_dir) {
                Some(resolve_existing_dir(Some(cs), data_dir)?)
            } else {
                None
            };

            let units = if let Some(ref raw_units) = p.units {
                if raw_units.is_empty() {
                    return Err(ApiError::bad_request("units list cannot be empty"));
                }
                let mut out = Vec::with_capacity(raw_units.len());
                for u in raw_units {
                    let parsed = LaxpcUnit::from_str(u.trim()).map_err(|_| {
                        ApiError::bad_request(format!("invalid LAXPC unit '{u}'"))
                    })?;
                    if !out.contains(&parsed) {
                        out.push(parsed);
                    }
                }
                out
            } else {
                LaxpcUnit::ALL.to_vec()
            };

            let tbin = p.tbin.unwrap_or(1.0);
            if !tbin.is_finite() || tbin <= 0.0 {
                return Err(ApiError::bad_request("tbin must be finite and > 0"));
            }
            let deadtime_sec = p.deadtime_sec.unwrap_or(43.0e-6);
            if !deadtime_sec.is_finite() || deadtime_sec < 0.0 {
                return Err(ApiError::bad_request("deadtime_sec must be finite and >= 0"));
            }
            let pha_min = p.pha_min.unwrap_or(0);
            let pha_max = p.pha_max.unwrap_or(1023);
            if pha_max > 1023 || pha_min > pha_max {
                return Err(ApiError::bad_request(format!(
                    "invalid PHA range [{pha_min}, {pha_max}]"
                )));
            }
            let energy_range_kev = if let Some(er) = non_empty(&p.energy_range) {
                let (lo, hi, _) = parse_and_canonicalize_laxpc_energy_range(er)?;
                Some((lo, hi))
            } else {
                None
            };
            if energy_range_kev.is_some() && caldb_dir.is_none() {
                return Err(ApiError::bad_request("energy_range requires caldb_dir"));
            }

            let mut gti_config = GtiConfig::default();
            if let Some(elv) = p.min_elv_deg {
                if !elv.is_finite() {
                    return Err(ApiError::bad_request("min_elv_deg must be finite"));
                }
                gti_config.min_elv_deg = elv;
            }
            if let Some(lon) = p.saa_entry_lon_deg {
                if !lon.is_finite() {
                    return Err(ApiError::bad_request("saa_entry_lon_deg must be finite"));
                }
                gti_config.saa_entry_lon_deg = lon;
            }
            if let Some(lon) = p.saa_exit_lon_deg {
                if !lon.is_finite() {
                    return Err(ApiError::bad_request("saa_exit_lon_deg must be finite"));
                }
                gti_config.saa_exit_lon_deg = lon;
            }
            if let Some(gap) = p.merge_gap_sec {
                if !gap.is_finite() || gap < 0.0 {
                    return Err(ApiError::bad_request(
                        "merge_gap_sec must be finite and >= 0",
                    ));
                }
                gti_config.merge_gap_sec = gap;
            }
            if let Some(dur) = p.min_duration_sec {
                if !dur.is_finite() || dur < 0.0 {
                    return Err(ApiError::bad_request(
                        "min_duration_sec must be finite and >= 0",
                    ));
                }
                gti_config.min_duration_sec = dur;
            }

            let anodes = if let Some(a) = non_empty(&p.anodes) {
                parse_laxpc_anodes(a)?
            } else if p.layer1_only.unwrap_or(false) {
                AnodeSelection::Layer1
            } else {
                AnodeSelection::AllMain
            };

            let user_gti_utc = if let Some(s) = non_empty(&p.user_gti_file) {
                let gti_path = resolve_existing_file(s, data_dir)?;
                Some(
                    load_user_gti(&gti_path).map_err(|e| {
                        ApiError::bad_request(format!(
                            "failed to load user_gti_file '{}': {e}",
                            gti_path.display()
                        ))
                    })?,
                )
            } else {
                None
            };
            let use_user_gti_only = p.user_gti_only.unwrap_or(false);
            if use_user_gti_only && user_gti_utc.is_none() {
                return Err(ApiError::bad_request(
                    "user_gti_only requires user_gti_file",
                ));
            }
            let background_model_path = if let Some(s) = non_empty(&p.background_model_path) {
                Some(resolve_existing_file(s, data_dir)?)
            } else {
                None
            };
            let background_pha_path = if let Some(s) = non_empty(&p.background_pha_path) {
                Some(resolve_existing_file(s, data_dir)?)
            } else {
                None
            };
            if let Some(ch) = p.background_peak_ch {
                if !ch.is_finite() || ch <= 0.0 {
                    return Err(ApiError::bad_request(
                        "background_peak_ch must be finite and > 0",
                    ));
                }
            }

            let filter = EventFilter {
                single_only: p.single_only.unwrap_or(false),
                anodes,
                pha_min,
                pha_max,
                energy_range_kev,
            };

            let mut cfg = LaxpcPipelineConfig {
                input_dir: Some(input_dir),
                output_dir,
                units,
                tbin,
                deadtime_sec,
                spectrum_channels: p.spectrum_channels,
                filter,
                gti_config,
                user_gti_utc,
                use_user_gti_only,
                caldb_dir,
                background_model_path,
                background_pha_path,
                background_peak_ch: p.background_peak_ch,
                write_format_a_combined: p.write_format_a_combined.unwrap_or(false),
                write_format_b_events: p.write_format_b_events.unwrap_or(true),
                write_orbit_file: p.write_orbit_file.unwrap_or(true),
                write_ascii_products: p.write_ascii_products.unwrap_or(false),
                ..LaxpcPipelineConfig::default()
            };
            if let Some(ch) = cfg.spectrum_channels {
                if !matches!(ch, 256 | 512 | 1024) {
                    return Err(ApiError::bad_request(
                        "spectrum_channels must be 256, 512, or 1024",
                    ));
                }
            }
            // Keep default orbits_by_unit empty for automatic discovery
            cfg.orbits_by_unit.clear();
            Ok(PreparedJob::Laxpc(Box::new(cfg)))
        }
        "bary" => {
            let p = req.bary.as_ref().ok_or_else(|| {
                ApiError::bad_request("missing 'bary' configuration object for bary job")
            })?;
            let in_s = non_empty(&p.input)
                .ok_or_else(|| ApiError::bad_request("missing bary input FITS file"))?;
            let input = resolve_existing_file(in_s, data_dir)?;
            let out_s = non_empty(&p.output)
                .ok_or_else(|| ApiError::bad_request("missing bary output FITS path"))?;
            let output = resolve_output_path(out_s, data_dir)?;
            let orbit = if let Some(s) = non_empty(&p.orbit) {
                Some(resolve_existing_file(s, data_dir)?)
            } else {
                None
            };
            let mkf = if let Some(s) = non_empty(&p.mkf) {
                Some(resolve_existing_file(s, data_dir)?)
            } else {
                None
            };
            if orbit.is_some() && mkf.is_some() {
                return Err(ApiError::bad_request(
                    "bary orbit and mkf are mutually exclusive; provide only one",
                ));
            }
            if let Some(ra) = p.ra_deg {
                if !ra.is_finite() || !(-360.0..=360.0).contains(&ra) {
                    return Err(ApiError::bad_request("ra_deg must be finite in [-360, 360]"));
                }
            }
            if let Some(dec) = p.dec_deg {
                if !dec.is_finite() || !(-90.0..=90.0).contains(&dec) {
                    return Err(ApiError::bad_request("dec_deg must be finite in [-90, 90]"));
                }
            }
            let frame = match non_empty(&p.frame) {
                None => CoordFrame::Icrs,
                Some(f) if f.eq_ignore_ascii_case("icrs") => CoordFrame::Icrs,
                Some(f) if f.eq_ignore_ascii_case("fk5") => CoordFrame::Fk5,
                Some(other) => {
                    return Err(ApiError::bad_request(format!(
                        "invalid frame '{other}' (expected 'icrs' or 'fk5')"
                    )))
                }
            };
            let ephemeris = match non_empty(&p.ephemeris) {
                None => EphemerisKind::Analytical,
                Some(e) if e.eq_ignore_ascii_case("analytical") => EphemerisKind::Analytical,
                Some(e) => {
                    let ep = resolve_existing_file(e, data_dir)?;
                    EphemerisKind::JplDe(ep)
                }
            };
            let opts = BaryFitsOptions {
                frame,
                ephemeris,
                write_checksums: p.write_checksums.unwrap_or(true),
            };
            Ok(PreparedJob::Bary {
                input,
                output,
                orbit,
                mkf,
                ra_deg: p.ra_deg,
                dec_deg: p.dec_deg,
                opts,
            })
        }
        "czti-addspec" => {
            let p = req.czti_addspec.as_ref().ok_or_else(|| {
                ApiError::bad_request(
                    "missing 'czti_addspec' configuration object for czti-addspec job",
                )
            })?;
            if p.pha.is_empty() {
                return Err(ApiError::bad_request(
                    "czti-addspec requires at least one input .pha file",
                ));
            }
            let mut pha_files = Vec::with_capacity(p.pha.len());
            for raw in &p.pha {
                let s = raw.trim();
                if s.is_empty() {
                    return Err(ApiError::bad_request(
                        "czti-addspec pha list contains an empty path",
                    ));
                }
                pha_files.push(resolve_existing_file(s, data_dir)?);
            }

            let mut rsp_files = Vec::new();
            if !p.rsp.is_empty() {
                if p.rsp.len() != pha_files.len() {
                    return Err(ApiError::bad_request(format!(
                        "czti-addspec rsp count ({}) must be 0 or match pha count ({})",
                        p.rsp.len(),
                        pha_files.len()
                    )));
                }
                rsp_files.reserve(p.rsp.len());
                for raw in &p.rsp {
                    let s = raw.trim();
                    if s.is_empty() {
                        return Err(ApiError::bad_request(
                            "czti-addspec rsp list contains an empty path",
                        ));
                    }
                    rsp_files.push(resolve_existing_file(s, data_dir)?);
                }
            }

            let out_pha_s = non_empty(&p.out_pha)
                .ok_or_else(|| ApiError::bad_request("missing czti-addspec out_pha path"))?;
            let out_pha = resolve_output_path(out_pha_s, data_dir)?;
            let out_rsp = if let Some(s) = non_empty(&p.out_rsp) {
                Some(resolve_output_path(s, data_dir)?)
            } else {
                None
            };
            if out_rsp.is_some() && rsp_files.is_empty() {
                return Err(ApiError::bad_request(
                    "out_rsp requires matching input rsp files",
                ));
            }

            Ok(PreparedJob::CztiAddspec {
                pha_files,
                rsp_files,
                out_pha,
                out_rsp,
                add_exposure: p.add_exposure.unwrap_or(false),
            })
        }
        "laxpc-backshift" => {
            let p = req.laxpc_backshift.as_ref().ok_or_else(|| {
                ApiError::bad_request(
                    "missing 'laxpc_backshift' configuration object for laxpc-backshift job",
                )
            })?;
            let unit_s = non_empty(&p.unit)
                .ok_or_else(|| ApiError::bad_request("missing laxpc-backshift unit"))?;
            let unit = LaxpcUnit::from_str(unit_s).map_err(|e| ApiError::bad_request(format!("{e}")))?;

            let caldb_s = non_empty(&p.caldb_dir)
                .ok_or_else(|| ApiError::bad_request("missing laxpc-backshift caldb_dir"))?;
            let caldb_dir = resolve_existing_dir(Some(caldb_s), data_dir)?;

            let fits_s = non_empty(&p.fits);
            if p.mission_day.is_some() && fits_s.is_some() {
                return Err(ApiError::bad_request(
                    "mission_day and fits are mutually exclusive; provide only one",
                ));
            }
            if p.mission_day.is_none() && fits_s.is_none() {
                return Err(ApiError::bad_request(
                    "provide either mission_day or fits for laxpc-backshift",
                ));
            }
            if let Some(md) = p.mission_day {
                if !md.is_finite() || md < 0.0 {
                    return Err(ApiError::bad_request(
                        "mission_day must be finite and >= 0",
                    ));
                }
            }
            let fits = if let Some(fs) = fits_s {
                Some(resolve_existing_file(fs, data_dir)?)
            } else {
                None
            };

            let bkg_pha = if let Some(s) = non_empty(&p.bkg_pha) {
                Some(resolve_existing_file(s, data_dir)?)
            } else {
                None
            };
            if let Some(ch) = p.bkg_peak_ch {
                if !ch.is_finite() || ch <= 0.0 {
                    return Err(ApiError::bad_request(
                        "bkg_peak_ch must be finite and > 0",
                    ));
                }
            }
            let output = if let Some(s) = non_empty(&p.output) {
                Some(resolve_output_path(s, data_dir)?)
            } else {
                None
            };
            if output.is_some() && bkg_pha.is_none() {
                return Err(ApiError::bad_request(
                    "output requires bkg_pha to write a shifted background spectrum",
                ));
            }

            let pha_min = p.pha_min.unwrap_or(0);
            let pha_max = p.pha_max.unwrap_or(1023);
            if pha_max > 1023 || pha_min > pha_max {
                return Err(ApiError::bad_request(format!(
                    "invalid PHA range [{pha_min}, {pha_max}]"
                )));
            }

            Ok(PreparedJob::LaxpcBackshift {
                unit,
                caldb_dir,
                mission_day: p.mission_day,
                fits,
                bkg_pha,
                bkg_peak_ch: p.bkg_peak_ch,
                output,
                layer1_only: p.layer1_only.unwrap_or(false),
                single_only: p.single_only.unwrap_or(false),
                pha_min,
                pha_max,
            })
        }
        other => Err(ApiError::bad_request(format!(
            "unsupported pipeline '{other}'"
        ))),
    }
}

fn execute_prepared_job(job: PreparedJob) -> Result<(Value, Vec<JobOutputFile>)> {
    match job {
        PreparedJob::Czti(cfg) => {
            let s = run_czti_pipeline(&cfg).context("CZTI pipeline failed")?;
            let mut outputs = Vec::new();
            let mut push = |label: &str, p: &Path| {
                if p.is_file() {
                    outputs.push(make_output_entry(label, p));
                }
            };
            if let Some(ref p) = s.clean_evt_path {
                push("Clean Events", p);
            }
            if let Some(ref p) = s.weighted_evt_path {
                push("Weighted Events", p);
            }
            if let Some(ref p) = s.bc_evt_path {
                push("Bunchclean Events", p);
            }
            if let Some(ref p) = s.gti_path {
                push("GTI", p);
            }
            if let Some(ref p) = s.dph_path {
                push("DPH", p);
            }
            if let Some(ref p) = s.dpi_path {
                push("DPI", p);
            }
            for p in &s.image_paths {
                push("Sky Image", p);
            }
            for p in &s.lightcurve_paths {
                push("Light Curve", p);
            }
            for p in &s.spectrum_paths {
                push("Spectrum", p);
            }
            for p in &s.response_paths {
                push("Response", p);
            }
            let summary = serde_json::to_value(&s)?;
            Ok((summary, outputs))
        }
        PreparedJob::Laxpc(cfg) => {
            let s = run_laxpc_pipeline(&cfg).context("LAXPC pipeline failed")?;
            let mut outputs = Vec::new();
            let mut push = |label: &str, p: &Path| {
                if p.is_file() {
                    outputs.push(make_output_entry(label, p));
                }
            };
            if let Some(ref p) = s.format_a_event_fits_path {
                push("Combined Format A Events", p);
            }
            if let Some(ref p) = s.orbit_fits_path {
                push("AstroSat Orbit", p);
            }
            for u in &s.unit_results {
                let tag = u.unit.prefix();
                if let Some(ref p) = u.event_fits_path {
                    push(&format!("{tag} Events"), p);
                }
                push(&format!("{tag} Light Curve"), &u.lightcurve_fits_path);
                push(&format!("{tag} Spectrum"), &u.spectrum_fits_path);
                push(&format!("{tag} GTI"), &u.gti_fits_path);
                if let Some(ref p) = u.bkg_spectrum_fits_path {
                    push(&format!("{tag} Background Spectrum"), p);
                }
                if let Some(ref p) = u.bkg_lightcurve_fits_path {
                    push(&format!("{tag} Background Light Curve"), p);
                }
            }
            let summary = json!({
                "units_processed": s.unit_results.len(),
                "format_a_total_events": s.format_a_total_events,
                "units": s.unit_results.iter().map(|u| json!({
                    "unit": u.unit.prefix(),
                    "num_orbits": u.num_orbits,
                    "total_valid_events": u.total_valid_events,
                    "gti_filtered_events": u.gti_filtered_events,
                    "exposure_sec": u.spectrum.exposure_sec,
                    "lightcurve_mean_rate_cps": u.lightcurve_mean_rate_cps,
                })).collect::<Vec<_>>()
            });
            Ok((summary, outputs))
        }
        PreparedJob::Bary {
            input,
            output,
            orbit,
            mkf,
            ra_deg,
            dec_deg,
            opts,
        } => {
            if let Some(parent) = output.parent() {
                if !parent.as_os_str().is_empty() {
                    std::fs::create_dir_all(parent)?;
                }
            }
            let (hdr_ra, hdr_dec) = read_ra_dec_from_fits(&input)?;
            let ra = ra_deg
                .or(hdr_ra)
                .context("RA not provided and not found in FITS headers")?;
            let dec = dec_deg
                .or(hdr_dec)
                .context("DEC not provided and not found in FITS headers")?;

            let mut synth_orb: Option<PathBuf> = None;
            let orbit_path = if let Some(orb) = orbit {
                orb
            } else if let Some(mkf_path) = mkf {
                let mkf_tbl = load_merged_mkf_files(std::slice::from_ref(&mkf_path))?;
                let mut meta = OgipMeta::default();
                if let Some((t0, t1)) = mkf_tbl.times().first().zip(mkf_tbl.times().last()) {
                    meta.tstart = *t0;
                    meta.tstop = *t1;
                }
                let tmp_orb = output.with_extension("tmp_synth.orb");
                write_orbit_from_mkf(&tmp_orb, &mkf_tbl, &meta)?;
                synth_orb = Some(tmp_orb.clone());
                tmp_orb
            } else if let Some(sibling) = input.parent().map(|d| d.join("AstroSat.orb")).filter(|p| p.is_file()) {
                sibling
            } else {
                bail!("missing orbit file (pass orbit or mkf, or place AstroSat.orb next to input)");
            };

            let res = barycenter_fits(&input, &output, &orbit_path, ra, dec, &opts);
            if let Some(ref tmp) = synth_orb {
                let _ = std::fs::remove_file(tmp);
            }
            let summary = res.context("barycentric correction failed")?;
            let outputs = if output.is_file() {
                vec![make_output_entry("Barycentered FITS", &output)]
            } else {
                Vec::new()
            };
            let summary_json = json!({
                "input": input.to_string_lossy(),
                "output": output.to_string_lossy(),
                "ra_deg": ra,
                "dec_deg": dec,
                "plephem": summary.plephem,
                "hdus_updated": summary.hdus_updated,
                "events_corrected": summary.events_corrected,
                "gti_rows_corrected": summary.gti_rows_corrected,
                "tstart_tdb": summary.tstart_tdb,
                "tstop_tdb": summary.tstop_tdb,
            });
            Ok((summary_json, outputs))
        }
        PreparedJob::CztiAddspec {
            pha_files,
            rsp_files,
            out_pha,
            out_rsp,
            add_exposure,
        } => {
            ensure_parent_dir(&out_pha)?;
            if let Some(ref rsp_p) = out_rsp {
                ensure_parent_dir(rsp_p)?;
            }
            let (spec, rsp) = run_addspec(
                &pha_files,
                &rsp_files,
                &out_pha,
                out_rsp.as_deref(),
                add_exposure,
            )
            .context("CZTI addspec failed")?;

            let mut outputs = Vec::new();
            if out_pha.is_file() {
                let mut e = make_output_entry("Co-added Spectrum", &out_pha);
                e.product = "spectrum".to_string();
                outputs.push(e);
            }
            if let (Some(ref p), Some(_)) = (&out_rsp, &rsp) {
                if p.is_file() {
                    let mut e = make_output_entry("Co-added Response", p);
                    e.product = "response".to_string();
                    outputs.push(e);
                }
            }

            let total_counts: f64 = spec.counts.iter().map(|&c| c as f64).sum();
            let summary = json!({
                "inputs": pha_files.len(),
                "rsp_inputs": rsp_files.len(),
                "out_pha": out_pha.to_string_lossy(),
                "out_rsp": if rsp.is_some() {
                    out_rsp.as_ref().map(|p| p.to_string_lossy().into_owned())
                } else {
                    None
                },
                "add_exposure": add_exposure,
                "exposure_sec": spec.exposure,
                "channels": spec.channel.len(),
                "total_counts": total_counts,
            });
            Ok((summary, outputs))
        }
        PreparedJob::LaxpcBackshift {
            unit,
            caldb_dir,
            mission_day,
            fits,
            bkg_pha,
            bkg_peak_ch,
            output,
            layer1_only,
            single_only,
            pha_min,
            pha_max,
        } => {
            let gain_path = caldb_dir.join("gaina8.dat");
            let gain_table = GainTable::from_file(&gain_path)
                .with_context(|| format!("failed to load {}", gain_path.display()))?;

            let hv_table = [caldb_dir.join("lxhv.dat"), caldb_dir.join("lxhv.data")]
                .into_iter()
                .find(|p| p.is_file())
                .map(|p| HvJumpTable::from_file(&p))
                .transpose()?;

            let resolved_day = if let Some(md) = mission_day {
                md
            } else if let Some(ref fp) = fits {
                let f = FitsFile::open(fp)?;
                let ts = f
                    .hdu(1)
                    .ok()
                    .and_then(|h| h.header().get_f64("TSTART"))
                    .or_else(|| f.hdu(0).ok().and_then(|h| h.header().get_f64("TSTART")))
                    .context("could not find TSTART in FITS file; pass mission_day")?;
                utc_to_mission_day(ts).max(0.0)
            } else {
                bail!("provide either mission_day or fits for laxpc-backshift");
            };

            let bkg_pha_loaded = bkg_pha
                .as_deref()
                .map(BackgroundPha::from_file)
                .transpose()?;
            let bkg_day = bkg_pha_loaded
                .as_ref()
                .map(|b| b.mission_day)
                .unwrap_or(resolved_day);
            let resolved_bkg_peak_ch = bkg_peak_ch.or_else(|| {
                if bkg_pha_loaded.is_some() {
                    let cat = BackgroundCatalog::load_or_builtin(&[
                        Some(caldb_dir.as_path()),
                        bkg_pha.as_deref().and_then(Path::parent),
                    ]);
                    cat.lookup_peak_30kev(unit, bkg_pha.as_deref(), bkg_day)
                } else {
                    None
                }
            });

            let gain_shift = compute_gain_shift(
                unit,
                resolved_day,
                bkg_day,
                resolved_bkg_peak_ch,
                &gain_table,
                hv_table.as_ref(),
                layer1_only,
                single_only,
            )?;

            let mut outputs = Vec::new();
            let mut lcbk_shift_fraction = None;
            let mut shifted_spectrum_fits: Option<PathBuf> = None;

            if let Some(pha) = &bkg_pha_loaded {
                let cmap = ChannelMap::from_caldb_for_day(&caldb_dir, unit, resolved_day)
                    .or_else(|| {
                        bkg_pha
                            .as_deref()
                            .and_then(Path::parent)
                            .and_then(|d| ChannelMap::from_caldb_for_day(d, unit, resolved_day))
                    })
                    .with_context(|| {
                        format!(
                            "failed to load channel map for {} (day {:.1}) from {}",
                            unit.detnam(),
                            resolved_day,
                            caldb_dir.display()
                        )
                    })?;

                let shifted = shift_background_spectrum(
                    unit,
                    gain_shift.delta_ch_30kev,
                    gain_shift.secular_scale,
                    resolved_day,
                    &cmap,
                    &pha.rates,
                )?;
                let filter = EventFilter {
                    single_only,
                    anodes: if layer1_only {
                        AnodeSelection::Layer1
                    } else {
                        AnodeSelection::AllMain
                    },
                    pha_min,
                    pha_max,
                    ..Default::default()
                };

                let nch = shifted.len().max(1);
                let ncf = (1024 / nch).max(1);
                let ch_lo = (pha_min as usize) / ncf;
                let ch_hi = (pha_max as usize) / ncf;
                let (_cdel, fshbk) = background_lightcurve_shift_fraction(
                    &pha.rates,
                    &shifted,
                    &filter.anodes,
                    ch_lo,
                    ch_hi,
                );
                lcbk_shift_fraction = Some(fshbk);

                if let Some(out_path) = &output {
                    let mut rates = Vec::with_capacity(shifted.len());
                    let mut errs = Vec::with_capacity(shifted.len());
                    for row in &shifted {
                        let s: f64 = (1..=10u8)
                            .filter(|&a| filter.anodes.matches(a))
                            .map(|a| row[(a - 1) as usize])
                            .sum();
                        rates.push(s);
                        errs.push((s / pha.exposure.max(1.0)).max(0.0).sqrt());
                    }
                    let mut meta = OgipMeta {
                        instrume: "LAXPC".into(),
                        detnam: Some(unit.detnam().into()),
                        object: Some("BACKGROUND".into()),
                        exposure: pha.exposure,
                        ..Default::default()
                    };
                    meta.tstart = (resolved_day + 2096.0) * 86_400.0;
                    meta.tstop = meta.tstart + pha.exposure;
                    ensure_parent_dir(out_path)?;
                    write_spectrum_fits(
                        out_path,
                        unit,
                        &rates,
                        &errs,
                        pha.exposure,
                        43e-6,
                        &filter,
                        &gain_shift.recommended_rmf,
                        "NONE",
                        &meta,
                    )?;
                    shifted_spectrum_fits = Some(out_path.clone());
                    if out_path.is_file() {
                        let mut e = make_output_entry("Shifted Background Spectrum", out_path);
                        e.product = "spectrum".to_string();
                        outputs.push(e);
                    }
                }
            }

            let summary = json!({
                "unit": unit.detnam(),
                "mission_day": resolved_day,
                "source_day": gain_shift.source_day,
                "background_day": gain_shift.background_day,
                "source_peak_ch": gain_shift.source_peak_ch,
                "background_peak_ch": gain_shift.background_peak_ch,
                "delta_ch_30kev": gain_shift.delta_ch_30kev,
                "gain_ratio": gain_shift.gain_ratio,
                "secular_scale": gain_shift.secular_scale,
                "crossed_hv_jumps": gain_shift.crossed_hv_jumps,
                "lxp20_gain_offset_kev": gain_shift.lxp20_gain_offset_kev,
                "recommended_rmf": gain_shift.recommended_rmf,
                "lcbk_shift_fraction": lcbk_shift_fraction,
                "shifted_spectrum_fits": shifted_spectrum_fits.as_ref().map(|p| p.to_string_lossy().into_owned()),
            });
            Ok((summary, outputs))
        }
    }
}

fn ensure_parent_dir(path: &Path) -> Result<()> {
    if let Some(parent) = path.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent)?;
        }
    }
    Ok(())
}

fn parse_laxpc_anodes(raw: &str) -> ApiResult<AnodeSelection> {
    let t = raw.trim().to_ascii_lowercase();
    if t == "all" || t == "a1-a7" {
        return Ok(AnodeSelection::AllMain);
    }
    if t == "layer1" || t == "l1" {
        return Ok(AnodeSelection::Layer1);
    }
    if t == "layer2" || t == "l2" {
        return Ok(AnodeSelection::Layer2);
    }
    let mut list = Vec::new();
    for part in t.split(',') {
        let p = part.trim().trim_start_matches('a');
        let id: u8 = p.parse().map_err(|_| {
            ApiError::bad_request(format!("invalid anode ID '{part}' in anodes"))
        })?;
        if !(1..=10).contains(&id) {
            return Err(ApiError::bad_request(format!(
                "anode ID must be in 1..=10, got {id}"
            )));
        }
        if !list.contains(&id) {
            list.push(id);
        }
    }
    if list.is_empty() {
        return Err(ApiError::bad_request("anodes selection cannot be empty"));
    }
    if list.len() == 1 {
        Ok(AnodeSelection::SingleAnode(list[0]))
    } else {
        Ok(AnodeSelection::Custom(list))
    }
}

fn load_user_gti(path: &Path) -> Result<GtiSet> {
    if let Ok(gti) = GtiSet::from_fits(path, None) {
        return Ok(gti);
    }
    let txt = std::fs::read_to_string(path)
        .with_context(|| format!("failed to read GTI file {}", path.display()))?;
    let mut intervals = Vec::new();
    for (idx, line) in txt.lines().enumerate() {
        let t = line.trim();
        if t.is_empty() || t.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = t.split_whitespace().collect();
        if cols.len() < 2 {
            bail!("GTI line {}: expected 'START STOP', got '{}'", idx + 1, t);
        }
        let a: f64 = cols[0]
            .parse()
            .with_context(|| format!("GTI line {} start", idx + 1))?;
        let b: f64 = cols[1]
            .parse()
            .with_context(|| format!("GTI line {} stop", idx + 1))?;
        intervals.push(GtiInterval::new(a, b));
    }
    if intervals.is_empty() {
        bail!("no GTI intervals found in {}", path.display());
    }
    Ok(GtiSet::new(intervals))
}

fn make_output_entry(label: &str, path: &Path) -> JobOutputFile {
    let name = path
        .file_name()
        .map(|s| s.to_string_lossy().to_ascii_lowercase())
        .unwrap_or_default();
    let ext = path
        .extension()
        .and_then(|s| s.to_str())
        .unwrap_or("")
        .to_ascii_lowercase();
    let size_bytes = std::fs::metadata(path).map(|m| m.len()).unwrap_or(0);
    JobOutputFile {
        label: label.to_string(),
        path: path.to_string_lossy().into_owned(),
        product: classify_product(&name, &ext).to_string(),
        size_bytes,
    }
}

fn read_ra_dec_from_fits(path: &Path) -> Result<(Option<f64>, Option<f64>)> {
    let fits = FitsFile::open(path)?;
    let mut ra = None;
    let mut dec = None;
    for h in fits.iter_hdus() {
        let hdr = h.header();
        if ra.is_none() {
            ra = hdr
                .get_f64("RA_OBJ")
                .or_else(|| hdr.get_f64("RA_NOM"))
                .or_else(|| hdr.get_f64("RA_PNT"));
        }
        if dec.is_none() {
            dec = hdr
                .get_f64("DEC_OBJ")
                .or_else(|| hdr.get_f64("DEC_NOM"))
                .or_else(|| hdr.get_f64("DEC_PNT"));
        }
        if ra.is_some() && dec.is_some() {
            break;
        }
    }
    Ok((ra, dec))
}

#[derive(Default)]
struct DiscoveredCztiL1 {
    m0: Option<PathBuf>,
    ss: Option<PathBuf>,
    tct: Option<PathBuf>,
    mkf: Option<PathBuf>,
    mkf_thresholds: Option<PathBuf>,
}

fn discover_czti_l1(root: &Path) -> DiscoveredCztiL1 {
    let mut found = DiscoveredCztiL1::default();
    let mut stack = vec![(root.to_path_buf(), 0usize)];
    while let Some((dir, depth)) = stack.pop() {
        if depth > 4 {
            continue;
        }
        let Ok(rd) = std::fs::read_dir(&dir) else {
            continue;
        };
        let mut entries: Vec<PathBuf> = rd.filter_map(|e| e.ok().map(|x| x.path())).collect();
        entries.sort();
        for p in entries {
            if p.is_dir() {
                stack.push((p, depth + 1));
            } else if p.is_file() {
                let name = p
                    .file_name()
                    .map(|s| s.to_string_lossy().to_ascii_lowercase())
                    .unwrap_or_default();
                if found.m0.is_none() && name.contains("cztm0") && name.ends_with(".fits") {
                    found.m0 = Some(p);
                } else if found.ss.is_none() && name.contains("cztss") && name.ends_with(".fits") {
                    found.ss = Some(p);
                } else if found.tct.is_none() && name.ends_with(".tct") {
                    found.tct = Some(p);
                } else if found.mkf.is_none() && name.ends_with(".mkf") {
                    found.mkf = Some(p);
                } else if found.mkf_thresholds.is_none() && name == "mkfthresholds.txt" {
                    found.mkf_thresholds = Some(p);
                }
            }
        }
    }
    found
}
