//! Lightweight, multi-threaded `tiny_http` server and REST router for `astrosat-ui`.

use std::net::SocketAddr;
use std::panic::{catch_unwind, AssertUnwindSafe};
use std::path::PathBuf;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, OnceLock};
use std::thread::JoinHandle;
use std::time::Duration;

use anyhow::{anyhow, Result};
use tiny_http::{Method, Request, Response, Server, StatusCode};

use crate::browse::{
    detect_presets, handle_browse, handle_health_with_presets, WorkspacePresets,
};
use crate::config_api::handle_config_parse;
use crate::fits_api::{handle_inspect, handle_validate};
use crate::http_util::{
    json_error, json_ok, json_status, parse_url, read_json_body, static_response, ApiError,
};
use crate::jobs::{build_cli_preview, JobManager, PipelineJobRequest};
use crate::pulse::{
    handle_pulse_multiband, handle_pulse_search, MultibandRequest, PulseSearchRequest,
};
use crate::science::{handle_image, handle_lightcurve, handle_spectrum};
use crate::UiConfig;

const ASSET_INDEX_HTML: &str = include_str!("../assets/index.html");
const ASSET_STYLES_CSS: &str = include_str!("../assets/styles.css");
const ASSET_CHARTS_JS: &str = include_str!("../assets/charts.js");
const ASSET_APP_JS: &str = include_str!("../assets/app.js");
const ASSET_STUDIO_INDEX_HTML: &str = include_str!("../assets/studio/index.html");
const ASSET_STUDIO_CSS: &str = include_str!("../assets/studio/studio.css");
const ASSET_STUDIO_JS: &str = include_str!("../assets/studio/studio.js");

struct AppState {
    data_dir: PathBuf,
    cwd: PathBuf,
    threads: usize,
    presets: OnceLock<WorkspacePresets>,
    jobs: Arc<JobManager>,
}

impl AppState {
    fn get_presets(&self) -> &WorkspacePresets {
        self.presets
            .get_or_init(|| detect_presets(&self.data_dir, &self.cwd))
    }
}

/// Running background HTTP server handle (primarily for integration tests and embedded callers).
pub struct UiServer {
    addr: SocketAddr,
    data_dir: PathBuf,
    threads: usize,
    shutdown: Arc<AtomicBool>,
    server: Arc<Server>,
    workers: Vec<JoinHandle<()>>,
}

impl UiServer {
    /// Start the HTTP server on background worker threads and return immediately with its bound address.
    pub fn spawn(cfg: UiConfig) -> Result<Self> {
        let threads = configure_threads(cfg.threads);
        let data_dir = resolve_initial_data_dir(&cfg.data_dir);
        let cwd = std::env::current_dir().unwrap_or_else(|_| data_dir.clone());
        let bind_addr = format!("{}:{}", cfg.host, cfg.port);
        let server = Server::http(&bind_addr)
            .map_err(|e| anyhow!("failed to bind astrosat-ui server on {bind_addr}: {e}"))?;
        let addr = server
            .server_addr()
            .to_ip()
            .ok_or_else(|| anyhow!("server did not bind to an IP socket address"))?;
        let server = Arc::new(server);
        let shutdown = Arc::new(AtomicBool::new(false));
        let state = Arc::new(AppState {
            data_dir: data_dir.clone(),
            cwd,
            threads,
            presets: OnceLock::new(),
            jobs: Arc::new(JobManager::new()),
        });

        // Pre-warm workspace presets asynchronously right after socket bind so startup stays
        // sub-50ms while `/api/health` serves cached results in < 1 ms.
        {
            let st_warm = Arc::clone(&state);
            let _ = std::thread::Builder::new()
                .name("astrosat-ui-presets".into())
                .spawn(move || {
                    let _ = st_warm.get_presets();
                });
        }

        // Use a modest HTTP worker pool so concurrent API requests (e.g. job polling + plotting)
        // never head-of-line block each other.
        let num_workers = std::thread::available_parallelism()
            .map(|n| n.get().clamp(4, 8))
            .unwrap_or(4);
        let mut workers = Vec::with_capacity(num_workers);
        for w in 0..num_workers {
            let srv = Arc::clone(&server);
            let sd = Arc::clone(&shutdown);
            let st = Arc::clone(&state);
            let handle = std::thread::Builder::new()
                .name(format!("astrosat-ui-http-{w}"))
                .spawn(move || worker_loop(srv, sd, st))
                .map_err(|e| anyhow!("failed to spawn HTTP worker thread {w}: {e}"))?;
            workers.push(handle);
        }

        Ok(Self {
            addr,
            data_dir,
            threads,
            shutdown,
            server,
            workers,
        })
    }

    pub fn addr(&self) -> SocketAddr {
        self.addr
    }

    /// Gracefully shut down the server and join worker threads.
    pub fn stop(mut self) {
        self.shutdown.store(true, Ordering::SeqCst);
        self.server.unblock();
        for w in self.workers.drain(..) {
            let _ = w.join();
        }
    }
}

impl Drop for UiServer {
    fn drop(&mut self) {
        self.shutdown.store(true, Ordering::SeqCst);
        self.server.unblock();
        for w in self.workers.drain(..) {
            let _ = w.join();
        }
    }
}

/// Start the local AstroSat web UI server and block serving requests until interrupted.
pub fn serve(cfg: UiConfig) -> Result<()> {
    let open_browser = cfg.open_browser;
    let display_host = if cfg.host == "0.0.0.0" {
        "127.0.0.1".to_string()
    } else {
        cfg.host.clone()
    };

    let handle = UiServer::spawn(cfg)?;
    let url = format!("http://{display_host}:{}", handle.addr().port());
    eprintln!(
        "astrosat ui listening on {url} (data_dir={}, threads={})",
        handle.data_dir.display(),
        handle.threads
    );
    eprintln!("Press Ctrl+C to stop.");

    if open_browser {
        try_open_browser(&url);
    }

    // Wait on workers indefinitely (until process receives SIGINT/SIGTERM).
    for w in handle.workers.iter() {
        let _ = w.thread().id();
    }
    loop {
        std::thread::sleep(Duration::from_secs(3600));
    }
}

fn configure_threads(req: Option<usize>) -> usize {
    if let Some(n) = req.filter(|&x| x > 0) {
        let _ = rayon::ThreadPoolBuilder::new().num_threads(n).build_global();
        return n;
    }
    // Avoid eagerly spawning the global Rayon pool at HTTP server startup;
    // Rayon initializes lazily on the first parallel compute request.
    std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1)
}

fn resolve_initial_data_dir(dir: &std::path::Path) -> PathBuf {
    if dir.as_os_str().is_empty() {
        std::env::current_dir().unwrap_or_else(|_| PathBuf::from("."))
    } else {
        std::fs::canonicalize(dir).unwrap_or_else(|_| dir.to_path_buf())
    }
}

fn worker_loop(server: Arc<Server>, shutdown: Arc<AtomicBool>, state: Arc<AppState>) {
    while !shutdown.load(Ordering::Relaxed) {
        match server.recv_timeout(Duration::from_millis(250)) {
            Ok(Some(req)) => handle_request_safe(req, &state),
            Ok(None) => continue,
            Err(_) => {
                if shutdown.load(Ordering::Relaxed) {
                    break;
                }
            }
        }
    }
}

fn handle_request_safe(mut req: Request, state: &AppState) {
    let resp = match catch_unwind(AssertUnwindSafe(|| route_request(&mut req, state))) {
        Ok(r) => r,
        Err(_) => json_error(&ApiError::internal(
            "internal server error: request handler panicked",
        )),
    };
    let _ = req.respond(resp);
}

fn route_request(req: &mut Request, state: &AppState) -> Response<std::io::Cursor<Vec<u8>>> {
    let method = req.method().clone();
    let (path, params) = parse_url(req.url());

    if method == Method::Options {
        return Response::from_data(Vec::new()).with_status_code(StatusCode(204));
    }

    match (method, path.as_str()) {
        // Static embedded assets (100% offline, zero CDN)
        (Method::Get, "/" | "/index.html" | "/studio" | "/studio/" | "/studio/index.html") => {
            static_response("text/html; charset=utf-8", ASSET_STUDIO_INDEX_HTML)
        }
        (Method::Get, "/studio/studio.css") => {
            static_response("text/css; charset=utf-8", ASSET_STUDIO_CSS)
        }
        (Method::Get, "/studio/studio.js") => {
            static_response("application/javascript; charset=utf-8", ASSET_STUDIO_JS)
        }
        (Method::Get, "/classic" | "/classic/" | "/classic/index.html") => {
            static_response("text/html; charset=utf-8", ASSET_INDEX_HTML)
        }
        (Method::Get, "/styles.css" | "/classic/styles.css") => {
            static_response("text/css; charset=utf-8", ASSET_STYLES_CSS)
        }
        (Method::Get, "/charts.js" | "/classic/charts.js") => {
            static_response("application/javascript; charset=utf-8", ASSET_CHARTS_JS)
        }
        (Method::Get, "/app.js" | "/classic/app.js") => {
            static_response("application/javascript; charset=utf-8", ASSET_APP_JS)
        }
        (Method::Get, "/favicon.ico") => {
            Response::from_data(Vec::new()).with_status_code(StatusCode(204))
        }

        // Health & workspace browser
        (Method::Get, "/api/health") => {
            let h = handle_health_with_presets(
                &state.data_dir,
                &state.cwd,
                state.threads,
                state.get_presets().clone(),
            );
            json_ok(&h)
        }
        (Method::Get, "/api/browse") => match handle_browse(&params, &state.data_dir) {
            Ok(v) => json_ok(&v),
            Err(e) => json_error(&e),
        },

        // FITS Inspector & Validator
        (Method::Get, "/api/inspect") => match handle_inspect(&params, &state.data_dir) {
            Ok(v) => json_ok(&v),
            Err(e) => json_error(&e),
        },
        (Method::Get, "/api/validate") => match handle_validate(&params, &state.data_dir) {
            Ok(v) => json_ok(&v),
            Err(e) => json_error(&e),
        },

        // Scientific data endpoints
        (Method::Get, "/api/lightcurve") => match handle_lightcurve(&params, &state.data_dir) {
            Ok(v) => json_ok(&v),
            Err(e) => json_error(&e),
        },
        (Method::Get, "/api/spectrum") => match handle_spectrum(&params, &state.data_dir) {
            Ok(v) => json_ok(&v),
            Err(e) => json_error(&e),
        },
        (Method::Get, "/api/image") => match handle_image(&params, &state.data_dir) {
            Ok(v) => json_ok(&v),
            Err(e) => json_error(&e),
        },

        // Multi-band phase folding & pulse search
        (Method::Post, "/api/pulse/search") => match read_json_body::<PulseSearchRequest>(req) {
            Ok(body) => match handle_pulse_search(&body, &state.data_dir) {
                Ok(v) => json_ok(&v),
                Err(e) => json_error(&e),
            },
            Err(e) => json_error(&e),
        },
        (Method::Post, "/api/pulse/multiband") => match read_json_body::<MultibandRequest>(req) {
            Ok(body) => match handle_pulse_multiband(&body, &state.data_dir) {
                Ok(v) => json_ok(&v),
                Err(e) => json_error(&e),
            },
            Err(e) => json_error(&e),
        },

        // Config import/parse
        (Method::Post, "/api/config/parse") => match handle_config_parse(req) {
            Ok(v) => json_ok(&v),
            Err(e) => json_error(&e),
        },

        // Reproducible CLI/TOML preview & background job runner
        (Method::Post, "/api/cli-preview") => match read_json_body::<PipelineJobRequest>(req) {
            Ok(body) => match build_cli_preview(&body) {
                Ok(v) => json_ok(&v),
                Err(e) => json_error(&e),
            },
            Err(e) => json_error(&e),
        },
        (Method::Post, "/api/jobs/run") => match read_json_body::<PipelineJobRequest>(req) {
            Ok(body) => match state.jobs.spawn_job(body, state.data_dir.clone()) {
                Ok(rec) => json_status(202, &rec),
                Err(e) => json_error(&e),
            },
            Err(e) => json_error(&e),
        },
        (Method::Get, "/api/jobs") => json_ok(&state.jobs.list()),
        (Method::Get, p) if p.starts_with("/api/jobs/") && p != "/api/jobs/run" => {
            let id = &p["/api/jobs/".len()..];
            if id.is_empty() || id.contains('/') {
                json_error(&ApiError::not_found("invalid job id"))
            } else if let Some(rec) = state.jobs.get(id) {
                json_ok(&rec)
            } else {
                json_error(&ApiError::not_found(format!("job '{id}' not found")))
            }
        }

        // Method mismatches on known API routes
        (_, p) if p.starts_with("/api/pulse/")
            || p == "/api/config/parse"
            || p == "/api/cli-preview"
            || p == "/api/jobs/run" =>
        {
            json_error(&ApiError::method_not_allowed(format!(
                "endpoint '{p}' requires POST"
            )))
        }
        (_, p) if p.starts_with("/api/") => {
            json_error(&ApiError::not_found(format!("unknown API route '{p}'")))
        }
        _ => json_error(&ApiError::not_found(format!("not found: '{path}'"))),
    }
}

fn try_open_browser(url: &str) {
    #[cfg(target_os = "macos")]
    let _ = std::process::Command::new("open").arg(url).spawn();

    #[cfg(target_os = "windows")]
    let _ = std::process::Command::new("cmd")
        .args(["/C", "start", "", url])
        .spawn();

    #[cfg(all(unix, not(target_os = "macos")))]
    {
        if std::env::var_os("DISPLAY").is_some() || std::env::var_os("WAYLAND_DISPLAY").is_some() {
            let _ = std::process::Command::new("xdg-open").arg(url).spawn();
        }
    }
}
