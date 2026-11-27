//! Workspace health/preset auto-detection (`GET /api/health`) and safe file browser (`GET /api/browse`).

use std::collections::HashMap;
use std::path::{Path, PathBuf};

use serde::Serialize;

use crate::http_util::{resolve_existing_dir, ApiError, ApiResult};

#[derive(Debug, Clone, Serialize)]
pub struct WorkspacePresets {
    pub czti_l1_dir: Option<String>,
    pub czti_l2_dir: Option<String>,
    pub czti_caldb: Option<String>,
    pub czti_orb: Option<String>,
    pub laxpc_l1_dir: Option<String>,
    pub laxpc_l2_dir: Option<String>,
    pub laxpc_caldb: Option<String>,
    pub laxpc_orb: Option<String>,
    pub jpleph_405: Option<String>,
    pub sample_lightcurves: Vec<String>,
    pub sample_spectra: Vec<String>,
    pub sample_images: Vec<String>,
    pub sample_events: Vec<String>,
}

#[derive(Debug, Clone, Serialize)]
pub struct HealthResponse {
    pub status: &'static str,
    pub version: &'static str,
    pub threads: usize,
    pub data_dir: String,
    pub cwd: String,
    pub presets: WorkspacePresets,
}

#[derive(Debug, Clone, Serialize)]
pub struct BrowseEntry {
    pub name: String,
    pub path: String,
    pub is_dir: bool,
    pub size_bytes: u64,
    pub ext: String,
    pub kind: String,
    pub product: String,
}

#[derive(Debug, Clone, Serialize)]
pub struct BrowseResponse {
    pub path: String,
    pub parent: Option<String>,
    pub entries: Vec<BrowseEntry>,
}

/// Build `/api/health` response from pre-computed/cached [`WorkspacePresets`].
pub fn handle_health_with_presets(
    data_dir: &Path,
    cwd: &Path,
    threads: usize,
    presets: WorkspacePresets,
) -> HealthResponse {
    HealthResponse {
        status: "ok",
        version: env!("CARGO_PKG_VERSION"),
        threads,
        data_dir: data_dir.to_string_lossy().into_owned(),
        cwd: cwd.to_string_lossy().into_owned(),
        presets,
    }
}

/// Build `/api/health` response with auto-detected sample datasets, CALDB, and JPL DE405 ephemeris.
pub fn handle_health(data_dir: &Path, threads: usize) -> HealthResponse {
    let cwd = std::env::current_dir().unwrap_or_else(|_| data_dir.to_path_buf());
    let presets = detect_presets(data_dir, &cwd);
    handle_health_with_presets(data_dir, &cwd, threads, presets)
}

/// Handle `GET /api/browse?path=...`.
pub fn handle_browse(
    params: &HashMap<String, String>,
    data_dir: &Path,
) -> ApiResult<BrowseResponse> {
    let raw_path = params.get("path").map(|s| s.as_str());
    let dir = resolve_existing_dir(raw_path, data_dir)?;
    let read_dir = std::fs::read_dir(&dir)
        .map_err(|e| ApiError::bad_request(format!("cannot read directory {}: {e}", dir.display())))?;

    let mut entries = Vec::new();
    for entry_res in read_dir {
        let Ok(entry) = entry_res else { continue };
        let name = entry.file_name().to_string_lossy().into_owned();
        if name.starts_with('.') {
            continue;
        }
        let path = entry.path();
        let Ok(meta) = entry.metadata() else { continue };
        let is_dir = meta.is_dir();
        let size_bytes = if is_dir { 0 } else { meta.len() };
        let ext = if is_dir {
            String::new()
        } else {
            path.extension()
                .and_then(|e| e.to_str())
                .unwrap_or("")
                .to_ascii_lowercase()
        };
        let (kind, product) = if is_dir {
            ("dir".to_string(), "other".to_string())
        } else {
            classify_astrosat_file(&path, &name, &ext)
        };
        entries.push(BrowseEntry {
            name,
            path: path.to_string_lossy().into_owned(),
            is_dir,
            size_bytes,
            ext,
            kind,
            product,
        });
    }

    entries.sort_by(|a, b| {
        b.is_dir
            .cmp(&a.is_dir)
            .then_with(|| a.name.to_ascii_lowercase().cmp(&b.name.to_ascii_lowercase()))
    });

    let parent = dir
        .parent()
        .map(|p| p.to_string_lossy().into_owned())
        .filter(|p| !p.is_empty());

    Ok(BrowseResponse {
        path: dir.to_string_lossy().into_owned(),
        parent,
        entries,
    })
}

/// Classify a file into an AstroSat `kind` and UI `product` hint.
pub fn classify_astrosat_file(path: &Path, name: &str, ext: &str) -> (String, String) {
    let lower_name = name.to_ascii_lowercase();
    let lower_path = path.to_string_lossy().to_ascii_lowercase();

    let product = classify_product(&lower_name, ext).to_string();

    if lower_name.starts_with("jpleph.") || lower_name.starts_with("lnxp") || lower_path.contains("/ephem/") || lower_path.contains("\\ephem\\") {
        return ("ephem".to_string(), "other".to_string());
    }
    if lower_path.contains("/caldb/") || lower_path.contains("\\caldb\\") || lower_name == "caldb.indx" || ext == "rmf" || ext == "channel" {
        return ("caldb".to_string(), product);
    }
    if ext == "orb" || lower_name.ends_with(".orb") {
        return ("orbit".to_string(), "orbit".to_string());
    }

    let is_czti = lower_name.contains("czt") || lower_path.contains("czti");
    let is_laxpc = lower_name.contains("lxp")
        || lower_name.contains("laxpc")
        || lower_path.contains("laxpc");
    let is_l1 = lower_name.contains("level1") || lower_path.contains("level1") || lower_path.contains("_l1");
    let is_l2 = lower_name.contains("level2")
        || lower_path.contains("level2")
        || lower_path.contains("_l2")
        || matches!(
            ext,
            "lc" | "lcurv" | "lcbk" | "pha" | "spec" | "dph" | "dpi" | "rsp" | "evn" | "evt" | "gtis" | "gti"
        )
        || lower_name.contains(".img_");

    let kind = if is_czti {
        if is_l2 && !is_l1 {
            "czti_l2"
        } else if is_l1 {
            "czti_l1"
        } else {
            "czti_l2"
        }
    } else if is_laxpc {
        if is_l2 && !is_l1 {
            "laxpc_l2"
        } else if is_l1 {
            "laxpc_l1"
        } else {
            "laxpc_l2"
        }
    } else if is_fits_like(ext, &lower_name) {
        "fits"
    } else {
        "other"
    };

    (kind.to_string(), product)
}

/// Infer the UI product role for a filename/extension.
pub fn classify_product(lower_name: &str, ext: &str) -> &'static str {
    if matches!(ext, "lc" | "lcurv" | "lcbk") || lower_name.ends_with("_lc.fits") {
        "lightcurve"
    } else if matches!(ext, "pha" | "spec") || lower_name.ends_with("_pha.fits") {
        "spectrum"
    } else if matches!(ext, "dph" | "dpi" | "img")
        || lower_name.contains(".img_")
        || lower_name.ends_with("_img.fits")
    {
        "image"
    } else if matches!(ext, "evt" | "evn" | "dblevt")
        || lower_name.contains("event")
        || lower_name.contains("_bary")
    {
        "events"
    } else if matches!(ext, "gti" | "gtis") {
        "gti"
    } else if ext == "orb" {
        "orbit"
    } else if matches!(ext, "rsp" | "rmf" | "arf") {
        "response"
    } else {
        "other"
    }
}

fn is_fits_like(ext: &str, lower_name: &str) -> bool {
    matches!(
        ext,
        "fits"
            | "fit"
            | "fts"
            | "evt"
            | "evn"
            | "dblevt"
            | "lc"
            | "lcurv"
            | "lcbk"
            | "pha"
            | "spec"
            | "dph"
            | "dpi"
            | "rsp"
            | "rmf"
            | "arf"
            | "gti"
            | "gtis"
            | "orb"
            | "mkf"
            | "tct"
            | "badpix"
            | "livetime"
    ) || lower_name.contains(".img_")
        || lower_name.ends_with(".fits.gz")
}

/// Auto-detect AstroSat sample workspaces, CALDBs, ephemeris, and sample science products.
pub fn detect_presets(data_dir: &Path, cwd: &Path) -> WorkspacePresets {
    let mut roots = vec![data_dir.to_path_buf()];
    let data_sub = data_dir.join("data");
    if data_sub != data_dir {
        roots.push(data_sub);
    }
    if cwd != data_dir {
        roots.push(cwd.to_path_buf());
        roots.push(cwd.join("data"));
    }
    if let Some(parent) = data_dir.parent() {
        if !roots.iter().any(|r| r == parent) {
            roots.push(parent.to_path_buf());
        }
    }

    let find_dir = |rels: &[&str]| -> Option<PathBuf> {
        for root in &roots {
            for rel in rels {
                let p = root.join(rel);
                if p.is_dir() {
                    return Some(p);
                }
            }
        }
        None
    };

    let find_file = |rels: &[&str]| -> Option<PathBuf> {
        for root in &roots {
            for rel in rels {
                let p = root.join(rel);
                if p.is_file() {
                    return Some(p);
                }
            }
        }
        None
    };

    // Run independent filesystem probes concurrently to overlap WSL/network filesystem stat latency.
    let (
        czti_l1,
        czti_l2,
        czti_caldb,
        czti_orb,
        laxpc_l1,
        laxpc_l2,
        laxpc_caldb,
        laxpc_orb,
        jpleph_405,
    ) = std::thread::scope(|s| {
        let h_czti = s.spawn(|| {
            let l1 = find_dir(&[
                "sample/czti/czti_sample_data/20160331_T01_112T01_9000000406_level1",
                "czti/czti_sample_data/20160331_T01_112T01_9000000406_level1",
            ]);
            let l2 = find_dir(&[
                "sample/czti/czti_sample_data/20160331_T01_112T01_9000000406_level2",
                "czti/czti_sample_data/20160331_T01_112T01_9000000406_level2",
            ]);
            let caldb = find_dir(&["caldb/czti"]);
            let orb = find_file(&[
                "sample/czti_orb/AS1T01_112T01_9000000406czt_level1.orb",
                "czti_orb/AS1T01_112T01_9000000406czt_level1.orb",
            ]);
            (l1, l2, caldb, orb)
        });
        let h_laxpc = s.spawn(|| {
            let l1 = find_dir(&[
                "sample/laxpc_l1/20160107_G02_010T01_9000000252_level1",
                "laxpc_l1/20160107_G02_010T01_9000000252_level1",
            ]);
            let l2 = find_dir(&[
                "sample/laxpc_l2/20160107_G02_010T01_9000000252_level2",
                "laxpc_l2/20160107_G02_010T01_9000000252_level2",
            ]);
            let caldb = find_dir(&["caldb/laxpc", "ref_src/laxpc"]);
            let orb = find_file(&[
                "sample/laxpc_l1/20160107_G02_010T01_9000000252_level1/laxpc/01509/aux/aux1/AS1G02_010T01_9000000252lxp_level1.orb",
                "laxpc_l1/20160107_G02_010T01_9000000252_level1/laxpc/01509/aux/aux1/AS1G02_010T01_9000000252lxp_level1.orb",
            ]);
            (l1, l2, caldb, orb)
        });
        let ephem = find_file(&["ephem/JPLEPH.405", "ephem/lnxp1600p2200.405"]);
        let (cz_l1, cz_l2, cz_caldb, cz_orb) = h_czti.join().unwrap_or((None, None, None, None));
        let (lx_l1, lx_l2, lx_caldb, lx_orb) = h_laxpc.join().unwrap_or((None, None, None, None));
        (cz_l1, cz_l2, cz_caldb, cz_orb, lx_l1, lx_l2, lx_caldb, lx_orb, ephem)
    });

    let (
        mut sample_lightcurves,
        mut sample_spectra,
        mut sample_images,
        mut sample_events,
    ) = (Vec::new(), Vec::new(), Vec::new(), Vec::new());

    // Fast targeted single-level scans in parallel across CZTI modeM0 and LAXPC modeEA.
    let (cz_samples, lx_samples) = std::thread::scope(|s| {
        let h_cz = s.spawn(|| {
            let (mut lcs, mut specs, mut imgs, mut evts) =
                (Vec::new(), Vec::new(), Vec::new(), Vec::new());
            if let Some(ref cz_dir) = czti_l2 {
                let mode_m0 = cz_dir.join("modeM0");
                let target = if mode_m0.is_dir() { &mode_m0 } else { cz_dir };
                collect_samples_flat(target, &mut lcs, &mut specs, &mut imgs, &mut evts);
            }
            (lcs, specs, imgs, evts)
        });
        let h_lx = s.spawn(|| {
            let (mut lcs, mut specs, mut imgs, mut evts) =
                (Vec::new(), Vec::new(), Vec::new(), Vec::new());
            if let Some(ref lx_dir) = laxpc_l2 {
                let lx_sub = lx_dir.join("laxpc");
                if lx_sub.is_dir() {
                    if let Ok(rd) = std::fs::read_dir(&lx_sub) {
                        let mut orbits: Vec<PathBuf> = rd
                            .filter_map(|e| e.ok().map(|x| x.path()))
                            .filter(|p| p.is_dir())
                            .collect();
                        orbits.sort();
                        if let Some(first_orbit) = orbits.first() {
                            for unit_dir in ["lxp1", "lxp2", "lxp3"] {
                                let ea = first_orbit.join(unit_dir).join("modeEA");
                                if ea.is_dir() {
                                    collect_samples_flat(
                                        &ea, &mut lcs, &mut specs, &mut imgs, &mut evts,
                                    );
                                }
                            }
                        }
                    }
                } else {
                    collect_samples_flat(lx_dir, &mut lcs, &mut specs, &mut imgs, &mut evts);
                }
            }
            (lcs, specs, imgs, evts)
        });
        (
            h_cz.join().unwrap_or_default(),
            h_lx.join().unwrap_or_default(),
        )
    });

    for (lcs, specs, imgs, evts) in [cz_samples, lx_samples] {
        for x in lcs {
            if sample_lightcurves.len() < 8 {
                sample_lightcurves.push(x);
            }
        }
        for x in specs {
            if sample_spectra.len() < 8 {
                sample_spectra.push(x);
            }
        }
        for x in imgs {
            if sample_images.len() < 12 {
                sample_images.push(x);
            }
        }
        for x in evts {
            if sample_events.len() < 8 {
                sample_events.push(x);
            }
        }
    }

    let pstr = |p: Option<PathBuf>| p.map(|x| x.to_string_lossy().into_owned());

    WorkspacePresets {
        czti_l1_dir: pstr(czti_l1),
        czti_l2_dir: pstr(czti_l2),
        czti_caldb: pstr(czti_caldb),
        czti_orb: pstr(czti_orb),
        laxpc_l1_dir: pstr(laxpc_l1),
        laxpc_l2_dir: pstr(laxpc_l2),
        laxpc_caldb: pstr(laxpc_caldb),
        laxpc_orb: pstr(laxpc_orb),
        jpleph_405: pstr(jpleph_405),
        sample_lightcurves,
        sample_spectra,
        sample_images,
        sample_events,
    }
}

fn collect_samples_flat(
    dir: &Path,
    lcs: &mut Vec<String>,
    specs: &mut Vec<String>,
    imgs: &mut Vec<String>,
    evts: &mut Vec<String>,
) {
    let Ok(rd) = std::fs::read_dir(dir) else {
        return;
    };
    let mut matched: Vec<(String, &'static str)> = Vec::new();
    for entry in rd.flatten() {
        let fname = entry.file_name();
        let name_str = fname.to_string_lossy();
        if name_str.starts_with('.') {
            continue;
        }
        let lower_name = name_str.to_ascii_lowercase();
        let ext = Path::new(&lower_name)
            .extension()
            .and_then(|s| s.to_str())
            .unwrap_or("");
        let prod = classify_product(&lower_name, ext);
        if !matches!(prod, "lightcurve" | "spectrum" | "image" | "events") {
            continue;
        }
        matched.push((entry.path().to_string_lossy().into_owned(), prod));
    }
    matched.sort_by(|a, b| a.0.cmp(&b.0));
    for (s, prod) in matched {
        match prod {
            "lightcurve" if lcs.len() < 8 => lcs.push(s),
            "spectrum" if specs.len() < 8 => specs.push(s),
            "image" if imgs.len() < 12 => imgs.push(s),
            "events" if evts.len() < 8 => evts.push(s),
            _ => {}
        }
    }
}
