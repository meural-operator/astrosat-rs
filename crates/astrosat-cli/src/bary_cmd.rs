//! `astrosat bary` — Solar System Barycentric timing correction for CZTI and LAXPC Format A & B.

use std::path::{Path, PathBuf};

use anyhow::{bail, Context, Result};
use astrosat_core::{barycenter_fits, BaryFitsOptions, CoordFrame, EphemerisKind};
use astrosat_fits::{FitsFile, OgipMeta};
use astrosat_laxpc::gti::load_merged_mkf_files;
use astrosat_laxpc::pipeline::write_orbit_from_mkf;
use clap::ArgMatches;
use serde::Serialize;
use tempfile::NamedTempFile;

use crate::cli::{is_cli_arg, resolve_flag, BaryArgs, FrameArg};
use crate::config::BaryFileConfig;

#[derive(Debug, Clone, Serialize)]
pub struct BaryReport {
    pub input: PathBuf,
    pub output: PathBuf,
    pub orbit: PathBuf,
    pub ra_deg: f64,
    pub dec_deg: f64,
    pub frame: String,
    pub plephem: String,
    pub hdus_updated: usize,
    pub events_corrected: usize,
    pub gti_rows_corrected: usize,
    pub tstart_tdb: Option<f64>,
    pub tstop_tdb: Option<f64>,
}

pub fn run_bary(args: &BaryArgs, cfg: &BaryFileConfig, matches: &ArgMatches) -> Result<BaryReport> {
    let input = args
        .input
        .clone()
        .or_else(|| args.input_flag.clone())
        .context("missing input FITS file (pass positional <INPUT_FITS> or --input <FITS>)")?;

    if !input.is_file() {
        bail!("input FITS file does not exist: {}", input.display());
    }

    let (hdr_ra, hdr_dec) = read_ra_dec_from_fits(&input)?;
    let ra_deg = args
        .ra
        .or(cfg.ra_deg)
        .or(hdr_ra)
        .context("could not determine target RA_OBJ/RA_NOM from input FITS; pass --ra <deg>")?;
    let dec_deg = args
        .dec
        .or(cfg.dec_deg)
        .or(hdr_dec)
        .context("could not determine target DEC_OBJ/DEC_NOM from input FITS; pass --dec <deg>")?;

    let frame = if is_cli_arg(matches, "frame") {
        match args.frame {
            FrameArg::Icrs => CoordFrame::Icrs,
            FrameArg::Fk5 => CoordFrame::Fk5,
        }
    } else if let Some(ref frame_str) = cfg.frame {
        if frame_str.eq_ignore_ascii_case("fk5") {
            CoordFrame::Fk5
        } else if frame_str.eq_ignore_ascii_case("icrs") {
            CoordFrame::Icrs
        } else {
            bail!(
                "invalid bary.frame '{}' in config (expected 'icrs' or 'fk5')",
                frame_str
            );
        }
    } else {
        CoordFrame::Icrs
    };

    let ephem_arg = if is_cli_arg(matches, "ephem") {
        args.ephem.as_str()
    } else {
        cfg.ephemeris.as_deref().unwrap_or(args.ephem.as_str())
    };
    let ephemeris = if ephem_arg.eq_ignore_ascii_case("analytical") {
        EphemerisKind::Analytical
    } else {
        let p = PathBuf::from(ephem_arg);
        if !p.is_file() {
            bail!("JPL DE ephemeris file not found: {}", p.display());
        }
        EphemerisKind::JplDe(p)
    };

    let write_checksums = resolve_flag(
        matches,
        "checksum",
        "no_checksum",
        cfg.write_checksums,
        true,
    );

    // Resolve orbit file: explicit --orbit -> explicit --mkf -> config -> sibling AstroSat.orb
    let mut temp_orb: Option<NamedTempFile> = None;
    let (orbit_path, orbit_display) = if let Some(ref orb) = args.orbit {
        (orb.clone(), orb.clone())
    } else if let Some(ref mkf_path) = args.mkf {
        let mkf = load_merged_mkf_files(std::slice::from_ref(mkf_path))
            .with_context(|| format!("failed to load MKF {}", mkf_path.display()))?;
        let tmp = NamedTempFile::new().context("failed to create temporary orbit FITS")?;
        let mut meta = OgipMeta::default();
        if let Some((t0, t1)) = mkf.times().first().zip(mkf.times().last()) {
            meta.tstart = *t0;
            meta.tstop = *t1;
        }
        write_orbit_from_mkf(tmp.path(), &mkf, &meta)
            .context("failed to synthesize orbit FITS from MKF")?;
        let path_owned = tmp.path().to_path_buf();
        temp_orb = Some(tmp);
        (path_owned, mkf_path.clone())
    } else if let Some(ref orb) = cfg.orbit_file {
        (orb.clone(), orb.clone())
    } else if let Some(sibling) = find_sibling_orbit(&input) {
        (sibling.clone(), sibling)
    } else {
        bail!(
            "missing AstroSat orbit file; pass --orbit <orb.fits> or --mkf <mkf.fits> (no sibling AstroSat.orb found near {})",
            input.display()
        );
    };

    if let Some(parent) = args.output.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent)?;
        }
    }

    let opts = BaryFitsOptions {
        frame,
        ephemeris,
        write_checksums,
    };

    let summary = barycenter_fits(&input, &args.output, &orbit_path, ra_deg, dec_deg, &opts)
        .with_context(|| format!("barycentric correction failed for {}", input.display()))?;
    drop(temp_orb);

    Ok(BaryReport {
        input,
        output: args.output.clone(),
        orbit: orbit_display,
        ra_deg,
        dec_deg,
        frame: match frame {
            CoordFrame::Icrs => "ICRS".to_string(),
            CoordFrame::Fk5 => "FK5".to_string(),
        },
        plephem: summary.plephem,
        hdus_updated: summary.hdus_updated,
        events_corrected: summary.events_corrected,
        gti_rows_corrected: summary.gti_rows_corrected,
        tstart_tdb: summary.tstart_tdb,
        tstop_tdb: summary.tstop_tdb,
    })
}

fn find_sibling_orbit(input: &Path) -> Option<PathBuf> {
    let dir = input.parent()?;
    let candidate = dir.join("AstroSat.orb");
    candidate.is_file().then_some(candidate)
}

fn read_ra_dec_from_fits(path: &Path) -> Result<(Option<f64>, Option<f64>)> {
    let fits = FitsFile::open(path)
        .with_context(|| format!("failed to open FITS file {}", path.display()))?;
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

pub fn print_bary_human(r: &BaryReport) {
    println!(
        "Barycentric Correction Complete: {} -> {}",
        r.input.display(),
        r.output.display()
    );
    println!(
        "  Target / Frame     : RA={:.6} deg, DEC={:.6} deg ({})",
        r.ra_deg, r.dec_deg, r.frame
    );
    println!(
        "  Orbit / Ephemeris  : {} ({})",
        r.orbit.display(),
        r.plephem
    );
    println!(
        "  Corrected          : {} events, {} GTI rows across {} HDUs",
        r.events_corrected, r.gti_rows_corrected, r.hdus_updated
    );
    if let (Some(ts), Some(te)) = (r.tstart_tdb, r.tstop_tdb) {
        println!("  TDB Window         : {ts:.6} .. {te:.6} s");
    }
}
