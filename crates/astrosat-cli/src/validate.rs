//! `astrosat validate <fits>` — built-in FITS 4.0 + CHECKSUM/DATASUM + OGIP schema validator.

use std::path::{Path, PathBuf};

use astrosat_fits::{validate_ogip, FitsFile, HduKind, IssueSeverity, FITS_BLOCK_SIZE};
use serde::Serialize;

#[derive(Debug, Clone, Serialize)]
pub struct ValidationFinding {
    pub hdu_index: usize,
    pub severity: String,
    pub message: String,
}

#[derive(Debug, Clone, Serialize)]
pub struct FileValidation {
    pub path: PathBuf,
    pub valid: bool,
    pub errors: usize,
    pub warnings: usize,
    pub num_hdus: usize,
    pub findings: Vec<ValidationFinding>,
}

#[derive(Debug, Clone, Serialize)]
pub struct ValidateReport {
    pub strict: bool,
    pub passed: bool,
    pub total_errors: usize,
    pub total_warnings: usize,
    pub files: Vec<FileValidation>,
}

pub fn run_validate(paths: &[PathBuf], strict: bool) -> ValidateReport {
    let mut files = Vec::with_capacity(paths.len());
    let mut total_errors = 0usize;
    let mut total_warnings = 0usize;

    for p in paths {
        let fv = validate_single_file(p, strict);
        total_errors += fv.errors;
        total_warnings += fv.warnings;
        files.push(fv);
    }

    let passed = total_errors == 0 && (!strict || total_warnings == 0);
    ValidateReport {
        strict,
        passed,
        total_errors,
        total_warnings,
        files,
    }
}

fn validate_single_file(path: &Path, strict: bool) -> FileValidation {
    let mut findings = Vec::new();

    match std::fs::metadata(path) {
        Ok(meta) => {
            if meta.len() == 0 {
                findings.push(ValidationFinding {
                    hdu_index: 0,
                    severity: "error".to_string(),
                    message: "FITS file is empty (0 bytes)".to_string(),
                });
            } else if !meta.len().is_multiple_of(FITS_BLOCK_SIZE as u64) {
                findings.push(ValidationFinding {
                    hdu_index: 0,
                    severity: "error".to_string(),
                    message: format!(
                        "FITS file size ({} bytes) is not a multiple of 2880 bytes",
                        meta.len()
                    ),
                });
            }
        }
        Err(e) => {
            findings.push(ValidationFinding {
                hdu_index: 0,
                severity: "error".to_string(),
                message: format!("cannot stat file: {e}"),
            });
            return finalize_validation(path, 0, findings, strict);
        }
    }

    let fits = match FitsFile::open(path) {
        Ok(f) => f,
        Err(e) => {
            findings.push(ValidationFinding {
                hdu_index: 0,
                severity: "error".to_string(),
                message: format!("FITS 4.0 parse failure: {e}"),
            });
            return finalize_validation(path, 0, findings, strict);
        }
    };

    // 1. OGIP + CHECKSUM/DATASUM validator from `astrosat-fits`.
    for issue in validate_ogip(&fits) {
        findings.push(ValidationFinding {
            hdu_index: issue.hdu_index,
            severity: match issue.severity {
                IssueSeverity::Error => "error".to_string(),
                IssueSeverity::Warning => "warning".to_string(),
            },
            message: issue.message,
        });
    }

    // 2. Additional FITS 4.0 and AstroSat schema checks.
    for hdu in fits.iter_hdus() {
        let idx = hdu.index();
        let hdr = hdu.header();

        if idx == 0 {
            if hdr.get_bool("SIMPLE") != Some(true) {
                findings.push(ValidationFinding {
                    hdu_index: 0,
                    severity: "error".into(),
                    message: "Primary HDU missing SIMPLE = T".into(),
                });
            }
            match hdr.get_i64("BITPIX") {
                Some(8 | 16 | 32 | 64 | -32 | -64) => {}
                other => findings.push(ValidationFinding {
                    hdu_index: 0,
                    severity: "error".into(),
                    message: format!("Primary HDU invalid BITPIX: {other:?}"),
                }),
            }
            if hdr.get_i64("NAXIS").is_none() {
                findings.push(ValidationFinding {
                    hdu_index: 0,
                    severity: "error".into(),
                    message: "Primary HDU missing NAXIS".into(),
                });
            }
            if fits.len() > 1 && hdr.get_bool("EXTEND") != Some(true) {
                findings.push(ValidationFinding {
                    hdu_index: 0,
                    severity: "warning".into(),
                    message: "File has extensions but Primary HDU lacks EXTEND = T".into(),
                });
            }
        } else if hdr.get_str("EXTNAME").is_none() {
            findings.push(ValidationFinding {
                hdu_index: idx,
                severity: "warning".into(),
                message: "Extension HDU missing EXTNAME keyword".into(),
            });
        }

        if hdr.get_str("CHECKSUM").is_none() || hdr.get_str("DATASUM").is_none() {
            findings.push(ValidationFinding {
                hdu_index: idx,
                severity: "warning".into(),
                message: "HDU missing CHECKSUM/DATASUM integrity keywords".into(),
            });
        }

        if let (Some(ts), Some(te)) = (hdr.get_f64("TSTART"), hdr.get_f64("TSTOP")) {
            if !ts.is_finite() || !te.is_finite() || te < ts {
                findings.push(ValidationFinding {
                    hdu_index: idx,
                    severity: "error".into(),
                    message: format!("Invalid time bounds TSTART={ts} > TSTOP={te}"),
                });
            }
        }

        if *hdu.kind() == HduKind::BinTable {
            if let Ok(bt) = hdu.as_bintable() {
                let ext = hdu
                    .extname()
                    .map(|s| s.trim().to_ascii_uppercase())
                    .unwrap_or_default();
                if ext.contains("GTI") {
                    if let (Ok(starts), Ok(stops)) =
                        (bt.read_col_f64("START"), bt.read_col_f64("STOP"))
                    {
                        let mut prev_stop = f64::NEG_INFINITY;
                        for (r, (&a, &b)) in starts.iter().zip(stops.iter()).enumerate() {
                            if !a.is_finite() || !b.is_finite() || b < a {
                                findings.push(ValidationFinding {
                                    hdu_index: idx,
                                    severity: "error".into(),
                                    message: format!("GTI row {r} has invalid interval [{a}, {b}]"),
                                });
                            } else {
                                if a < prev_stop {
                                    findings.push(ValidationFinding {
                                        hdu_index: idx,
                                        severity: "warning".into(),
                                        message: format!(
                                            "GTI row {r} overlaps or is unsorted ({a} < {prev_stop})"
                                        ),
                                    });
                                }
                                prev_stop = b;
                            }
                        }
                    }
                } else if ext == "ORBIT" {
                    let has_pos = bt.column_index("X").is_some() || bt.column_index("POS_X").is_some();
                    let has_vel = bt.column_index("VX").is_some() || bt.column_index("VEL_X").is_some();
                    if bt.column_index("TIME").is_none() || !has_pos || !has_vel {
                        findings.push(ValidationFinding {
                            hdu_index: idx,
                            severity: "error".into(),
                            message: "ORBIT HDU missing TIME or position/velocity columns".into(),
                        });
                    }
                }
            }
        }
    }

    finalize_validation(path, fits.len(), findings, strict)
}

fn finalize_validation(
    path: &Path,
    num_hdus: usize,
    findings: Vec<ValidationFinding>,
    strict: bool,
) -> FileValidation {
    let errors = findings.iter().filter(|f| f.severity == "error").count();
    let warnings = findings.iter().filter(|f| f.severity == "warning").count();
    let valid = errors == 0 && (!strict || warnings == 0);
    FileValidation {
        path: path.to_path_buf(),
        valid,
        errors,
        warnings,
        num_hdus,
        findings,
    }
}

pub fn print_validate_human(report: &ValidateReport) {
    for f in &report.files {
        let status = if f.valid { "PASS" } else { "FAIL" };
        println!(
            "[{status}] {} ({} HDUs, {} errors, {} warnings)",
            f.path.display(),
            f.num_hdus,
            f.errors,
            f.warnings
        );
        for item in &f.findings {
            println!(
                "  HDU {:>2} [{:<7}] {}",
                item.hdu_index, item.severity, item.message
            );
        }
    }
}
