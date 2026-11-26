//! OGIP / AstroSat FITS product builders (`GTI`, `SPECTRUM`, `LIGHTCURVE`, `RESPONSE`, `ORBIT`, WCS)
//! and structural/OGIP conformance validator.

use crate::bintable::BinTableBuilder;
use crate::card::Header;
use crate::error::Result;
use crate::io::{FitsFile, Hdu, HduKind};

/// Common AstroSat / OGIP header metadata applied across Level-2 science products.
#[derive(Debug, Clone)]
pub struct OgipMeta {
    pub telescop: String,
    pub instrume: String,
    pub detnam: Option<String>,
    pub obs_id: Option<String>,
    pub object: Option<String>,
    pub ra_obj: Option<f64>,
    pub dec_obj: Option<f64>,
    pub tstart: f64,
    pub tstop: f64,
    pub exposure: f64,
    pub ontime: Option<f64>,
    pub livetime: Option<f64>,
    pub mjdrefi: i64,
    pub mjdreff: f64,
    pub timesys: String,
    pub timeref: String,
    pub timeunit: String,
    pub date_obs: Option<String>,
    pub date_end: Option<String>,
    pub creator: String,
    pub origin: String,
}

impl Default for OgipMeta {
    fn default() -> Self {
        Self {
            telescop: "ASTROSAT".to_string(),
            instrume: "CZTI".to_string(),
            detnam: None,
            obs_id: None,
            object: None,
            ra_obj: None,
            dec_obj: None,
            tstart: 0.0,
            tstop: 0.0,
            exposure: 0.0,
            ontime: None,
            livetime: None,
            mjdrefi: 55197,
            mjdreff: 0.0,
            timesys: "UTC".to_string(),
            timeref: "LOCAL".to_string(),
            timeunit: "s".to_string(),
            date_obs: None,
            date_end: None,
            creator: "astrosat-rs 0.1.0".to_string(),
            origin: "astrosat-rs".to_string(),
        }
    }
}

/// Apply standard AstroSat/OGIP timing & observation keywords to a header.
pub fn apply_ogip_meta(hdr: &mut Header, meta: &OgipMeta) {
    hdr.set_str("TELESCOP", &meta.telescop, "mission/satellite name");
    hdr.set_str("INSTRUME", &meta.instrume, "instrument/detector name");
    if let Some(det) = &meta.detnam {
        hdr.set_str("DETNAM", det, "specific detector name");
    }
    if let Some(obs) = &meta.obs_id {
        hdr.set_str("OBS_ID", obs, "observation identifier");
    }
    if let Some(obj) = &meta.object {
        hdr.set_str("OBJECT", obj, "name of observed object");
    }
    if let Some(ra) = meta.ra_obj {
        hdr.set_f64("RA_OBJ", ra, "[deg] Right Ascension of target");
    }
    if let Some(dec) = meta.dec_obj {
        hdr.set_f64("DEC_OBJ", dec, "[deg] Declination of target");
    }
    hdr.set_f64("TSTART", meta.tstart, "start time");
    hdr.set_f64("TSTOP", meta.tstop, "stop time");
    hdr.set_f64("EXPOSURE", meta.exposure, "[s] exposure time");
    if let Some(on) = meta.ontime {
        hdr.set_f64("ONTIME", on, "[s] total good time");
    }
    if let Some(lt) = meta.livetime {
        hdr.set_f64("LIVETIME", lt, "[s] deadtime-corrected exposure");
    }
    hdr.set_i64("MJDREFI", meta.mjdrefi, "integer part of MJD reference");
    hdr.set_f64("MJDREFF", meta.mjdreff, "fractional part of MJD reference");
    hdr.set_str("TIMESYS", &meta.timesys, "time system");
    hdr.set_str("TIMEREF", &meta.timeref, "reference frame for times");
    hdr.set_str("TIMEUNIT", &meta.timeunit, "unit for time keywords");
    hdr.set_f64("TIMEZERO", 0.0, "time zero offset");
    if let Some(d) = &meta.date_obs {
        hdr.set_str("DATE-OBS", d, "start date/time of observation");
    }
    if let Some(d) = &meta.date_end {
        hdr.set_str("DATE-END", d, "end date/time of observation");
    }
    hdr.set_str("CREATOR", &meta.creator, "software that created this file");
    hdr.set_str("ORIGIN", &meta.origin, "origin of FITS file");
}

/// Build an OGIP-compliant `GTI` (Good Time Interval) BINTABLE HDU.
pub fn build_gti_hdu(
    extname: &str,
    starts: &[f64],
    stops: &[f64],
    meta: &OgipMeta,
) -> Result<Hdu> {
    let mut hdu = BinTableBuilder::new(extname)
        .add_f64("START", Some("s"), starts.to_vec())
        .add_f64("STOP", Some("s"), stops.to_vec())
        .header_str("HDUCLASS", "OGIP", "format conforms to OGIP standard")
        .header_str("HDUCLAS1", "GTI", "table contains Good Time Intervals")
        .header_str("HDUCLAS2", "STANDARD", "standard GTI")
        .build()?;
    apply_ogip_meta(&mut hdu.header, meta);
    Ok(hdu)
}

/// Spectrum payload: either integer `COUNTS` or floating-point `RATE` (counts/s).
#[derive(Debug, Clone)]
pub enum SpectrumData {
    Counts(Vec<i32>),
    Rate(Vec<f32>),
}

/// OGIP 1.2.0 `SPECTRUM` extension metadata.
#[derive(Debug, Clone)]
pub struct SpectrumMeta {
    pub chantype: String,
    pub detchans: i64,
    pub tlmin1: i64,
    pub tlmax1: i64,
    pub poisserr: bool,
    pub areascal: f64,
    pub backscal: f64,
    pub corrscal: f64,
    pub sys_err: f64,
    pub backfile: String,
    pub respfile: String,
    pub ancrfile: String,
    pub corrfile: String,
}

impl Default for SpectrumMeta {
    fn default() -> Self {
        Self {
            chantype: "PI".to_string(),
            detchans: 512,
            tlmin1: 0,
            tlmax1: 511,
            poisserr: false,
            areascal: 1.0,
            backscal: 1.0,
            corrscal: 0.0,
            sys_err: 0.0,
            backfile: "NONE".to_string(),
            respfile: "NONE".to_string(),
            ancrfile: "NONE".to_string(),
            corrfile: "NONE".to_string(),
        }
    }
}

/// Build an OGIP 1.2.0 Type-I `SPECTRUM` BINTABLE HDU.
#[allow(clippy::too_many_arguments)]
pub fn build_spectrum_hdu(
    extname: &str,
    channels: &[i32],
    data: SpectrumData,
    stat_err: Option<&[f32]>,
    quality: Option<&[i16]>,
    grouping: Option<&[i16]>,
    meta: &OgipMeta,
    spec: &SpectrumMeta,
) -> Result<Hdu> {
    let mut b = BinTableBuilder::new(extname)
        .add_i32("CHANNEL", None, channels.to_vec())
        .with_last_col_range(spec.tlmin1, spec.tlmax1);

    match data {
        SpectrumData::Counts(c) => {
            b = b.add_i32("COUNTS", Some("count"), c);
        }
        SpectrumData::Rate(r) => {
            b = b.add_f32("RATE", Some("count/s"), r);
        }
    }

    if let Some(err) = stat_err {
        b = b.add_f32("STAT_ERR", Some("count/s"), err.to_vec());
    }
    if let Some(q) = quality {
        b = b.add_i16("QUALITY", None, q.to_vec());
    }
    if let Some(g) = grouping {
        b = b.add_i16("GROUPING", None, g.to_vec());
    }

    let mut hdu = b
        .header_str("HDUCLASS", "OGIP", "format conforms to OGIP standard")
        .header_str("HDUCLAS1", "SPECTRUM", "PHA dataset (OGIP memo OGIP-92-007)")
        .header_str("HDUVERS", "1.2.0", "version of format (OGIP memo OGIP-92-007)")
        .header_str("CHANTYPE", &spec.chantype, "channel type (PHA or PI)")
        .header_i64("DETCHANS", spec.detchans, "total number of detector channels")
        .header_bool("POISSERR", spec.poisserr, "Poissonian errors applicable")
        .header_f64("AREASCAL", spec.areascal, "area scaling factor")
        .header_f64("BACKSCAL", spec.backscal, "background scaling factor")
        .header_f64("CORRSCAL", spec.corrscal, "correction scaling factor")
        .header_f64("SYS_ERR", spec.sys_err, "fractional systematic error")
        .header_str("BACKFILE", &spec.backfile, "associated background filename")
        .header_str("RESPFILE", &spec.respfile, "associated redistribution matrix filename")
        .header_str("ANCRFILE", &spec.ancrfile, "associated ancillary response filename")
        .header_str("CORRFILE", &spec.corrfile, "associated correction filename")
        .build()?;

    apply_ogip_meta(&mut hdu.header, meta);
    Ok(hdu)
}

/// Build an OGIP `RATE` (`LIGHTCURVE`) BINTABLE HDU.
pub fn build_lightcurve_hdu(
    extname: &str,
    time: &[f64],
    rate: &[f32],
    error: &[f32],
    fracexp: Option<&[f32]>,
    timedel: f64,
    meta: &OgipMeta,
) -> Result<Hdu> {
    let mut b = BinTableBuilder::new(extname)
        .add_f64("TIME", Some("s"), time.to_vec())
        .add_f32("RATE", Some("count/s"), rate.to_vec())
        .add_f32("ERROR", Some("count/s"), error.to_vec());

    if let Some(fe) = fracexp {
        b = b.add_f32("FRACEXP", None, fe.to_vec());
    }

    let mut hdu = b
        .header_str("HDUCLASS", "OGIP", "format conforms to OGIP standard")
        .header_str("HDUCLAS1", "LIGHTCURVE", "dataset contains a light curve")
        .header_str("HDUCLAS2", "TOTAL", "total counts/rate")
        .header_str("HDUCLAS3", "RATE", "data stored as count rate")
        .header_f64("TIMEDEL", timedel, "[s] time bin width")
        .header_f64("TIMEPIXR", 0.5, "bin time is bin center")
        .build()?;

    apply_ogip_meta(&mut hdu.header, meta);
    Ok(hdu)
}

/// Sparse OGIP 1.3.0 Redistribution Matrix (`MATRIX` + `EBOUNDS`) inputs.
#[derive(Debug, Clone)]
pub struct ResponseMatrixData {
    pub chantype: String,
    pub detchans: i64,
    pub tlmin4: i64,
    pub energ_lo: Vec<f32>,
    pub energ_hi: Vec<f32>,
    pub n_grp: Vec<i16>,
    pub f_chan: Vec<Vec<i16>>,
    pub n_chan: Vec<Vec<i16>>,
    pub matrix: Vec<Vec<f32>>,
    pub channels: Vec<i16>,
    pub e_min: Vec<f32>,
    pub e_max: Vec<f32>,
    pub lo_thresh: f64,
}

/// Build OGIP 1.3.0 `MATRIX` and `EBOUNDS` HDUs for `.rmf` / `.rsp` files.
pub fn build_response_hdus(
    rmf: &ResponseMatrixData,
    meta: &OgipMeta,
) -> Result<(Hdu, Hdu)> {
    let mut matrix_hdu = BinTableBuilder::new("MATRIX")
        .add_f32("ENERG_LO", Some("keV"), rmf.energ_lo.clone())
        .add_f32("ENERG_HI", Some("keV"), rmf.energ_hi.clone())
        .add_i16("N_GRP", None, rmf.n_grp.clone())
        .add_var_i16("F_CHAN", None, rmf.f_chan.clone())
        .with_last_col_range(rmf.tlmin4, rmf.tlmin4 + rmf.detchans - 1)
        .add_var_i16("N_CHAN", None, rmf.n_chan.clone())
        .add_var_f32("MATRIX", None, rmf.matrix.clone())
        .header_str("HDUCLASS", "OGIP", "format conforms to OGIP standard")
        .header_str("HDUCLAS1", "RESPONSE", "dataset relates to spectral response")
        .header_str("HDUCLAS2", "RSP_MATRIX", "redistribution matrix")
        .header_str("HDUVERS", "1.3.0", "OGIP RMF format version")
        .header_str("RMFVERSN", "1992a", "OGIP classification of FITS format")
        .header_str("CHANTYPE", &rmf.chantype, "channel type")
        .header_i64("DETCHANS", rmf.detchans, "total number of detector channels")
        .header_i64("NUMGRP", rmf.n_grp.iter().map(|&g| g as i64).sum(), "total number of channel subsets")
        .header_i64("NUMELT", rmf.matrix.iter().map(|m| m.len() as i64).sum(), "total number of response elements")
        .header_f64("LO_THRES", rmf.lo_thresh, "lower probability threshold")
        .build()?;
    apply_ogip_meta(&mut matrix_hdu.header, meta);

    let mut ebounds_hdu = BinTableBuilder::new("EBOUNDS")
        .add_i16("CHANNEL", None, rmf.channels.clone())
        .add_f32("E_MIN", Some("keV"), rmf.e_min.clone())
        .add_f32("E_MAX", Some("keV"), rmf.e_max.clone())
        .header_str("HDUCLASS", "OGIP", "format conforms to OGIP standard")
        .header_str("HDUCLAS1", "RESPONSE", "dataset relates to spectral response")
        .header_str("HDUCLAS2", "EBOUNDS", "nominal energy boundaries of channels")
        .header_str("HDUVERS", "1.2.0", "OGIP EBOUNDS format version")
        .header_str("CHANTYPE", &rmf.chantype, "channel type")
        .header_i64("DETCHANS", rmf.detchans, "total number of detector channels")
        .build()?;
    apply_ogip_meta(&mut ebounds_hdu.header, meta);

    Ok((matrix_hdu, ebounds_hdu))
}

/// Build an AstroSat `ORBIT` BINTABLE HDU (`TIME` [s], `X,Y,Z` [km], `VX,VY,VZ` [km/s]).
pub fn build_orbit_hdu(
    time: &[f64],
    pos_km: &[[f64; 3]],
    vel_kms: &[[f64; 3]],
    meta: &OgipMeta,
) -> Result<Hdu> {
    let x: Vec<f64> = pos_km.iter().map(|p| p[0]).collect();
    let y: Vec<f64> = pos_km.iter().map(|p| p[1]).collect();
    let z: Vec<f64> = pos_km.iter().map(|p| p[2]).collect();
    let vx: Vec<f64> = vel_kms.iter().map(|v| v[0]).collect();
    let vy: Vec<f64> = vel_kms.iter().map(|v| v[1]).collect();
    let vz: Vec<f64> = vel_kms.iter().map(|v| v[2]).collect();

    let mut hdu = BinTableBuilder::new("ORBIT")
        .add_f64("TIME", Some("s"), time.to_vec())
        .add_f64("X", Some("km"), x)
        .add_f64("Y", Some("km"), y)
        .add_f64("Z", Some("km"), z)
        .add_f64("VX", Some("km/s"), vx)
        .add_f64("VY", Some("km/s"), vy)
        .add_f64("VZ", Some("km/s"), vz)
        .build()?;
    apply_ogip_meta(&mut hdu.header, meta);
    Ok(hdu)
}

/// Celestial 2D TAN/GNOMONIC WCS parameters for FITS `IMAGE` HDUs.
#[derive(Debug, Clone)]
pub struct WcsParams {
    pub ctype1: String,
    pub ctype2: String,
    pub crpix1: f64,
    pub crpix2: f64,
    pub crval1: f64,
    pub crval2: f64,
    pub cdelt1: f64,
    pub cdelt2: f64,
    pub crota2: f64,
    pub cunit1: String,
    pub cunit2: String,
    pub radesys: String,
    pub equinox: f64,
}

impl Default for WcsParams {
    fn default() -> Self {
        Self {
            ctype1: "RA---TAN".to_string(),
            ctype2: "DEC--TAN".to_string(),
            crpix1: 1.0,
            crpix2: 1.0,
            crval1: 0.0,
            crval2: 0.0,
            cdelt1: -0.05,
            cdelt2: 0.05,
            crota2: 0.0,
            cunit1: "deg".to_string(),
            cunit2: "deg".to_string(),
            radesys: "ICRS".to_string(),
            equinox: 2000.0,
        }
    }
}

/// Write WCS keywords onto an image header.
pub fn apply_wcs(hdr: &mut Header, wcs: &WcsParams) {
    hdr.set_str("CTYPE1", &wcs.ctype1, "coordinate type for axis 1");
    hdr.set_str("CTYPE2", &wcs.ctype2, "coordinate type for axis 2");
    hdr.set_f64("CRPIX1", wcs.crpix1, "reference pixel along axis 1");
    hdr.set_f64("CRPIX2", wcs.crpix2, "reference pixel along axis 2");
    hdr.set_f64("CRVAL1", wcs.crval1, "[deg] coordinate value at reference pixel 1");
    hdr.set_f64("CRVAL2", wcs.crval2, "[deg] coordinate value at reference pixel 2");
    hdr.set_f64("CDELT1", wcs.cdelt1, "[deg/pixel] coordinate increment along axis 1");
    hdr.set_f64("CDELT2", wcs.cdelt2, "[deg/pixel] coordinate increment along axis 2");
    hdr.set_f64("CROTA2", wcs.crota2, "[deg] rotation angle");
    hdr.set_str("CUNIT1", &wcs.cunit1, "units of coordinate increment axis 1");
    hdr.set_str("CUNIT2", &wcs.cunit2, "units of coordinate increment axis 2");
    hdr.set_str("RADESYS", &wcs.radesys, "celestial reference system");
    hdr.set_f64("EQUINOX", wcs.equinox, "equinox of celestial coord system");
}

/// Severity of a FITS/OGIP validation finding.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IssueSeverity {
    Warning,
    Error,
}

/// A single validation issue discovered in a FITS file.
#[derive(Debug, Clone)]
pub struct ValidationIssue {
    pub hdu_index: usize,
    pub severity: IssueSeverity,
    pub message: String,
}

/// Validate FITS 4.0 structural rules, `DATASUM`/`CHECKSUM` integrity, and OGIP conventions.
pub fn validate_ogip(fits: &FitsFile) -> Vec<ValidationIssue> {
    let mut issues = Vec::new();

    for hdu in fits.iter_hdus() {
        let idx = hdu.index();
        let hdr = hdu.header();

        // Checksum verification if present.
        let (d_ok, c_ok) = hdu.verify_checksum();
        if !d_ok {
            issues.push(ValidationIssue {
                hdu_index: idx,
                severity: IssueSeverity::Error,
                message: "DATASUM mismatch".into(),
            });
        }
        if !c_ok {
            issues.push(ValidationIssue {
                hdu_index: idx,
                severity: IssueSeverity::Error,
                message: "CHECKSUM mismatch".into(),
            });
        }

        if *hdu.kind() == HduKind::BinTable {
            if let Err(e) = hdu.as_bintable() {
                issues.push(ValidationIssue {
                    hdu_index: idx,
                    severity: IssueSeverity::Error,
                    message: format!("invalid BINTABLE: {e}"),
                });
                continue;
            }
            let bt = hdu.as_bintable().unwrap();
            let hduclas1 = hdr
                .get_str("HDUCLAS1")
                .map(|s| s.trim().to_ascii_uppercase());

            match hduclas1.as_deref() {
                Some("SPECTRUM") => {
                    for req in ["TELESCOP", "INSTRUME", "EXPOSURE", "CHANTYPE", "DETCHANS", "HDUVERS"] {
                        if hdr.get(req).is_none() {
                            issues.push(ValidationIssue {
                                hdu_index: idx,
                                severity: IssueSeverity::Error,
                                message: format!("OGIP SPECTRUM missing required keyword {req}"),
                            });
                        }
                    }
                    if bt.column_index("CHANNEL").is_none() {
                        issues.push(ValidationIssue {
                            hdu_index: idx,
                            severity: IssueSeverity::Error,
                            message: "OGIP SPECTRUM missing CHANNEL column".into(),
                        });
                    }
                    if bt.column_index("COUNTS").is_none() && bt.column_index("RATE").is_none() {
                        issues.push(ValidationIssue {
                            hdu_index: idx,
                            severity: IssueSeverity::Error,
                            message: "OGIP SPECTRUM missing COUNTS or RATE column".into(),
                        });
                    }
                }
                Some("LIGHTCURVE") => {
                    for req in ["TELESCOP", "INSTRUME", "TIMEDEL"] {
                        if hdr.get(req).is_none() {
                            issues.push(ValidationIssue {
                                hdu_index: idx,
                                severity: IssueSeverity::Error,
                                message: format!("OGIP LIGHTCURVE missing required keyword {req}"),
                            });
                        }
                    }
                    if bt.column_index("TIME").is_none() {
                        issues.push(ValidationIssue {
                            hdu_index: idx,
                            severity: IssueSeverity::Error,
                            message: "OGIP LIGHTCURVE missing TIME column".into(),
                        });
                    }
                }
                Some("GTI") => {
                    if bt.column_index("START").is_none() || bt.column_index("STOP").is_none() {
                        issues.push(ValidationIssue {
                            hdu_index: idx,
                            severity: IssueSeverity::Error,
                            message: "OGIP GTI missing START or STOP column".into(),
                        });
                    }
                }
                Some("RESPONSE") => {
                    let hduclas2 = hdr
                        .get_str("HDUCLAS2")
                        .map(|s| s.trim().to_ascii_uppercase());
                    if hduclas2.as_deref() == Some("RSP_MATRIX") {
                        for col in ["ENERG_LO", "ENERG_HI", "N_GRP", "F_CHAN", "N_CHAN", "MATRIX"] {
                            if bt.column_index(col).is_none() {
                                issues.push(ValidationIssue {
                                    hdu_index: idx,
                                    severity: IssueSeverity::Error,
                                    message: format!("OGIP MATRIX missing {col} column"),
                                });
                            }
                        }
                    } else if hduclas2.as_deref() == Some("EBOUNDS") {
                        for col in ["CHANNEL", "E_MIN", "E_MAX"] {
                            if bt.column_index(col).is_none() {
                                issues.push(ValidationIssue {
                                    hdu_index: idx,
                                    severity: IssueSeverity::Error,
                                    message: format!("OGIP EBOUNDS missing {col} column"),
                                });
                            }
                        }
                    }
                }
                _ => {}
            }
        }
    }

    issues
}
