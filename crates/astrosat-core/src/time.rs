//! AstroSat mission time system (`MJDREFI = 55197`), embedded `tai-utc.dat` leap-second table,
//! UTC -> TAI -> TT -> TDB conversions, and Fairhead & Bretagnon (1990) Einstein term (`TDB - TT`).

use chrono::{Duration, NaiveDate, NaiveDateTime};

/// AstroSat reference MJD integer part (2010-01-01T00:00:00.000 UTC).
pub const ASTROSAT_MJDREFI: i64 = 55197;
/// AstroSat reference MJD fractional part.
pub const ASTROSAT_MJDREFF: f64 = 0.0;
/// Exact offset `TT - TAI` in seconds.
pub const TT_MINUS_TAI: f64 = 32.184;
/// Exact seconds per SI day.
pub const SECS_PER_DAY: f64 = 86_400.0;
/// `JD - MJD`.
pub const JD_MINUS_MJD: f64 = 2_400_000.5;
/// J2000.0 Julian Date (`2000-01-01T12:00:00 TT`).
pub const JD_J2000: f64 = 2_451_545.0;

/// Embedded IERS `tai-utc.dat` leap-second table: `(mjd_utc_start, tai_minus_utc_seconds)`.
/// Covers all integer leap seconds from 1972-01-01 through 2026+.
pub const LEAP_SECONDS_TABLE: &[(f64, f64)] = &[
    (41317.0, 10.0), // 1972-01-01
    (41499.0, 11.0), // 1972-07-01
    (41683.0, 12.0), // 1973-01-01
    (42048.0, 13.0), // 1974-01-01
    (42413.0, 14.0), // 1975-01-01
    (42778.0, 15.0), // 1976-01-01
    (43144.0, 16.0), // 1977-01-01
    (43509.0, 17.0), // 1978-01-01
    (43874.0, 18.0), // 1979-01-01
    (44239.0, 19.0), // 1980-01-01
    (44786.0, 20.0), // 1981-07-01
    (45151.0, 21.0), // 1982-07-01
    (45516.0, 22.0), // 1983-07-01
    (46247.0, 23.0), // 1985-07-01
    (47161.0, 24.0), // 1988-01-01
    (47892.0, 25.0), // 1990-01-01
    (48257.0, 26.0), // 1991-01-01
    (48804.0, 27.0), // 1992-07-01
    (49169.0, 28.0), // 1993-07-01
    (49534.0, 29.0), // 1994-07-01
    (50083.0, 30.0), // 1996-01-01
    (50630.0, 31.0), // 1997-07-01
    (51179.0, 32.0), // 1999-01-01
    (53736.0, 33.0), // 2006-01-01
    (54832.0, 34.0), // 2009-01-01 (active at AstroSat epoch MJD 55197)
    (56109.0, 35.0), // 2012-07-01
    (57204.0, 36.0), // 2015-07-01 (active at AstroSat launch 2015-09-28)
    (57754.0, 37.0), // 2017-01-01 (active through 2026+)
];

/// Return `TAI - UTC` in seconds at a given `mjd_utc`.
pub fn tai_minus_utc_mjd(mjd_utc: f64) -> f64 {
    let mut ls = LEAP_SECONDS_TABLE[0].1;
    for &(mjd_start, val) in LEAP_SECONDS_TABLE {
        if mjd_utc >= mjd_start {
            ls = val;
        } else {
            break;
        }
    }
    ls
}

/// Return `TAI - UTC` in seconds for an AstroSat UTC timestamp (seconds since MJD 55197 UTC).
pub fn tai_minus_utc_astrosat(astrosat_utc_sec: f64) -> f64 {
    tai_minus_utc_mjd(astrosat_utc_to_mjd_utc(astrosat_utc_sec))
}

/// Number of leap seconds inserted between AstroSat epoch (2010-01-01, MJD 55197, TAI-UTC = 34 s)
/// and `astrosat_utc_sec`.
pub fn leap_seconds_since_epoch(astrosat_utc_sec: f64) -> f64 {
    tai_minus_utc_astrosat(astrosat_utc_sec) - 34.0
}

/// Convert AstroSat UTC seconds (since MJD 55197.0 UTC) to `MJD (UTC)`.
#[inline]
pub fn astrosat_utc_to_mjd_utc(t_utc_sec: f64) -> f64 {
    (ASTROSAT_MJDREFI as f64) + ASTROSAT_MJDREFF + t_utc_sec / SECS_PER_DAY
}

/// Convert `MJD (UTC)` to AstroSat UTC seconds (since MJD 55197.0 UTC).
#[inline]
pub fn mjd_utc_to_astrosat_utc(mjd_utc: f64) -> f64 {
    (mjd_utc - (ASTROSAT_MJDREFI as f64) - ASTROSAT_MJDREFF) * SECS_PER_DAY
}

/// Convert AstroSat UTC seconds to `TT` seconds since `MJD 55197.0 TT`.
///
/// `TT = UTC + (TAI - UTC) + 32.184 s`.
#[inline]
pub fn astrosat_utc_to_tt_sec(t_utc_sec: f64) -> f64 {
    let tai_utc = tai_minus_utc_astrosat(t_utc_sec);
    t_utc_sec + tai_utc + TT_MINUS_TAI
}

/// Convert `TT` seconds since `MJD 55197.0 TT` back to AstroSat UTC seconds.
pub fn tt_sec_to_astrosat_utc(t_tt_sec: f64) -> f64 {
    // Two-step fixed-point iteration handles leap-second boundaries cleanly.
    let mut t_utc = t_tt_sec - 36.0 - TT_MINUS_TAI;
    for _ in 0..3 {
        let tai_utc = tai_minus_utc_astrosat(t_utc);
        t_utc = t_tt_sec - tai_utc - TT_MINUS_TAI;
    }
    t_utc
}

/// Convert AstroSat UTC seconds to `MJD (TT)`.
#[inline]
pub fn astrosat_utc_to_mjd_tt(t_utc_sec: f64) -> f64 {
    (ASTROSAT_MJDREFI as f64) + astrosat_utc_to_tt_sec(t_utc_sec) / SECS_PER_DAY
}

/// Convert AstroSat UTC seconds to high-precision split `(jd_hi, jd_lo)` in TT.
///
/// Keeping integer and fractional day parts separate preserves sub-nanosecond precision
/// over multi-decade spans.
#[inline]
pub fn astrosat_utc_to_jd_tt(t_utc_sec: f64) -> (f64, f64) {
    let t_tt = astrosat_utc_to_tt_sec(t_utc_sec);
    let days = (t_tt / SECS_PER_DAY).floor();
    let rem_sec = t_tt - days * SECS_PER_DAY;
    let jd_hi = (ASTROSAT_MJDREFI as f64) + JD_MINUS_MJD + days;
    let jd_lo = rem_sec / SECS_PER_DAY;
    (jd_hi, jd_lo)
}

/// Compute the relativistic Einstein time dilation `TDB - TT` in seconds
/// using the full 127-term Fairhead & Bretagnon (1990) analytical series
/// (exact 1:1 match to `as1bary/ctatv.c` lines 67–203: 93 T0 + 27 T1 + 6 T2 + 1 T3).
///
/// Guaranteed accuracy `< 100 ns` across 1900–2100.
#[allow(clippy::inconsistent_digit_grouping)]
pub fn tdb_minus_tt(jd_tt_hi: f64, jd_tt_lo: f64) -> f64 {
    // Julian millennia since J2000.0 (365250 days).
    let t = ((jd_tt_hi - JD_J2000) + jd_tt_lo) / 365_250.0;
    let tt = t * t;

    // Exact coefficients from Fairhead & Bretagnon (1990) / as1bary/ctatv.c (amplitudes in us,
    // frequencies in rad/millennium, phases in rad).
    const TERMS_T0: &[(f64, f64, f64)] = &[
        // t1 (20 terms)
        (1656.674564, 6283.075943033, 6.240054195),
        (22.417471, 5753.384970095, 4.296977442),
        (13.839792, 12566.151886066, 6.196904410),
        (4.770086, 529.690965095, 0.444401603),
        (4.676740, 6069.776754553, 4.021195093),
        (2.256707, 213.299095438, 5.543113262),
        (1.694205, -3.523118349, 5.025132748),
        (1.554905, 77713.772618729, 5.198467090),
        (1.276839, 7860.419392439, 5.988822341),
        (1.193379, 5223.693919802, 3.649823730),
        (1.115322, 3930.209696220, 1.422745069),
        (0.794185, 11506.769769794, 2.322313077),
        (0.600309, 1577.343542448, 2.678271909),
        (0.496817, 6208.294251424, 5.696701824),
        (0.486306, 5884.926846583, 0.520007179),
        (0.468597, 6244.942814354, 5.866398759),
        (0.447061, 26.298319800, 3.615796498),
        (0.435206, -398.149003408, 4.349338347),
        (0.432392, 74.781598567, 2.435898309),
        (0.375510, 5507.553238667, 4.103476804),
        // t2 (20 terms)
        (0.243085, -775.522611324, 3.651837925),
        (0.230685, 5856.477659115, 4.773852582),
        (0.203747, 12036.460734888, 4.333987818),
        (0.173435, 18849.227549974, 6.153743485),
        (0.159080, 10977.078804699, 1.890075226),
        (0.143935, -796.298006816, 5.957517795),
        (0.137927, 11790.629088659, 1.135934669),
        (0.119979, 38.133035638, 4.551585768),
        (0.118971, 5486.777843175, 1.914547226),
        (0.116120, 1059.381930189, 0.873504123),
        (0.101868, -5573.142801634, 5.984503847),
        (0.098358, 2544.314419883, 0.092793886),
        (0.080164, 206.185548437, 2.095377709),
        (0.079645, 4694.002954708, 2.949233637),
        (0.075019, 2942.463423292, 4.980931759),
        (0.064397, 5746.271337896, 1.280308748),
        (0.063814, 5760.498431898, 4.167901731),
        (0.062617, 20.775395492, 2.654394814),
        (0.058844, 426.598190876, 4.839650148),
        (0.054139, 17260.154654690, 3.411091093),
        // t3 (20 terms)
        (0.048373, 155.420399434, 2.251573730),
        (0.048042, 2146.165416475, 1.495846011),
        (0.046551, -0.980321068, 0.921573539),
        (0.042732, 632.783739313, 5.720622217),
        (0.042560, 161000.685737473, 1.270837679),
        (0.042411, 6275.962302991, 2.869567043),
        (0.040759, 12352.852604545, 3.981496998),
        (0.040480, 15720.838784878, 2.546610123),
        (0.040184, -7.113547001, 3.565975565),
        (0.036955, 3154.687084896, 5.071801441),
        (0.036564, 5088.628839767, 3.324679049),
        (0.036507, 801.820931124, 6.248866009),
        (0.034867, 522.577418094, 5.210064075),
        (0.033529, 9437.762934887, 2.404714239),
        (0.033477, 6062.663207553, 4.144987272),
        (0.032438, 6076.890301554, 0.749317412),
        (0.032423, 8827.390269875, 5.541473556),
        (0.030215, 7084.896781115, 3.389610345),
        (0.029862, 12139.553509107, 1.770181024),
        (0.029247, -71430.695617928, 4.183178762),
        // t4 (20 terms)
        (0.028244, -6286.598968340, 5.069663519),
        (0.027567, 6279.552731642, 5.040846034),
        (0.025196, 1748.016413067, 2.901883301),
        (0.024816, -1194.447010225, 1.087136918),
        (0.022567, 6133.512652857, 3.307984806),
        (0.022509, 10447.387839604, 1.460726241),
        (0.021691, 14143.495242431, 5.952658009),
        (0.020937, 8429.241266467, 0.652303414),
        (0.020322, 419.484643875, 3.735430632),
        (0.017673, 6812.766815086, 3.186129845),
        (0.017806, 73.297125859, 3.475975097),
        (0.016155, 10213.285546211, 1.331103168),
        (0.015974, -2352.866153772, 6.145309371),
        (0.015949, -220.412642439, 4.005298270),
        (0.015078, 19651.048481098, 3.969480770),
        (0.014751, 1349.867409659, 4.308933301),
        (0.014318, 16730.463689596, 3.016058075),
        (0.014223, 17789.845619785, 2.104551349),
        (0.013671, -536.804512095, 5.971672571),
        (0.012462, 103.092774219, 1.737438797),
        // t5 (13 terms)
        (0.012420, 4690.479836359, 4.734090399),
        (0.011942, 8031.092263058, 2.053414715),
        (0.011847, 5643.178563677, 5.489005403),
        (0.011707, -4705.732307544, 2.654125618),
        (0.011622, 5120.601145584, 4.863931876),
        (0.010962, 3.590428652, 2.196567739),
        (0.010825, 553.569402842, 0.842715011),
        (0.010396, 951.718406251, 5.717799605),
        (0.010453, 5863.591206116, 1.913704550),
        (0.010099, 283.859318865, 1.942176992),
        (0.009858, 6309.374169791, 1.061816410),
        (0.009963, 149.563197135, 4.870690598),
        (0.009370, 149854.400135205, 0.673880395),
    ];

    const TERMS_T1: &[(f64, f64, f64)] = &[
        // t24 (6 terms)
        (102.156724, 6283.075849991, 4.249032005),
        (1.706807, 12566.151699983, 4.205904248),
        (0.269668, 213.299095438, 3.400290479),
        (0.265919, 529.690965095, 5.836047367),
        (0.210568, -3.523118349, 6.262738348),
        (0.077996, 5223.693919802, 4.670344204),
        // t25 (21 terms)
        (0.059146, 26.298319800, 1.083044735),
        (0.054764, 1577.343542448, 4.534800170),
        (0.034420, -398.149003408, 5.980077351),
        (0.033595, 5507.553238667, 5.980162321),
        (0.032088, 18849.227549974, 4.162913471),
        (0.029198, 5856.477659115, 0.623811863),
        (0.027764, 155.420399434, 3.745318113),
        (0.025190, 5746.271337896, 2.980330535),
        (0.024976, 5760.498431898, 2.467913690),
        (0.022997, -796.298006816, 1.174411803),
        (0.021774, 206.185548437, 3.854787540),
        (0.017925, -775.522611324, 1.092065955),
        (0.013794, 426.598190876, 2.699831988),
        (0.013276, 6062.663207553, 5.845801920),
        (0.012869, 6076.890301554, 5.333425680),
        (0.012152, 1059.381930189, 6.222874454),
        (0.011774, 12036.460734888, 2.292832062),
        (0.011081, -7.113547001, 5.154724984),
        (0.010143, 4694.002954708, 4.044013795),
        (0.010084, 522.577418094, 0.749320262),
        (0.009357, 5486.777843175, 3.416081409),
    ];

    const TERMS_T2: &[(f64, f64, f64)] = &[
        // t29 (6 terms)
        (0.370115, 0.0, 4.712388980),
        (4.322990, 6283.075849991, 2.642893748),
        (0.122605, 12566.151699983, 2.438140634),
        (0.019476, 213.299095438, 1.642186981),
        (0.016916, 529.690965095, 4.510959344),
        (0.013374, -3.523118349, 1.502210314),
    ];

    let mut s0 = 0.0f64;
    for &(amp, freq, phase) in TERMS_T0 {
        s0 += amp * (freq * t + phase).sin();
    }
    let mut s1 = 0.0f64;
    for &(amp, freq, phase) in TERMS_T1 {
        s1 += amp * (freq * t + phase).sin();
    }
    let mut s2 = 0.0f64;
    for &(amp, freq, phase) in TERMS_T2 {
        s2 += amp * (freq * t + phase).sin();
    }
    let s3 = 0.143388 * (6283.075849991 * t + 1.131453581).sin();

    (s0 + t * s1 + tt * s2 + t * tt * s3) * 1.0e-6
}

/// Convert AstroSat UTC seconds to `MJD (TDB)` (geocentric TDB, before Roemer/Shapiro).
pub fn astrosat_utc_to_mjd_tdb(t_utc_sec: f64) -> f64 {
    let (jd_hi, jd_lo) = astrosat_utc_to_jd_tt(t_utc_sec);
    let einstein = tdb_minus_tt(jd_hi, jd_lo);
    (jd_hi - JD_MINUS_MJD) + (jd_lo + einstein / SECS_PER_DAY)
}

/// Format an AstroSat UTC timestamp (seconds since 2010-01-01T00:00:00 UTC) as ISO-8601
/// (`YYYY-MM-DDTHH:MM:SS.fff`) suitable for FITS `DATE-OBS` / `DATE-END`.
///
/// Never panics on non-finite or extreme timestamps (returns `"0000-00-00T00:00:00.000"`).
pub fn astrosat_utc_to_iso8601(t_utc_sec: f64) -> String {
    if !t_utc_sec.is_finite() {
        return "0000-00-00T00:00:00.000".to_string();
    }
    let floored = t_utc_sec.floor();
    // Guard against f64 -> i64 saturation and chrono::NaiveDateTime year bounds.
    if !(-6.0e12..=6.0e12).contains(&floored) {
        return "0000-00-00T00:00:00.000".to_string();
    }
    let epoch = NaiveDate::from_ymd_opt(2010, 1, 1)
        .unwrap()
        .and_hms_opt(0, 0, 0)
        .unwrap();
    let whole_secs = floored as i64;
    let nanos = ((t_utc_sec - floored) * 1e9)
        .round()
        .clamp(0.0, 999_999_999.0) as i64;
    let Some(dur_s) = Duration::try_seconds(whole_secs) else {
        return "0000-00-00T00:00:00.000".to_string();
    };
    let Some(dt1) = epoch.checked_add_signed(dur_s) else {
        return "0000-00-00T00:00:00.000".to_string();
    };
    let Some(dt) = dt1.checked_add_signed(Duration::nanoseconds(nanos)) else {
        return "0000-00-00T00:00:00.000".to_string();
    };
    dt.format("%Y-%m-%dT%H:%M:%S%.3f").to_string()
}

/// Parse an ISO-8601 date (`YYYY-MM-DD` or `YYYY-MM-DDTHH:MM:SS`) into MJD.
pub fn parse_date_to_mjd(date_str: &str) -> Option<f64> {
    let s = date_str.trim();
    if s.is_empty() {
        return None;
    }
    let epoch = NaiveDate::from_ymd_opt(1858, 11, 17)?.and_hms_opt(0, 0, 0)?;
    if let Ok(dt) = NaiveDateTime::parse_from_str(s, "%Y-%m-%dT%H:%M:%S%.f") {
        let diff = dt.signed_duration_since(epoch);
        return Some(diff.num_milliseconds() as f64 / 86_400_000.0);
    }
    if let Ok(dt) = NaiveDateTime::parse_from_str(s, "%Y-%m-%d %H:%M:%S%.f") {
        let diff = dt.signed_duration_since(epoch);
        return Some(diff.num_milliseconds() as f64 / 86_400_000.0);
    }
    if let Ok(d) = NaiveDate::parse_from_str(s, "%Y-%m-%d") {
        let dt = d.and_hms_opt(0, 0, 0)?;
        let diff = dt.signed_duration_since(epoch);
        return Some(diff.num_seconds() as f64 / 86_400.0);
    }
    None
}
