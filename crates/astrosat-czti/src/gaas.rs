//! Stage 2b: `cztgaas` — Generate AstroSat CZTI Average Aspect (`RA`, `DEC`, `TWIST`) per quadrant.
//!
//! Transforms CZTI detector body axes `(0, 0, 1)` (pointing) and `(0, 1, 0)` (twist reference)
//! through the CALDB `TELDEF` alignment matrix (`ALIGNM11..33`) and MKF Roll/Pitch/Yaw inertial
//! direction cosines into the Inertial Coordinate System (ICS).

use crate::types::QuadAspect;

/// 3x3 alignment matrix from CALDB `TELDEF` (`AS1czt{0,1,2,3}*.teldef`), mapping CZTI detector
/// frame `(detX, detY, detZ)` -> AstroSat spacecraft frame `(satX, satY, satZ)`.
#[derive(Debug, Clone, Copy)]
pub struct TeldefMatrix {
    pub m: [[f64; 3]; 3],
}

impl Default for TeldefMatrix {
    fn default() -> Self {
        Self {
            m: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        }
    }
}

/// Attitude sample from MKF containing Roll, Pitch, and Yaw axis `(RA, DEC)` in radians.
#[derive(Debug, Clone, Copy, Default)]
pub struct RpySampleRad {
    pub time: f64,
    pub roll_ra: f64,
    pub roll_dec: f64,
    pub pitch_ra: f64,
    pub pitch_dec: f64,
    pub yaw_ra: f64,
    pub yaw_dec: f64,
}

/// Transform a unit vector `(det_x, det_y, det_z)` in CZTI detector frame into normalized ICS
/// `(nX, nY, nZ)` using `teldef` and MKF Roll/Pitch/Yaw (`coordinateTrasformation.cpp::get_inertial_vector`).
#[inline]
pub fn det_to_inertial(
    rpy: &RpySampleRad,
    teldef: &TeldefMatrix,
    det: [f64; 3],
) -> [f64; 3] {
    // Detector frame -> satellite reference frame
    let sat_x = teldef.m[0][0] * det[0] + teldef.m[0][1] * det[1] + teldef.m[0][2] * det[2];
    let sat_y = teldef.m[1][0] * det[0] + teldef.m[1][1] * det[1] + teldef.m[1][2] * det[2];
    let sat_z = teldef.m[2][0] * det[0] + teldef.m[2][1] * det[1] + teldef.m[2][2] * det[2];

    // Direction cosines of Yaw, Roll, Pitch axes in ICS
    let yaw_l = rpy.yaw_ra.cos() * rpy.yaw_dec.cos();
    let yaw_m = rpy.yaw_ra.sin() * rpy.yaw_dec.cos();
    let yaw_n = rpy.yaw_dec.sin();

    let roll_l = rpy.roll_ra.cos() * rpy.roll_dec.cos();
    let roll_m = rpy.roll_ra.sin() * rpy.roll_dec.cos();
    let roll_n = rpy.roll_dec.sin();

    let pitch_l = rpy.pitch_ra.cos() * rpy.pitch_dec.cos();
    let pitch_m = rpy.pitch_ra.sin() * rpy.pitch_dec.cos();
    let pitch_n = rpy.pitch_dec.sin();

    // Satellite frame -> ICS
    let ics_x = yaw_l * sat_x + roll_l * sat_y + pitch_l * sat_z;
    let ics_y = yaw_m * sat_x + roll_m * sat_y + pitch_m * sat_z;
    let ics_z = yaw_n * sat_x + roll_n * sat_y + pitch_n * sat_z;

    let norm = (ics_x * ics_x + ics_y * ics_y + ics_z * ics_z).sqrt();
    if norm > 0.0 {
        [ics_x / norm, ics_y / norm, ics_z / norm]
    } else {
        [0.0, 0.0, 1.0]
    }
}

/// Compute average `(RA_deg, DEC_deg, TWIST_deg)` from ICS pointing `(nX, nY, nZ)` and twist
/// reference `(nXt, nYt, nZt)` series (`coordinateTrasformation.cpp::calculate_RA_DEC_Twist`).
pub fn calculate_ra_dec_twist_deg(
    nx: &[f32],
    ny: &[f32],
    nz: &[f32],
    nxt: &[f32],
    nyt: &[f32],
    nzt: &[f32],
) -> (f64, f64, f64) {
    if nx.is_empty() {
        return (0.0, 0.0, 0.0);
    }
    let inv_n = 1.0 / (nx.len() as f64);
    let ax: f64 = nx.iter().map(|&v| v as f64).sum::<f64>() * inv_n;
    let ay: f64 = ny.iter().map(|&v| v as f64).sum::<f64>() * inv_n;
    let az: f64 = nz.iter().map(|&v| v as f64).sum::<f64>() * inv_n;
    let tx: f64 = nxt.iter().map(|&v| v as f64).sum::<f64>() * inv_n;
    let ty: f64 = nyt.iter().map(|&v| v as f64).sum::<f64>() * inv_n;
    let tz: f64 = nzt.iter().map(|&v| v as f64).sum::<f64>() * inv_n;

    let dec = az.clamp(-1.0, 1.0).asin();
    let mut ra = ay.atan2(ax);
    if ra < 0.0 {
        ra += std::f64::consts::TAU;
    }

    let temp_ct = tz * dec.cos() - tx * dec.sin() * ra.cos() - ty * dec.sin() * ra.sin();
    let temp_st = ty * ra.cos() - tx * ra.sin();
    let twist = temp_st.atan2(temp_ct);

    (ra.to_degrees(), dec.to_degrees(), twist.to_degrees())
}

/// Interpolate MKF `RpySampleRad` at time `t` via binary search + linear interpolation.
pub fn interpolate_rpy(samples: &[RpySampleRad], t: f64) -> RpySampleRad {
    let n = samples.len();
    if n == 0 {
        return RpySampleRad::default();
    }
    if n == 1 || t <= samples[0].time {
        return samples[0];
    }
    if t >= samples[n - 1].time {
        return samples[n - 1];
    }
    let hi = samples.partition_point(|s| s.time < t);
    if hi == 0 {
        return samples[0];
    }
    let lo = hi - 1;
    let a = &samples[lo];
    let b = &samples[hi];
    let dt = b.time - a.time;
    if dt <= 0.0 {
        return *a;
    }
    let f = (t - a.time) / dt;
    let lerp_angle = |u: f64, v: f64| -> f64 {
        let mut d = v - u;
        if d > std::f64::consts::PI {
            d -= std::f64::consts::TAU;
        } else if d < -std::f64::consts::PI {
            d += std::f64::consts::TAU;
        }
        u + f * d
    };
    RpySampleRad {
        time: t,
        roll_ra: lerp_angle(a.roll_ra, b.roll_ra),
        roll_dec: a.roll_dec + f * (b.roll_dec - a.roll_dec),
        pitch_ra: lerp_angle(a.pitch_ra, b.pitch_ra),
        pitch_dec: a.pitch_dec + f * (b.pitch_dec - a.pitch_dec),
        yaw_ra: lerp_angle(a.yaw_ra, b.yaw_ra),
        yaw_dec: a.yaw_dec + f * (b.yaw_dec - a.yaw_dec),
    }
}

/// Compute a quadrant's aspect time series at 1 s cadence over `[tstart, tstop]` (or restricted
/// to `gti_intervals` if non-empty), matching `cztgaas_v2.cpp`.
pub fn compute_quad_aspect(
    quad_id: u8,
    teldef_name: String,
    teldef: &TeldefMatrix,
    mkf_rpy: &[RpySampleRad],
    tstart: f64,
    tstop: f64,
    gti_intervals: Option<&[(f64, f64)]>,
) -> QuadAspect {
    let num = ((tstop - tstart).max(0.0)) as usize;
    let mut time = Vec::with_capacity(num.max(1));
    let mut nx = Vec::with_capacity(num.max(1));
    let mut ny = Vec::with_capacity(num.max(1));
    let mut nz = Vec::with_capacity(num.max(1));
    let mut nxt = Vec::with_capacity(num.max(1));
    let mut nyt = Vec::with_capacity(num.max(1));
    let mut nzt = Vec::with_capacity(num.max(1));

    let in_gti = |t: f64| -> bool {
        match gti_intervals {
            None | Some([]) => true,
            Some(g) => g.iter().any(|&(a, b)| t >= a && t <= b),
        }
    };

    let steps = num.max(1);
    for j in 0..steps {
        let t = tstart + (j as f64);
        if !in_gti(t) {
            continue;
        }
        let rpy = interpolate_rpy(mkf_rpy, t);
        let p = det_to_inertial(&rpy, teldef, [0.0, 0.0, 1.0]);
        let tw = det_to_inertial(&rpy, teldef, [0.0, 1.0, 0.0]);
        time.push(t);
        nx.push(p[0] as f32);
        ny.push(p[1] as f32);
        nz.push(p[2] as f32);
        nxt.push(tw[0] as f32);
        nyt.push(tw[1] as f32);
        nzt.push(tw[2] as f32);
    }

    // If GTI filtering excluded every integer step, fall back to midpoint sample
    if time.is_empty() && !mkf_rpy.is_empty() {
        let t = 0.5 * (tstart + tstop);
        let rpy = interpolate_rpy(mkf_rpy, t);
        let p = det_to_inertial(&rpy, teldef, [0.0, 0.0, 1.0]);
        let tw = det_to_inertial(&rpy, teldef, [0.0, 1.0, 0.0]);
        time.push(t);
        nx.push(p[0] as f32);
        ny.push(p[1] as f32);
        nz.push(p[2] as f32);
        nxt.push(tw[0] as f32);
        nyt.push(tw[1] as f32);
        nzt.push(tw[2] as f32);
    }

    let (ra_deg, dec_deg, twist_deg) =
        calculate_ra_dec_twist_deg(&nx, &ny, &nz, &nxt, &nyt, &nzt);

    QuadAspect {
        quad_id,
        teldef_name,
        time,
        nx,
        ny,
        nz,
        nxt,
        nyt,
        nzt,
        ra_deg,
        dec_deg,
        twist_deg,
    }
}
