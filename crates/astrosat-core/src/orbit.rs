//! AstroSat spacecraft orbit (`.orb`) reader and 3-point centered-difference
//! 2nd-order Taylor interpolator (matching `as1bary/scorbit.c`).

use std::path::Path;

use astrosat_fits::{FitsFile, HduKind};

use crate::error::{CoreError, Result};

/// Geocentric inertial (J2000/ICRS) spacecraft position (m) and velocity (m/s).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct OrbitState {
    pub pos_m: [f64; 3],
    pub vel_ms: [f64; 3],
}

/// AstroSat spacecraft orbit ephemeris table (stored internally in SI meters and m/s).
#[derive(Debug, Clone)]
pub struct OrbitTable {
    time: Vec<f64>,
    pos_m: Vec<[f64; 3]>,
    vel_ms: Vec<[f64; 3]>,
}

impl OrbitTable {
    /// Construct from arrays in km and km/s (standard AstroSat `.orb` units).
    pub fn from_arrays(
        time: Vec<f64>,
        pos_km: Vec<[f64; 3]>,
        vel_kms: Vec<[f64; 3]>,
    ) -> Result<Self> {
        if time.len() < 2 || time.len() != pos_km.len() || time.len() != vel_kms.len() {
            return Err(CoreError::InvalidOrbit(format!(
                "orbit arrays must match and have >= 2 rows (t={}, pos={}, vel={})",
                time.len(),
                pos_km.len(),
                vel_kms.len()
            )));
        }
        let pos_m: Vec<[f64; 3]> = pos_km
            .into_iter()
            .map(|p| [p[0] * 1000.0, p[1] * 1000.0, p[2] * 1000.0])
            .collect();
        let vel_ms: Vec<[f64; 3]> = vel_kms
            .into_iter()
            .map(|v| [v[0] * 1000.0, v[1] * 1000.0, v[2] * 1000.0])
            .collect();
        Self::from_si(time, pos_m, vel_ms)
    }

    /// Construct directly from SI meters and m/s arrays.
    pub fn from_si(
        time: Vec<f64>,
        pos_m: Vec<[f64; 3]>,
        vel_ms: Vec<[f64; 3]>,
    ) -> Result<Self> {
        if time.len() < 2 || time.len() != pos_m.len() || time.len() != vel_ms.len() {
            return Err(CoreError::InvalidOrbit("invalid orbit array lengths".into()));
        }
        // Sort by time if needed and deduplicate identical timestamps.
        let mut rows: Vec<(f64, [f64; 3], [f64; 3])> = time
            .into_iter()
            .zip(pos_m)
            .zip(vel_ms)
            .filter(|((t, _, ), _)| t.is_finite())
            .map(|((t, p), v)| (t, p, v))
            .collect();
        rows.sort_by(|a, b| a.0.total_cmp(&b.0));
        rows.dedup_by(|a, b| (a.0 - b.0).abs() < 1e-9);
        if rows.len() < 2 {
            return Err(CoreError::InvalidOrbit(
                "fewer than 2 distinct orbit timestamps".into(),
            ));
        }

        let mut t_out = Vec::with_capacity(rows.len());
        let mut p_out = Vec::with_capacity(rows.len());
        let mut v_out = Vec::with_capacity(rows.len());
        for (t, p, v) in rows {
            t_out.push(t);
            p_out.push(p);
            v_out.push(v);
        }

        Ok(Self {
            time: t_out,
            pos_m: p_out,
            vel_ms: v_out,
        })
    }

    /// Load an AstroSat `.orb` FITS file (`ORBIT` HDU, converting km -> m automatically).
    pub fn from_fits(path: impl AsRef<Path>) -> Result<Self> {
        let fits = FitsFile::open(path)?;
        let hdu = fits
            .hdu_by_name("ORBIT")
            .or_else(|_| {
                fits.iter_hdus()
                    .find(|h| *h.kind() == HduKind::BinTable)
                    .ok_or_else(|| astrosat_fits::FitsError::HduNotFound("ORBIT".into()))
            })?;
        let bt = hdu.as_bintable()?;

        let time = bt.read_col_f64("TIME")?;

        // Support either scalar X,Y,Z / VX,VY,VZ or POS_X.. / VEL_X..
        let (x_name, y_name, z_name) = if bt.column_index("X").is_some() {
            ("X", "Y", "Z")
        } else if bt.column_index("POS_X").is_some() {
            ("POS_X", "POS_Y", "POS_Z")
        } else {
            return Err(CoreError::InvalidOrbit(
                "ORBIT HDU missing X/Y/Z or POS_X/POS_Y/POS_Z columns".into(),
            ));
        };
        let (vx_name, vy_name, vz_name) = if bt.column_index("VX").is_some() {
            ("VX", "VY", "VZ")
        } else if bt.column_index("VEL_X").is_some() {
            ("VEL_X", "VEL_Y", "VEL_Z")
        } else {
            return Err(CoreError::InvalidOrbit(
                "ORBIT HDU missing VX/VY/VZ or VEL_X/VEL_Y/VEL_Z columns".into(),
            ));
        };

        let xs = bt.read_col_f64(x_name)?;
        let ys = bt.read_col_f64(y_name)?;
        let zs = bt.read_col_f64(z_name)?;
        let vxs = bt.read_col_f64(vx_name)?;
        let vys = bt.read_col_f64(vy_name)?;
        let vzs = bt.read_col_f64(vz_name)?;

        if time.len() < 2
            || xs.len() != time.len()
            || ys.len() != time.len()
            || zs.len() != time.len()
            || vxs.len() != time.len()
            || vys.len() != time.len()
            || vzs.len() != time.len()
        {
            return Err(CoreError::InvalidOrbit(format!(
                "ORBIT table requires >= 2 rows (got {})",
                time.len()
            )));
        }

        // Detect unit from TUNIT or orbital radius magnitude (~7000 km vs ~7e6 m).
        let pos_unit_m = bt
            .column_spec(x_name)
            .and_then(|s| s.unit.as_deref())
            .map(|u| u.trim().eq_ignore_ascii_case("m"))
            .unwrap_or_else(|| {
                let r0 = (xs[0] * xs[0] + ys[0] * ys[0] + zs[0] * zs[0]).sqrt();
                r0 > 1.0e5
            });
        let vel_unit_ms = bt
            .column_spec(vx_name)
            .and_then(|s| s.unit.as_deref())
            .map(|u| u.trim().eq_ignore_ascii_case("m/s"))
            .unwrap_or(pos_unit_m);

        let p_scale = if pos_unit_m { 1.0 } else { 1000.0 };
        let v_scale = if vel_unit_ms { 1.0 } else { 1000.0 };

        let pos_m: Vec<[f64; 3]> = (0..time.len())
            .map(|i| [xs[i] * p_scale, ys[i] * p_scale, zs[i] * p_scale])
            .collect();
        let vel_ms: Vec<[f64; 3]> = (0..time.len())
            .map(|i| [vxs[i] * v_scale, vys[i] * v_scale, vzs[i] * v_scale])
            .collect();

        Self::from_si(time, pos_m, vel_ms)
    }

    pub fn len(&self) -> usize {
        self.time.len()
    }

    pub fn is_empty(&self) -> bool {
        self.time.is_empty()
    }

    pub fn tstart(&self) -> f64 {
        self.time[0]
    }

    pub fn tstop(&self) -> f64 {
        *self.time.last().unwrap()
    }

    pub fn times(&self) -> &[f64] {
        &self.time
    }

    /// Interpolate geocentric spacecraft position (m) and velocity (m/s) at time `t_sec`
    /// using 3-point centered-difference 2nd-order Taylor expansion (matching `as1bary/scorbit.c`).
    pub fn interpolate(&self, t_sec: f64) -> Result<OrbitState> {
        if !t_sec.is_finite() {
            return Err(CoreError::TimeOutOfRange(t_sec));
        }
        let n = self.time.len();
        // Allow small edge extrapolation (up to 120 s) cleanly, matching as1bary tolerance.
        let tol = 120.0;
        if t_sec < self.time[0] - tol || t_sec > self.time[n - 1] + tol {
            return Err(CoreError::TimeOutOfRange(t_sec));
        }

        // Find nearest sample index `i`.
        let pos = self.time.partition_point(|&t| t <= t_sec);
        let i = if pos == 0 {
            0
        } else if pos >= n {
            n - 1
        } else if (t_sec - self.time[pos - 1]).abs() <= (self.time[pos] - t_sec).abs() {
            pos - 1
        } else {
            pos
        };

        // Compute centered acceleration `a = dv/dt` at index `i`.
        let (i_lo, i_hi) = if n >= 3 && i > 0 && i < n - 1 {
            (i - 1, i + 1)
        } else if i == 0 {
            (0, 1)
        } else {
            (n - 2, n - 1)
        };
        let dt_span = self.time[i_hi] - self.time[i_lo];
        let mut acc = [0.0f64; 3];
        if dt_span > 0.0 {
            for k in 0..3 {
                acc[k] = (self.vel_ms[i_hi][k] - self.vel_ms[i_lo][k]) / dt_span;
            }
        }

        let dt = t_sec - self.time[i];
        let half_dt2 = 0.5 * dt * dt;
        let mut pos_m = [0.0f64; 3];
        let mut vel_ms = [0.0f64; 3];
        for k in 0..3 {
            pos_m[k] = self.pos_m[i][k] + self.vel_ms[i][k] * dt + acc[k] * half_dt2;
            vel_ms[k] = self.vel_ms[i][k] + acc[k] * dt;
        }

        Ok(OrbitState { pos_m, vel_ms })
    }

    /// Interpolate spacecraft position (km) and velocity (km/s).
    pub fn interpolate_km(&self, t_sec: f64) -> Result<([f64; 3], [f64; 3])> {
        let s = self.interpolate(t_sec)?;
        Ok((
            [s.pos_m[0] * 1e-3, s.pos_m[1] * 1e-3, s.pos_m[2] * 1e-3],
            [s.vel_ms[0] * 1e-3, s.vel_ms[1] * 1e-3, s.vel_ms[2] * 1e-3],
        ))
    }
}
