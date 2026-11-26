//! Core domain types and configuration structures for AstroSat LAXPC.

use crate::error::{LaxpcError, Result};
use serde::{Deserialize, Serialize};
use std::fmt;
use std::str::FromStr;

/// AstroSat LAXPC detector unit (LAXPC10, LAXPC20, LAXPC30).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize, Deserialize)]
pub enum LaxpcUnit {
    /// LAXPC10 (`lxp1`, detector ID 1).
    Lxp10 = 1,
    /// LAXPC20 (`lxp2`, detector ID 2).
    Lxp20 = 2,
    /// LAXPC30 (`lxp3`, detector ID 3).
    Lxp30 = 3,
}

impl LaxpcUnit {
    /// All three LAXPC units in canonical order.
    pub const ALL: [LaxpcUnit; 3] = [LaxpcUnit::Lxp10, LaxpcUnit::Lxp20, LaxpcUnit::Lxp30];

    /// Numeric detector ID (1, 2, or 3) as encoded in telemetry frames and FITS `DETECTOR_ID`.
    #[inline]
    pub const fn id(self) -> u8 {
        self as u8
    }

    /// Zero-based index (0, 1, or 2) into per-unit arrays.
    #[inline]
    pub const fn index(self) -> usize {
        (self as usize) - 1
    }

    /// Construct from numeric detector ID (1, 2, or 3).
    pub fn from_id(id: u8) -> Result<Self> {
        match id {
            1 => Ok(Self::Lxp10),
            2 => Ok(Self::Lxp20),
            3 => Ok(Self::Lxp30),
            other => Err(LaxpcError::InvalidUnit(other.to_string())),
        }
    }

    /// File prefix used in Level-1 and Format B Level-2 filenames (`lxp1`, `lxp2`, `lxp3`).
    #[inline]
    pub const fn prefix(self) -> &'static str {
        match self {
            Self::Lxp10 => "lxp1",
            Self::Lxp20 => "lxp2",
            Self::Lxp30 => "lxp3",
        }
    }

    /// Response/channel prefix (`lx10`, `lx20`, `lx30`).
    #[inline]
    pub const fn short_name(self) -> &'static str {
        match self {
            Self::Lxp10 => "lx10",
            Self::Lxp20 => "lx20",
            Self::Lxp30 => "lx30",
        }
    }

    /// FITS `DETNAM` keyword (`LX10`, `LX20`, `LX30`).
    #[inline]
    pub const fn detnam(self) -> &'static str {
        match self {
            Self::Lxp10 => "LX10",
            Self::Lxp20 => "LX20",
            Self::Lxp30 => "LX30",
        }
    }

    /// Default number of spectrum channels (`512` for LAXPC10/30, `256` for LAXPC20)
    /// to avoid even-odd ADC binning fluctuations, matching `laxpcl1.f`.
    #[inline]
    pub const fn default_spectrum_channels(self) -> usize {
        match self {
            Self::Lxp10 | Self::Lxp30 => 512,
            Self::Lxp20 => 256,
        }
    }

    /// Lower and upper ULD rate bounds (counts/s) used for 4-quartile ULD binning (`laxpcl1.f:98`).
    #[inline]
    pub const fn uld_quartile_bounds(self) -> (f64, f64) {
        match self {
            Self::Lxp10 => (750.0, 1050.0),
            Self::Lxp20 => (780.0, 1100.0),
            Self::Lxp30 => (740.0, 1090.0),
        }
    }

    /// Default minimum all-anode count rate threshold (cps) for active-detector GTI detection
    /// (`laxpcl1.f:1210,1224-1225`).
    #[inline]
    pub fn count_rate_threshold(self, t_utc_start: f64) -> f64 {
        match self {
            Self::Lxp10 | Self::Lxp20 => 80.0,
            Self::Lxp30 => {
                if t_utc_start > 2.37e8 {
                    20.0
                } else {
                    50.0
                }
            }
        }
    }

    /// Reference 30 keV calibration peak channel in anode A8 for nominal RMFs (`backshiftv3.f:44,417-421`).
    #[inline]
    pub fn reference_30kev_channel(self, ndays: f64) -> f64 {
        match self {
            Self::Lxp10 => {
                if ndays > 2944.0 {
                    35.0
                } else if ndays > 1112.0 {
                    27.0
                } else if ndays > 1036.0 {
                    42.0
                } else if ndays > 991.0 {
                    63.6
                } else if ndays > 917.0 {
                    99.4
                } else {
                    239.1
                }
            }
            Self::Lxp20 => 224.1,
            Self::Lxp30 => 248.0,
        }
    }

    /// Quadratic energy non-linearity parameter `e2i` used in `backshiftv3.f:gainsh`.
    #[inline]
    pub const fn gain_shift_e2(self) -> f64 {
        match self {
            Self::Lxp10 => -1.98e-3,
            Self::Lxp20 => -2.67e-3,
            Self::Lxp30 => -1.97e-3,
        }
    }

    /// Nominal ULD cutoff channel (1024-channel scale) used in `backshiftv3.f:gainsh`.
    #[inline]
    pub const fn uld_cutoff_channel_1024(self) -> usize {
        match self {
            Self::Lxp10 => 631,
            Self::Lxp20 => 556,
            Self::Lxp30 => 629,
        }
    }
}

impl fmt::Display for LaxpcUnit {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.detnam())
    }
}

impl FromStr for LaxpcUnit {
    type Err = LaxpcError;

    fn from_str(s: &str) -> Result<Self> {
        match s.trim().to_ascii_lowercase().as_str() {
            "1" | "10" | "lxp1" | "lxp10" | "lx10" | "laxpc10" => Ok(Self::Lxp10),
            "2" | "20" | "lxp2" | "lxp20" | "lx20" | "laxpc20" => Ok(Self::Lxp20),
            "3" | "30" | "lxp3" | "lxp30" | "lx30" | "laxpc30" => Ok(Self::Lxp30),
            _ => Err(LaxpcError::InvalidUnit(s.to_string())),
        }
    }
}

/// Anode / layer selection matching `laxpcl1.f` parameter `ian` plus custom sets.
///
/// LAXPC has 5 anode layers (7 main anodes):
/// - Layer 1 (top): A1, A2
/// - Layer 2: A3, A4
/// - Layer 3: A5
/// - Layer 4: A6
/// - Layer 5: A7
/// - Veto / Calibration: A8 (Am-241 calibration source), A9, A10 (side/bottom veto)
#[derive(Debug, Clone, Default, PartialEq, Eq, Serialize, Deserialize)]
pub enum AnodeSelection {
    /// All main anodes A1–A7 (`ian = 0`).
    #[default]
    AllMain,
    /// Top layer A1 + A2 (`ian = 1`).
    Layer1,
    /// Second layer A3 + A4 (`ian = 3`).
    Layer2,
    /// Single specific anode 1..=10 (`ian = 2, 4..=10`, or explicit single anode).
    SingleAnode(u8),
    /// Custom subset of anode IDs (1..=10).
    Custom(Vec<u8>),
}

impl AnodeSelection {
    /// Construct from legacy `laxpcl1.f` `ian` integer (`0` = A1-A7, `1` = Layer 1 (A1+A2),
    /// `3` = Layer 2 (A3+A4), `2`/`4..=10` = individual anode).
    pub fn from_ian(ian: u8) -> Result<Self> {
        match ian {
            0 => Ok(Self::AllMain),
            1 => Ok(Self::Layer1),
            3 => Ok(Self::Layer2),
            2 | 4..=10 => Ok(Self::SingleAnode(ian)),
            _ => Err(LaxpcError::Config(format!(
                "invalid anode selector ian={ian} (expected 0..=10)"
            ))),
        }
    }

    /// Map an anode ID (1..=10) to physical LAXPC layer number (1..=5 for main anodes, 6..=8 for veto/cal).
    #[inline]
    pub const fn anode_to_layer(anode: u8) -> u8 {
        match anode {
            1 | 2 => 1,
            3 | 4 => 2,
            5 => 3,
            6 => 4,
            7 => 5,
            8 => 6,
            9 => 7,
            10 => 8,
            _ => 0,
        }
    }

    /// Check whether an event with `anode` (1..=10) passes this selection.
    #[inline]
    pub fn matches(&self, anode: u8) -> bool {
        match self {
            Self::AllMain => (1..=7).contains(&anode),
            Self::Layer1 => anode == 1 || anode == 2,
            Self::Layer2 => anode == 3 || anode == 4,
            Self::SingleAnode(a) => anode == *a,
            Self::Custom(list) => list.contains(&anode),
        }
    }

    /// Whether this selection represents top layer only (used for RMF selection).
    #[inline]
    pub fn is_layer1_only(&self) -> bool {
        matches!(self, Self::Layer1)
    }

    /// Return the FITS `FILTER` header string matching `laxpcl1.f`.
    pub fn filter_keyword(&self, single_only: bool) -> String {
        let ev = if single_only {
            "Single Events"
        } else {
            "All Events"
        };
        match self {
            Self::AllMain => format!("A1-A7 {ev}"),
            Self::Layer1 => format!("A1 {ev}"),
            Self::Layer2 => format!("A3 {ev}"),
            Self::SingleAnode(a) => format!("A{a} {ev}"),
            Self::Custom(list) => {
                let s: Vec<String> = list.iter().map(|a| format!("A{a}")).collect();
                format!("{} {ev}", s.join(","))
            }
        }
    }
}

/// Event filtering parameters for EA mode extraction.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EventFilter {
    /// If `true`, accept only single-anode events (`ic2 == 0`); if `false`, also include
    /// Xe K-escape double events (`laxpcl1.f` default `abs(iev) = 2`).
    pub single_only: bool,
    /// Anode / layer selection (default: `AllMain`, A1–A7).
    pub anodes: AnodeSelection,
    /// Inclusive lower 1024-scale channel bound for light curves and filtered events (`0..=1023`).
    pub pha_min: u16,
    /// Inclusive upper 1024-scale channel bound for light curves and filtered events (`0..=1023`).
    pub pha_max: u16,
    /// Optional energy range in keV `(e_min, e_max)`; if provided alongside a `ChannelMap`,
    /// overrides `(pha_min, pha_max)`.
    pub energy_range_kev: Option<(f64, f64)>,
}

impl Default for EventFilter {
    fn default() -> Self {
        Self {
            single_only: false,
            anodes: AnodeSelection::AllMain,
            pha_min: 0,
            pha_max: 1023,
            energy_range_kev: None,
        }
    }
}

/// Convert AstroSat UTC seconds (since MJD 55197 = 2010-01-01T00:00:00 UTC) to
/// days since AstroSat launch (2015-09-28 = MJD 57293, `t / 86400 - 2096.0`, `laxpcl1.f:869`).
#[inline]
pub fn utc_to_mission_day(t_utc: f64) -> f64 {
    t_utc / 86_400.0 - 2096.0
}

/// Default deadtime constant (43 microseconds, `laxpcl1.f:354`).
pub const DEFAULT_DEADTIME_SEC: f64 = 43.0e-6;

/// Singular rate guard for non-paralyzable deadtime (`23000 cps`, `laxpcl1.f:1303,1660`).
pub const DEADTIME_RATE_GUARD_CPS: f64 = 23_000.0;

/// Standard LAXPC Level-1 frame size in bytes (both EA and BB modes).
pub const L1_FRAME_SIZE: usize = 2242;

/// LAXPC 100 kHz STBG clock ticks per second.
pub const STBG_CLOCK_HZ: f64 = 100_000.0;

/// BBC coarse accumulation bin size (`0.128 * 256 = 32.768 s`, `laxpcl1.f:350`).
pub const BBC_BIN_SEC: f64 = 0.128 * 256.0;
