"""Real numerical head-to-head verification: astrosat bary vs official legacy NASA/IUCAA as1bary C oracle.

Runs across Scenarios B1 (real CZTI 9000000406), B2 (real LAXPC Crab 9000000252 Format B & Format A +
Crab pulsar timing), B3 (2016-12-31 -> 2017-01-01 leap-second crossing), and B4 (extreme geometry &
edge cases: solar conjunction Shapiro delay, FK5->ICRS c200to405, TIMEZERO != 0, ecliptic pole vs plane,
and exact jd == jd_end ephemeris boundary).

All outputs and temp artifacts are written strictly under D:\\ISRO\\.tmp\\bary_verify (WSL /mnt/d/ISRO/.tmp/bary_verify).
"""
from __future__ import annotations

import argparse
import json
import math
import os
import struct
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fits_reader import read_fits
from check_crab_timing import search_crab

ISRO_WIN = Path(r"d:\ISRO")
ISRO_WSL = "/mnt/d/ISRO"
TMP_WIN = ISRO_WIN / ".tmp" / "bary_verify"
TMP_WSL = f"{ISRO_WSL}/.tmp/bary_verify"

ASTROSAT_BIN_WSL = f"{ISRO_WSL}/target/release/astrosat"
AS1BARY_BIN_WSL = f"{ISRO_WSL}/.tmp/legacy_bins/as1bary/as1bary"
TIMING_DIR_WSL = f"{ISRO_WSL}/.tmp/legacy_bins/as1bary/timing_dir"

DE405_BIN_WSL = f"{ISRO_WSL}/data/ephem/lnxp1600p2200.405"
JPLEPH_405_WSL = f"{ISRO_WSL}/data/ephem/JPLEPH.405"
JPLEPH_200_WSL = f"{ISRO_WSL}/data/ephem/JPLEPH.200"

CZTI_EVT_WSL = f"{ISRO_WSL}/data/sample/czti/czti_sample_data/20160331_T01_112T01_9000000406_level2/modeM0/AS1T01_112T01_9000000406czt_level2_quad_clean.evt"
CZTI_ORB_AUX1_WSL = f"{ISRO_WSL}/data/sample/czti/czti_sample_data/20160331_T01_112T01_9000000406_level1/aux/aux1/AS1T01_112T01_9000000406czt_level1.orb"
CZTI_ORB_ARCH_WSL = f"{ISRO_WSL}/data/sample/czti_orb/AS1T01_112T01_9000000406czt_level1.orb"

LAXPC_FMTB_WSL = f"{ISRO_WSL}/.tmp/verify_run/laxpc/lxp1level2.evn"
LAXPC_FMTA_WSL = f"{ISRO_WSL}/.tmp/verify_run/laxpc/level2.event.fits"
LAXPC_ORB_MKF_WSL = f"{ISRO_WSL}/.tmp/verify_run/laxpc/AstroSat.orb"
LAXPC_ORB_L2_WSL = f"{ISRO_WSL}/data/sample/laxpc_l2/20160107_G02_010T01_9000000252_level2/laxpc/01509/aux/aux1/AS1G02_010T01_9000000252lxp_level2.orb"


def win_to_wsl(p: Path) -> str:
    s = str(p.resolve())
    if len(s) >= 2 and s[1] == ":":
        drive = s[0].lower()
        rest = s[2:].replace("\\", "/")
        return f"/mnt/{drive}{rest}"
    return s.replace("\\", "/")


def wsl_to_win(p: str) -> Path:
    if p.startswith("/mnt/") and len(p) > 6 and p[6] == "/":
        return Path(f"{p[5].upper()}:\\" + p[7:].replace("/", "\\"))
    return Path(p)


def run_wsl(cmd: str, check: bool = True) -> subprocess.CompletedProcess:
    full = f"source {ISRO_WSL}/scripts/env.sh && {cmd}"
    res = subprocess.run(
        ["wsl", "bash", "-lc", full],
        capture_output=True,
        text=True,
    )
    if check and res.returncode != 0:
        raise RuntimeError(
            f"WSL command failed (rc={res.returncode}):\nCMD: {cmd}\nSTDOUT:\n{res.stdout}\nSTDERR:\n{res.stderr}"
        )
    return res


def run_astrosat_bary(
    in_wsl: str,
    out_wsl: str,
    orb_wsl: str,
    ephem: str,
    ra: Optional[float] = None,
    dec: Optional[float] = None,
    frame: str = "icrs",
    check: bool = True,
) -> Tuple[subprocess.CompletedProcess, float]:
    args = [
        ASTROSAT_BIN_WSL,
        "bary",
        in_wsl,
        "--orbit",
        orb_wsl,
        "--ephem",
        ephem,
        "--frame",
        frame,
        "-o",
        out_wsl,
        "--json",
    ]
    if ra is not None:
        args += ["--ra", f"{ra:.10f}"]
    if dec is not None:
        args += ["--dec", f"{dec:.10f}"]
    t0 = time.perf_counter()
    res = run_wsl(" ".join(args), check=check)
    return res, time.perf_counter() - t0


def run_legacy_as1bary(
    in_wsl: str,
    out_wsl: str,
    orb_wsl: str,
    ra: Optional[float] = None,
    dec: Optional[float] = None,
    ref_frame: Optional[str] = "ICRS",
    check: bool = True,
) -> Tuple[subprocess.CompletedProcess, float]:
    args = [
        f"TIMING_DIR={TIMING_DIR_WSL}",
        AS1BARY_BIN_WSL,
        "-i",
        orb_wsl,
        "-f",
        in_wsl,
        "-o",
        out_wsl,
    ]
    if ra is not None and dec is not None:
        args += ["-ra", f"{ra:.10f}", "-dec", f"{dec:.10f}"]
    if ref_frame is not None:
        args += ["-ref", ref_frame]
    t0 = time.perf_counter()
    res = run_wsl(" ".join(args), check=check)
    return res, time.perf_counter() - t0


def residual_stats_us(t_a: np.ndarray, t_b: np.ndarray) -> Dict[str, Any]:
    diff_us = (np.asarray(t_a, dtype=np.float64) - np.asarray(t_b, dtype=np.float64)) * 1e6
    abs_us = np.abs(diff_us)
    return {
        "n": int(len(diff_us)),
        "max_abs_us": float(np.max(abs_us)) if len(diff_us) else 0.0,
        "rms_us": float(np.sqrt(np.mean(diff_us * diff_us))) if len(diff_us) else 0.0,
        "p99_abs_us": float(np.percentile(abs_us, 99)) if len(diff_us) else 0.0,
        "median_us": float(np.median(diff_us)) if len(diff_us) else 0.0,
        "min_us": float(np.min(diff_us)) if len(diff_us) else 0.0,
        "max_us": float(np.max(diff_us)) if len(diff_us) else 0.0,
    }


# ---------------------------------------------------------------------------
# Minimal Pure-Python FITS Writer for Synthetic B3 / B4 Test Fixtures
# ---------------------------------------------------------------------------
def _card_str(kw: str, val: str, comment: str = "") -> bytes:
    padded = val if len(val) >= 8 else val.ljust(8)
    s = f"{kw:<8}= '{padded}'"
    if comment:
        s = f"{s:<30} / {comment}"
    return s[:80].ljust(80).encode("ascii")


def _card_int(kw: str, val: int, comment: str = "") -> bytes:
    s = f"{kw:<8}= {val:>20d}"
    if comment:
        s = f"{s} / {comment}"
    return s[:80].ljust(80).encode("ascii")


def _card_float(kw: str, val: float, comment: str = "") -> bytes:
    s = f"{kw:<8}= {val:>20.15G}"
    if comment:
        s = f"{s} / {comment}"
    return s[:80].ljust(80).encode("ascii")


def _card_bool(kw: str, val: bool, comment: str = "") -> bytes:
    s = f"{kw:<8}= {'T' if val else 'F':>20}"
    if comment:
        s = f"{s} / {comment}"
    return s[:80].ljust(80).encode("ascii")


def _pad_block(data: bytes, fill: bytes = b" ") -> bytes:
    rem = len(data) % 2880
    return data if rem == 0 else data + fill * (2880 - rem)


def _build_header(cards: List[bytes]) -> bytes:
    raw = b"".join(cards) + b"END".ljust(80)
    return _pad_block(raw, b" ")


def write_synthetic_orbit_fits(
    out_path: Path,
    t_start: float,
    n_rows: int,
    deltat: float = 1.0,
    r_km: float = 7028.137,  # ~650 km altitude LEO
    inc_deg: float = 6.0,
) -> None:
    """Write a physically exact circular Keplerian AstroSat orbit FITS (km, km/s)."""
    mu_earth = 398600.4418  # km^3 / s^2
    omega = math.sqrt(mu_earth / (r_km ** 3))
    inc = math.radians(inc_deg)
    ci, si = math.cos(inc), math.sin(inc)

    ts = t_start + np.arange(n_rows, dtype=np.float64) * deltat
    theta = omega * (ts - t_start)
    ct, st = np.cos(theta), np.sin(theta)

    x = r_km * ct
    y = r_km * st * ci
    z = r_km * st * si
    vx = -r_km * omega * st
    vy = r_km * omega * ct * ci
    vz = r_km * omega * ct * si

    t_end = float(ts[-1])
    prim = _build_header([
        _card_bool("SIMPLE", True),
        _card_int("BITPIX", 8),
        _card_int("NAXIS", 0),
        _card_bool("EXTEND", True),
        _card_str("TELESCOP", "ASTROSAT"),
        _card_str("MISSION", "ASTROSAT"),
        _card_str("TIMESYS", "UTC"),
        _card_int("MJDREFI", 55197),
        _card_float("MJDREFF", 0.0),
        _card_float("TSTART", t_start),
        _card_float("TSTOP", t_end),
    ])

    cols = [("TIME", "s"), ("X", "km"), ("Y", "km"), ("Z", "km"), ("VX", "km/s"), ("VY", "km/s"), ("VZ", "km/s")]
    hcards = [
        _card_str("XTENSION", "BINTABLE"),
        _card_int("BITPIX", 8),
        _card_int("NAXIS", 2),
        _card_int("NAXIS1", 7 * 8),
        _card_int("NAXIS2", n_rows),
        _card_int("PCOUNT", 0),
        _card_int("GCOUNT", 1),
        _card_int("TFIELDS", 7),
        _card_str("EXTNAME", "ORBIT"),
        _card_str("TELESCOP", "ASTROSAT"),
        _card_str("MISSION", "ASTROSAT"),
        _card_str("TIMESYS", "UTC"),
        _card_str("TIMEUNIT", "s"),
        _card_int("MJDREFI", 55197),
        _card_float("MJDREFF", 0.0),
        _card_float("TSTART", t_start),
        _card_float("TSTOP", t_end),
        _card_float("DELTAT", deltat),
        _card_float("TIMEDEL", deltat),
    ]
    for idx, (cname, unit) in enumerate(cols, start=1):
        hcards.append(_card_str(f"TTYPE{idx}", cname))
        hcards.append(_card_str(f"TFORM{idx}", "1D"))
        hcards.append(_card_str(f"TUNIT{idx}", unit))

    tbl = np.column_stack([ts, x, y, z, vx, vy, vz]).astype(">f8").tobytes()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "wb") as fh:
        fh.write(prim)
        fh.write(_build_header(hcards))
        fh.write(_pad_block(tbl, b"\x00"))


def write_synthetic_event_fits(
    out_path: Path,
    times: np.ndarray,
    ra_deg: float,
    dec_deg: float,
    radecsys: str = "ICRS",
    timezero: float = 0.0,
    gti_intervals: Optional[List[Tuple[float, float]]] = None,
    tstart_hdr: Optional[float] = None,
    tstop_hdr: Optional[float] = None,
) -> None:
    """Write a minimal AstroSat event (+ optional GTI) FITS readable by both as1bary and astrosat bary."""
    times = np.asarray(times, dtype=np.float64)
    ts_val = float(times[0]) if tstart_hdr is None else tstart_hdr
    te_val = float(times[-1]) if tstop_hdr is None else tstop_hdr

    common_meta = [
        _card_str("TELESCOP", "ASTROSAT"),
        _card_str("MISSION", "ASTROSAT"),
        _card_str("INSTRUME", "LAXPC"),
        _card_str("TIMESYS", "UTC"),
        _card_str("TIMEUNIT", "s"),
        _card_str("TIMEREF", "LOCAL"),
        _card_int("MJDREFI", 55197),
        _card_float("MJDREFF", 0.0),
        _card_float("TIMEZERO", timezero),
        _card_float("TSTART", ts_val),
        _card_float("TSTOP", te_val),
        _card_float("RA_OBJ", ra_deg),
        _card_float("DEC_OBJ", dec_deg),
        _card_float("RA_NOM", ra_deg),
        _card_float("DEC_NOM", dec_deg),
        _card_float("RA_PNT", ra_deg),
        _card_float("DEC_PNT", dec_deg),
        _card_float("EQUINOX", 2000.0),
        _card_str("RADECSYS", radecsys),
    ]

    prim = _build_header([
        _card_bool("SIMPLE", True),
        _card_int("BITPIX", 8),
        _card_int("NAXIS", 0),
        _card_bool("EXTEND", True),
    ] + common_meta)

    ev_hdr = _build_header([
        _card_str("XTENSION", "BINTABLE"),
        _card_int("BITPIX", 8),
        _card_int("NAXIS", 2),
        _card_int("NAXIS1", 8),
        _card_int("NAXIS2", len(times)),
        _card_int("PCOUNT", 0),
        _card_int("GCOUNT", 1),
        _card_int("TFIELDS", 1),
        _card_str("EXTNAME", "EVENTS"),
        _card_str("TTYPE1", "TIME"),
        _card_str("TFORM1", "1D"),
        _card_str("TUNIT1", "s"),
    ] + common_meta)
    ev_data = _pad_block(times.astype(">f8").tobytes(), b"\x00")

    chunks = [prim, ev_hdr, ev_data]
    if gti_intervals:
        gti_arr = np.asarray(gti_intervals, dtype=">f8")
        gti_hdr = _build_header([
            _card_str("XTENSION", "BINTABLE"),
            _card_int("BITPIX", 8),
            _card_int("NAXIS", 2),
            _card_int("NAXIS1", 16),
            _card_int("NAXIS2", len(gti_intervals)),
            _card_int("PCOUNT", 0),
            _card_int("GCOUNT", 1),
            _card_int("TFIELDS", 2),
            _card_str("EXTNAME", "GTI"),
            _card_str("TTYPE1", "START"),
            _card_str("TFORM1", "1D"),
            _card_str("TUNIT1", "s"),
            _card_str("TTYPE2", "STOP"),
            _card_str("TFORM2", "1D"),
            _card_str("TUNIT2", "s"),
        ] + common_meta)
        chunks.extend([gti_hdr, _pad_block(gti_arr.tobytes(), b"\x00")])

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "wb") as fh:
        for ch in chunks:
            fh.write(ch)


# ---------------------------------------------------------------------------
# Scenario B1: Real CZTI 9000000406
# ---------------------------------------------------------------------------
def run_b1(test_jpleph_fits: bool = False) -> Dict[str, Any]:
    print("\n=== SCENARIO B1: Real CZTI 9000000406 ===")
    out: Dict[str, Any] = {}
    ra, dec = 83.63308, 22.01446

    # B1a: Live 0.128 s orbit (aux1)
    rust_de405_a = f"{TMP_WSL}/b1a_czti_rust_de405.evt"
    rust_anal_a = f"{TMP_WSL}/b1a_czti_rust_anal.evt"
    leg_a = f"{TMP_WSL}/b1a_czti_legacy.evt"
    _, t_r_de = run_astrosat_bary(CZTI_EVT_WSL, rust_de405_a, CZTI_ORB_AUX1_WSL, DE405_BIN_WSL, ra, dec)
    _, t_r_an = run_astrosat_bary(CZTI_EVT_WSL, rust_anal_a, CZTI_ORB_AUX1_WSL, "analytical", ra, dec)
    _, t_leg = run_legacy_as1bary(CZTI_EVT_WSL, leg_a, CZTI_ORB_AUX1_WSL, ra, dec, "ICRS")

    r_de = read_fits(str(wsl_to_win(rust_de405_a)))
    r_an = read_fits(str(wsl_to_win(rust_anal_a)))
    lg = read_fits(str(wsl_to_win(leg_a)))

    all_de, all_an, all_lg = [], [], []
    gti_de, gti_an, gti_lg = [], [], []
    per_hdu = {}
    for idx in [1, 2, 3, 4]:
        t1 = np.asarray(r_de[idx].read_column("TIME"), dtype=np.float64)
        t2 = np.asarray(r_an[idx].read_column("TIME"), dtype=np.float64)
        t3 = np.asarray(lg[idx].read_column("TIME"), dtype=np.float64)
        all_de.append(t1); all_an.append(t2); all_lg.append(t3)
        per_hdu[r_de[idx].extname] = {
            "de405_vs_legacy": residual_stats_us(t1, t3),
            "analytical_vs_legacy": residual_stats_us(t2, t3),
        }
    for idx in [8, 9, 10, 11]:
        for col in ("START", "STOP"):
            gti_de.append(np.asarray(r_de[idx].read_column(col), dtype=np.float64))
            gti_an.append(np.asarray(r_an[idx].read_column(col), dtype=np.float64))
            gti_lg.append(np.asarray(lg[idx].read_column(col), dtype=np.float64))

    all_de_v = np.concatenate(all_de); all_an_v = np.concatenate(all_an); all_lg_v = np.concatenate(all_lg)
    gti_de_v = np.concatenate(gti_de); gti_an_v = np.concatenate(gti_an); gti_lg_v = np.concatenate(gti_lg)

    b1a = {
        "orbit": "aux1 (live 0.128 s cadence)",
        "runtime_s": {"rust_de405": t_r_de, "rust_analytical": t_r_an, "legacy_as1bary": t_leg},
        "events_de405_vs_legacy": residual_stats_us(all_de_v, all_lg_v),
        "events_analytical_vs_legacy": residual_stats_us(all_an_v, all_lg_v),
        "gti_de405_vs_legacy": residual_stats_us(gti_de_v, gti_lg_v),
        "gti_analytical_vs_legacy": residual_stats_us(gti_an_v, gti_lg_v),
        "per_quadrant": per_hdu,
        "header_diff_q0": {
            "rust_tstart": r_de[1].header.get("TSTART"),
            "legacy_tstart": lg[1].header.get("TSTART"),
            "rust_q0gti_tstart": r_de[8].header.get("TSTART"),
            "legacy_q0gti_tstart": lg[8].header.get("TSTART"),
            "q0gti_tstart_diff_us": (float(r_de[8].header.get("TSTART", 0.0)) - float(lg[8].header.get("TSTART", 0.0))) * 1e6,
            "rust_primary_tstop": r_de[0].header.get("TSTOP"),
            "legacy_primary_tstop": lg[0].header.get("TSTOP"),
            "rust_date_obs": r_de[1].header.get("DATE-OBS"),
            "legacy_date_obs": lg[1].header.get("DATE-OBS"),
        },
    }
    print(f"  B1a (live aux1 .orb): events DE405 max|d|={b1a['events_de405_vs_legacy']['max_abs_us']:.6f} us, "
          f"RMS={b1a['events_de405_vs_legacy']['rms_us']:.6f} us | "
          f"analytical RMS={b1a['events_analytical_vs_legacy']['rms_us']:.2f} us")
    out["B1a_live_aux1_orb"] = b1a

    if test_jpleph_fits:
        rust_fits405 = f"{TMP_WSL}/b1a_czti_rust_jpleph405.evt"
        res_f, t_f = run_astrosat_bary(CZTI_EVT_WSL, rust_fits405, CZTI_ORB_AUX1_WSL, JPLEPH_405_WSL, ra, dec, check=False)
        if res_f.returncode == 0:
            rf = read_fits(str(wsl_to_win(rust_fits405)))
            all_rf = np.concatenate([np.asarray(rf[i].read_column("TIME"), dtype=np.float64) for i in [1, 2, 3, 4]])
            st_rf = residual_stats_us(all_rf, all_lg_v)
            st_bin = residual_stats_us(all_rf, all_de_v)
            out["B1a_fits_jpleph405"] = {
                "supported": True,
                "runtime_s": t_f,
                "vs_legacy": st_rf,
                "vs_binary_de405": st_bin,
            }
            print(f"  B1a (FITS JPLEPH.405): vs_legacy max|d|={st_rf['max_abs_us']:.6f} us, vs_bin_de405 max|d|={st_bin['max_abs_us']:.6f} us")
        else:
            out["B1a_fits_jpleph405"] = {"supported": False, "error": res_f.stderr.strip()}

    # B1b: Archive 1 s orbit (contains frozen GPS rows in obs window)
    rust_de405_b = f"{TMP_WSL}/b1b_czti_arch_rust_de405.evt"
    leg_b = f"{TMP_WSL}/b1b_czti_arch_legacy.evt"
    run_astrosat_bary(CZTI_EVT_WSL, rust_de405_b, CZTI_ORB_ARCH_WSL, DE405_BIN_WSL, ra, dec)
    run_legacy_as1bary(CZTI_EVT_WSL, leg_b, CZTI_ORB_ARCH_WSL, ra, dec, "ICRS")

    orig = read_fits(str(wsl_to_win(CZTI_EVT_WSL)))
    r_b = read_fits(str(wsl_to_win(rust_de405_b)))
    l_b = read_fits(str(wsl_to_win(leg_b)))
    ot_b = np.concatenate([np.asarray(orig[i].read_column("TIME"), dtype=np.float64) for i in [1, 2, 3, 4]])
    rt_b = np.concatenate([np.asarray(r_b[i].read_column("TIME"), dtype=np.float64) for i in [1, 2, 3, 4]])
    lt_b = np.concatenate([np.asarray(l_b[i].read_column("TIME"), dtype=np.float64) for i in [1, 2, 3, 4]])
    frac_b = ot_b % 1.0
    m_hi = frac_b >= 0.5
    m_lo = ~m_hi
    out["B1b_archive_frozen_orb"] = {
        "orbit": "czti_orb (1 s cadence; frozen/stale GPS X/VX during obs window)",
        "all_events": residual_stats_us(rt_b, lt_b),
        "frac_ge_0_5_same_row": residual_stats_us(rt_b[m_hi], lt_b[m_hi]),
        "frac_lt_0_5_scorbit_off_by_1": residual_stats_us(rt_b[m_lo], lt_b[m_lo]),
    }
    print(f"  B1b (archive frozen .orb): frac>=0.5 max|d|={out['B1b_archive_frozen_orb']['frac_ge_0_5_same_row']['max_abs_us']:.6f} us, "
          f"frac<0.5 max|d|={out['B1b_archive_frozen_orb']['frac_lt_0_5_scorbit_off_by_1']['max_abs_us']:.3f} us (scorbit.c:185 +0.5 offset on frozen GPS rows)")
    return out


# ---------------------------------------------------------------------------
# Scenario B2: Real LAXPC Crab 9000000252 (Format B + Format A + Crab timing)
# ---------------------------------------------------------------------------
def run_b2(reuse_existing: bool = True, test_jpleph_fits: bool = False) -> Dict[str, Any]:
    print("\n=== SCENARIO B2: Real LAXPC Crab 9000000252 ===")
    out: Dict[str, Any] = {}
    ra, dec = 83.633083, 22.0145

    # B2a: Format B (lxp1level2.evn, 8,939,817 events) with AstroSat.orb
    rust_de_b = f"{TMP_WSL}/laxpc_fmtB_rust_de405.evn"
    rust_an_b = f"{TMP_WSL}/laxpc_fmtB_rust_anal.evn"
    leg_b = f"{TMP_WSL}/laxpc_fmtB_legacy.evn"

    if not (reuse_existing and wsl_to_win(rust_de_b).exists()):
        _, t_r_de = run_astrosat_bary(LAXPC_FMTB_WSL, rust_de_b, LAXPC_ORB_MKF_WSL, DE405_BIN_WSL, ra, dec)
    else:
        t_r_de = -1.0
    _, t_r_an = run_astrosat_bary(LAXPC_FMTB_WSL, rust_an_b, LAXPC_ORB_MKF_WSL, "analytical", ra, dec)
    if not (reuse_existing and wsl_to_win(leg_b).exists()):
        _, t_leg = run_legacy_as1bary(LAXPC_FMTB_WSL, leg_b, LAXPC_ORB_MKF_WSL, ra, dec, "ICRS")
    else:
        t_leg = -1.0

    r_de = read_fits(str(wsl_to_win(rust_de_b)))
    r_an = read_fits(str(wsl_to_win(rust_an_b)))
    lg = read_fits(str(wsl_to_win(leg_b)))

    t_de = np.asarray(r_de[1].read_column("TIME"), dtype=np.float64)
    t_an = np.asarray(r_an[1].read_column("TIME"), dtype=np.float64)
    t_lg = np.asarray(lg[1].read_column("TIME"), dtype=np.float64)

    g_de = np.concatenate([r_de[2].read_column("START"), r_de[2].read_column("STOP")])
    g_an = np.concatenate([r_an[2].read_column("START"), r_an[2].read_column("STOP")])
    g_lg = np.concatenate([lg[2].read_column("START"), lg[2].read_column("STOP")])

    st_de = residual_stats_us(t_de, t_lg)
    st_an = residual_stats_us(t_an, t_lg)
    print(f"  B2a Format B + AstroSat.orb (N={len(t_de)}): DE405 vs legacy max|d|={st_de['max_abs_us']:.3f} us, "
          f"RMS={st_de['rms_us']:.3f} us, p99={st_de['p99_abs_us']:.3f} us | "
          f"Analytical RMS={st_an['rms_us']:.2f} us")

    # Also split Format B AstroSat.orb residuals by frac(t) >= 0.5 vs < 0.5 to isolate MKF dx/dt != vx + scorbit.c:185
    orig_b = read_fits(str(wsl_to_win(LAXPC_FMTB_WSL)))
    ot_b_all = np.asarray(orig_b[1].read_column("TIME"), dtype=np.float64)
    frac_b = ot_b_all % 1.0
    m_b_hi = frac_b >= 0.5
    m_b_lo = ~m_b_hi
    st_de_hi = residual_stats_us(t_de[m_b_hi], t_lg[m_b_hi])
    st_de_lo = residual_stats_us(t_de[m_b_lo], t_lg[m_b_lo])
    print(f"    AstroSat.orb split: frac>=0.5 max|d|={st_de_hi['max_abs_us']:.6f} us (RMS={st_de_hi['rms_us']:.6f} us); "
          f"frac<0.5 max|d|={st_de_lo['max_abs_us']:.3f} us (RMS={st_de_lo['rms_us']:.3f} us)")

    # Test Format B slice (99.3k strided events) with live 0.128 s LAXPC L2 .orb:
    # 1) Raw official L2 .orb (lacks TSTART/TSTOP/DELTAT -> legacy as1bary silently skips, Rust succeeds!)
    # 2) Patched-header copy of L2 .orb (adds TSTART/TSTOP/DELTAT=0.128 so legacy as1bary can run apples-to-apples!)
    sub_times = ot_b_all[::90]  # ~99,332 events spanning entire observation
    sub_in_win = TMP_WIN / "b2a_sub_liveorb_in.fits"
    write_synthetic_event_fits(sub_in_win, sub_times, ra, dec, "ICRS")
    sub_in_wsl = win_to_wsl(sub_in_win)
    sub_r_wsl = f"{TMP_WSL}/b2a_sub_liveorb_rust.fits"
    sub_l_raw_wsl = f"{TMP_WSL}/b2a_sub_liveorb_leg_raw.fits"
    run_astrosat_bary(sub_in_wsl, sub_r_wsl, LAXPC_ORB_L2_WSL, DE405_BIN_WSL, ra, dec)
    run_legacy_as1bary(sub_in_wsl, sub_l_raw_wsl, LAXPC_ORB_L2_WSL, ra, dec, "ICRS")
    sr_h = read_fits(str(wsl_to_win(sub_r_wsl)))
    sl_raw_h = read_fits(str(wsl_to_win(sub_l_raw_wsl)))
    sr = np.asarray(sr_h[1].read_column("TIME"), dtype=np.float64)
    legacy_raw_timesys = str(sl_raw_h[1].header.get("TIMESYS", "")).strip()

    # Build a clean header-patched copy of live 0.128s LAXPC L2 .orb (uniform 0.128s rows + TSTART/TSTOP/DELTAT)
    l2_orb_h = read_fits(str(wsl_to_win(LAXPC_ORB_L2_WSL)))[1]
    l2_ot = np.asarray(l2_orb_h.read_column("TIME"), dtype=np.float64)
    l2_cols = [np.asarray(l2_orb_h.read_column(c), dtype=np.float64) for c in ("X", "Y", "Z", "VX", "VY", "VZ")]
    patched_orb_win = TMP_WIN / "b2a_laxpc_l2_patched.orb"
    p_prim = _build_header([
        _card_bool("SIMPLE", True), _card_int("BITPIX", 8), _card_int("NAXIS", 0), _card_bool("EXTEND", True),
        _card_str("TELESCOP", "ASTROSAT"), _card_str("MISSION", "ASTROSAT"),
    ])
    p_hcards = [
        _card_str("XTENSION", "BINTABLE"), _card_int("BITPIX", 8), _card_int("NAXIS", 2),
        _card_int("NAXIS1", 56), _card_int("NAXIS2", len(l2_ot)), _card_int("PCOUNT", 0), _card_int("GCOUNT", 1),
        _card_int("TFIELDS", 7), _card_str("EXTNAME", "ORBIT"),
        _card_float("TSTART", float(l2_ot[0])), _card_float("TSTOP", float(l2_ot[-1])),
        _card_float("DELTAT", 0.128), _card_float("TIMEDEL", 0.128),
    ]
    for ci, (cn, cu) in enumerate([("TIME","s"),("X","km"),("Y","km"),("Z","km"),("VX","km/s"),("VY","km/s"),("VZ","km/s")], 1):
        p_hcards += [_card_str(f"TTYPE{ci}", cn), _card_str(f"TFORM{ci}", "1D"), _card_str(f"TUNIT{ci}", cu)]
    p_tbl = np.column_stack([l2_ot] + l2_cols).astype(">f8").tobytes()
    with open(patched_orb_win, "wb") as fh:
        fh.write(p_prim); fh.write(_build_header(p_hcards)); fh.write(_pad_block(p_tbl, b"\x00"))
    sub_l_patched_wsl = f"{TMP_WSL}/b2a_sub_liveorb_leg_patched.fits"
    run_legacy_as1bary(sub_in_wsl, sub_l_patched_wsl, win_to_wsl(patched_orb_win), ra, dec, "ICRS")
    sl = np.asarray(read_fits(str(wsl_to_win(sub_l_patched_wsl)))[1].read_column("TIME"), dtype=np.float64)
    st_live = residual_stats_us(sr, sl)
    print(f"  B2a-live (99.3k real Crab timestamps + live 0.128s L2 .orb): raw legacy TIMESYS={legacy_raw_timesys!r} (skipped); "
          f"header-patched max|d|={st_live['max_abs_us']:.6f} us, RMS={st_live['rms_us']:.6f} us")

    # B2c: Crab pulsar frequency & pulse-phase comparison (DE405 vs Analytical vs Legacy)
    print("  B2c: Running Crab pulsar timing (~29.66 Hz) on all 8.94M events (DE405 / Analytical / Legacy)...")
    crab_de = search_crab(np.sort(t_de))
    crab_an = search_crab(np.sort(t_an))
    crab_lg = search_crab(np.sort(t_lg))

    # Measure exact phase offset at shared reference epoch & frequency (crab_lg['best_freq_hz'])
    f_ref = crab_lg["best_freq_hz"]
    fdot = -3.693e-10
    t0_ref = float(np.min(t_lg))
    def harmonic1_phase(times: np.ndarray) -> float:
        dt = times - t0_ref
        ph = np.mod(f_ref * dt + 0.5 * fdot * dt * dt, 1.0)
        ang = 2.0 * math.pi * ph
        return math.atan2(float(np.sum(np.sin(ang))), float(np.sum(np.cos(ang)))) / (2.0 * math.pi)

    ph_lg = harmonic1_phase(t_lg)
    ph_de = harmonic1_phase(t_de)
    ph_an = harmonic1_phase(t_an)
    def wrap_phase_diff(a: float, b: float) -> float:
        d = (a - b + 0.5) % 1.0 - 0.5
        return d

    dphi_de_lg = wrap_phase_diff(ph_de, ph_lg)
    dphi_an_lg = wrap_phase_diff(ph_an, ph_lg)
    dt_phase_de_us = (dphi_de_lg / f_ref) * 1e6
    dt_phase_an_us = (dphi_an_lg / f_ref) * 1e6
    print(f"  B2c Crab: legacy f={crab_lg['best_freq_hz']:.7f} Hz (Z2^2={crab_lg['z2_2']:.1f}), "
          f"DE405 f={crab_de['best_freq_hz']:.7f} Hz (Z2^2={crab_de['z2_2']:.1f}, dphi={dphi_de_lg:+.6e} cyc = {dt_phase_de_us:+.3f} us), "
          f"Anal f={crab_an['best_freq_hz']:.7f} Hz (Z2^2={crab_an['z2_2']:.1f}, dphi={dphi_an_lg:+.5f} cyc = {dt_phase_an_us:+.1f} us)")

    # B2b: Format A (level2.event.fits, 25.8M events)
    # Run astrosat bary on full 25.8M Format A file, AND run legacy as1bary on a 50k Format A slice
    # (and test if legacy as1bary even accepts Format A schema directly!)
    rust_fmta_wsl = f"{TMP_WSL}/laxpc_fmtA_rust_de405.fits"
    _, t_fmta_rust = run_astrosat_bary(LAXPC_FMTA_WSL, rust_fmta_wsl, LAXPC_ORB_MKF_WSL, DE405_BIN_WSL, ra, dec)
    # Verify Format A slice against legacy as1bary using a 50k-row Format A FITS preserving exact Format A structure
    fa_orig = read_fits(str(wsl_to_win(LAXPC_FMTA_WSL)))
    fa_rust = read_fits(str(wsl_to_win(rust_fmta_wsl)))
    fa_ot = np.asarray(fa_orig[1].read_column("TIME"), dtype=np.float64)
    fa_rt = np.asarray(fa_rust[1].read_column("TIME"), dtype=np.float64)
    stride = max(1, len(fa_ot) // 50_000)
    fa_sub_idx = np.arange(0, len(fa_ot), stride)
    fa_sub_in_win = TMP_WIN / "b2b_fmtA_50k_in.fits"
    write_synthetic_event_fits(fa_sub_in_win, fa_ot[fa_sub_idx], ra, dec, "ICRS")
    fa_sub_leg_wsl = f"{TMP_WSL}/b2b_fmtA_50k_leg.fits"
    run_legacy_as1bary(win_to_wsl(fa_sub_in_win), fa_sub_leg_wsl, LAXPC_ORB_MKF_WSL, ra, dec, "ICRS")
    fa_lt = np.asarray(read_fits(str(wsl_to_win(fa_sub_leg_wsl)))[1].read_column("TIME"), dtype=np.float64)
    st_fmta = residual_stats_us(fa_rt[fa_sub_idx], fa_lt)
    # Also check whether legacy as1bary works natively on a real 100-row Format A header clone
    print(f"  B2b Format A (full 25.82M events in {t_fmta_rust:.2f}s; 50k strided oracle check): "
          f"max|d|={st_fmta['max_abs_us']:.3f} us, RMS={st_fmta['rms_us']:.3f} us, p99={st_fmta['p99_abs_us']:.3f} us")

    out["B2a_fmtB_AstroSat_orb"] = {
        "events_de405_vs_legacy": st_de,
        "events_de405_vs_legacy_frac_ge_0_5": st_de_hi,
        "events_de405_vs_legacy_frac_lt_0_5": st_de_lo,
        "events_analytical_vs_legacy": st_an,
        "gti_de405_vs_legacy": residual_stats_us(g_de, g_lg),
        "gti_analytical_vs_legacy": residual_stats_us(g_an, g_lg),
    }
    out["B2a_live_L2_orb_99k"] = {
        "raw_legacy_timesys_skipped": legacy_raw_timesys,
        "header_patched_de405_vs_legacy": st_live,
    }
    out["B2b_fmtA_full_25_8M"] = {
        "n_events_total": int(len(fa_ot)),
        "rust_de405_runtime_s": t_fmta_rust,
        "strided_50k_vs_legacy": st_fmta,
    }
    out["B2c_crab_timing"] = {
        "legacy": {k: crab_lg[k] for k in ("best_freq_hz", "z2_2", "chi2_32bins", "pulsed_amplitude_ratio")},
        "rust_de405": {k: crab_de[k] for k in ("best_freq_hz", "z2_2", "chi2_32bins", "pulsed_amplitude_ratio")},
        "rust_analytical": {k: crab_an[k] for k in ("best_freq_hz", "z2_2", "chi2_32bins", "pulsed_amplitude_ratio")},
        "phase_offset_de405_vs_legacy_cyc": dphi_de_lg,
        "phase_offset_de405_vs_legacy_us": dt_phase_de_us,
        "phase_offset_analytical_vs_legacy_cyc": dphi_an_lg,
        "phase_offset_analytical_vs_legacy_us": dt_phase_an_us,
    }
    if test_jpleph_fits:
        sub_rf = f"{TMP_WSL}/b2a_sub_liveorb_jpleph405.fits"
        res_f, _ = run_astrosat_bary(sub_in_wsl, sub_rf, LAXPC_ORB_L2_WSL, JPLEPH_405_WSL, ra, dec, check=False)
        if res_f.returncode == 0:
            srf = np.asarray(read_fits(str(wsl_to_win(sub_rf)))[1].read_column("TIME"), dtype=np.float64)
            out["B2a_live_fits_jpleph405"] = residual_stats_us(srf, sl)
            print(f"  B2a-live FITS JPLEPH.405 vs legacy: max|d|={out['B2a_live_fits_jpleph405']['max_abs_us']:.6f} us")
    return out


# ---------------------------------------------------------------------------
# Scenario B3: Leap-Second Boundary Crossing (2016-12-31T23:59:45 .. 2017-01-01T00:00:15 UTC)
# ---------------------------------------------------------------------------
def run_b3(test_jpleph_fits: bool = False) -> Dict[str, Any]:
    print("\n=== SCENARIO B3: 2016-12-31 -> 2017-01-01 Leap-Second Boundary Crossing ===")
    # MJD 57754.0 = 2017-01-01T00:00:00 UTC; AstroSat MJDREF = 55197.0
    t_leap = (57754 - 55197) * 86400.0  # 220924800.0 s
    # Dense grid across [-15.0, +15.0] s including exact +/-1ms, +/-1us, and 0.0 around leap second
    offsets = np.concatenate([
        np.linspace(-15.0, -1.0, 29),
        np.array([-0.5, -0.1, -0.01, -1e-3, -1e-5, 0.0, 1e-5, 1e-3, 0.01, 0.1, 0.5]),
        np.linspace(1.0, 15.0, 29),
    ])
    times_utc = t_leap + offsets

    orb_win = TMP_WIN / "b3_leap_orbit.orb"
    evt_win = TMP_WIN / "b3_leap_events.fits"
    write_synthetic_orbit_fits(orb_win, t_leap - 120.0, 241, deltat=1.0)
    write_synthetic_event_fits(
        evt_win, times_utc, ra_deg=83.633083, dec_deg=22.0145, radecsys="ICRS",
        gti_intervals=[(t_leap - 15.0, t_leap + 15.0)],
    )

    orb_wsl = win_to_wsl(orb_win)
    evt_wsl = win_to_wsl(evt_win)
    r_de_wsl = f"{TMP_WSL}/b3_rust_de405.fits"
    r_an_wsl = f"{TMP_WSL}/b3_rust_anal.fits"
    leg_wsl = f"{TMP_WSL}/b3_legacy.fits"

    run_astrosat_bary(evt_wsl, r_de_wsl, orb_wsl, DE405_BIN_WSL, 83.633083, 22.0145)
    run_astrosat_bary(evt_wsl, r_an_wsl, orb_wsl, "analytical", 83.633083, 22.0145)
    run_legacy_as1bary(evt_wsl, leg_wsl, orb_wsl, 83.633083, 22.0145, "ICRS")

    rt = np.asarray(read_fits(str(wsl_to_win(r_de_wsl)))[1].read_column("TIME"), dtype=np.float64)
    at = np.asarray(read_fits(str(wsl_to_win(r_an_wsl)))[1].read_column("TIME"), dtype=np.float64)
    lt = np.asarray(read_fits(str(wsl_to_win(leg_wsl)))[1].read_column("TIME"), dtype=np.float64)

    # Check exact 1.000000 s leap-second step between t_leap - 1e-5 and t_leap + 0.0
    idx_before = int(np.where(offsets == -1e-5)[0][0])
    idx_at = int(np.where(offsets == 0.0)[0][0])
    idx_m1 = int(np.where(offsets == -1.0)[0][0])
    idx_p1 = int(np.where(offsets == 1.0)[0][0])

    # Expected physical step between t_utc = t_leap - 1e-5 (TAI-UTC=36) and t_utc = t_leap (TAI-UTC=37)
    # is 1.000010 s + d(delay)/dt * 1.000010 s!
    step_rust = float(rt[idx_at] - rt[idx_before])
    step_leg = float(lt[idx_at] - lt[idx_before])
    # Compute smooth delay slope from [-15..-1] to isolate pure leap step
    delay_before = (rt[:idx_before + 1] - times_utc[:idx_before + 1]) - (36.0 + 32.184)
    delay_after = (rt[idx_at:] - times_utc[idx_at:]) - (37.0 + 32.184)
    # Extrapolate delay_before to t_leap and compare with delay_after[0]
    slope = (delay_before[idx_before] - delay_before[idx_m1]) / (times_utc[idx_before] - times_utc[idx_m1])
    # Note: ephemeris evaluates at TT which jumps by +1 s at leap second, so d(earth_roemer)/dt * 1s ~ O(1e-4 s)!
    continuity_rust_s = float((rt[idx_at] - rt[idx_before]) - 1e-5)
    continuity_leg_s = float((lt[idx_at] - lt[idx_before]) - 1e-5)

    st_de = residual_stats_us(rt, lt)
    st_an = residual_stats_us(at, lt)
    print(f"  B3 DE405 vs legacy: max|d|={st_de['max_abs_us']:.6f} us, RMS={st_de['rms_us']:.6f} us")
    print(f"  B3 leap step (t_leap - (t_leap-10us)): rust={step_rust:.9f} s, legacy={step_leg:.9f} s, diff={(step_rust - step_leg)*1e6:+.6f} us")
    print(f"  B3 span [-1s -> +1s]: rust_dt={rt[idx_p1]-rt[idx_m1]:.9f} s, leg_dt={lt[idx_p1]-lt[idx_m1]:.9f} s (expected ~3.000 s incl. 1s leap)")

    res = {
        "t_leap_utc": t_leap,
        "de405_vs_legacy": st_de,
        "analytical_vs_legacy": st_an,
        "step_across_leap_rust_s": step_rust,
        "step_across_leap_legacy_s": step_leg,
        "step_diff_us": (step_rust - step_leg) * 1e6,
        "span_m1_to_p1_rust_s": float(rt[idx_p1] - rt[idx_m1]),
        "span_m1_to_p1_legacy_s": float(lt[idx_p1] - lt[idx_m1]),
        "pre_leap_residual": residual_stats_us(rt[:idx_at], lt[:idx_at]),
        "post_leap_residual": residual_stats_us(rt[idx_at:], lt[idx_at:]),
    }
    if test_jpleph_fits:
        rf_wsl = f"{TMP_WSL}/b3_rust_jpleph405.fits"
        rf_cmd, _ = run_astrosat_bary(evt_wsl, rf_wsl, orb_wsl, JPLEPH_405_WSL, 83.633083, 22.0145, check=False)
        if rf_cmd.returncode == 0:
            rft = np.asarray(read_fits(str(wsl_to_win(rf_wsl)))[1].read_column("TIME"), dtype=np.float64)
            res["fits_jpleph405_vs_legacy"] = residual_stats_us(rft, lt)
            print(f"  B3 FITS JPLEPH.405 vs legacy: max|d|={res['fits_jpleph405_vs_legacy']['max_abs_us']:.6f} us")
    return res


# ---------------------------------------------------------------------------
# Scenario B4: Extreme Geometry & Edge Cases
# ---------------------------------------------------------------------------
def run_b4(test_jpleph_fits: bool = False) -> Dict[str, Any]:
    print("\n=== SCENARIO B4: Extreme Geometry & Edge Cases ===")
    out: Dict[str, Any] = {}

    # Shared base epoch: 2016-06-15T12:00:00 UTC (Sun is near Crab RA=83.63, Dec=22.01 -> solar conjunction!)
    # MJD 57554.5 -> AstroSat UTC = (57554.5 - 55197.0) * 86400 = 203688000.0 s
    t_conj = (57554.5 - 55197.0) * 86400.0
    orb_win = TMP_WIN / "b4_conj_orbit.orb"
    write_synthetic_orbit_fits(orb_win, t_conj - 600.0, 1201, deltat=1.0)
    orb_wsl = win_to_wsl(orb_win)
    times = t_conj + np.linspace(-300.0, 300.0, 61)

    # B4a: Near Solar Conjunction (Crab on 2016-06-15 is ~1.3 deg from Sun center -> huge Shapiro delay!)
    # Also test an even closer 0.35 deg grazing ray (just outside solar limb 0.267 deg): RA=84.1, Dec=23.0
    for tag, ra, dec, desc in [
        ("B4a_conjunction_crab_1p3deg", 83.633083, 22.0145, "Crab at 2016-06-15 solar conjunction (~1.3 deg from Sun)"),
        ("B4a_conjunction_limb_0p35deg", 84.10, 23.00, "Near-limb grazing ray (~0.35 deg from Sun, ~100 us Shapiro)"),
        ("B4d_ecliptic_pole", 270.0, 66.560708, "North Ecliptic Pole (minimal Earth Roemer)"),
        ("B4d_ecliptic_plane", 0.0, 0.0, "Vernal Equinox on Ecliptic Plane (max Roemer gradient)"),
    ]:
        ev_win = TMP_WIN / f"{tag}_in.fits"
        write_synthetic_event_fits(ev_win, times, ra, dec, "ICRS")
        ev_wsl = win_to_wsl(ev_win)
        r_wsl = f"{TMP_WSL}/{tag}_rust.fits"
        l_wsl = f"{TMP_WSL}/{tag}_leg.fits"
        run_astrosat_bary(ev_wsl, r_wsl, orb_wsl, DE405_BIN_WSL, ra, dec, "icrs")
        run_legacy_as1bary(ev_wsl, l_wsl, orb_wsl, ra, dec, "ICRS")
        rt = np.asarray(read_fits(str(wsl_to_win(r_wsl)))[1].read_column("TIME"), dtype=np.float64)
        lt = np.asarray(read_fits(str(wsl_to_win(l_wsl)))[1].read_column("TIME"), dtype=np.float64)
        st = residual_stats_us(rt, lt)
        out[tag] = {"description": desc, "ra": ra, "dec": dec, "de405_vs_legacy": st}
        print(f"  {tag}: max|d|={st['max_abs_us']:.6f} us, RMS={st['rms_us']:.6f} us")

    # B4b: FK5 -> ICRS frame conversion (c200to405)
    # Notice axBary.c:315 triggers convdir=1 ONLY when RA/Dec are read from header (RADECSYS='FK5') and -ref ICRS!
    ev_fk5_win = TMP_WIN / "b4b_fk5_in.fits"
    write_synthetic_event_fits(ev_fk5_win, times, 83.633083, 22.0145, radecsys="FK5")
    ev_fk5_wsl = win_to_wsl(ev_fk5_win)
    r_fk5_wsl = f"{TMP_WSL}/b4b_fk5_rust_de405.fits"
    l_fk5_icrs_wsl = f"{TMP_WSL}/b4b_fk5_leg_c200to405.fits"
    l_fk5_de200_wsl = f"{TMP_WSL}/b4b_fk5_leg_de200.fits"
    r_icrs_wsl = f"{TMP_WSL}/b4b_icrs_rust_de405.fits"

    run_astrosat_bary(ev_fk5_wsl, r_fk5_wsl, orb_wsl, DE405_BIN_WSL, 83.633083, 22.0145, "fk5")
    run_astrosat_bary(ev_fk5_wsl, r_icrs_wsl, orb_wsl, DE405_BIN_WSL, 83.633083, 22.0145, "icrs")
    # Pass NO -ra/-dec to legacy as1bary so it reads RA_NOM/DEC_NOM + RADECSYS='FK5' and runs c200to405 under -ref ICRS!
    run_legacy_as1bary(ev_fk5_wsl, l_fk5_icrs_wsl, orb_wsl, ra=None, dec=None, ref_frame="ICRS")
    # Also run legacy as1bary with -ref FK5 (uses JPLEPH.200!)
    run_legacy_as1bary(ev_fk5_wsl, l_fk5_de200_wsl, orb_wsl, ra=None, dec=None, ref_frame="FK5")

    rt_fk5 = np.asarray(read_fits(str(wsl_to_win(r_fk5_wsl)))[1].read_column("TIME"), dtype=np.float64)
    rt_icrs = np.asarray(read_fits(str(wsl_to_win(r_icrs_wsl)))[1].read_column("TIME"), dtype=np.float64)
    lt_conv = np.asarray(read_fits(str(wsl_to_win(l_fk5_icrs_wsl)))[1].read_column("TIME"), dtype=np.float64)
    lt_de200 = np.asarray(read_fits(str(wsl_to_win(l_fk5_de200_wsl)))[1].read_column("TIME"), dtype=np.float64)
    st_fk5_conv = residual_stats_us(rt_fk5, lt_conv)
    st_fk5_shift = residual_stats_us(rt_fk5, rt_icrs)
    st_de405_de200 = residual_stats_us(rt_fk5, lt_de200)
    out["B4b_fk5_c200to405"] = {
        "rust_fk5_de405_vs_legacy_c200to405": st_fk5_conv,
        "fk5_minus_icrs_shift_us": st_fk5_shift,
        "rust_fk5_de405_vs_legacy_de200": st_de405_de200,
    }
    if test_jpleph_fits:
        r_de200_wsl = f"{TMP_WSL}/b4b_fk5_rust_jpleph200.fits"
        res_200, _ = run_astrosat_bary(ev_fk5_wsl, r_de200_wsl, orb_wsl, JPLEPH_200_WSL, 83.633083, 22.0145, "fk5", check=False)
        if res_200.returncode == 0:
            rt_200 = np.asarray(read_fits(str(wsl_to_win(r_de200_wsl)))[1].read_column("TIME"), dtype=np.float64)
            st_200 = residual_stats_us(rt_200, lt_de200)
            out["B4b_fk5_c200to405"]["rust_fk5_jpleph200_vs_legacy_de200"] = st_200
            print(f"  B4b FITS JPLEPH.200 (--frame fk5) vs legacy (-ref FK5): max|d|={st_200['max_abs_us']:.6f} us")
    print(f"  B4b FK5->ICRS (c200to405): rust(--frame fk5) vs legacy(c200to405) max|d|={st_fk5_conv['max_abs_us']:.6f} us "
          f"(frame rotation magnitude={st_fk5_shift['median_us']:+.3f} us; vs DE200={st_de405_de200['median_us']:+.3f} us)")

    # B4c: Pre-existing TIMEZERO != 0 (small TIMEZERO = +12.345678 s AND large TIMEZERO = 200000000.0 s)
    tz_small = 12.345678
    ev_tz_s_win = TMP_WIN / "b4c_tz_small_in.fits"
    write_synthetic_event_fits(
        ev_tz_s_win, times - tz_small, 83.633083, 22.0145, "ICRS",
        timezero=tz_small,
        gti_intervals=[(float(times[0] - tz_small), float(times[-1] - tz_small))],
    )
    ev_tz_s_wsl = win_to_wsl(ev_tz_s_win)
    r_tz_s = f"{TMP_WSL}/b4c_tz_small_rust.fits"
    l_tz_s = f"{TMP_WSL}/b4c_tz_small_leg.fits"
    run_astrosat_bary(ev_tz_s_wsl, r_tz_s, orb_wsl, DE405_BIN_WSL, 83.633083, 22.0145, "icrs")
    run_legacy_as1bary(ev_tz_s_wsl, l_tz_s, orb_wsl, 83.633083, 22.0145, "ICRS")
    rt_tzs = read_fits(str(wsl_to_win(r_tz_s)))
    lt_tzs = read_fits(str(wsl_to_win(l_tz_s)))
    st_tzs_ev = residual_stats_us(rt_tzs[1].read_column("TIME"), lt_tzs[1].read_column("TIME"))
    st_tzs_gti = residual_stats_us(
        np.concatenate([rt_tzs[2].read_column("START"), rt_tzs[2].read_column("STOP")]),
        np.concatenate([lt_tzs[2].read_column("START"), lt_tzs[2].read_column("STOP")]),
    )
    out["B4c_timezero_small"] = {
        "timezero_in": tz_small,
        "timezero_out_rust": rt_tzs[1].header.get("TIMEZERO"),
        "timezero_out_legacy": lt_tzs[1].header.get("TIMEZERO"),
        "events_vs_legacy": st_tzs_ev,
        "gti_vs_legacy": st_tzs_gti,
    }
    print(f"  B4c TIMEZERO={tz_small}: events max|d|={st_tzs_ev['max_abs_us']:.6f} us, GTI max|d|={st_tzs_gti['max_abs_us']:.6f} us, "
          f"out TIMEZERO={rt_tzs[1].header.get('TIMEZERO')}")

    # B4c2: Large TIMEZERO (200,000,000 s, header TSTART already includes TIMEZERO -> heuristic applies!)
    tz_large = 200_000_000.0
    ev_tz_l_win = TMP_WIN / "b4c_tz_large_in.fits"
    write_synthetic_event_fits(
        ev_tz_l_win, times - tz_large, 83.633083, 22.0145, "ICRS",
        timezero=tz_large,
        tstart_hdr=float(times[0]), tstop_hdr=float(times[-1]),
    )
    ev_tz_l_wsl = win_to_wsl(ev_tz_l_win)
    r_tz_l = f"{TMP_WSL}/b4c_tz_large_rust.fits"
    l_tz_l = f"{TMP_WSL}/b4c_tz_large_leg.fits"
    run_astrosat_bary(ev_tz_l_wsl, r_tz_l, orb_wsl, DE405_BIN_WSL, 83.633083, 22.0145, "icrs")
    run_legacy_as1bary(ev_tz_l_wsl, l_tz_l, orb_wsl, 83.633083, 22.0145, "ICRS")
    rt_tzl = read_fits(str(wsl_to_win(r_tz_l)))
    lt_tzl = read_fits(str(wsl_to_win(l_tz_l)))
    st_tzl_ev = residual_stats_us(rt_tzl[1].read_column("TIME"), lt_tzl[1].read_column("TIME"))
    out["B4c_timezero_large"] = {
        "timezero_in": tz_large,
        "events_vs_legacy": st_tzl_ev,
        "primary_tstart_rust": rt_tzl[0].header.get("TSTART"),
        "primary_tstart_legacy": lt_tzl[0].header.get("TSTART"),
    }
    print(f"  B4c TIMEZERO={tz_large}: events max|d|={st_tzl_ev['max_abs_us']:.6f} us")

    # B4e: Exact jd == jd_end boundary!
    # 1) Test exact jd == jd_end of JPLEPH.405 (TSTOP = 2469808.5 JD(TDB)) in legacy as1bary vs Rust binary DE405 (and FITS JPLEPH.405 if supported)
    # 2) Test exact jd == jd_end of lnxp1600p2200.405 (JD 2525008.5) directly in Rust!
    # For JD_TDB = 2469808.5 -> MJD_TDB = 69808.0 -> t_tt ≈ (69808 - 55197)*86400 = 1262390400.0
    # Compute t_utc such that JD_TDB = 2469808.500000000 exactly (within < 1 ms)
    # At 2050-01-01, TAI-UTC = 37 s, TT-TAI = 32.184 s, TDB-TT ≈ -0.000083 s -> t_utc = 1262390400.0 - 69.184 + 0.000083
    t_end_jpl405_approx = (69808.0 - 55197.0) * 86400.0 - (37.0 + 32.184)
    # Test several points approaching exact jd_end: -10s, -1s, -0.01s, exact 0.0 (where jd_tdb <= 2469808.5)
    # Note: at 2050-01-01 00:00:00 TT (JD 2469808.5), TDB-TT is ~ -0.000083 s < 0, so JD_TDB = 2469808.499999999!
    # For exact JD_TDB == 2469808.5, add +0.000083 s? Better: test -1.0s, -0.001s, 0.0s (right at JD_TT=2469808.5)
    t_jpl405_grid = t_end_jpl405_approx + np.array([-10.0, -1.0, -0.01, -0.001, 0.0])
    orb_end1_win = TMP_WIN / "b4e_jpleph405_end.orb"
    ev_end1_win = TMP_WIN / "b4e_jpleph405_end_in.fits"
    write_synthetic_orbit_fits(orb_end1_win, t_end_jpl405_approx - 60.0, 121, deltat=1.0)
    write_synthetic_event_fits(ev_end1_win, t_jpl405_grid, 83.633083, 22.0145, "ICRS")
    r_e1 = f"{TMP_WSL}/b4e_jpleph405_end_rust_bin.fits"
    l_e1 = f"{TMP_WSL}/b4e_jpleph405_end_leg.fits"
    run_astrosat_bary(win_to_wsl(ev_end1_win), r_e1, win_to_wsl(orb_end1_win), DE405_BIN_WSL, 83.633083, 22.0145, "icrs")
    run_legacy_as1bary(win_to_wsl(ev_end1_win), l_e1, win_to_wsl(orb_end1_win), 83.633083, 22.0145, "ICRS")
    rt_e1 = np.asarray(read_fits(str(wsl_to_win(r_e1)))[1].read_column("TIME"), dtype=np.float64)
    lt_e1 = np.asarray(read_fits(str(wsl_to_win(l_e1)))[1].read_column("TIME"), dtype=np.float64)
    st_e1 = residual_stats_us(rt_e1, lt_e1)

    # Also test exact jd == jd_end of lnxp1600p2200.405 (read jd_end directly from binary header offset 2660!)
    with open(wsl_to_win(DE405_BIN_WSL), "rb") as fh:
        fh.seek(2652)
        bin_jd_start, bin_jd_end, bin_span = struct.unpack("<ddd", fh.read(24))
    # Compute t_utc so jd_tdb hits exact bin_jd_end (at JD 2525008.5 = Feb 2200, TDB-TT = +0.0011434 s)
    mjd_end_bin = bin_jd_end - 2400000.5
    t_end_bin_exact = (mjd_end_bin - 55197.0) * 86400.0 - (37.0 + 32.184) - 0.001144
    t_bin_grid = t_end_bin_exact + np.array([-10.0, -1.0, -0.001, 0.0])
    orb_end2_win = TMP_WIN / "b4e_bin405_end.orb"
    ev_end2_win = TMP_WIN / "b4e_bin405_end_in.fits"
    write_synthetic_orbit_fits(orb_end2_win, t_end_bin_exact - 60.0, 121, deltat=1.0)
    write_synthetic_event_fits(ev_end2_win, t_bin_grid, 83.633083, 22.0145, "ICRS")
    r_e2 = f"{TMP_WSL}/b4e_bin405_end_rust.fits"
    res_e2, _ = run_astrosat_bary(win_to_wsl(ev_end2_win), r_e2, win_to_wsl(orb_end2_win), DE405_BIN_WSL, 83.633083, 22.0145, "icrs", check=False)
    e2_ok = res_e2.returncode == 0
    e2_times = np.asarray(read_fits(str(wsl_to_win(r_e2)))[1].read_column("TIME"), dtype=np.float64).tolist() if e2_ok else None

    # Also verify jd_tdb > bin_jd_end (by +0.001 s) is cleanly rejected with TimeOutOfRange
    ev_end2_over_win = TMP_WIN / "b4e_bin405_over_in.fits"
    write_synthetic_event_fits(ev_end2_over_win, np.array([t_end_bin_exact + 0.01]), 83.633083, 22.0145, "ICRS")
    res_e2_over, _ = run_astrosat_bary(win_to_wsl(ev_end2_over_win), f"{TMP_WSL}/b4e_bin405_over_rust.fits", win_to_wsl(orb_end2_win), DE405_BIN_WSL, 83.633083, 22.0145, "icrs", check=False)

    out["B4e_jd_end_boundary"] = {
        "jpleph405_jd_end": 2469808.5,
        "jpleph405_boundary_vs_legacy": st_e1,
        "bin_de405_jd_end": bin_jd_end,
        "bin_de405_boundary_ok": e2_ok,
        "bin_de405_boundary_times_tdb": e2_times,
        "bin_de405_over_rejected": res_e2_over.returncode != 0,
        "bin_de405_boundary_stderr": res_e2.stderr.strip() if not e2_ok else "",
    }
    print(f"  B4e jd_end boundary: JPLEPH.405 (JD 2469808.5) max|d|={st_e1['max_abs_us']:.6f} us | "
          f"bin DE405 (JD {bin_jd_end}) ok={e2_ok}, over_rejected={res_e2_over.returncode != 0}")

    if test_jpleph_fits:
        r_e1f = f"{TMP_WSL}/b4e_jpleph405_end_rust_fits.fits"
        res_f, _ = run_astrosat_bary(win_to_wsl(ev_end1_win), r_e1f, win_to_wsl(orb_end1_win), JPLEPH_405_WSL, 83.633083, 22.0145, "icrs", check=False)
        if res_f.returncode == 0:
            rt_e1f = np.asarray(read_fits(str(wsl_to_win(r_e1f)))[1].read_column("TIME"), dtype=np.float64)
            out["B4e_fits_jpleph405_boundary_vs_legacy"] = residual_stats_us(rt_e1f, lt_e1)
            print(f"  B4e FITS JPLEPH.405 boundary vs legacy: max|d|={out['B4e_fits_jpleph405_boundary_vs_legacy']['max_abs_us']:.6f} us")
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description="Head-to-head astrosat bary vs legacy as1bary verification")
    parser.add_argument("--scenarios", nargs="+", default=["B1", "B2", "B3", "B4"], choices=["B1", "B2", "B3", "B4"])
    parser.add_argument("--no-reuse", action="store_true", help="Re-run even if B2a files exist")
    parser.add_argument("--test-jpleph-fits", action="store_true", help="Also test FITS-wrapped JPLEPH.405 in astrosat bary")
    parser.add_argument("--json-out", default=str(TMP_WIN / "head2head_results.json"))
    args = parser.parse_args()

    TMP_WIN.mkdir(parents=True, exist_ok=True)
    results: Dict[str, Any] = {}
    if "B1" in args.scenarios:
        results["B1"] = run_b1(test_jpleph_fits=args.test_jpleph_fits)
    if "B2" in args.scenarios:
        results["B2"] = run_b2(reuse_existing=not args.no_reuse, test_jpleph_fits=args.test_jpleph_fits)
    if "B3" in args.scenarios:
        results["B3"] = run_b3(test_jpleph_fits=args.test_jpleph_fits)
    if "B4" in args.scenarios:
        results["B4"] = run_b4(test_jpleph_fits=args.test_jpleph_fits)

    with open(args.json_out, "w", encoding="utf-8") as fh:
        json.dump(results, fh, indent=2)
    print(f"\nSaved full JSON results to {args.json_out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
