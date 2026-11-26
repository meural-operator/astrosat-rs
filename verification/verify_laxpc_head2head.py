#!/usr/bin/env python3
"""Head-to-head numerical verification of `astrosat laxpc` against compiled legacy
LaxpcSoft v3.4.6 (`laxpcl1.e`, `backshiftv3.e`), official AstroSat L2 reference products,
and synthetic ground-truth / adversarial edge scenarios.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import re
import shutil
import struct
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fits_reader import read_fits

ASTROSAT_BIN = "/mnt/d/ISRO/target/release/astrosat"
LAXPCL1_BIN = "/mnt/d/ISRO/.tmp/legacy_bins/laxpc/laxpcl1.e"
BACKSHIFT_BIN = "/mnt/d/ISRO/.tmp/legacy_bins/laxpc/backshiftv3.e"
FITSVERIFY_BIN = "/mnt/d/ISRO/.toolchain/bin/fitsverify"

REF_SRC = "/mnt/d/ISRO/data/ref_src/laxpc"
CALDB_DIR = "/mnt/d/ISRO/data/caldb/laxpc"
CRAB_L1_DIR = "/mnt/d/ISRO/data/sample/laxpc_l1/20160107_G02_010T01_9000000252_level1"
CRAB_L1_ORB = f"{CRAB_L1_DIR}/laxpc/01509"
CRAB_L2_REF = "/mnt/d/ISRO/data/sample/laxpc_l2/20160107_G02_010T01_9000000252_level2/laxpc/01509"

WORK_ROOT = "/mnt/d/ISRO/.tmp/legacy_bins/laxpc/head2head"

UNITS = [
    {"id": 1, "prefix": "lxp1", "det": "lxp10", "def_ch": 512},
    {"id": 2, "prefix": "lxp2", "det": "lxp20", "def_ch": 256},
    {"id": 3, "prefix": "lxp3", "det": "lxp30", "def_ch": 512},
]


def ensure_clean_dir(path: str) -> str:
    if os.path.exists(path):
        shutil.rmtree(path)
    os.makedirs(path, exist_ok=True)
    return path


def run_cmd(cmd: list[str], cwd: str | None = None, stdin_text: str | None = None, timeout: int = 180) -> subprocess.CompletedProcess:
    return subprocess.run(
        cmd,
        cwd=cwd,
        input=stdin_text,
        text=True,
        capture_output=True,
        timeout=timeout,
    )


def pearson_r(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    if len(a) != len(b) or len(a) == 0:
        return 0.0
    da = a - np.mean(a)
    db = b - np.mean(b)
    denom = float(np.sqrt(np.sum(da * da) * np.sum(db * db)))
    return float(np.sum(da * db) / denom) if denom > 0.0 else (1.0 if np.allclose(a, b) else 0.0)


def compare_arrays(a: np.ndarray, b: np.ndarray) -> dict[str, float]:
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    n = min(len(a), len(b))
    if n == 0:
        return {"len_a": float(len(a)), "len_b": float(len(b)), "pearson_r": 0.0}
    a = a[:n]
    b = b[:n]
    diff = a - b
    abs_diff = np.abs(diff)
    denom = np.maximum(np.abs(a), np.abs(b))
    rel = np.divide(abs_diff, denom, out=np.zeros_like(abs_diff), where=denom > 1e-12)
    return {
        "n": int(n),
        "sum_a": float(np.sum(a)),
        "sum_b": float(np.sum(b)),
        "max_abs": float(np.max(abs_diff)),
        "rms": float(np.sqrt(np.mean(diff * diff))),
        "max_rel": float(np.max(rel)),
        "mean_rel": float(np.mean(rel)),
        "pearson_r": pearson_r(a, b),
    }


def fitsverify_check(path: str) -> tuple[bool, str]:
    if not os.path.isfile(FITSVERIFY_BIN) or not os.path.isfile(path):
        return False, "missing binary or file"
    p = run_cmd([FITSVERIFY_BIN, "-q", path])
    ok = (p.returncode == 0) and ("0 errors" in p.stdout or "verification OK" in p.stdout)
    return ok, (p.stdout + p.stderr).strip()


def write_laxpcl1_inputs(
    work_dir: str,
    orbits: list[tuple[str, str, str, str]],
    gti_ins_intervals: list[tuple[float, float]] | None = None,
    bkg_model_lxp10: str = f"{REF_SRC}/backfit.lxp10mar16",
    bkg_pha_lxp10: str = f"{REF_SRC}/backlxp10jan16.pha",
) -> None:
    with open(os.path.join(work_dir, "laxpcl1.inp"), "w") as f:
        for mkf, tct, bb, ea in orbits:
            f.write(f"{mkf}\n{tct}\n{bb}\n{ea}\n")
    with open(os.path.join(work_dir, "back4.inp"), "w") as f:
        f.write(f"{bkg_model_lxp10}\n")
        for k in range(1, 5):
            f.write(f"{REF_SRC}/backlxp10mar16ul{k}.pha\n")
        f.write(f"{bkg_pha_lxp10}\n")
        f.write(f"{REF_SRC}/backfit3.lxp10mar16\n")
    with open(os.path.join(work_dir, "gti.inp"), "w") as f:
        if not gti_ins_intervals:
            f.write("1\n0.0 1.0e12\n")
        else:
            f.write(f"{len(gti_ins_intervals)}\n")
            for a, b in gti_ins_intervals:
                f.write(f"{a:.14e} {b:.14e}\n")


def run_laxpcl1(
    work_dir: str,
    lxp_id: int,
    tbin: float,
    ian: int,
    icl: int,
    icu: int,
    nca: int,
    nul: int,
    je: int,
    iev: int,
) -> subprocess.CompletedProcess:
    stdin_text = f"{lxp_id} {tbin} {ian}\n{icl} {icu}\n{nca}\n{nul} {je} {iev}\n"
    return run_cmd([LAXPCL1_BIN], cwd=work_dir, stdin_text=stdin_text, timeout=180)


def parse_laxpcl1_frame_file(path: str) -> dict[str, Any]:
    out = {}
    with open(path, "r") as f:
        for line in f:
            if line.strip().startswith("Frame sum:"):
                parts = line.split()
                # Frame sum: lxp1  idobs  norb  np+nloss  nloss  ngap
                if len(parts) >= 8:
                    out["unit"] = parts[2]
                    out["total_frames_incl_loss"] = int(parts[5])
                    out["nloss"] = int(parts[6])
                    out["ngap"] = int(parts[7])
    return out


def parse_laxpcl1_gti(path: str) -> list[tuple[float, float, float, float]]:
    rows = []
    with open(path, "r") as f:
        lines = [ln.strip() for ln in f if ln.strip()]
    if not lines:
        return rows
    n = int(lines[0].split()[0])
    for ln in lines[1 : 1 + n]:
        vals = [float(x) for x in ln.split()]
        if len(vals) >= 4:
            rows.append((vals[0], vals[1], vals[2], vals[3]))
    return rows


def parse_laxpcl1_stdout_metrics(stdout: str) -> dict[str, Any]:
    res: dict[str, Any] = {}
    for line in stdout.splitlines():
        if "tobs=" in line:
            # tobs=   2487.0000000000000        0.0000000000000000
            m = re.findall(r"[-+]?\d*\.\d+(?:[eEdD][-+]?\d+)?|[-+]?\d+", line)
            if len(m) >= 2:
                res["tobs_gross"] = float(m[0].replace("D", "E"))
                res["tgap"] = float(m[1].replace("D", "E"))
        if line.strip().startswith("count "):
            m = re.findall(r"[-+]?\d*\.\d+(?:[eEdD][-+]?\d+)?|[-+]?\d+", line)
            if len(m) >= 2:
                res["spec_scount"] = float(m[0].replace("D", "E"))
                res["lc_icount"] = float(m[1].replace("D", "E"))
    return res


def parse_ascii_pha(path: str) -> dict[str, Any]:
    with open(path, "r") as f:
        lines = [ln.strip() for ln in f if ln.strip()]
    hdr = lines[0]
    m_exp = re.search(r"Exposure\s*=\s*([-+0-9.eEdD]+)", hdr)
    m_rate = re.search(r"rate:\s*([-+0-9.eEdD]+)", hdr)
    m_day = re.search(r"day:\s*([0-9]+)", hdr)
    rows = []
    for ln in lines[1:]:
        if ln.startswith("#"):
            continue
        parts = [float(x.replace("D", "E")) for x in ln.split()]
        if len(parts) >= 12:
            rows.append(parts[1:12])  # 10 anodes + total
    arr = np.array(rows, dtype=np.float64)
    return {
        "header": hdr,
        "exposure": float(m_exp.group(1)) if m_exp else 0.0,
        "tot_corr_rate": float(m_rate.group(1)) if m_rate else 0.0,
        "day": int(m_day.group(1)) if m_day else 0,
        "rates": arr,
    }


def read_spec_fits(path: str) -> dict[str, Any]:
    hdus = read_fits(path)
    sh = next(h for h in hdus if h.extname.upper() == "SPECTRUM")
    cols = {c[0].upper(): c[0] for c in sh.bintable_columns()}
    vec = np.asarray(sh.read_column(cols["RATE"] if "RATE" in cols else cols["COUNTS"]), dtype=np.float64)
    err = np.asarray(sh.read_column(cols["STAT_ERR"]), dtype=np.float64) if "STAT_ERR" in cols else np.zeros_like(vec)
    return {
        "exposure": float(sh.header.get("EXPOSURE", 0.0)),
        "deadtime": float(sh.header.get("DEADTIME", 0.0)),
        "respfile": str(sh.header.get("RESPFILE", "")),
        "backfile": str(sh.header.get("BACKFILE", "")),
        "rate": vec,
        "err": err,
        "is_counts": "COUNTS" in cols and "RATE" not in cols,
    }


def read_lc_fits(path: str) -> dict[str, Any]:
    hdus = read_fits(path)
    rh = next(h for h in hdus if h.xtension == "BINTABLE" and h.extname.upper() != "GTI")
    cols = {c[0].strip().upper(): c[0] for c in rh.bintable_columns()}
    t = np.asarray(rh.read_column(cols["TIME"]), dtype=np.float64)
    r = np.asarray(rh.read_column(cols["RATE"]), dtype=np.float64)
    err_key = next((cols[k] for k in ("ERROR", "STAT_ERR", "ERR") if k in cols), None)
    e = np.asarray(rh.read_column(err_key), dtype=np.float64) if err_key else np.zeros_like(r)
    f = np.asarray(rh.read_column(cols["FRACEXP"]), dtype=np.float64) if "FRACEXP" in cols else np.ones_like(r)
    return {
        "tstart": float(rh.header.get("TSTART", 0.0)),
        "tstop": float(rh.header.get("TSTOP", 0.0)),
        "timedel": float(rh.header.get("TIMEDEL", 0.0)),
        "exposure": float(rh.header.get("EXPOSURE", 0.0)),
        "backapp": rh.header.get("BACKAPP"),
        "time": t,
        "rate": r,
        "error": e,
        "fracexp": f,
    }


def read_evn_fits(path: str) -> dict[str, Any]:
    hdus = read_fits(path)
    eh = next(h for h in hdus if h.xtension == "BINTABLE" and h.extname.upper() not in ("GTI", "BTI"))
    cols = {c[0].upper(): c[0] for c in eh.bintable_columns()}
    t = np.asarray(eh.read_column(cols["TIME"]), dtype=np.float64)
    ch_col = "PHA" if "PHA" in cols else "CHANNEL"
    an_col = next((cols[k] for k in ("ANODE_ID", "ANODE ID", "LAYER") if k in cols), None)
    det_col = next((cols[k] for k in ("DETECTOR_ID", "DETECTOR  ID", "LAXPC_NO") if k in cols), None)
    gti_h = next((h for h in hdus if h.extname.upper() == "GTI"), None)
    gti_exp = 0.0
    if gti_h is not None and int(gti_h.header.get("NAXIS2", 0)) > 0:
        gti_exp = float(np.sum(gti_h.read_column("STOP") - gti_h.read_column("START")))
    return {
        "nrows": len(t),
        "time": t,
        "ch": np.asarray(eh.read_column(cols[ch_col]), dtype=np.int64),
        "an": np.asarray(eh.read_column(an_col), dtype=np.int64) if an_col else np.zeros(len(t), dtype=np.int64),
        "det": np.asarray(eh.read_column(det_col), dtype=np.int64) if det_col else np.zeros(len(t), dtype=np.int64),
        "gti_exp": gti_exp,
    }


def read_crab_tct() -> tuple[np.ndarray, np.ndarray]:
    tct_path = f"{CRAB_L1_ORB}/aux/AS1G02_010T01_9000000252lxp_level1.tct"
    th = read_fits(tct_path)[1]
    cols = [c[0] for c in th.bintable_columns()]
    tsps = np.asarray(th.read_column(cols[0]), dtype=np.float64)
    tins = np.asarray(th.read_column(cols[2]), dtype=np.float64)
    # Apply same TCT glitch filter as laxpcl1.f:692-703 / tct.rs
    kept_s = [tsps[0]]
    kept_i = [tins[0]]
    for k in range(1, len(tsps)):
        if tsps[k] <= 1e5:
            continue
        dt_s = tsps[k] - kept_s[-1]
        dt_i = tins[k] - kept_i[-1]
        if dt_i <= 0:
            continue
        ratio = dt_s / dt_i
        if 1.000 < ratio <= 1.001:
            kept_s.append(tsps[k])
            kept_i.append(tins[k])
    return np.array(kept_s), np.array(kept_i)


def crab_orbits() -> list[tuple[str, str, str, str]]:
    mkf = f"{CRAB_L1_ORB}/AS1G02_010T01_9000000252lxp_level1.mkf"
    tct = f"{CRAB_L1_ORB}/aux/AS1G02_010T01_9000000252lxp_level1.tct"
    bb = f"{CRAB_L1_ORB}/lxp1/modeBB/AS1G02_010T01_9000000252lxp1BB_level1.fits"
    ea = f"{CRAB_L1_ORB}/lxp1/modeEA/AS1G02_010T01_9000000252lxp1EA_level1.fits"
    return [(mkf, tct, bb, ea)]


# ---------------------------------------------------------------------------
# Scenario L1: Real Crab 9000000252 (All 3 Units LX10, LX20, LX30)
# ---------------------------------------------------------------------------
def run_scenario_l1() -> dict[str, Any]:
    print("=== Running Scenario L1: Crab 9000000252 (LX10, LX20, LX30) ===")
    l1_dir = ensure_clean_dir(f"{WORK_ROOT}/L1")
    orbits = crab_orbits()
    tct_utc, tct_ins = read_crab_tct()

    # 1. Run Rust astrosat laxpc run (auto-GTI, default channels, Format A + B + Orbit + ASCII + Bkg)
    rust_auto_dir = ensure_clean_dir(f"{l1_dir}/rust_auto_1s")
    p_rust = run_cmd([
        ASTROSAT_BIN, "laxpc", "run",
        "--input", CRAB_L1_DIR,
        "--output", rust_auto_dir,
        "--units", "all",
        "--tbin", "1.0",
        "--caldb", REF_SRC,
        "--bkg-model", f"{REF_SRC}/backfit.lxp10mar16",
        "--bkg-pha", f"{REF_SRC}/backlxp10jan16.pha",
        "--format-a", "--format-b", "--orbit", "--ascii",
        "--json",
    ])
    if p_rust.returncode != 0:
        raise RuntimeError(f"Rust L1 auto failed: {p_rust.stderr}\n{p_rust.stdout}")
    rust_auto_json = json.loads(p_rust.stdout)

    # Run Rust high-res tbin = 0.01s
    rust_hi_dir = ensure_clean_dir(f"{l1_dir}/rust_auto_0p01s")
    p_rust_hi = run_cmd([
        ASTROSAT_BIN, "laxpc", "run",
        "--input", CRAB_L1_DIR,
        "--output", rust_hi_dir,
        "--units", "all",
        "--tbin", "0.01",
        "--caldb", REF_SRC,
        "--bkg-model", f"{REF_SRC}/backfit.lxp10mar16",
        "--bkg-pha", f"{REF_SRC}/backlxp10jan16.pha",
        "--no-format-b", "--no-orbit",
        "--json",
    ])
    if p_rust_hi.returncode != 0:
        raise RuntimeError(f"Rust L1 0.01s failed: {p_rust_hi.stderr}")
    rust_hi_json = json.loads(p_rust_hi.stdout)

    results: dict[str, Any] = {"units": {}, "orbit": {}, "format_a": {}, "fitsverify": {}}

    # Validate FITS with fitsverify
    for fn in os.listdir(rust_auto_dir):
        if fn.endswith((".spec", ".lcurv", ".lcbk", ".evn", ".gtis", ".orb", ".fits")):
            ok, msg = fitsverify_check(os.path.join(rust_auto_dir, fn))
            results["fitsverify"][fn] = {"ok": ok, "msg": msg[:200]}

    for u in UNITS:
        uid = u["id"]
        pref = u["prefix"]
        def_ch = u["def_ch"]
        u_res: dict[str, Any] = {}

        # A. Legacy Pass 1 (auto-GTI generation)
        leg_p1 = ensure_clean_dir(f"{l1_dir}/leg_{pref}_p1")
        write_laxpcl1_inputs(leg_p1, orbits, None)
        p1 = run_laxpcl1(leg_p1, uid, 1.0, 0, 0, 1023, def_ch, -1, 1, 2)
        if p1.returncode != 0:
            raise RuntimeError(f"laxpcl1 pass1 {pref} failed: {p1.stderr}")
        auto_gti_rows = parse_laxpcl1_gti(f"{leg_p1}/{pref}level2.gti")
        ins_gtis = [(r[0], r[1]) for r in auto_gti_rows]
        utc_gtis_from_leg = [(r[2], r[3]) for r in auto_gti_rows]
        u_res["legacy_auto_gti"] = utc_gtis_from_leg

        # B. Legacy Pass 2 (auto-GTI screened, 1.0s, with .evn, .spec, .lcurv, .lcbk, AstroSat.orb)
        leg_p2 = ensure_clean_dir(f"{l1_dir}/leg_{pref}_p2_1s")
        write_laxpcl1_inputs(leg_p2, orbits, ins_gtis)
        p2 = run_laxpcl1(leg_p2, uid, 1.0, 0, 0, 1023, def_ch, -1, 1, -2)
        if p2.returncode != 0:
            raise RuntimeError(f"laxpcl1 pass2 {pref} failed: {p2.stderr}")
        leg_frame = parse_laxpcl1_frame_file(f"{leg_p2}/{pref}level2.frame")
        leg_out_m = parse_laxpcl1_stdout_metrics(p2.stdout)

        # C. Legacy Pass 2 high-res tbin=0.01s
        leg_hi = ensure_clean_dir(f"{l1_dir}/leg_{pref}_p2_0p01s")
        write_laxpcl1_inputs(leg_hi, orbits, ins_gtis)
        p_hi = run_laxpcl1(leg_hi, uid, 0.01, 0, 0, 1023, def_ch, -1, 1, 2)

        # Compare auto-GTI intervals & frames & exposures
        rust_u = next(x for x in rust_auto_json["unit_results"] if x["unit"].lower().endswith(f"{uid}0"))
        rust_gtis_hdu = read_fits(f"{rust_auto_dir}/{pref}level2.gtis")[1]
        rust_gtis = list(zip(
            rust_gtis_hdu.read_column("START").tolist(),
            rust_gtis_hdu.read_column("STOP").tolist(),
        ))
        u_res["rust_auto_gti"] = rust_gtis
        leg_gross_gti = sum(b - a for a, b in utc_gtis_from_leg)
        rust_gross_gti = rust_u["gti"]["gross_exposure_sec"]
        u_res["frames"] = {
            "rust_total_raw": rust_u["ea_frame_diagnostics"]["total_raw_frames"],
            "rust_retained": rust_u["ea_frame_diagnostics"]["retained_frames"],
            "rust_missing": rust_u["ea_frame_diagnostics"]["total_missing_frames"],
            "leg_total_incl_loss": leg_frame.get("total_frames_incl_loss"),
            "leg_nloss": leg_frame.get("nloss"),
            "leg_ngap": leg_frame.get("ngap"),
        }
        u_res["exposure_auto"] = {
            "leg_gross_ins": leg_out_m.get("tobs_gross"),
            "leg_tgap_ins": leg_out_m.get("tgap"),
            "leg_net_ins": (leg_out_m.get("tobs_gross", 0) - leg_out_m.get("tgap", 0)),
            "leg_gross_utc": leg_gross_gti,
            "rust_gross": rust_gross_gti,
            "rust_gap_loss": rust_u["gti"]["frame_gap_loss_sec"],
            "rust_net": rust_u["gti"]["net_exposure_sec"],
        }

        # D. Compare auto-GTI .evn, .spec, .lcurv, .lcbk
        leg_evn = read_evn_fits(f"{leg_p2}/{pref}level2.evn")
        rust_evn = read_evn_fits(f"{rust_auto_dir}/{pref}level2.evn")
        u_res["events_auto"] = {
            "leg_count": leg_evn["nrows"],
            "rust_count": rust_evn["nrows"],
            "rel_diff_pct": 100.0 * abs(rust_evn["nrows"] - leg_evn["nrows"]) / max(leg_evn["nrows"], 1),
        }

        leg_spec = read_spec_fits(f"{leg_p2}/{pref}level2.spec")
        rust_spec = read_spec_fits(f"{rust_auto_dir}/{pref}level2.spec")
        leg_pha_asc = parse_ascii_pha(f"{leg_p2}/{pref}level2.pha")
        u_res["spec_auto_defch"] = {
            "channels": def_ch,
            "leg_tot_corr_rate": leg_pha_asc["tot_corr_rate"],
            "rust_tot_corr_rate": rust_u["spectrum"]["total_corrected_rate_cps"],
            "rust_fcd": rust_u["spectrum"]["deadtime_factor"],
            "cmp": compare_arrays(rust_spec["rate"], leg_spec["rate"]),
        }

        # E. Apples-to-apples run: Feed Legacy's EXACT auto-GTI (`utc_gtis_from_leg`) into Rust!
        rust_matched_dir = ensure_clean_dir(f"{l1_dir}/rust_matched_gti_{pref}")
        gti_txt = f"{rust_matched_dir}/leg_gti.txt"
        with open(gti_txt, "w") as gf:
            for a, b in utc_gtis_from_leg:
                gf.write(f"{a:.9f} {b:.9f}\n")
        p_rm = run_cmd([
            ASTROSAT_BIN, "laxpc", "run",
            "--input", CRAB_L1_DIR,
            "--output", rust_matched_dir,
            "--units", str(uid * 10),
            "--tbin", "1.0",
            "--user-gti", gti_txt, "--user-gti-only",
            "--caldb", REF_SRC,
            "--bkg-model", f"{REF_SRC}/backfit.lxp10mar16",
            "--bkg-pha", f"{REF_SRC}/backlxp10jan16.pha",
            "--format-b", "--no-orbit", "--ascii", "--json",
        ])
        if p_rm.returncode != 0:
            raise RuntimeError(f"Rust matched-GTI {pref} failed: {p_rm.stderr}")
        rm_json = json.loads(p_rm.stdout)["unit_results"][0]
        rm_evn = read_evn_fits(f"{rust_matched_dir}/{pref}level2.evn")
        rm_spec = read_spec_fits(f"{rust_matched_dir}/{pref}level2.spec")
        rm_lc = read_lc_fits(f"{rust_matched_dir}/{pref}level2.lcurv")
        leg_lc = read_lc_fits(f"{leg_p2}/{pref}level2.lcurv")
        rm_lcbk = read_lc_fits(f"{rust_matched_dir}/{pref}level2.lcbk")
        leg_lcbk = read_lc_fits(f"{leg_p2}/{pref}level2.lcbk")

        # Event-by-event timestamp comparison on matched GTI
        n_ev_min = min(leg_evn["nrows"], rm_evn["nrows"])
        ev_time_diff = float(np.max(np.abs(rm_evn["time"][:n_ev_min] - leg_evn["time"][:n_ev_min]))) if (leg_evn["nrows"] == rm_evn["nrows"] and n_ev_min > 0) else -1.0
        ev_ch_mismatch = int(np.sum(rm_evn["ch"][:n_ev_min] != leg_evn["ch"][:n_ev_min])) if leg_evn["nrows"] == rm_evn["nrows"] else -1

        # Match lightcurve bins by nearest timestamp within 0.05s
        t_diff_lc = np.abs(rm_lc["time"][:, None] - leg_lc["time"][None, :]) if len(rm_lc["time"]) < 5000 and len(leg_lc["time"]) < 5000 else None
        lc_cmp = {}
        lc_err_cmp = {}
        lcbk_cmp = {}
        if t_diff_lc is not None and t_diff_lc.size > 0:
            idx_leg = np.argmin(t_diff_lc, axis=1)
            mask = t_diff_lc[np.arange(len(rm_lc["time"])), idx_leg] < 0.05
            lc_cmp = compare_arrays(rm_lc["rate"][mask], leg_lc["rate"][idx_leg[mask]])
            lc_err_cmp = compare_arrays(rm_lc["error"][mask], leg_lc["error"][idx_leg[mask]])
            lc_cmp["matched_bins"] = int(np.sum(mask))
            lc_cmp["rust_bins"] = int(len(rm_lc["time"]))
            lc_cmp["leg_bins"] = int(len(leg_lc["time"]))
            lcbk_cmp = compare_arrays(rm_lcbk["rate"][mask], leg_lcbk["rate"][idx_leg[mask]])
            lcbk_cmp["rust_backapp"] = rm_lcbk["backapp"]
            lcbk_cmp["leg_backapp"] = leg_lcbk["backapp"]

        # High-res 0.01s matched-GTI comparison (testing 1-s running-mean deadtime smoothing)
        rust_matched_hi = ensure_clean_dir(f"{l1_dir}/rust_matched_0p01s_{pref}")
        run_cmd([
            ASTROSAT_BIN, "laxpc", "run",
            "--input", CRAB_L1_DIR,
            "--output", rust_matched_hi,
            "--units", str(uid * 10),
            "--tbin", "0.01",
            "--user-gti", gti_txt, "--user-gti-only",
            "--caldb", REF_SRC,
            "--bkg-model", f"{REF_SRC}/backfit.lxp10mar16",
            "--bkg-pha", f"{REF_SRC}/backlxp10jan16.pha",
            "--no-format-b", "--no-orbit", "--json",
        ])
        rm_hi_lc = read_lc_fits(f"{rust_matched_hi}/{pref}level2.lcurv")
        leg_hi_lc = read_lc_fits(f"{leg_hi}/{pref}level2.lcurv")
        # Align high-res LC by nearest timestamp within 0.002s via searchsorted
        order_l = np.argsort(leg_hi_lc["time"])
        lt_sorted = leg_hi_lc["time"][order_l]
        lr_sorted = leg_hi_lc["rate"][order_l]
        pos = np.searchsorted(lt_sorted, rm_hi_lc["time"])
        best_r = []
        best_l = []
        for idx_r, p in enumerate(pos):
            cand = []
            if p > 0:
                cand.append(p - 1)
            if p < len(lt_sorted):
                cand.append(p)
            if cand:
                j_best = min(cand, key=lambda c: abs(lt_sorted[c] - rm_hi_lc["time"][idx_r]))
                if abs(lt_sorted[j_best] - rm_hi_lc["time"][idx_r]) < 0.002:
                    best_r.append(rm_hi_lc["rate"][idx_r])
                    best_l.append(lr_sorted[j_best])
        hi_cmp = compare_arrays(np.array(best_r), np.array(best_l))
        hi_cmp["rust_bins"] = len(rm_hi_lc["time"])
        hi_cmp["leg_bins"] = len(leg_hi_lc["time"])
        hi_cmp["rust_mean_cps"] = float(np.mean(rm_hi_lc["rate"]))
        hi_cmp["leg_mean_cps"] = float(np.mean(leg_hi_lc["rate"]))

        # F. Multi-channel spectra comparison (1024, 512, 256) on matched GTI
        # Plus Official L2 reference comparison using Official L2 GTI!
        off_gti_p = f"{CRAB_L2_REF}/{pref}/modeEA/AS1G02_010T01_9000000252{pref}EA_level2.gti"
        off_pha_p = f"{CRAB_L2_REF}/{pref}/modeEA/AS1G02_010T01_9000000252{pref}EA_level2.pha"
        off_lc_p = f"{CRAB_L2_REF}/{pref}/modeEA/AS1G02_010T01_9000000252{pref}EA_level2.lc"
        off_gti_hdu = read_fits(off_gti_p)[1]
        off_gti_utc = list(zip(
            off_gti_hdu.read_column("START").tolist(),
            off_gti_hdu.read_column("STOP").tolist(),
        ))
        off_gti_ins = []
        for a, b in off_gti_utc:
            if a <= tct_utc[0]:
                ia = float(tct_ins[0] + (a - tct_utc[0]) * (tct_ins[1] - tct_ins[0]) / (tct_utc[1] - tct_utc[0]))
            elif a >= tct_utc[-1]:
                ia = float(tct_ins[-1] + (a - tct_utc[-1]) * (tct_ins[-1] - tct_ins[-2]) / (tct_utc[-1] - tct_utc[-2]))
            else:
                ia = float(np.interp(a, tct_utc, tct_ins))
            if b <= tct_utc[0]:
                ib = float(tct_ins[0] + (b - tct_utc[0]) * (tct_ins[1] - tct_ins[0]) / (tct_utc[1] - tct_utc[0]))
            elif b >= tct_utc[-1]:
                ib = float(tct_ins[-1] + (b - tct_utc[-1]) * (tct_ins[-1] - tct_ins[-2]) / (tct_utc[-1] - tct_utc[-2]))
            else:
                ib = float(np.interp(b, tct_utc, tct_ins))
            off_gti_ins.append((ia, ib))
        off_gti_txt = f"{rust_matched_dir}/off_gti.txt"
        with open(off_gti_txt, "w") as gf:
            for a, b in off_gti_utc:
                gf.write(f"{a:.9f} {b:.9f}\n")

        off_pha = read_spec_fits(off_pha_p)
        off_lc = read_lc_fits(off_lc_p)

        ch_parity = {}
        for nch in (1024, 512, 256):
            # Legacy with nul=1 stops cleanly after writing .spec/.pha/.lcurv/.gti (no bkg pha length clash)
            leg_ch_dir = ensure_clean_dir(f"{l1_dir}/leg_{pref}_offgti_{nch}")
            write_laxpcl1_inputs(leg_ch_dir, orbits, off_gti_ins)
            run_laxpcl1(leg_ch_dir, uid, 1.0, 0, 0, 1023, nch, 1, 1, 2)
            ls_ch = read_spec_fits(f"{leg_ch_dir}/{pref}level2.spec")
            l_asc = parse_ascii_pha(f"{leg_ch_dir}/{pref}level2.pha")

            rust_ch_dir = ensure_clean_dir(f"{l1_dir}/rust_{pref}_offgti_{nch}")
            run_cmd([
                ASTROSAT_BIN, "laxpc", "run",
                "--input", CRAB_L1_DIR,
                "--output", rust_ch_dir,
                "--units", str(uid * 10),
                "--channels", str(nch),
                "--tbin", "1.0",
                "--user-gti", off_gti_txt, "--user-gti-only",
                "--no-format-b", "--no-orbit", "--json",
            ])
            rs_ch = read_spec_fits(f"{rust_ch_dir}/{pref}level2.spec")
            entry = {
                "rust_vs_leg": compare_arrays(rs_ch["rate"], ls_ch["rate"]),
                "leg_exp": l_asc["exposure"],
                "rust_exp": rs_ch["exposure"],
                "leg_corr_rate": l_asc["tot_corr_rate"],
            }
            if nch == 1024:
                # Official L2 .pha stores raw COUNTS (sum over all anodes or A1-7, deadtime-uncorrected)
                off_counts = off_pha["rate"]
                entry["rust_vs_off_pearson_r"] = pearson_r(rs_ch["rate"], off_counts)
                entry["leg_vs_off_pearson_r"] = pearson_r(ls_ch["rate"], off_counts)
                entry["off_exposure"] = off_pha["exposure"]
                entry["off_counts_sum"] = float(np.sum(off_counts))
            ch_parity[str(nch)] = entry

        u_res["matched_gti"] = {
            "leg_events": leg_evn["nrows"],
            "rust_events": rm_evn["nrows"],
            "event_max_dt_sec": ev_time_diff,
            "event_ch_mismatches": ev_ch_mismatch,
            "leg_net_exp": leg_out_m.get("tobs_gross", 0) - leg_out_m.get("tgap", 0),
            "rust_net_exp": rm_json["gti"]["net_exposure_sec"],
            "spec_cmp": compare_arrays(rm_spec["rate"], leg_spec["rate"]),
            "spec_err_cmp": compare_arrays(rm_spec["err"], leg_spec["err"]),
            "lc_1s_cmp": lc_cmp,
            "lc_1s_err_cmp": lc_err_cmp,
            "lcbk_1s_cmp": lcbk_cmp,
            "lc_0p01s_cmp": hi_cmp,
            "leg_gti_mean_lc": float(np.mean(leg_lc["rate"])),
            "rust_gti_mean_lc": float(np.mean(rm_lc["rate"])),
            "off_l2_mean_lc": float(np.mean(off_lc["rate"])),
        }
        u_res["channels_off_gti"] = ch_parity

        # G. Orbit comparison (use leg_p2 AstroSat.orb vs rust_auto AstroSat.orb)
        if uid == 1:
            leg_orb = read_fits(f"{leg_p2}/AstroSat.orb")[1]
            rust_orb = read_fits(f"{rust_auto_dir}/AstroSat.orb")[1]
            lt_o = np.asarray(leg_orb.read_column("TIME"), dtype=np.float64)
            rt_o = np.asarray(rust_orb.read_column("TIME"), dtype=np.float64)
            common_t, il, ir = np.intersect1d(np.round(lt_o).astype(np.int64), np.round(rt_o).astype(np.int64), return_indices=True)
            orb_cmp = {"common_rows": int(len(common_t)), "leg_rows": int(len(lt_o)), "rust_rows": int(len(rt_o))}
            for col in ("X", "Y", "Z", "VX", "VY", "VZ"):
                la = np.asarray(leg_orb.read_column(col), dtype=np.float64)[il]
                ra = np.asarray(rust_orb.read_column(col), dtype=np.float64)[ir]
                orb_cmp[col] = compare_arrays(ra, la)
            results["orbit"] = orb_cmp

        results["units"][pref] = u_res

    # Check Format A level2.event.fits
    fa_evn = read_evn_fits(f"{rust_auto_dir}/level2.event.fits")
    results["format_a"] = {
        "total_events": fa_evn["nrows"],
        "ch_min": int(np.min(fa_evn["ch"])),
        "ch_max": int(np.max(fa_evn["ch"])),
        "units_present": sorted(set(fa_evn["det"].tolist())),
        "sum_format_b_events": sum(results["units"][u["prefix"]]["events_auto"]["rust_count"] for u in UNITS),
    }

    return results


# ---------------------------------------------------------------------------
# Scenario L2: Detector Filtering Modes (Single-only, Layer-1, Anodes, Energy)
# ---------------------------------------------------------------------------
def run_scenario_l2() -> dict[str, Any]:
    print("=== Running Scenario L2: Detector Filtering Modes ===")
    l2_dir = ensure_clean_dir(f"{WORK_ROOT}/L2")
    orbits = crab_orbits()

    # First obtain LXP10 & LXP20 legacy Pass-1 auto-GTI so all filter comparisons use identical GTIs
    res: dict[str, Any] = {}
    for u in UNITS:
        pref = u["prefix"]
        uid = u["id"]
        def_ch = u["def_ch"]
        p1_dir = ensure_clean_dir(f"{l2_dir}/p1_{pref}")
        write_laxpcl1_inputs(p1_dir, orbits, None)
        run_laxpcl1(p1_dir, uid, 1.0, 0, 0, 1023, def_ch, -1, 1, 2)
        gti_rows = parse_laxpcl1_gti(f"{p1_dir}/{pref}level2.gti")
        ins_gtis = [(r[0], r[1]) for r in gti_rows]
        utc_gtis = [(r[2], r[3]) for r in gti_rows]
        gti_txt = f"{p1_dir}/gti_utc.txt"
        with open(gti_txt, "w") as gf:
            for a, b in utc_gtis:
                gf.write(f"{a:.9f} {b:.9f}\n")

        # Load channel map for energy-range -> PHA mapping check
        cmap_path = f"{REF_SRC}/{u['det']}.channel"
        e_lo, e_hi = [], []
        with open(cmap_path, "r") as cf:
            for ln in cf:
                ps = ln.split()
                if len(ps) >= 2:
                    e_lo.append(float(ps[0]))
                    e_hi.append(float(ps[1]))
        scale = 1024 // len(e_lo)

        def map_erange(emin: float, emax: float) -> tuple[int, int]:
            lo_i = min(i for i in range(len(e_lo)) if e_hi[i] >= emin)
            hi_i = max(i for i in range(len(e_lo)) if e_lo[i] <= emax)
            return lo_i * scale, (hi_i + 1) * scale - 1

        pha_3_20 = map_erange(3.0, 20.0)
        pha_20_80 = map_erange(20.0, 80.0)

        modes = [
            ("single_only", 0, 0, 1023, -1, ["--single-only"]),
            ("layer1_only", 1, 0, 1023, -2, ["--layer1-only"]),
            ("anode2_only", 2, 0, 1023, -2, ["--anodes", "2"]),
            ("layer2_a34", 3, 0, 1023, -2, ["--anodes", "3,4"]),
            # laxpcl1.f:363-364 shifts icl/icu by +1 internally, so pass pha-1 to match 0-based [pha_min, pha_max]
            ("band_3_20kev", 0, max(0, pha_3_20[0] - 1), max(0, pha_3_20[1] - 1), -2, ["--energy-range", "3:20"]),
            ("band_20_80kev", 0, max(0, pha_20_80[0] - 1), max(0, pha_20_80[1] - 1), -2, ["--energy-range", "20:80"]),
        ]

        u_modes = {}
        for mname, ian, icl, icu, iev, rust_flags in modes:
            ld = ensure_clean_dir(f"{l2_dir}/leg_{pref}_{mname}")
            write_laxpcl1_inputs(ld, orbits, ins_gtis)
            p_l = run_laxpcl1(ld, uid, 1.0, ian, icl, icu, def_ch, -1, 1, iev)
            lm = parse_laxpcl1_stdout_metrics(p_l.stdout)
            ls = read_spec_fits(f"{ld}/{pref}level2.spec")
            ll = read_lc_fits(f"{ld}/{pref}level2.lcurv")
            le = read_evn_fits(f"{ld}/{pref}level2.evn")

            rd = ensure_clean_dir(f"{l2_dir}/rust_{pref}_{mname}")
            cmd = [
                ASTROSAT_BIN, "laxpc", "run",
                "--input", CRAB_L1_DIR,
                "--output", rd,
                "--units", str(uid * 10),
                "--tbin", "1.0",
                "--user-gti", gti_txt, "--user-gti-only",
                "--caldb", REF_SRC,
                "--format-b", "--no-orbit", "--json",
            ] + rust_flags
            p_r = run_cmd(cmd)
            if p_r.returncode != 0:
                raise RuntimeError(f"Rust L2 {pref} {mname} failed: {p_r.stderr}")
            rj = json.loads(p_r.stdout)["unit_results"][0]
            rs = read_spec_fits(f"{rd}/{pref}level2.spec")
            rl = read_lc_fits(f"{rd}/{pref}level2.lcurv")
            re_evn = read_evn_fits(f"{rd}/{pref}level2.evn")

            # Align LC bins
            t_diff = np.abs(rl["time"][:, None] - ll["time"][None, :])
            idx_l = np.argmin(t_diff, axis=1)
            mask = t_diff[np.arange(len(rl["time"])), idx_l] < 0.05

            u_modes[mname] = {
                "icl_icu": [icl, icu],
                "leg_evn_rows": le["nrows"],
                "rust_evn_rows": re_evn["nrows"],
                "leg_lc_icount": lm.get("lc_icount"),
                "rust_lc_mean_cps": float(np.mean(rl["rate"])),
                "leg_lc_mean_cps": float(np.mean(ll["rate"])),
                "rust_spec_rate_sum": float(np.sum(rs["rate"])),
                "leg_spec_rate_sum": float(np.sum(ls["rate"])),
                "rust_fcd": rj["spectrum"]["deadtime_factor"],
                "spec_cmp": compare_arrays(rs["rate"], ls["rate"]),
                "lc_cmp": compare_arrays(rl["rate"][mask], ll["rate"][idx_l[mask]]),
            }

        # Also verify standalone `astrosat laxpc lightcurve` & `spectrum` match `astrosat laxpc run`
        sa_dir = ensure_clean_dir(f"{l2_dir}/rust_standalone_{pref}")
        ea_f = f"{CRAB_L1_ORB}/{pref}/modeEA/AS1G02_010T01_9000000252{pref}EA_level1.fits"
        bb_f = f"{CRAB_L1_ORB}/{pref}/modeBB/AS1G02_010T01_9000000252{pref}BB_level1.fits"
        mkf_f = f"{CRAB_L1_ORB}/AS1G02_010T01_9000000252lxp_level1.mkf"
        tct_f = f"{CRAB_L1_ORB}/aux/AS1G02_010T01_9000000252lxp_level1.tct"
        run_cmd([
            ASTROSAT_BIN, "laxpc", "lightcurve",
            "--ea", ea_f, "--mkf", mkf_f, "--tct", tct_f, "--bb", bb_f,
            "--gti", gti_txt, "--unit", str(uid * 10), "--tbin", "1.0",
            "--caldb", REF_SRC, "--energy-range", "3:20", "--single-only",
            "-o", f"{sa_dir}/{pref}_3_20_se.lcurv", "--json",
        ])
        p_sp = run_cmd([
            ASTROSAT_BIN, "laxpc", "spectrum",
            "--ea", ea_f, "--mkf", mkf_f, "--tct", tct_f, "--bb", bb_f,
            "--gti", gti_txt, "--unit", str(uid * 10), "--channels", str(def_ch),
            "--energy-range", "3:20", "--single-only", "--caldb", REF_SRC,
            "-o", f"{sa_dir}/{pref}_3_20_se.spec", "--json",
        ])
        sa_lc = read_lc_fits(f"{sa_dir}/{pref}_3_20_se.lcurv")
        sa_sp = read_spec_fits(f"{sa_dir}/{pref}_3_20_se.spec") if p_sp.returncode == 0 else None
        # Compare against `run --energy-range 3:20 --single-only`
        rd_comb = ensure_clean_dir(f"{l2_dir}/rust_run_{pref}_3_20_se")
        run_cmd([
            ASTROSAT_BIN, "laxpc", "run",
            "--input", CRAB_L1_DIR, "--output", rd_comb, "--units", str(uid * 10),
            "--tbin", "1.0", "--user-gti", gti_txt, "--user-gti-only",
            "--caldb", REF_SRC, "--energy-range", "3:20", "--single-only",
            "--no-format-b", "--no-orbit", "--json",
        ])
        run_lc = read_lc_fits(f"{rd_comb}/{pref}level2.lcurv")
        run_sp = read_spec_fits(f"{rd_comb}/{pref}level2.spec")
        u_modes["standalone_lightcurve_vs_run"] = compare_arrays(sa_lc["rate"], run_lc["rate"])
        if sa_sp is not None:
            u_modes["standalone_spectrum_vs_run"] = compare_arrays(sa_sp["rate"], run_sp["rate"])
        res[pref] = u_modes
    return res


# ---------------------------------------------------------------------------
# Scenario L3: Background Modeling & Gain Shift Across Epochs (backshiftv3.e)
# ---------------------------------------------------------------------------
def build_backshiftv3_inp(dest_dir: str) -> str:
    """Create a backshiftv3.inp pointing to absolute files in REF_SRC."""
    inp_p = os.path.join(dest_dir, "backshiftv3.inp")
    with open(f"{REF_SRC}/backshiftv3.inp", "r") as src, open(inp_p, "w") as dst:
        lines = src.readlines()
        for idx, ln in enumerate(lines):
            if idx < 5:
                fn = ln.strip()
                dst.write(f"{REF_SRC}/{fn}\n")
            elif idx < 9:
                dst.write(ln)
            else:
                # Format: a5,7x,3i3,i5,3f7.1,2x,a200
                prefix_part = ln[:49]
                fname = ln[49:].strip()
                if fname:
                    dst.write(f"{prefix_part}  {REF_SRC}/{fname}\n")
    return inp_p


def run_backshift_case(
    case_dir: str,
    lxp_id: int,
    ian: int,
    ndays: int,
    bkg_pha_lxp10_name: str,
    bkg_day: int,
    pha_min: int = 0,
    pha_max: int = 1023,
    nul: int = -1,
    make_lcbk: bool = True,
) -> dict[str, Any]:
    ensure_clean_dir(case_dir)
    build_backshiftv3_inp(case_dir)
    pref = f"lxp{lxp_id}"
    det = f"lxp{lxp_id}0"
    bkg_pha_unit_path = f"{REF_SRC}/{bkg_pha_lxp10_name.replace('lxp10', det)}"
    bkg_spec_unit_path = bkg_pha_unit_path.replace(".pha", ".spec")

    # Parse bkg_pha to build synthetic input `lxp?level2back.pha` with exact (ndayb, ndays) header
    parsed_bk = parse_ascii_pha(bkg_pha_unit_path)
    obs_pha = os.path.join(case_dir, f"{pref}level2back.pha")
    obs_spec = os.path.join(case_dir, f"{pref}level2back.spec")
    nch = parsed_bk["rates"].shape[0]
    with open(obs_pha, "w") as f:
        f.write(f"#  Exposure ={parsed_bk['exposure']:15.6E}  lxp{lxp_id}  Channels:{0:5d}{nch-1:5d}  day:{bkg_day:5d}{ndays:5d}  rate:{np.sum(parsed_bk['rates'][:, 10]):15.6E}\n")
        for i in range(nch):
            vals = " ".join(f"{parsed_bk['rates'][i, j]:12.4E}" for j in range(11))
            f.write(f"{i:5d} {vals}\n")
    shutil.copyfile(bkg_spec_unit_path, obs_spec)

    # Optionally create dummy .lc, .lcbk, .lcurv so backshiftv3.e prints fshbk / shifts .lcbk!
    if make_lcbk:
        # Take Crab L1 .lcbk/.lcurv/.lc if present or build minimal FITS .lcbk with constant 100.0 cps
        l1_sample = f"{WORK_ROOT}/L1/leg_{pref}_p2_1s"
        if os.path.isfile(f"{l1_sample}/{pref}level2.lcbk"):
            shutil.copyfile(f"{l1_sample}/{pref}level2.lc", os.path.join(case_dir, f"{pref}level2.lc"))
            shutil.copyfile(f"{l1_sample}/{pref}level2.lcbk", os.path.join(case_dir, f"{pref}level2.lcbk"))
            shutil.copyfile(f"{l1_sample}/{pref}level2.lcurv", os.path.join(case_dir, f"{pref}level2.lcurv"))

    stdin_text = f"{lxp_id} {ian} 0\n1.0 {pha_min} {pha_max} {nul}\n{pref}level2back.pha\n"
    p_leg = run_cmd([BACKSHIFT_BIN], cwd=case_dir, stdin_text=stdin_text, timeout=60)

    # Parse backshiftv3 stdout
    leg_info: dict[str, Any] = {"returncode": p_leg.returncode, "stdout": p_leg.stdout}
    for ln in p_leg.stdout.splitlines():
        if "Source-ch=" in ln:
            nums = re.findall(r"[-+]?\d*\.\d+|[-+]?\d+", ln)
            if len(nums) >= 3:
                leg_info["src_ch"] = float(nums[0])
                leg_info["bkg_ch"] = float(nums[1])
                leg_info["dch"] = float(nums[2])
        if "Mean count rate:" in ln:
            nums = re.findall(r"[-+]?\d*\.\d+|[-+]?\d+", ln)
            if len(nums) >= 3:
                leg_info["secular_scale"] = float(nums[2])
        if "use the response file" in ln:
            leg_info["rmf"] = ln.split("use the response file")[-1].strip()
        if "change in count rate=" in ln:
            nums = re.findall(r"[-+]?\d*\.\d+(?:[eEdD][-+]?\d+)?", ln)
            if nums:
                leg_info["cdel"] = float(nums[-1].replace("D", "E"))

    leg_shifted_spec = None
    if os.path.isfile(os.path.join(case_dir, f"{pref}level2back_shifted.spec")):
        try:
            leg_shifted_spec = read_spec_fits(os.path.join(case_dir, f"{pref}level2back_shifted.spec"))
        except Exception:
            leg_shifted_spec = None

    leg_lcbk_ratio = None
    if os.path.isfile(os.path.join(case_dir, f"{pref}level2_shifted.lcbk")) and os.path.isfile(os.path.join(case_dir, f"{pref}level2.lcbk")):
        b0 = read_lc_fits(os.path.join(case_dir, f"{pref}level2.lcbk"))["rate"]
        b1 = read_lc_fits(os.path.join(case_dir, f"{pref}level2_shifted.lcbk"))["rate"]
        valid = b0 > 0
        if np.any(valid):
            leg_lcbk_ratio = float(np.median(b1[valid] / b0[valid])) - 1.0

    # Run Rust `astrosat laxpc backshift`
    rust_spec_out = os.path.join(case_dir, f"rust_{pref}_shifted.spec")
    cmd_r = [
        ASTROSAT_BIN, "laxpc", "backshift",
        "--unit", str(lxp_id * 10),
        "--caldb", REF_SRC,
        "--mission-day", str(ndays),
        "--bkg-pha", bkg_pha_unit_path,
        "--pha-min", str(pha_min),
        "--pha-max", str(pha_max),
        "-o", rust_spec_out,
        "--json",
    ]
    if ian == 1:
        cmd_r.append("--layer1-only")
    p_r = run_cmd(cmd_r, cwd=case_dir, timeout=60)
    rust_info: dict[str, Any] = {"returncode": p_r.returncode, "stderr": p_r.stderr}
    if p_r.returncode == 0:
        rj = json.loads(p_r.stdout)
        rust_info["report"] = rj
        rust_shifted_spec = read_spec_fits(rust_spec_out)
        if leg_shifted_spec is not None:
            rust_info["spec_cmp"] = compare_arrays(rust_shifted_spec["rate"], leg_shifted_spec["rate"])
            rust_info["spec_err_cmp"] = compare_arrays(rust_shifted_spec["err"], leg_shifted_spec["err"])
    return {
        "unit": pref,
        "ndays": ndays,
        "bkg_pha": bkg_pha_lxp10_name,
        "bkg_day": bkg_day,
        "ian": ian,
        "pha_range": [pha_min, pha_max],
        "leg": {k: v for k, v in leg_info.items() if k != "stdout"},
        "leg_lcbk_fshbk": leg_lcbk_ratio,
        "rust": rust_info,
    }


def run_scenario_l3() -> dict[str, Any]:
    print("=== Running Scenario L3: Multi-Epoch Backshift & Gain Shift ===")
    l3_dir = ensure_clean_dir(f"{WORK_ROOT}/L3")
    cases = [
        # Epoch 1: Crab Jan 2016 (day 101) vs jan16 (day 100) and mar16 (day 178)
        ("ep1_lxp1_jan16", 1, 0, 101, "backlxp10jan16.pha", 100, 0, 1023),
        ("ep1_lxp2_jan16", 2, 0, 101, "backlxp10jan16.pha", 100, 0, 1023),
        ("ep1_lxp3_jan16", 3, 0, 101, "backlxp10jan16.pha", 100, 0, 1023),
        ("ep1_lxp1_mar16", 1, 0, 101, "backlxp10mar16.pha", 178, 0, 1023),
        ("ep1_lxp2_mar16_band", 2, 1, 101, "backlxp10mar16.pha", 178, 20, 300),
        # Epoch 2: Mid-mission day 700 (post-day 550 LAXPC20 linear drift regime)
        ("ep2_lxp1_d700", 1, 0, 700, "backlxp10aug17.pha", 694, 0, 1023),
        ("ep2_lxp2_d700", 2, 0, 700, "backlxp10aug17.pha", 694, 0, 1023),
        ("ep2_lxp3_d700", 3, 0, 700, "backlxp10aug17.pha", 694, 0, 1023),
        # Epoch 3: Late-mission day 3100 (post-day 2977 LAXPC20 HV jump regime) + cross-regime
        ("ep3_lxp2_d3100_paired", 2, 0, 3100, "backlxp10jan24.pha", 3025, 0, 1023),
        ("ep3_lxp2_d3100_unpaired", 2, 0, 3100, "backlxp10sep23.pha", 2909, 0, 1023),
        # LAXPC10 post-917 cshm regimes & LAXPC30 density/off regimes
        ("lxp1_d950_cshm", 1, 0, 950, "backlxp10apr18.pha", 917, 0, 1023),
        ("lxp1_d1010_cshm", 1, 0, 1010, "backlxp10jun18.pha", 982, 0, 1023),
        ("lxp1_d1060_cshm", 1, 0, 1060, "backlxp10aug18.pha", 1039, 0, 1023),
        ("lxp1_d1150_cshm", 1, 0, 1150, "backlxp10nov18.pha", 1151, 0, 1023),
        ("lxp3_d888_exp", 3, 0, 888, "backlxp10feb18.pha", 870, 0, 1023),
        ("lxp3_d895_off", 3, 0, 895, "backlxp10feb18.pha", 870, 0, 1023),
    ]
    out = {}
    for name, uid, ian, nd, bk_f, bk_d, plo, phi in cases:
        out[name] = run_backshift_case(
            f"{l3_dir}/{name}", uid, ian, nd, bk_f, bk_d, plo, phi,
            make_lcbk=True,
        )
    return out


# ---------------------------------------------------------------------------
# Scenario L4: Synthetic Ground-Truth & Adversarial Edge Scenarios
# ---------------------------------------------------------------------------
def _fits_card(key: str, val_str: str, comment: str = "") -> bytes:
    if comment:
        s = f"{key:<8}= {val_str:<20} / {comment}"
    else:
        s = f"{key:<8}= {val_str}"
    return s[:80].ljust(80).encode("ascii")


def _pad2880(data: bytes) -> bytes:
    rem = len(data) % 2880
    if rem == 0:
        return data
    return data + b"\x00" * (2880 - rem)


def _pad_hdr(cards: list[bytes]) -> bytes:
    cards = list(cards) + [b"END".ljust(80)]
    raw = b"".join(cards)
    rem = len(raw) % 2880
    if rem != 0:
        raw += b" " * (2880 - rem)
    return raw


def write_simple_bintable_fits(
    path: str,
    primary_cards: list[bytes],
    ext_cards: list[bytes],
    rows_bytes: bytes,
) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "wb") as f:
        f.write(_pad_hdr(primary_cards))
        f.write(_pad_hdr(ext_cards))
        f.write(_pad2880(rows_bytes))


def make_synthetic_tct(path: str, t_ins_vals: list[float], t_utc_vals: list[float], tstart: float, tstop: float) -> None:
    p_cards = [
        _fits_card("SIMPLE", "T"),
        _fits_card("BITPIX", "8"),
        _fits_card("NAXIS", "0"),
        _fits_card("EXTEND", "T"),
        _fits_card("TSTART", f"{tstart:.6f}"),
        _fits_card("TSTOP", f"{tstop:.6f}"),
        _fits_card("TSTARTI", str(int(tstart))),
        _fits_card("TSTARTF", f"{tstart - int(tstart):.6f}"),
        _fits_card("TSTOPI", str(int(tstop))),
        _fits_card("TSTOPF", f"{tstop - int(tstop):.6f}"),
    ]
    n = len(t_ins_vals)
    e_cards = [
        _fits_card("XTENSION", "'BINTABLE'"),
        _fits_card("BITPIX", "8"),
        _fits_card("NAXIS", "2"),
        _fits_card("NAXIS1", "24"),
        _fits_card("NAXIS2", str(n)),
        _fits_card("PCOUNT", "0"),
        _fits_card("GCOUNT", "1"),
        _fits_card("TFIELDS", "3"),
        _fits_card("TTYPE1", "'SPS_TIME'"),
        _fits_card("TFORM1", "'1D'"),
        _fits_card("TTYPE2", "'SPS_OBT'"),
        _fits_card("TFORM2", "'1D'"),
        _fits_card("TTYPE3", "'INSTRUMENT_TIME'"),
        _fits_card("TFORM3", "'1D'"),
        _fits_card("EXTNAME", "'TCT'"),
    ]
    buf = bytearray()
    for ti, tu in zip(t_ins_vals, t_utc_vals):
        buf += struct.pack(">ddd", tu, ti, ti)
    write_simple_bintable_fits(path, p_cards, e_cards, bytes(buf))


def make_synthetic_mkf(path: str, times_utc: list[float], elv: float = 35.0, lon: float = 120.0, lat: float = 5.0, rate_alt: float = 550.0) -> None:
    t0, t1 = times_utc[0], times_utc[-1]
    p_cards = [
        _fits_card("SIMPLE", "T"),
        _fits_card("BITPIX", "8"),
        _fits_card("NAXIS", "0"),
        _fits_card("EXTEND", "T"),
        _fits_card("TELESCOP", "'ASTROSAT'"),
        _fits_card("INSTRUME", "'LAXPC'"),
        _fits_card("OBJECT", "'SYNTH_TEST'"),
        _fits_card("OBSERVER", "'ORACLE'"),
        _fits_card("OBS_ID", "'9000000999'"),
        _fits_card("RA_NOM", "83.6332"),
        _fits_card("DEC_NOM", "22.0145"),
        _fits_card("RA_REF", "83.6332"),
        _fits_card("DEC_REF", "22.0145"),
        _fits_card("DATE-OBS", "'2016-01-07T00:00:00'"),
        _fits_card("TIME-OBS", "'00:00:00'"),
        _fits_card("DATE-END", "'2016-01-07T01:00:00'"),
        _fits_card("TIME-END", "'01:00:00'"),
        _fits_card("DATE", "'2026-11-26T00:00:00'"),
        _fits_card("TSTART", f"{t0:.6f}"),
        _fits_card("TSTOP", f"{t1:.6f}"),
        _fits_card("TSTARTI", str(int(t0))),
        _fits_card("TSTARTF", f"{t0 - int(t0):.6f}"),
        _fits_card("TSTOPI", str(int(t1))),
        _fits_card("TSTOPF", f"{t1 - int(t1):.6f}"),
        _fits_card("TIMEDEL", "1.0"),
        _fits_card("TIMEUNIT", "'s'"),
        _fits_card("ORIGIN", "'ASTROSAT-ORACLE'"),
        _fits_card("CREATOR", "'SYNTH_MKF'"),
        _fits_card("RADESSYS", "'ICRS'"),
        _fits_card("CHECKSUM", "'0000000000000000'"),
    ]
    # Use the exact 233-byte layout from converk in laxpcl1.f:2499-2568
    col_defs = [
        ("TIME", "1D", 0, "d"),
        ("ROLL_RA", "1E", 40, "f"),
        ("ROLL_DEC", "1E", 44, "f"),
        ("PITCH_RA", "1E", 56, "f"),
        ("PITCH_DEC", "1E", 60, "f"),
        ("YAW_RA", "1E", 64, "f"),
        ("POSX", "1E", 68, "f"),
        ("POSY", "1E", 72, "f"),
        ("POSZ", "1E", 76, "f"),
        ("VELX", "1E", 80, "f"),
        ("VELY", "1E", 84, "f"),
        ("VELZ", "1E", 88, "f"),
        ("EARTHLAT", "1E", 92, "f"),
        ("EARTHLON", "1E", 96, "f"),
        ("ALTITUDE", "1E", 100, "f"),
        ("ELV", "1E", 148, "f"),
        ("ANG_OFFSET", "1D", 152, "d"),
        ("CPM_RATE", "1J", 196, "i"),
    ]
    # Build contiguous 233-byte row by declaring dummy filler columns so NAXIS1 == 233 and FITS column offsets match converk
    # Easier: copy column header cards directly from real Crab MKF and zero-fill 233 bytes then pack fields!
    real_mkf_hdu = read_fits(f"{CRAB_L1_ORB}/AS1G02_010T01_9000000252lxp_level1.mkf")[1]
    e_cards = []
    for c in real_mkf_hdu.cards:
        k = c[:8].strip()
        if k == "END":
            continue
        if k == "NAXIS2":
            e_cards.append(_fits_card("NAXIS2", str(len(times_utc))))
        elif k in ("TSTART", "TSTOP"):
            e_cards.append(_fits_card(k, f"{t0 if k=='TSTART' else t1:.6f}"))
        else:
            e_cards.append(c.encode("ascii"))
    buf = bytearray()
    for t in times_utc:
        row = bytearray(233)
        struct.pack_into(">d", row, 0, t)
        struct.pack_into(">ff", row, 40, 83.6332, 22.0145)
        struct.pack_into(">ffffff", row, 68, 7000.0, 100.0, -50.0, 1.0, 7.5, 0.2)
        struct.pack_into(">fff", row, 92, lat, lon, rate_alt)
        struct.pack_into(">f", row, 148, elv)
        struct.pack_into(">d", row, 152, 0.01)
        struct.pack_into(">i", row, 196, 10)
        buf += row
    write_simple_bintable_fits(path, p_cards, e_cards, bytes(buf))


def build_ea_frame(
    l1_time: float,
    unit_id: int,
    frame_counter: int,
    events_spec: list[tuple[int, int, int]],  # (ticks_10us_rel_to_coarse, pha, anode)
    stbg_ticks_start: int,
) -> bytes:
    """Construct a valid 2242-byte EA frame matching laxpcl1.f:conver1 and ea.rs."""
    frame = bytearray(2242)
    struct.pack_into(">d", frame, 0, l1_time)
    struct.pack_into(">q", frame, 36, frame_counter)
    frame[44] = 0xDE
    frame[45] = unit_id
    frame[46] = 0xEA
    # 56-bit STBG at bytes 47..54 (40-bit at 49..54)
    stbg_bytes = stbg_ticks_start.to_bytes(7, "big", signed=False)
    frame[47:54] = stbg_bytes
    frame[57] = 0
    struct.pack_into(">H", frame, 58, frame_counter & 0xFFFF)
    # 5-byte records starting at byte 60 (Fortran 1-based 61) up to byte 2090
    pos = 60
    last_coarse_b234 = (-1, -1, -1)
    for abs_tick, pha, anode in events_spec:
        b1 = (abs_tick >> 8) & 0xFF
        b2 = (abs_tick >> 16) & 0xFF
        b3 = (abs_tick >> 24) & 0xFF
        if (b1, b2, b3) != last_coarse_b234:
            if pos + 10 > 2090:
                break
            # Emit 0xEF coarse marker in big-endian: [0xEF, b3, b2, b1, 0] (`laxpcl1.f:2231`)
            frame[pos : pos + 5] = bytes([0xEF, b3, b2, b1, 0])
            pos += 5
            last_coarse_b234 = (b1, b2, b3)
        fine = ((abs_tick & 0xFF) + 2) & 0xFF
        # Pack single event: ic2=0, ic1=pha, ia1=anode
        b0 = anode & 0x0F
        b2_ev = (pha >> 2) & 0xFF
        l1 = pha & 0x03
        b4_ev = (l1 << 2) & 0x0C
        frame[pos : pos + 5] = bytes([b0, fine, b2_ev, 0, b4_ev])
        pos += 5
    if pos + 5 <= 2090 and len(events_spec) > 0:
        # Append 1 dummy event at end before 0xEE because laxpcl1.f:2330-2338 & ea.rs:382
        # drop the final record before 0xEE (`nad = cut - 1`)
        last_tick = events_spec[-1][0] + 10
        fine = ((last_tick & 0xFF) + 2) & 0xFF
        frame[pos : pos + 5] = bytes([1, fine, 25, 0, 0])
        pos += 5
    if pos < 2090:
        frame[pos] = 0xEE
    # ie1, ie2, ie3 = 0 at 2092..2098 (already 0)
    return bytes(frame)


def build_bb_frame(unit_id: int, frame_counter: int, stbg_ticks_start: int, uld_per_sub: int = 10, jan_per_sub: int = 10) -> bytes:
    frame = bytearray(2242)
    frame[44] = 0xDE
    frame[45] = unit_id
    frame[46] = 0xBB
    frame[57] = 0
    struct.pack_into(">H", frame, 58, frame_counter & 0xFFFF)
    for j in range(1, 32):
        base = j * 64
        t_tick = stbg_ticks_start + (j - 1) * 6400  # 64 ms steps -> 128 ms per 2 subframes
        frame[base + 47 : base + 54] = t_tick.to_bytes(7, "big")
        struct.pack_into(">H", frame, base + 56, 10)
        struct.pack_into(">H", frame, base + 58, jan_per_sub)
        struct.pack_into(">H", frame, base + 60, uld_per_sub)
        struct.pack_into(">H", frame, base + 64, 5)
    return bytes(frame)


def write_l1_frame_fits(path: str, frames: list[bytes], tstart: float, tstop: float, obs_id: str = "9000000999") -> None:
    p_cards = [
        _fits_card("SIMPLE", "T"),
        _fits_card("BITPIX", "8"),
        _fits_card("NAXIS", "0"),
        _fits_card("EXTEND", "T"),
        _fits_card("TSTART", f"{tstart:.6f}"),
        _fits_card("TSTOP", f"{tstop:.6f}"),
        _fits_card("TSTARTI", str(int(tstart))),
        _fits_card("TSTARTF", f"{tstart - int(tstart):.6f}"),
        _fits_card("TSTOPI", str(int(tstop))),
        _fits_card("TSTOPF", f"{tstop - int(tstop):.6f}"),
        _fits_card("OBS_ID", f"'{obs_id}'"),
        _fits_card("CREATOR", "'SYNTH_L1'"),
        _fits_card("ORIGIN", "'ORACLE'"),
    ]
    e_cards = [
        _fits_card("XTENSION", "'BINTABLE'"),
        _fits_card("BITPIX", "8"),
        _fits_card("NAXIS", "2"),
        _fits_card("NAXIS1", "2242"),
        _fits_card("NAXIS2", str(len(frames))),
        _fits_card("PCOUNT", "0"),
        _fits_card("GCOUNT", "1"),
        _fits_card("TFIELDS", "1"),
        _fits_card("TTYPE1", "'FRAME'"),
        _fits_card("TFORM1", "'2242B'"),
        _fits_card("EXTNAME", "'LAXPC_DATA'"),
        _fits_card("TSTARTI", str(int(tstart))),
        _fits_card("TSTARTF", f"{tstart - int(tstart):.6f}"),
        _fits_card("TSTOPI", str(int(tstop))),
        _fits_card("TSTOPF", f"{tstop - int(tstop):.6f}"),
        _fits_card("CREATOR", "'SYNTH_L1'"),
    ]
    write_simple_bintable_fits(path, p_cards, e_cards, b"".join(frames))


def run_scenario_l4() -> dict[str, Any]:
    print("=== Running Scenario L4: Synthetic Ground-Truth & Adversarial Edge Scenarios ===")
    l4_dir = ensure_clean_dir(f"{WORK_ROOT}/L4")
    results: dict[str, Any] = {}

    # Constants from ea.rs / laxpcl1.f:2191
    offset = 1.8121823219164e8
    base_utc = 190_000_000.0
    base_ins = 8_781_767.80836  # so base_ins + offset == base_utc exactly (n4 = 0)
    base_tick = int(round(base_ins * 1e5))

    # ------------------------------------------------------------------
    # L4a: Multi-orbit merging, out-of-order, duplicates, stray, 16-bit rollover, missing-frame gap
    # ------------------------------------------------------------------
    l4a_dir = ensure_clean_dir(f"{l4_dir}/L4a_frames_gaps")
    orb1_dir = f"{l4a_dir}/l1/laxpc/02001"
    orb2_dir = f"{l4a_dir}/l1/laxpc/02002"
    os.makedirs(f"{orb1_dir}/aux", exist_ok=True)
    os.makedirs(f"{orb1_dir}/lxp1/modeEA", exist_ok=True)
    os.makedirs(f"{orb1_dir}/lxp1/modeBB", exist_ok=True)
    os.makedirs(f"{orb2_dir}/aux", exist_ok=True)
    os.makedirs(f"{orb2_dir}/lxp1/modeEA", exist_ok=True)
    os.makedirs(f"{orb2_dir}/lxp1/modeBB", exist_ok=True)

    # Orbit 1: frames 65530..65535 -> rollover -> 0..5 (12 frames, 0.2s apart), plus 1 stray (counter 30000), 1 duplicate of 65532, out-of-order 65534 before 65533
    # Orbit 2: overlaps frame 4..5, then jumps by missing gap of 0.50 s (10 missing frames -> counter 6..15 missing? No: counter 6, then gap 0.50s -> counter 10), ends at 14
    t_utc_grid = [base_utc - 5.0 + i for i in range(40)]
    t_ins_grid = [base_ins - 5.0 + i * (1.0 / 1.0001) for i in range(40)]
    # Exact 1:1.0001 linear TCT
    tct_1 = f"{orb1_dir}/aux/AS1SYNTH_0000099999lxp_level1.tct"
    tct_2 = f"{orb2_dir}/aux/AS1SYNTH_0000099998lxp_level1.tct"
    make_synthetic_tct(tct_1, t_ins_grid[:25], [base_utc - 5.0 + i for i in range(25)], base_utc - 5.0, base_utc + 20.0)
    make_synthetic_tct(tct_2, t_ins_grid[10:], [base_utc - 5.0 + i for i in range(10, 40)], base_utc + 5.0, base_utc + 35.0)
    mkf_1 = f"{orb1_dir}/AS1SYNTH_0000099999lxp_level1.mkf"
    mkf_2 = f"{orb2_dir}/AS1SYNTH_0000099998lxp_level1.mkf"
    make_synthetic_mkf(mkf_1, [base_utc - 5.0 + i for i in range(25)])
    make_synthetic_mkf(mkf_2, [base_utc - 5.0 + i for i in range(10, 40)])

    def make_ev_frame(fc: int, dt_sec: float, n_ev: int = 10) -> bytes:
        tk0 = base_tick + int(round(dt_sec * 1e5))
        evs = [(tk0 + k * 200, 100 + k * 5, (k % 7) + 1) for k in range(n_ev)]
        return build_ea_frame(base_utc + dt_sec, 1, fc, evs, tk0)

    # Build Orbit 1 frames (counters 65530..65535, 0..5 at dt=0.0, 0.2, ..., 2.2 s)
    f65530 = make_ev_frame(65530, 0.0, 10)
    f65531 = make_ev_frame(65531, 0.2, 10)
    f65532 = make_ev_frame(65532, 0.4, 10)
    f65533 = make_ev_frame(65533, 0.6, 10)
    f65534 = make_ev_frame(65534, 0.8, 10)
    f65535 = make_ev_frame(65535, 1.0, 10)
    f0 = make_ev_frame(0, 1.2, 10)
    f1 = make_ev_frame(1, 1.4, 10)
    f2 = make_ev_frame(2, 1.6, 10)
    f3 = make_ev_frame(3, 1.8, 10)
    f4 = make_ev_frame(4, 2.0, 10)
    f5 = make_ev_frame(5, 2.2, 10)
    f_stray = make_ev_frame(30000, 0.5, 5)
    f65532_dup = make_ev_frame(65532, 0.4, 10)

    # Sequence with out-of-order (65534 before 65533), stray (30000), duplicate (65532), and 16-bit rollover (65535 -> 0)
    orb1_frames = [f65530, f65531, f65532, f_stray, f65534, f65533, f65532_dup, f65535, f0, f1, f2, f3, f4, f5]
    ea_1 = f"{orb1_dir}/lxp1/modeEA/AS1SYNTH_0000099999lxp1EA_level1.fits"
    bb_1 = f"{orb1_dir}/lxp1/modeBB/AS1SYNTH_0000099999lxp1BB_level1.fits"
    write_l1_frame_fits(ea_1, orb1_frames, base_utc, base_utc + 3.0)
    write_l1_frame_fits(bb_1, [build_bb_frame(1, k, base_tick + k * 200000) for k in range(3)], base_utc, base_utc + 3.0)

    # Orbit 2: overlaps f4 (dt=2.0), f5 (dt=2.2), then gap of 0.60 s (counters 6,7,8 missing -> f9 at dt=3.0..3.8)
    orb2_frames = [
        make_ev_frame(4, 2.0, 10),
        make_ev_frame(5, 2.2, 10),
        make_ev_frame(9, 3.0, 10),
        make_ev_frame(10, 3.2, 10),
        make_ev_frame(11, 3.4, 10),
        make_ev_frame(12, 3.6, 10),
    ]
    ea_2 = f"{orb2_dir}/lxp1/modeEA/AS1SYNTH_0000099998lxp1EA_level1.fits"
    bb_2 = f"{orb2_dir}/lxp1/modeBB/AS1SYNTH_0000099998lxp1BB_level1.fits"
    write_l1_frame_fits(ea_2, orb2_frames, base_utc + 2.0, base_utc + 4.0)
    write_l1_frame_fits(bb_2, [build_bb_frame(1, k, base_tick + 200000 + k * 200000) for k in range(3)], base_utc + 2.0, base_utc + 4.0)

    # Ground truth:
    # Orbit 1 retained: 12 unique valid frames (65530..65535, 0..5), 120 events.
    # Orbit 2: f4 & f5 overlap <= t00 (dt <= 2.218 s) skipped, f9..f12 retained (40 events).
    # Total ground-truth unique events = 160!
    # Gap between last event of f5 (dt = 2.2 + 9*0.002 = 2.218 s) and first event of f9 (dt = 3.0 s) is 0.782 s (>= 50 ms).
    gti_ins_a = [(base_ins - 0.1, base_ins + 4.0)]
    # Map exact UTC GTI via TCT (dt_s / dt_i = 1.0001 -> t_utc = base_utc - 5.0 + (t_ins - (base_ins - 5.0)) * 1.0001)
    def ins2utc(ti: float) -> float:
        return (base_utc - 5.0) + (ti - (base_ins - 5.0)) * 1.0001

    gti_utc_a = [(ins2utc(a), ins2utc(b)) for a, b in gti_ins_a]
    gti_a_txt = f"{l4a_dir}/gti_utc.txt"
    with open(gti_a_txt, "w") as gf:
        for a, b in gti_utc_a:
            gf.write(f"{a:.9f} {b:.9f}\n")

    leg_l4a = ensure_clean_dir(f"{l4a_dir}/leg_out")
    write_laxpcl1_inputs(leg_l4a, [(mkf_1, tct_1, bb_1, ea_1), (mkf_2, tct_2, bb_2, ea_2)], gti_ins_a)
    p_leg4a = run_laxpcl1(leg_l4a, 1, 1.0, 0, 0, 1023, 512, 1, 1, -2)
    leg4a_m = parse_laxpcl1_stdout_metrics(p_leg4a.stdout)
    leg4a_evn = read_evn_fits(f"{leg_l4a}/lxp1level2.evn")

    rust_l4a = ensure_clean_dir(f"{l4a_dir}/rust_out")
    p_r4a = run_cmd([
        ASTROSAT_BIN, "laxpc", "run",
        "--input", f"{l4a_dir}/l1",
        "--output", rust_l4a,
        "--units", "10",
        "--tbin", "1.0",
        "--user-gti", gti_a_txt, "--user-gti-only",
        "--format-b", "--no-orbit", "--json",
    ])
    if p_r4a.returncode != 0:
        raise RuntimeError(f"Rust L4a failed: {p_r4a.stderr}")
    r4a_j = json.loads(p_r4a.stdout)["unit_results"][0]
    r4a_evn = read_evn_fits(f"{rust_l4a}/lxp1level2.evn")

    expected_gap = 3.0 - (2.2 + 9 * 0.002)  # 0.782 s
    results["L4a_multi_orbit_rollover_gap"] = {
        "expected_events": 160,
        "rust_events": r4a_evn["nrows"],
        "leg_events": leg4a_evn["nrows"],
        "expected_gap_sec": expected_gap,
        "rust_gap_loss_sec": r4a_j["gti"]["frame_gap_loss_sec"],
        "leg_tgap_sec": leg4a_m.get("tgap"),
        "rust_gross_exp": r4a_j["gti"]["gross_exposure_sec"],
        "rust_net_exp": r4a_j["gti"]["net_exposure_sec"],
        "rust_diag_orb": r4a_j["ea_frame_diagnostics"],
    }

    # ------------------------------------------------------------------
    # L4b: High-count-rate deadtime test (R_all = 1000, 10000, 22500, 23500 cps)
    # ------------------------------------------------------------------
    l4b_res = {}
    for target_rate in (1000, 10000, 22500, 23500):
        case_d = ensure_clean_dir(f"{l4_dir}/L4b_rate_{target_rate}")
        os.makedirs(f"{case_d}/l1/laxpc/03001/aux", exist_ok=True)
        os.makedirs(f"{case_d}/l1/laxpc/03001/lxp1/modeEA", exist_ok=True)
        os.makedirs(f"{case_d}/l1/laxpc/03001/lxp1/modeBB", exist_ok=True)
        dur = 4.0
        n_total = int(target_rate * dur)
        # Pack up to 400 events per frame (spacing 10us..40us)
        # Keep inter-frame gap < 50ms by emitting frames every 10 ms!
        frames = []
        ev_per_ms10 = target_rate // 100  # e.g. 23500/100 = 235 events per 10 ms frame (< 400 max!)
        step_ticks = max(1, int(1000 // max(ev_per_ms10, 1)))
        for fi in range(400):  # 400 * 10 ms = 4.0 s
            dt_s = fi * 0.01
            tk0 = base_tick + int(round(dt_s * 1e5))
            evs = [(tk0 + k * step_ticks, 200, 1) for k in range(ev_per_ms10)]
            frames.append(build_ea_frame(base_utc + dt_s, 1, fi + 1, evs, tk0))
        ea_p = f"{case_d}/l1/laxpc/03001/lxp1/modeEA/AS1SYNTH_0000099901lxp1EA_level1.fits"
        bb_p = f"{case_d}/l1/laxpc/03001/lxp1/modeBB/AS1SYNTH_0000099901lxp1BB_level1.fits"
        tct_p = f"{case_d}/l1/laxpc/03001/aux/AS1SYNTH_0000099901lxp_level1.tct"
        mkf_p = f"{case_d}/l1/laxpc/03001/AS1SYNTH_0000099901lxp_level1.mkf"
        make_synthetic_tct(tct_p, t_ins_grid[:25], [base_utc - 5.0 + i for i in range(25)], base_utc - 5.0, base_utc + 20.0)
        make_synthetic_mkf(mkf_p, [base_utc - 5.0 + i for i in range(25)])
        write_l1_frame_fits(ea_p, frames, base_utc, base_utc + 4.0)
        write_l1_frame_fits(bb_p, [build_bb_frame(1, k, base_tick + k * 200000) for k in range(3)], base_utc, base_utc + 4.0)

        gti_ins_b = [(base_ins + 0.5, base_ins + 3.5)]  # 3.0 s strictly inside [0, 4.0]
        gti_utc_b = [(ins2utc(a), ins2utc(b)) for a, b in gti_ins_b]
        gti_b_txt = f"{case_d}/gti_utc.txt"
        with open(gti_b_txt, "w") as gf:
            for a, b in gti_utc_b:
                gf.write(f"{a:.9f} {b:.9f}\n")

        ld = ensure_clean_dir(f"{case_d}/leg")
        write_laxpcl1_inputs(ld, [(mkf_p, tct_p, bb_p, ea_p)], gti_ins_b)
        p_l = run_laxpcl1(ld, 1, 1.0, 0, 0, 1023, 512, -1, 1, 2)
        lm = parse_laxpcl1_stdout_metrics(p_l.stdout)
        ls = read_spec_fits(f"{ld}/lxp1level2.spec")
        ll = read_lc_fits(f"{ld}/lxp1level2.lcurv")

        rd = ensure_clean_dir(f"{case_d}/rust")
        p_r = run_cmd([
            ASTROSAT_BIN, "laxpc", "run",
            "--input", f"{case_d}/l1", "--output", rd, "--units", "10",
            "--tbin", "1.0", "--user-gti", gti_b_txt, "--user-gti-only",
            "--no-format-b", "--no-orbit", "--json",
        ])
        rj = json.loads(p_r.stdout)["unit_results"][0]
        rs = read_spec_fits(f"{rd}/lxp1level2.spec")
        rl = read_lc_fits(f"{rd}/lxp1level2.lcurv")

        r_obs = ev_per_ms10 * 100.0
        analytic_corr = r_obs / (1.0 - 43.0e-6 * r_obs) if r_obs < 23000.0 else r_obs
        l4b_res[str(target_rate)] = {
            "r_obs_cps": r_obs,
            "analytic_corr_cps": analytic_corr,
            "rust_spec_corr_cps": rj["spectrum"]["total_corrected_rate_cps"],
            "rust_fcd": rj["spectrum"]["deadtime_factor"],
            "leg_spec_corr_cps": float(np.sum(ls["rate"])),
            "rust_lc_mean_cps": float(np.mean(rl["rate"])),
            "leg_lc_mean_cps": float(np.mean(ll["rate"])),
            "rust_all_finite": bool(np.all(np.isfinite(rs["rate"])) and np.all(np.isfinite(rl["rate"]))),
            "leg_all_finite": bool(np.all(np.isfinite(ls["rate"])) and np.all(np.isfinite(ll["rate"]))),
        }
    results["L4b_high_rate_deadtime"] = l4b_res

    # ------------------------------------------------------------------
    # L4c: Format A vs Format B self-consistency & Barycenter readiness
    # ------------------------------------------------------------------
    l4c_dir = ensure_clean_dir(f"{l4_dir}/L4c_fmt_a_vs_b")
    rust_l1_dir = f"{WORK_ROOT}/L1/rust_auto_1s"
    if os.path.isfile(f"{rust_l1_dir}/level2.event.fits"):
        fa = read_evn_fits(f"{rust_l1_dir}/level2.event.fits")
        per_unit_cmp = {}
        for u in UNITS:
            uid = u["id"]
            pref = u["prefix"]
            fb = read_evn_fits(f"{rust_l1_dir}/{pref}level2.evn")
            mask = fa["det"] == uid
            fa_t = fa["time"][mask]
            fa_ch = fa["ch"][mask]
            # Format A default is 1024-ch; Format B default is 512 (Lxp10/30) or 256 (Lxp20)
            scale = 1024 // u["def_ch"]
            # Check exact time equality (sort both since Format A sorts globally by time_utc)
            ord_a = np.argsort(fa_t, kind="mergesort")
            ord_b = np.argsort(fb["time"], kind="mergesort")
            same_len = len(fa_t) == len(fb["time"])
            max_dt = float(np.max(np.abs(fa_t[ord_a] - fb["time"][ord_b]))) if same_len else -1.0
            ch_scaled_match = bool(np.array_equal(fa_ch[ord_a] // scale, fb["ch"][ord_b])) if same_len else False
            per_unit_cmp[pref] = {
                "fmt_a_rows": int(np.sum(mask)),
                "fmt_b_rows": fb["nrows"],
                "exact_count_match": same_len,
                "max_dt_sec": max_dt,
                "channel_scale_factor": scale,
                "scaled_channel_exact_match": ch_scaled_match,
            }

        # Run `astrosat bary` on BOTH Format B lxp1level2.evn and Format A level2.event.fits
        bary_fb_out = f"{l4c_dir}/lxp1_bary.evn"
        bary_fa_out = f"{l4c_dir}/level2_bary.event.fits"
        p_bfb = run_cmd([
            ASTROSAT_BIN, "bary",
            "-i", f"{rust_l1_dir}/lxp1level2.evn",
            "-o", bary_fb_out,
            "--orbit", f"{rust_l1_dir}/AstroSat.orb",
            "--ra", "83.633208", "--dec", "22.014472", "--json",
        ])
        p_bfa = run_cmd([
            ASTROSAT_BIN, "bary",
            "-i", f"{rust_l1_dir}/level2.event.fits",
            "-o", bary_fa_out,
            "--orbit", f"{rust_l1_dir}/AstroSat.orb",
            "--ra", "83.633208", "--dec", "22.014472", "--json",
        ])
        bary_cmp = {"fb_rc": p_bfb.returncode, "fa_rc": p_bfa.returncode}
        if p_bfb.returncode == 0 and p_bfa.returncode == 0:
            bfb_ev = read_evn_fits(bary_fb_out)
            bfa_ev = read_evn_fits(bary_fa_out)
            mask1 = bfa_ev["det"] == 1
            bary_cmp["lxp1_bary_max_dt_sec"] = float(np.max(np.abs(np.sort(bfa_ev["time"][mask1]) - np.sort(bfb_ev["time"]))))
        results["L4c_fmt_a_vs_b"] = {"units": per_unit_cmp, "bary": bary_cmp}

    return results


def main() -> int:
    parser = argparse.ArgumentParser(description="LAXPC Head-to-Head Legacy Oracle & Scenario Verifier")
    parser.add_argument("--scenarios", nargs="+", default=["L1", "L2", "L3", "L4"], choices=["L1", "L2", "L3", "L4"])
    parser.add_argument("--json-out", default=f"{WORK_ROOT}/laxpc_head2head_results.json")
    args = parser.parse_args()

    os.makedirs(WORK_ROOT, exist_ok=True)
    report: dict[str, Any] = {"timestamp": time.strftime("%Y-%m-%d %H:%M:%S")}

    if "L1" in args.scenarios:
        report["L1"] = run_scenario_l1()
    if "L2" in args.scenarios:
        report["L2"] = run_scenario_l2()
    if "L3" in args.scenarios:
        report["L3"] = run_scenario_l3()
    if "L4" in args.scenarios:
        report["L4"] = run_scenario_l4()

    with open(args.json_out, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2)
    print(f"Wrote results to {args.json_out}")
    print(json.dumps(report, indent=2)[:4000])
    return 0


if __name__ == "__main__":
    sys.exit(main())
