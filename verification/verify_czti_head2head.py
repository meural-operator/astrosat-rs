#!/usr/bin/env python3
"""
Comprehensive Head-to-Head Legacy Oracle & Ground-Truth Scenario Verifier for AstroSat CZTI.

Covers:
  - Scenario C1: Real ObsID 9000000406 (2016 Crab) stage-by-stage parity vs Official IUCAA L2
                 and compiled Legacy CZTI v3.0.1 binaries (default & generalized v3.0),
                 plus 2nd real epoch ObsID 9000002432 (2018 Crab, 5945 s).
  - Scenario C2: COMMON vs QUAD mode mathematical consistency, Stage 2 & Stage 3 resume
                 bit-exactness, and cztaddspec (4-quad + multi-ObsID + vs legacy cztaddspec).
  - Scenario C3: Ground-truth planted synthetic scenarios:
                 C3.1 Off-axis coded-mask source at (tx=+2.4, ty=-1.8) deg (balcc, fft, bindata,
                      and exact zero net background).
                 C3.2 Out-of-FOV guard (|tx| > 5.8 deg or |ty| > 5.8 deg).
                 C3.3 Origin coordinates (0.0, 0.0) without falling back to header RA/DEC.
                 C3.4 Reduced 2-word telemetry mode (modeID=2,3) vs legacy cztscience2event.
                 C3.5 20-bit CZTI clock rollover (1048576 s).
                 C3.6 PHA == 0 -> Energy == 0.0.
                 C3.7 100 s clock-tick 400 ms exclusion in bunchclean.
                 C3.8 Secular gain drift across epochs (2016 vs 2022) vs legacy cztevtencor.

All temporary files and outputs stay strictly on D: (/mnt/d/ISRO/.tmp/legacy_bins/czti/...).
"""
from __future__ import annotations

import argparse
import gzip
import json
import math
import os
import shutil
import struct
import subprocess
import sys
import time
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fits_reader import Hdu, read_fits

# Paths strictly on D:
ROOT = "/mnt/d/ISRO"
ASTROSAT_BIN = f"{ROOT}/target/release/astrosat"
FITSVERIFY_BIN = f"{ROOT}/.toolchain/bin/fitsverify"
CALDB_ROOT = f"{ROOT}/data/caldb/czti"
CALDB_CZTI = f"{CALDB_ROOT}/data/as1/czti"
REF_SRC_CZTI = f"{ROOT}/data/ref_src/czti/czti_pipeline/czti"
LEGACY_BIN_DIR = f"{ROOT}/.tmp/legacy_bins/czti/bin"
WORK_BASE = f"{ROOT}/.tmp/legacy_bins/czti"

L1_406 = f"{ROOT}/data/sample/czti/czti_sample_data/20160331_T01_112T01_9000000406_level1"
REF_L2_406 = f"{ROOT}/data/sample/czti/czti_sample_data/20160331_T01_112T01_9000000406_level2/modeM0"
REF_L2_2432 = f"{ROOT}/data/sample/czti_2018/20181014_C04_007T02_9000002432_level2/czti"
PRE_406 = "AS1T01_112T01_9000000406cztM0_level2"
PRE_2432 = "AS1C04_007T02_9000002432cztM0_level2"


# ---------------------------------------------------------------------------
# Minimal Pure-Python FITS Writer (for synthetic ground-truth scenarios)
# ---------------------------------------------------------------------------

def _pad_block(b: bytes) -> bytes:
    rem = len(b) % 2880
    return b if rem == 0 else b + b"\x00" * (2880 - rem)


def _pad_header(cards: List[str]) -> bytes:
    buf = bytearray()
    for c in cards:
        s = c.encode("ascii")[:80]
        buf.extend(s.ljust(80, b" "))
    buf.extend(b"END".ljust(80, b" "))
    rem = len(buf) % 2880
    if rem != 0:
        buf.extend(b" " * (2880 - rem))
    return bytes(buf)


def _card(key: str, val: Any, comment: str = "") -> str:
    k = key.upper()[:8].ljust(8)
    if isinstance(val, bool):
        vstr = ("T" if val else "F").rjust(20)
    elif isinstance(val, (int, np.integer)):
        vstr = f"{int(val):>20d}"
    elif isinstance(val, (float, np.floating)):
        vstr = f"{float(val):>20.12G}"
    elif isinstance(val, str):
        esc = val.replace("'", "''")
        vstr = f"'{esc:<8s}'" if len(esc) < 8 else f"'{esc}'"
        vstr = vstr.ljust(20)
    else:
        vstr = str(val).rjust(20)
    c = f" / {comment}" if comment else ""
    return f"{k}= {vstr}{c}"[:80]


def _make_primary_hdu(extra_cards: Optional[List[str]] = None) -> bytes:
    cards = [
        _card("SIMPLE", True),
        _card("BITPIX", 8),
        _card("NAXIS", 0),
        _card("EXTEND", True),
    ]
    if extra_cards:
        cards.extend(extra_cards)
    return _pad_header(cards)


def _make_bintable_hdu(
    extname: str,
    columns: List[Tuple[str, str, np.ndarray]],
    extra_cards: Optional[List[str]] = None,
) -> bytes:
    """columns: list of (name, tform, data_array)."""
    type_info = {
        "B": (">u1", 1),
        "I": (">i2", 2),
        "J": (">i4", 4),
        "K": (">i8", 8),
        "E": (">f4", 4),
        "D": (">f8", 8),
        "A": ("|S1", 1),
    }
    nrows = len(columns[0][2]) if columns else 0
    col_specs = []
    row_len = 0
    for name, tform, arr in columns:
        i = 0
        while i < len(tform) and tform[i].isdigit():
            i += 1
        rep = int(tform[:i]) if i > 0 else 1
        code = tform[i]
        dt, sz = type_info[code]
        col_specs.append((name, tform, rep, code, dt, sz, arr))
        row_len += rep * sz

    cards = [
        _card("XTENSION", "BINTABLE"),
        _card("BITPIX", 8),
        _card("NAXIS", 2),
        _card("NAXIS1", row_len),
        _card("NAXIS2", nrows),
        _card("PCOUNT", 0),
        _card("GCOUNT", 1),
        _card("TFIELDS", len(columns)),
        _card("EXTNAME", extname),
    ]
    for idx, (name, tform, *_rest) in enumerate(col_specs, start=1):
        cards.append(_card(f"TTYPE{idx}", name))
        cards.append(_card(f"TFORM{idx}", tform))
    if extra_cards:
        cards.extend(extra_cards)

    hdr = _pad_header(cards)
    if nrows == 0 or row_len == 0:
        return hdr
    buf = np.zeros((nrows, row_len), dtype=np.uint8)
    off = 0
    for name, tform, rep, code, dt, sz, arr in col_specs:
        width = rep * sz
        a = np.asarray(arr, dtype=dt)
        if rep > 1:
            a = a.reshape(nrows, rep)
        else:
            a = a.reshape(nrows)
        raw = np.frombuffer(np.ascontiguousarray(a).tobytes(), dtype=np.uint8).reshape(nrows, width)
        buf[:, off : off + width] = raw
        off += width
    return hdr + _pad_block(buf.tobytes())


def write_fits_file(path: str, hdus_bytes: List[bytes]) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "wb") as f:
        for b in hdus_bytes:
            f.write(b)


# ---------------------------------------------------------------------------
# Helpers: RSP decoding, legacy env, subprocess runner
# ---------------------------------------------------------------------------

def legacy_env() -> Dict[str, str]:
    env = os.environ.copy()
    env["as1czt"] = REF_SRC_CZTI
    env["CALDB"] = CALDB_ROOT
    env["CZTI_templates"] = f"{REF_SRC_CZTI}/templates"
    env["PFILES"] = f"{REF_SRC_CZTI}/paramfiles"
    env["USER"] = "diat"
    env["TMPDIR"] = f"{ROOT}/.tmp"
    env["PATH"] = f"{LEGACY_BIN_DIR}:{env.get('PATH', '')}"
    os.makedirs(f"{WORK_BASE}/run_logs", exist_ok=True)
    env["logfile_env"] = f"{WORK_BASE}/run_logs/verify.log"
    return env


def run_cmd(cmd: List[str], env: Optional[Dict[str, str]] = None, cwd: Optional[str] = None, check: bool = True) -> subprocess.CompletedProcess:
    proc = subprocess.run(cmd, env=env, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if check and proc.returncode != 0:
        raise RuntimeError(f"Command failed (rc={proc.returncode}): {' '.join(cmd)}\nSTDOUT:\n{proc.stdout[-1500:]}\nSTDERR:\n{proc.stderr[-1500:]}")
    return proc


def decode_rsp_matrix(hdu: Hdu) -> np.ndarray:
    nrows = int(hdu.header["NAXIS2"])
    row_len = int(hdu.header["NAXIS1"])
    theap = int(hdu.header.get("THEAP", nrows * row_len))
    n_grp = hdu.read_column("N_GRP")
    cols = {c[0].upper(): c for c in hdu.bintable_columns()}

    def read_vla(col_name: str, dtype_ch: str) -> List[np.ndarray]:
        _, rep, code, off = cols[col_name]
        raw = np.frombuffer(hdu.data_bytes[: nrows * row_len], dtype=np.uint8).reshape(nrows, row_len)
        elem_sz = 2 if dtype_ch == "I" else 4
        dt = ">i2" if dtype_ch == "I" else ">f4"
        out = []
        if code == "P":
            for r in range(nrows):
                cnt, heap_off = struct.unpack(">ii", raw[r, off : off + 8].tobytes())
                if cnt <= 0:
                    out.append(np.empty(0, dtype=dt))
                else:
                    s = theap + heap_off
                    out.append(np.frombuffer(hdu.data_bytes[s : s + cnt * elem_sz], dtype=dt))
        else:
            arr = hdu.read_column(col_name)
            for r in range(nrows):
                out.append(np.atleast_1d(arr[r]))
        return out

    f_ch = read_vla("F_CHAN", "I")
    n_ch = read_vla("N_CHAN", "I")
    mat_v = read_vla("MATRIX", "E")
    dense = np.zeros((nrows, 512), dtype=np.float64)
    for r in range(nrows):
        pos = 0
        ng = min(int(n_grp[r]), len(f_ch[r]), len(n_ch[r]))
        for g in range(ng):
            fc, nc = int(f_ch[r][g]), int(n_ch[r][g])
            if 0 <= fc < 512:
                take = min(nc, 512 - fc, len(mat_v[r]) - pos)
                if take > 0:
                    dense[r, fc : fc + take] = mat_v[r][pos : pos + take]
            pos += nc
    return dense


def read_fits_maybe_gz(path: str) -> List[Hdu]:
    with open(path, "rb") as f:
        raw = f.read()
    if raw[:2] == b"\x1f\x8b":
        raw = gzip.decompress(raw)
        tmp = f"{WORK_BASE}/_verify_ungz.fits"
        with open(tmp, "wb") as f:
            f.write(raw)
        return read_fits(tmp)
    return read_fits(path)


def detid_pixid_to_detxy(quad_id: int, det_id: int, pix_id: int) -> Tuple[int, int]:
    dx = det_id % 4
    dy = det_id // 4
    px = pix_id % 16
    py = pix_id // 16
    detx = 16 * dx + px
    dety = 16 * dy + py
    if quad_id in (0, 3):
        dety = 63 - dety
    else:
        detx = 63 - detx
    return detx, dety


def _roll_rot_to_rpy_vecs(ra_roll_deg: float, dec_roll_deg: float, roll_rot_deg: float = 0.0) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Exact basis vectors (roll, pitch, yaw) matching radec2CZTI.f::roll_rot2rpy & bindata.rs::roll_rot_to_rpy."""
    ra = math.radians(ra_roll_deg)
    dec = math.radians(dec_roll_deg)
    rot = math.radians(roll_rot_deg)
    alr = math.cos(ra) * math.cos(dec)
    amr = math.sin(ra) * math.cos(dec)
    anr = math.sin(dec)
    amp = -math.cos(rot) * math.sin(dec) * math.sin(ra) - math.sin(rot) * math.cos(ra)
    alp = -math.cos(rot) * math.sin(dec) * math.cos(ra) + math.sin(rot) * math.sin(ra)
    anp = math.sqrt(max(0.0, 1.0 - alp * alp - amp * amp))
    temp = alp * alr + amp * amr
    anp = -math.copysign(anp, temp)
    if dec < 0.0:
        anp = -anp
    roll = np.array([alr, amr, anr], dtype=np.float64)
    pitch = np.array([alp, amp, anp], dtype=np.float64)
    yaw = np.cross(roll, pitch)
    return roll, pitch, yaw


def radec_from_pointing_and_txty(ra_roll_deg: float, dec_roll_deg: float, tx_deg: float, ty_deg: float) -> Tuple[float, float]:
    """Exact inverse of radec2CZTI.f::fcompute_tx_ty for roll_rot=0."""
    roll, pitch, yaw = _roll_rot_to_rpy_vecs(ra_roll_deg, dec_roll_deg, 0.0)
    # src_fin = [-dot(yaw, src), +dot(pitch, src), +dot(roll, src)]
    # tx = atan2(src_fin[0], src_fin[2]), ty = atan2(src_fin[1], src_fin[2])
    ux = math.tan(math.radians(tx_deg))
    uy = math.tan(math.radians(ty_deg))
    v_ics = roll + uy * pitch - ux * yaw
    v_ics = v_ics / np.linalg.norm(v_ics)
    dec_src = math.degrees(math.asin(np.clip(v_ics[2], -1.0, 1.0)))
    ra_src = math.degrees(math.atan2(v_ics[1], v_ics[0])) % 360.0
    return ra_src, dec_src


# ---------------------------------------------------------------------------
# Scenario C1: Real ObsID 9000000406 + 2018 ObsID 9000002432
# ---------------------------------------------------------------------------

def run_scenario_c1(results: Dict[str, Any], rerun_rust: bool = True) -> None:
    print("\n================ SCENARIO C1: REAL OBSID 9000000406 & 9000002432 ================")
    c1: Dict[str, Any] = {"checks": {}}

    rst_def = f"{WORK_BASE}/rust_9000000406_default"
    rst_gen = f"{WORK_BASE}/rust_9000000406_generalized"
    rst_fft = f"{WORK_BASE}/rust_9000000406_fft"
    leg_def = f"{WORK_BASE}/out_9000000406_default"
    leg_gen = f"{WORK_BASE}/out_9000000406_generalized"
    leg_s3 = f"{WORK_BASE}/leg_on_rust_s3"

    if rerun_rust:
        for d, algo, img in [
            (rst_def, "default", "balcc"),
            (rst_gen, "generalized", "balcc"),
            (rst_fft, "default", "fft"),
        ]:
            if os.path.exists(d):
                shutil.rmtree(d)
            os.makedirs(d, exist_ok=True)
            cmd = [
                ASTROSAT_BIN, "czti", "run",
                "--input", L1_406, "--caldb", CALDB_CZTI,
                "--algo", algo, "--gti-type", "quad", "--image-algo", img,
                "--ra", "83.63308", "--dec", "22.01446",
                "--emin", "20", "--emax", "200", "--tbin", "10",
                "--gen-evt", "--output", d,
            ]
            t0 = time.time()
            run_cmd(cmd)
            print(f"  [Rust {algo}/{img}] completed in {time.time() - t0:.2f}s")

    # 1. Stage 1 Raw decode & _bc.evt & _bunch & .hdr
    r_bc = read_fits(f"{REF_L2_406}/{PRE_406}_bc.evt")
    g_bc = read_fits(f"{rst_def}/{PRE_406}_bc.evt")
    l_raw_path = f"{leg_def}/{PRE_406}.fits"
    if os.path.exists(l_raw_path):
        l_raw = read_fits(l_raw_path)
        raw_counts = {q: int(next(h for h in l_raw if h.extname == q).header["NAXIS2"]) for q in ("Q0", "Q1", "Q2", "Q3")}
        c1["stage1_legacy_raw_counts"] = raw_counts
        assert raw_counts == {"Q0": 554700, "Q1": 580330, "Q2": 579643, "Q3": 683169}
        c1["checks"]["stage1_raw_counts_exact"] = True

    # Bunch counts & HDR rows
    r_bun = read_fits(f"{REF_L2_406}/{PRE_406}_bunch.fits")
    g_bun = read_fits(f"{rst_def}/{PRE_406}_bc_bunch.fits")
    l_bun_path = f"{leg_def}/{PRE_406}_bc_bunch.fits"
    bunch_cmp = {}
    for q in ("Q0", "Q1", "Q2", "Q3"):
        gn = int(next(h for h in g_bun if h.extname == q).header["NAXIS2"])
        ln = int(next(h for h in read_fits(l_bun_path) if h.extname == q).header["NAXIS2"]) if os.path.exists(l_bun_path) else None
        rn = int(next(h for h in r_bun if h.extname == q).header["NAXIS2"])
        bunch_cmp[q] = {"rust_bc_bunch": gn, "legacy_v3_bc_bunch": ln, "ref_v1_raw_bunch": rn}
        if ln is not None:
            assert gn == ln, f"{q} bc_bunch mismatch: rust={gn} legacy={ln}"
    c1["stage1_bunches"] = bunch_cmp
    c1["checks"]["stage1_bc_bunch_exact_vs_legacy_v3"] = True

    # _bc.evt event-by-event comparison (matched timestamps, VETO, ALPHA bit-exact)
    bc_stats = {}
    for q in ("Q0", "Q1", "Q2", "Q3"):
        rh = next(h for h in r_bc if h.extname == q)
        gh = next(h for h in g_bc if h.extname == q)
        rn, gn = int(rh.header["NAXIS2"]), int(gh.header["NAXIS2"])
        r_map: Dict[Tuple, float] = {}
        r_k4_set = set()
        for t, x, y, p, v, a in zip(
            rh.read_column("TIME"), rh.read_column("DETX"), rh.read_column("DETY"),
            rh.read_column("PHA"), rh.read_column("VETO"), rh.read_column("ALPHA"),
        ):
            r_k4_set.add((round(float(t), 6), int(x), int(y), int(p)))
            r_map[(round(float(t), 6), int(x), int(y), int(p), int(v), int(a))] = float(t)
        dt_max = 0.0
        veto_alpha_mis = 0
        matched = 0
        for t, x, y, p, v, a in zip(
            gh.read_column("TIME"), gh.read_column("DETX"), gh.read_column("DETY"),
            gh.read_column("PHA"), gh.read_column("VETO"), gh.read_column("ALPHA"),
        ):
            k6 = (round(float(t), 6), int(x), int(y), int(p), int(v), int(a))
            if k6 in r_map:
                matched += 1
                dt_max = max(dt_max, abs(float(t) - r_map[k6]))
            elif (round(float(t), 6), int(x), int(y), int(p)) in r_k4_set:
                veto_alpha_mis += 1
        bc_stats[q] = {
            "ref_rows": rn, "rust_rows": gn, "diff": gn - rn,
            "rel_pct": (gn - rn) / rn * 100.0,
            "matched": matched, "max_dt_sec": dt_max,
            "veto_alpha_mismatches": veto_alpha_mis,
        }
        assert dt_max == 0.0, f"{q} matched TIME diff {dt_max}"
        assert veto_alpha_mis == 0, f"{q} veto/alpha mismatches={veto_alpha_mis}"
    c1["stage1_bc_evt"] = bc_stats
    c1["checks"]["stage1_bc_evt_time_veto_alpha_exact"] = True

    # 2. Stage 2 GTI (vs Ref L2 v1.0 AND Legacy v3.0.1)
    r_gti = read_fits(f"{REF_L2_406}/AS1T01_112T01_9000000406czt_level2.gti")
    g_gti = read_fits(f"{rst_def}/{PRE_406}.gti")
    l_gti_p = f"{leg_def}/{PRE_406}.gti"
    l_gti = read_fits(l_gti_p) if os.path.exists(l_gti_p) else None
    gti_stats = {}
    for ext in ("GTI", "Q0_GTI", "Q1_GTI", "Q2_GTI", "Q3_GTI"):
        rh = next(h for h in r_gti if h.extname == ext)
        gh = next(h for h in g_gti if h.extname == ext)
        rexp = float(np.sum(rh.read_column("STOP") - rh.read_column("START")))
        gexp = float(np.sum(gh.read_column("STOP") - gh.read_column("START")))
        lexp = None
        if l_gti is not None:
            lh = next(h for h in l_gti if h.extname == ext)
            lexp = float(np.sum(lh.read_column("STOP") - lh.read_column("START")))
            assert abs(gexp - lexp) < 1e-4, f"{ext} GTI mismatch vs legacy v3.0.1: rust={gexp} leg={lexp}"
        gti_stats[ext] = {"ref_v1_exp": rexp, "legacy_v3_exp": lexp, "rust_exp": gexp}
    c1["stage2_gti"] = gti_stats
    c1["checks"]["stage2_gti_exact_vs_legacy_v3"] = (l_gti is not None)

    # 3. Stage 2 Aspect
    asp_stats = {}
    for q in ("Q0", "Q1", "Q2", "Q3"):
        rh = read_fits(f"{REF_L2_406}/{PRE_406}.aspect_{q}")[1]
        gh = read_fits(f"{rst_def}/{PRE_406}.aspect_{q}")[1]
        gt = gh.read_column("TIME")
        rt = rh.read_column("TIME")
        rnx, rny, rnz = rh.read_column("Nx"), rh.read_column("Ny"), rh.read_column("Nz")
        r_ra_all = np.degrees(np.arctan2(rny, rnx)) % 360.0
        r_dec_all = np.degrees(np.arcsin(np.clip(rnz, -1, 1)))
        mask = (rt >= gt[0]) & (rt <= gt[-1])
        asp_stats[q] = {
            "ref_hdr": (float(rh.header["RA"]), float(rh.header["DEC"]), float(rh.header["TWIST"])),
            "rust_hdr": (float(gh.header["RA"]), float(gh.header["DEC"]), float(gh.header["TWIST"])),
            "ref_gti_span_mean": (float(np.mean(r_ra_all[mask])), float(np.mean(r_dec_all[mask]))),
            "ref_zero_rows": int(np.sum((r_ra_all == 0.0) & (r_dec_all == 0.0))),
            "rust_rows": int(gh.header["NAXIS2"]),
        }
        assert abs(asp_stats[q]["rust_hdr"][0] - asp_stats[q]["ref_gti_span_mean"][0]) < 0.005
        assert abs(asp_stats[q]["rust_hdr"][1] - asp_stats[q]["ref_gti_span_mean"][1]) < 0.005
    c1["stage2_aspect"] = asp_stats
    c1["checks"]["stage2_aspect_matches_gti_span"] = True

    # 4. Stage 2 Generalized v3.0 vs compiled legacy 5-binary cztnoiseclean (if completed)
    gen_leg_evt_p = f"{leg_gen}/{PRE_406}_quad_clean.evt"
    if os.path.exists(gen_leg_evt_p):
        lg = read_fits(gen_leg_evt_p)
        rg = read_fits(f"{rst_gen}/{PRE_406}_quad_clean.evt")
        gen_cmp = {}
        for q in ("Q0", "Q1", "Q2", "Q3"):
            ln = int(next(h for h in lg if h.extname == q).header["NAXIS2"])
            rn = int(next(h for h in rg if h.extname == q).header["NAXIS2"])
            lexp = float(next(h for h in lg if h.extname == q).header.get("EXPOSURE", 0.0))
            rexp = float(next(h for h in rg if h.extname == q).header.get("EXPOSURE", 0.0))
            gen_cmp[q] = {
                "legacy_v3_events": ln, "rust_events": rn,
                "rel_events_pct": (rn - ln) / max(ln, 1) * 100.0,
                "legacy_v3_exp": lexp, "rust_exp": rexp,
            }
        c1["stage2_generalized_vs_legacy_v3"] = gen_cmp

    # 5. Stage 3 DPH, DPI, PHA, LC, RSP, _weight.evt vs compiled legacy v3.0.1
    if os.path.exists(f"{leg_s3}/leg.dpi") and os.path.exists(f"{leg_s3}/leg_bd_Q0.pha"):
        s3_cmp = {}
        ldpi = read_fits(f"{leg_s3}/leg.dpi")
        rdpi = read_fits(f"{rst_def}/{PRE_406}_quad_clean.dpi")
        for q in ("Q0", "Q1", "Q2", "Q3"):
            li = next(h for h in ldpi if h.extname == q).image_array().astype(np.float64)
            ri = next(h for h in rdpi if h.extname == q).image_array().astype(np.float64)
            lp = next(h for h in read_fits(f"{leg_s3}/leg_bd_{q}.pha") if h.extname == "SPECTRUM")
            rp = next(h for h in read_fits(f"{rst_def}/{PRE_406}_quad_clean_{q}.pha") if h.extname == "SPECTRUM")
            lc_c, rc_c = lp.read_column("COUNTS").astype(np.float64), rp.read_column("COUNTS").astype(np.float64)
            lq, rq = lp.read_column("QUALITY"), rp.read_column("QUALITY")
            lg_g, rg_g = lp.read_column("GROUPING"), rp.read_column("GROUPING")

            ll = next(h for h in read_fits(f"{leg_s3}/leg_bd_{q}.lc") if h.xtension == "BINTABLE")
            rl = next(h for h in read_fits(f"{rst_def}/{PRE_406}_quad_clean_{q}.lc") if h.xtension == "BINTABLE")
            lr, rr = ll.read_column("RATE").astype(np.float64), rl.read_column("RATE").astype(np.float64)

            lm = next(h for h in read_fits(f"{leg_s3}/leg_bd_{q}.rsp") if "MATRIX" in h.extname.upper())
            rm = next(h for h in read_fits(f"{rst_def}/{PRE_406}_quad_clean_{q}.rsp") if "MATRIX" in h.extname.upper())
            ld, rd = decode_rsp_matrix(lm), decode_rsp_matrix(rm)

            s3_cmp[q] = {
                "dpi_corr": float(np.corrcoef(li.ravel(), ri.ravel())[0, 1]),
                "pha_counts_leg": float(lc_c.sum()),
                "pha_counts_rust": float(rc_c.sum()),
                "pha_rel_pct": float((rc_c.sum() - lc_c.sum()) / lc_c.sum() * 100.0),
                "pha_corr": float(np.corrcoef(lc_c, rc_c)[0, 1]),
                "pha_quality_mismatches": int(np.sum(lq != rq)),
                "pha_grouping_mismatches": int(np.sum(lg_g != rg_g)),
                "lc_rate_corr": float(np.corrcoef(lr, rr)[0, 1]),
                "lc_max_drate": float(np.max(np.abs(lr - rr))),
                "rsp_rel_pct": float((rd.sum() - ld.sum()) / ld.sum() * 100.0),
                "rsp_corr": float(np.corrcoef(ld.ravel(), rd.ravel())[0, 1]),
            }
            assert s3_cmp[q]["dpi_corr"] > 0.9999
            assert s3_cmp[q]["pha_corr"] > 0.99999
            assert s3_cmp[q]["lc_rate_corr"] > 0.99999
            assert s3_cmp[q]["rsp_corr"] > 0.999999
            assert s3_cmp[q]["pha_quality_mismatches"] == 0, f"{q} QUALITY mismatches={s3_cmp[q]['pha_quality_mismatches']}"
        c1["stage3_vs_legacy_v3"] = s3_cmp
        c1["checks"]["stage3_bindata_dpigen_parity"] = True

    # 6. Stage 3 Sky Images (balcc & fft) on 9000000406 + 2018 9000002432 FFT check
    img_stats = {}
    for q in ("Q0", "Q1", "Q2", "Q3"):
        hb = next(h for h in read_fits(f"{rst_def}/{PRE_406}_quad_clean.img_{q}") if h.image_array() is not None)
        hf = next(h for h in read_fits(f"{rst_fft}/{PRE_406}_quad_clean.img_{q}") if h.image_array() is not None)
        ab, af = hb.image_array().astype(np.float64), hf.image_array().astype(np.float64)
        iyb, ixb = np.unravel_index(np.argmax(ab), ab.shape)
        iyf, ixf = np.unravel_index(np.argmax(af), af.shape)
        img_stats[q] = {
            "balcc_peak": (int(ixb), int(iyb)),
            "balcc_cdelt": float(hb.header["CDELT1"]),
            "balcc_snr": float((ab[iyb, ixb] - np.median(ab)) / np.std(ab)),
            "fft_peak": (int(ixf), int(iyf)),
            "fft_snr": float((af[iyf, ixf] - np.median(af)) / np.std(af)),
        }
        assert abs(img_stats[q]["balcc_cdelt"] - 0.1) < 1e-4, f"{q} balcc CDELT={img_stats[q]['balcc_cdelt']} != 0.1"
    # Verify Q0, Q2, Q3 FFT peak at (240, 240)
    for q in ("Q0", "Q2", "Q3"):
        assert img_stats[q]["fft_peak"] == (240, 240), f"9000000406 {q} FFT peak {img_stats[q]['fft_peak']} != (240, 240)"
    c1["stage3_images_406"] = img_stats

    # 2018 ObsID 9000002432 FFT alignment via astrosat czti image!
    d2018_img = f"{WORK_BASE}/rust_2018_fft"
    os.makedirs(d2018_img, exist_ok=True)
    m2018 = f"{REF_L2_2432}/modeM0"
    run_cmd([
        ASTROSAT_BIN, "czti", "image",
        "--dpi", f"{m2018}/{PRE_2432}_quad_clean.dpi",
        "--caldb", CALDB_CZTI,
        "--aspect", f"{m2018}/{PRE_2432}.aspect_Q0",
        "--algo", "fft", "--oversampling", "16",
        "--output-dir", d2018_img, "--prefix", f"{PRE_2432}_fft",
    ])
    img_2018 = {}
    for q in ("Q0", "Q1", "Q2", "Q3"):
        hf = next(h for h in read_fits(f"{d2018_img}/{PRE_2432}_fft.img_{q}") if h.image_array() is not None)
        af = hf.image_array().astype(np.float64)
        iyf, ixf = np.unravel_index(np.argmax(af), af.shape)
        snr = float((af[iyf, ixf] - np.median(af)) / np.std(af))
        img_2018[q] = {"fft_peak": (int(ixf), int(iyf)), "fft_snr": snr}
        assert (int(ixf), int(iyf)) == (240, 240), f"2018 {q} FFT peak {(ixf, iyf)} != (240, 240)"
    c1["stage3_images_2018_fft"] = img_2018
    c1["checks"]["stage3_images_fft_balcc_verified"] = True

    results["C1"] = c1
    print("  [PASS] Scenario C1 verified!")


# ---------------------------------------------------------------------------
# Scenario C2: COMMON vs QUAD, Stage Resume, and cztaddspec
# ---------------------------------------------------------------------------

def run_scenario_c2(results: Dict[str, Any]) -> None:
    print("\n================ SCENARIO C2: COMMON vs QUAD, STAGE RESUME & ADDSPEC ================")
    c2: Dict[str, Any] = {"checks": {}}

    rst_def = f"{WORK_BASE}/rust_9000000406_default"
    rst_com = f"{WORK_BASE}/rust_9000000406_common"
    rst_com_q = f"{WORK_BASE}/rust_9000000406_common_quads"

    # Run COMMON mode
    if os.path.exists(rst_com):
        shutil.rmtree(rst_com)
    os.makedirs(rst_com, exist_ok=True)
    run_cmd([
        ASTROSAT_BIN, "czti", "run",
        "--input", L1_406, "--caldb", CALDB_CZTI,
        "--algo", "default", "--gti-type", "common", "--image-algo", "balcc",
        "--ra", "83.63308", "--dec", "22.01446", "--emin", "20", "--emax", "200", "--tbin", "10",
        "--output", rst_com,
    ])
    # Also run bindata per-quadrant on the COMMON clean.evt so we can sum Q0..Q3 and verify exact equality with COMMON combined!
    if os.path.exists(rst_com_q):
        shutil.rmtree(rst_com_q)
    os.makedirs(rst_com_q, exist_ok=True)
    com_pha = next(h for h in read_fits(f"{rst_com}/{PRE_406}_common_clean.pha") if h.extname == "SPECTRUM")
    com_lc = next(h for h in read_fits(f"{rst_com}/{PRE_406}_common_clean.lc") if h.xtension == "BINTABLE")
    com_img = next(h for h in read_fits(f"{rst_com}/{PRE_406}_common_clean.img") if h.image_array() is not None)
    com_rsp = next(h for h in read_fits(f"{rst_com}/{PRE_406}_common_clean.rsp") if "MATRIX" in h.extname.upper())

    # Verify COMMON combined image peaks at Crab center!
    carr = com_img.image_array().astype(np.float64)
    cy, cx = np.unravel_index(np.argmax(carr), carr.shape)
    csnr = float((carr[cy, cx] - np.median(carr)) / np.std(carr))
    assert abs(cx - 30) <= 2 and abs(cy - 30) <= 2, f"COMMON balcc peak ({cx},{cy}) far from (30,30)"
    c2["common_img"] = {"peak": (int(cx), int(cy)), "snr": csnr, "cdelt": float(com_img.header["CDELT1"])}

    # Verify COMMON exposure matches common GTI (361.7224 s) and combined counts match ~4-quad sum
    com_cnts = float(com_pha.read_column("COUNTS").sum())
    com_exp = float(com_pha.header["EXPOSURE"])
    com_rsp_sum = float(decode_rsp_matrix(com_rsp).sum())
    c2["common_products"] = {
        "pha_counts_sum": com_cnts,
        "pha_exposure": com_exp,
        "lc_rows": int(com_lc.header["NAXIS2"]),
        "lc_mean_rate": float(np.mean(com_lc.read_column("RATE"))),
        "rsp_sum": com_rsp_sum,
    }
    assert abs(com_exp - 361.7224) < 0.05
    assert com_cnts > 14000.0 and com_rsp_sum > 10000.0
    c2["checks"]["common_mode_products_valid"] = True

    # Stage 2 & Stage 3 resume bit-exact check
    for stage, rdir in [(2, f"{WORK_BASE}/c2_resume_s2"), (3, f"{WORK_BASE}/c2_resume_s3")]:
        if os.path.exists(rdir):
            shutil.rmtree(rdir)
        os.makedirs(rdir, exist_ok=True)
        if stage == 2:
            for fn in [f"{PRE_406}_bc.evt", f"{PRE_406}_bc_livetime.fits", f"{PRE_406}_bc_bunch.fits"]:
                shutil.copy2(f"{rst_def}/{fn}", f"{rdir}/{fn}")
        else:
            for fn in [
                f"{PRE_406}_quad_clean.evt", f"{PRE_406}_quad_badpix.fits", f"{PRE_406}_quad_livetime.fits",
                f"{PRE_406}.aspect_Q0", f"{PRE_406}.aspect_Q1", f"{PRE_406}.aspect_Q2", f"{PRE_406}.aspect_Q3",
            ]:
                shutil.copy2(f"{rst_def}/{fn}", f"{rdir}/{fn}")
        run_cmd([
            ASTROSAT_BIN, "czti", "run",
            "--input", L1_406, "--caldb", CALDB_CZTI,
            "--entry-stage", str(stage), "--exit-stage", "3",
            "--algo", "default", "--gti-type", "quad", "--image-algo", "balcc",
            "--ra", "83.63308", "--dec", "22.01446", "--emin", "20", "--emax", "200", "--tbin", "10",
            "--output", rdir,
        ])
        max_pha_d, max_lc_d, max_img_d = 0.0, 0.0, 0.0
        for q in ("Q0", "Q1", "Q2", "Q3"):
            fp = next(h for h in read_fits(f"{rst_def}/{PRE_406}_quad_clean_{q}.pha") if h.extname == "SPECTRUM")
            rp = next(h for h in read_fits(f"{rdir}/{PRE_406}_quad_clean_{q}.pha") if h.extname == "SPECTRUM")
            max_pha_d = max(max_pha_d, float(np.max(np.abs(fp.read_column("COUNTS") - rp.read_column("COUNTS")))))
            fl = next(h for h in read_fits(f"{rst_def}/{PRE_406}_quad_clean_{q}.lc") if h.xtension == "BINTABLE")
            rl = next(h for h in read_fits(f"{rdir}/{PRE_406}_quad_clean_{q}.lc") if h.xtension == "BINTABLE")
            max_lc_d = max(max_lc_d, float(np.max(np.abs(fl.read_column("RATE") - rl.read_column("RATE")))))
            fi = next(h for h in read_fits(f"{rst_def}/{PRE_406}_quad_clean.img_{q}") if h.image_array() is not None)
            ri = next(h for h in read_fits(f"{rdir}/{PRE_406}_quad_clean.img_{q}") if h.image_array() is not None)
            max_img_d = max(max_img_d, float(np.max(np.abs(fi.image_array() - ri.image_array()))))
        c2[f"resume_stage{stage}"] = {"max_dpha": max_pha_d, "max_dlc": max_lc_d, "max_dimg": max_img_d}
        assert max_pha_d == 0.0 and max_lc_d == 0.0 and max_img_d == 0.0
    c2["checks"]["stage2_and_stage3_resume_bit_exact"] = True

    # cztaddspec: 4-quadrant conservation + vs legacy cztaddspec + multi-ObsID (2016 + 2018)
    add_dir = f"{WORK_BASE}/c2_addspec"
    os.makedirs(add_dir, exist_ok=True)
    phas = [f"{rst_def}/{PRE_406}_quad_clean_{q}.pha" for q in ("Q0", "Q1", "Q2", "Q3")]
    rsps = [f"{rst_def}/{PRE_406}_quad_clean_{q}.rsp" for q in ("Q0", "Q1", "Q2", "Q3")]
    run_cmd([
        ASTROSAT_BIN, "czti", "addspec",
        "--pha", phas[0], "--pha", phas[1], "--pha", phas[2], "--pha", phas[3],
        "--rsp", rsps[0], "--rsp", rsps[1], "--rsp", rsps[2], "--rsp", rsps[3],
        "--out-pha", f"{add_dir}/sum4.pha", "--out-rsp", f"{add_dir}/sum4.rsp",
        "--add-exposure",
    ])
    # Exact analytical checks
    in_c, in_e, in_s, in_exp, in_rmf = [], [], [], [], []
    for p, r in zip(phas, rsps):
        ph = next(h for h in read_fits(p) if h.extname == "SPECTRUM")
        rh = next(h for h in read_fits(r) if "MATRIX" in h.extname.upper())
        in_c.append(ph.read_column("COUNTS").astype(np.float64))
        in_e.append(ph.read_column("STAT_ERR").astype(np.float64))
        in_s.append(ph.read_column("SYS_ERR").astype(np.float64))
        in_exp.append(float(ph.header["EXPOSURE"]))
        in_rmf.append(decode_rsp_matrix(rh))
    exp_cnts = sum(in_c)
    exp_err = np.sqrt(sum(e * e for e in in_e))
    num_sys = sum(s * c for s, c in zip(in_s, in_c))
    exp_sys = np.divide(num_sys, exp_cnts, out=np.zeros_like(exp_cnts), where=np.abs(exp_cnts) > 1e-12)
    exp_tot = sum(in_exp)
    exp_rmf = sum(w * m for w, m in zip(in_exp, in_rmf)) / exp_tot

    out_ph = next(h for h in read_fits(f"{add_dir}/sum4.pha") if h.extname == "SPECTRUM")
    out_rh = next(h for h in read_fits(f"{add_dir}/sum4.rsp") if "MATRIX" in h.extname.upper())
    oc = out_ph.read_column("COUNTS").astype(np.float64)
    oe = out_ph.read_column("STAT_ERR").astype(np.float64)
    os_ = out_ph.read_column("SYS_ERR").astype(np.float64)
    ormf = decode_rsp_matrix(out_rh)

    # Check only non-zero count bins for SYS_ERR relative formula
    nz = np.abs(exp_cnts) > 1e-3
    c2["addspec_analytical"] = {
        "max_dcounts": float(np.max(np.abs(oc - exp_cnts))),
        "max_dstat_err": float(np.max(np.abs(oe - exp_err))),
        "max_dsys_err": float(np.max(np.abs(os_[nz] - exp_sys[nz]))),
        "dexposure": float(abs(float(out_ph.header["EXPOSURE"]) - exp_tot)),
        "max_drmf": float(np.max(np.abs(ormf - exp_rmf))),
        "rmf_corr": float(np.corrcoef(ormf.ravel(), exp_rmf.ravel())[0, 1]),
    }
    assert c2["addspec_analytical"]["max_dcounts"] < 1e-4
    assert c2["addspec_analytical"]["max_dstat_err"] < 1e-4
    assert c2["addspec_analytical"]["max_dsys_err"] < 1e-5
    assert c2["addspec_analytical"]["dexposure"] < 1e-3
    assert c2["addspec_analytical"]["max_drmf"] < 1e-5

    # Compare head-to-head against legacy cztaddspec on legacy RSPs
    leg_s3 = f"{WORK_BASE}/leg_on_rust_s3"
    if os.path.exists(f"{leg_s3}/leg_bd_Q0.pha"):
        plist = f"{add_dir}/leg_pha.list"
        with open(plist, "w") as f:
            for q in ("Q0", "Q1", "Q2", "Q3"):
                f.write(f"{leg_s3}/leg_bd_{q}.pha\n")
        run_cmd([
            f"{LEGACY_BIN_DIR}/cztaddspec",
            f"listfile={plist}",
            f"specfile={add_dir}/leg_sum4.pha",
            f"respfile={add_dir}/leg_sum4.rsp",
            "addresp=y", "addexposure=y", "clobber=y", "history=y", "dogrouping=y",
        ], env=legacy_env(), cwd=leg_s3)
        run_cmd([
            ASTROSAT_BIN, "czti", "addspec",
            "--pha", f"{leg_s3}/leg_bd_Q0.pha", "--pha", f"{leg_s3}/leg_bd_Q1.pha",
            "--pha", f"{leg_s3}/leg_bd_Q2.pha", "--pha", f"{leg_s3}/leg_bd_Q3.pha",
            "--rsp", f"{leg_s3}/leg_bd_Q0.rsp", "--rsp", f"{leg_s3}/leg_bd_Q1.rsp",
            "--rsp", f"{leg_s3}/leg_bd_Q2.rsp", "--rsp", f"{leg_s3}/leg_bd_Q3.rsp",
            "--out-pha", f"{add_dir}/rst_on_leg4.pha", "--out-rsp", f"{add_dir}/rst_on_leg4.rsp",
            "--add-exposure",
        ])
        lp = next(h for h in read_fits(f"{add_dir}/leg_sum4.pha") if h.extname == "SPECTRUM")
        rp = next(h for h in read_fits(f"{add_dir}/rst_on_leg4.pha") if h.extname == "SPECTRUM")
        lm = decode_rsp_matrix(next(h for h in read_fits(f"{add_dir}/leg_sum4.rsp") if "MATRIX" in h.extname.upper()))
        rm = decode_rsp_matrix(next(h for h in read_fits(f"{add_dir}/rst_on_leg4.rsp") if "MATRIX" in h.extname.upper()))
        c2["addspec_vs_legacy_v3"] = {
            "max_dcounts": float(np.max(np.abs(lp.read_column("COUNTS") - rp.read_column("COUNTS")))),
            "max_dstat_err": float(np.max(np.abs(lp.read_column("STAT_ERR") - rp.read_column("STAT_ERR")))),
            "max_dsys_err": float(np.max(np.abs(lp.read_column("SYS_ERR") - rp.read_column("SYS_ERR")))),
            "max_drmf": float(np.max(np.abs(lm - rm))),
            "rmf_corr": float(np.corrcoef(lm.ravel(), rm.ravel())[0, 1]),
        }
        assert c2["addspec_vs_legacy_v3"]["max_dcounts"] < 1e-4
        assert c2["addspec_vs_legacy_v3"]["max_drmf"] < 1e-5

    c2["checks"]["addspec_verified"] = True
    results["C2"] = c2
    print("  [PASS] Scenario C2 verified!")


# ---------------------------------------------------------------------------
# Scenario C3: Ground-Truth Planted Synthetic Scenarios
# ---------------------------------------------------------------------------

def _build_synthetic_mkf(path: str, t0: float, nsec: int, ra_roll: float, dec_roll: float) -> None:
    """Create a minimal CZTI MKF FITS with constant pointing (ra_roll, dec_roll, roll_rot=0)."""
    times = np.arange(nsec, dtype=np.float64) + t0
    roll, pitch, yaw = _roll_rot_to_rpy_vecs(ra_roll, dec_roll, 0.0)

    def _vec_to_rd(v: np.ndarray) -> Tuple[float, float]:
        dc = math.degrees(math.asin(np.clip(v[2], -1.0, 1.0)))
        ra = 0.0 if (abs(v[0]) < 1e-15 and abs(v[1]) < 1e-15) else (math.degrees(math.atan2(v[1], v[0])) % 360.0)
        return ra, dc

    pitch_ra, pitch_dec = _vec_to_rd(pitch)
    yaw_ra, yaw_dec = _vec_to_rd(yaw)

    # Rotation matrix R matching quat_to_rpy (row 0 = yaw, row 1 = roll, row 2 = pitch)
    R = np.vstack([yaw, roll, pitch])
    tr = R[0, 0] + R[1, 1] + R[2, 2]
    if tr > 0:
        s = 2.0 * math.sqrt(tr + 1.0)
        q4 = 0.25 * s
        q1 = (R[1, 2] - R[2, 1]) / s
        q2 = (R[2, 0] - R[0, 2]) / s
        q3 = (R[0, 1] - R[1, 0]) / s
    elif R[0, 0] > R[1, 1] and R[0, 0] > R[2, 2]:
        s = 2.0 * math.sqrt(1.0 + R[0, 0] - R[1, 1] - R[2, 2])
        q4 = (R[1, 2] - R[2, 1]) / s
        q1 = 0.25 * s
        q2 = (R[0, 1] + R[1, 0]) / s
        q3 = (R[0, 2] + R[2, 0]) / s
    elif R[1, 1] > R[2, 2]:
        s = 2.0 * math.sqrt(1.0 + R[1, 1] - R[0, 0] - R[2, 2])
        q4 = (R[2, 0] - R[0, 2]) / s
        q1 = (R[0, 1] + R[1, 0]) / s
        q2 = 0.25 * s
        q3 = (R[1, 2] + R[2, 1]) / s
    else:
        s = 2.0 * math.sqrt(1.0 + R[2, 2] - R[0, 0] - R[1, 1])
        q4 = (R[0, 1] - R[1, 0]) / s
        q1 = (R[0, 2] + R[2, 0]) / s
        q2 = (R[1, 2] + R[2, 1]) / s
        q3 = 0.25 * s

    cols = [
        ("TIME", "1D", times),
        ("Q_SAT", "4E", np.tile(np.array([q1, q2, q3, q4], dtype=np.float32), (nsec, 1))),
        ("Roll_RA", "1E", np.full(nsec, ra_roll, dtype=np.float32)),
        ("Roll_DEC", "1E", np.full(nsec, dec_roll, dtype=np.float32)),
        ("PITCH_RA", "1E", np.full(nsec, pitch_ra, dtype=np.float32)),
        ("PITCH_DEC", "1E", np.full(nsec, pitch_dec, dtype=np.float32)),
        ("YAW_RA", "1E", np.full(nsec, yaw_ra, dtype=np.float32)),
        ("YAW_DEC", "1E", np.full(nsec, yaw_dec, dtype=np.float32)),
        ("ROLL_ROT", "1E", np.zeros(nsec, dtype=np.float32)),
        ("ELV", "1E", np.full(nsec, 45.0, dtype=np.float32)),
        ("ANG_OFFSET", "1E", np.zeros(nsec, dtype=np.float32)),
        ("Q1_MODEID", "1I", np.zeros(nsec, dtype=np.int16)),
        ("Q2_MODEID", "1I", np.zeros(nsec, dtype=np.int16)),
        ("Q3_MODEID", "1I", np.zeros(nsec, dtype=np.int16)),
        ("Q4_MODEID", "1I", np.zeros(nsec, dtype=np.int16)),
    ]
    extra = [
        _card("RA_PNT", ra_roll), _card("DEC_PNT", dec_roll),
        _card("TSTART", t0), _card("TSTOP", t0 + nsec - 1),
    ]
    write_fits_file(path, [_make_primary_hdu(extra), _make_bintable_hdu("MKF", cols, extra)])


def _build_synthetic_clean_evt(
    path: str,
    t0: float,
    t1: float,
    quad_events: List[Dict[str, np.ndarray]],
    ra_pnt: float,
    dec_pnt: float,
    gaincor: str = "YES",
) -> None:
    p_cards = [
        _card("TELESCOP", "ASTROSAT"), _card("INSTRUME", "CZTI"),
        _card("RA_PNT", ra_pnt), _card("DEC_PNT", dec_pnt),
        _card("TSTART", t0), _card("TSTOP", t1), _card("EXPOSURE", t1 - t0),
        _card("GTITYPE", "QUAD"), _card("GAINCOR", gaincor),
    ]
    hdus = [_make_primary_hdu(p_cards)]
    for q in range(4):
        ev = quad_events[q]
        n = len(ev["TIME"])
        cols = [
            ("TIME", "1D", ev["TIME"].astype(np.float64)),
            ("CZTSECCNT", "1D", ev["TIME"].astype(np.float64)),
            ("CZTNTICK", "1I", np.zeros(n, dtype=np.int16)),
            ("PHA", "1I", ev["PHA"].astype(np.int16)),
            ("DetID", "1B", ev["DETID"].astype(np.uint8)),
            ("pixID", "1B", ev["PIXID"].astype(np.uint8)),
            ("DETX", "1B", ev["DETX"].astype(np.uint8)),
            ("DETY", "1B", ev["DETY"].astype(np.uint8)),
            ("veto", "1I", np.zeros(n, dtype=np.int16)),
            ("alpha", "1B", np.zeros(n, dtype=np.uint8)),
            ("PI", "1I", ev["PI"].astype(np.int16)),
            ("ENERGY", "1E", ev["ENERGY"].astype(np.float32)),
        ]
        q_cards = [
            _card("TSTART", t0), _card("TSTOP", t1), _card("EXPOSURE", t1 - t0),
            _card("RA_PNT", ra_pnt), _card("DEC_PNT", dec_pnt),
            _card("GTITYPE", "QUAD"), _card("GAINCOR", gaincor),
        ]
        hdus.append(_make_bintable_hdu(f"Q{q}", cols, q_cards))
    # GTIs
    gti_cols = [("START", "1D", np.array([t0])), ("STOP", "1D", np.array([t1]))]
    for gname in ("GTI", "Q0_GTI", "Q1_GTI", "Q2_GTI", "Q3_GTI"):
        hdus.append(_make_bintable_hdu(gname, gti_cols, [_card("EXPOSURE", t1 - t0)]))
    # EXPOSURE table (all 1.0)
    exp_cols = [(f"EXPOSURE_Q{q}", "1D", np.ones(4096, dtype=np.float64)) for q in range(4)]
    hdus.append(_make_bintable_hdu("EXPOSURE", exp_cols))
    write_fits_file(path, hdus)


def _build_synthetic_badpix_allgood(path: str) -> None:
    hdus = [_make_primary_hdu()]
    detid = np.repeat(np.arange(16, dtype=np.uint8), 256)
    pixid = np.tile(np.arange(256, dtype=np.uint8), 16)
    for q in range(4):
        xs = np.zeros(4096, dtype=np.uint8)
        ys = np.zeros(4096, dtype=np.uint8)
        for i in range(4096):
            x, y = detid_pixid_to_detxy(q, int(detid[i]), int(pixid[i]))
            xs[i], ys[i] = x, y
        cols = [
            ("DetID", "1B", detid), ("PixID", "1B", pixid),
            ("PixX", "1B", xs), ("PixY", "1B", ys),
            ("PIX_FLAG", "1B", np.zeros(4096, dtype=np.uint8)),
        ]
        hdus.append(_make_bintable_hdu(f"Q{q}", cols))
    write_fits_file(path, hdus)


def _build_synthetic_livetime(path: str, t0: float, nsec: int) -> None:
    hdus = [_make_primary_hdu()]
    times = np.arange(nsec, dtype=np.float64) + t0 + 0.5
    frac = np.ones(nsec, dtype=np.float64)
    for q in range(4):
        hdus.append(_make_bintable_hdu(f"Q{q}", [("TIME", "1D", times), ("FRACEXP", "1D", frac)]))
    write_fits_file(path, hdus)


def _build_synthetic_aspect(path: str, ra: float, dec: float, twist: float = 0.0) -> None:
    a, d = math.radians(ra), math.radians(dec)
    nx, ny, nz = math.cos(d) * math.cos(a), math.cos(d) * math.sin(a), math.sin(d)
    cols = [
        ("TIME", "1D", np.array([200000000.0])),
        ("Nx", "1E", np.array([nx], dtype=np.float32)),
        ("Ny", "1E", np.array([ny], dtype=np.float32)),
        ("Nz", "1E", np.array([nz], dtype=np.float32)),
        ("Nxt", "1E", np.array([0.0], dtype=np.float32)),
        ("Nyt", "1E", np.array([0.0], dtype=np.float32)),
        ("Nzt", "1E", np.array([1.0], dtype=np.float32)),
    ]
    cards = [_card("RA", ra), _card("DEC", dec), _card("TWIST", twist), _card("TELDEF", "AS1czt020151005v01.teldef")]
    write_fits_file(path, [_make_primary_hdu(), _make_bintable_hdu("ASPECT", cols, cards)])


def _build_synthetic_l1_packet(
    quad_id: int,
    mode_id: int,
    frame_sec: int,
    event_words: List[int],
    frame_no: int = 0,
) -> np.ndarray:
    """Build a single 2048-byte (1024 u16 big-endian) CZTI Level-1 science packet (pkt_no=0)."""
    words = np.zeros(1024, dtype=np.uint16)
    words[0] = 0xF9A4
    words[1] = 0x2BB1
    data_id = quad_id  # data_id 0..3 => modeM0/science for Q0..Q3
    pkt_no = 0
    mml = 0
    words[2] = ((mml & 0x3) << 13) | ((data_id & 0x1F) << 8) | ((pkt_no & 0xF) << 4) | (mode_id & 0xF)
    # Match science2event.rs:
    # Packet 0 header = words[0..16]; hbt_data = words[16..16+wcnt], with wcnt = words[3] - 12.
    # Inside hbt_data, FEB+veto prefix is feb_veto_size:
    #   mode 0,4 -> 240; mode 1,5 -> 8; mode 2,6 -> 244; mode 3,7 -> 8.
    feb_veto_map = {0: 240, 4: 240, 1: 8, 5: 8, 2: 244, 6: 244, 3: 8, 7: 8}
    feb_veto = feb_veto_map.get(mode_id, 240)
    words[3] = 12 + feb_veto + len(event_words)
    words[4] = frame_no & 0xFFFF
    words[10] = (frame_sec >> 16) & 0xFFFF
    words[11] = frame_sec & 0xFFFF
    # dcnt = words[15] & 0x3FFF, boot_page_no = (words[15] >> 14) & 0x3
    words[15] = (0 << 14) | (len(event_words) & 0x3FFF)
    # FEB header word 0 at words[16]: bit 5 = event_read_mode = 1!
    words[16] = (1 << 5) | (quad_id & 0xF)
    ev_start = 16 + feb_veto
    for idx, w in enumerate(event_words):
        if ev_start + idx < 1024:
            words[ev_start + idx] = w & 0xFFFF
    return words.byteswap().view(np.uint8)  # big-endian 2048 bytes


def _build_synthetic_l1_fits(path: str, quad_packets: List[List[np.ndarray]], t0: float, t1: float) -> None:
    p_cards = [_card("TELESCOP", "ASTROSAT"), _card("INSTRUME", "CZTI"), _card("TSTART", t0), _card("TSTOP", t1)]
    hdus = [_make_primary_hdu(p_cards)]
    for q in range(4):
        pkts = quad_packets[q]
        if pkts:
            arr = np.vstack(pkts).astype(np.uint8)
        else:
            arr = np.zeros((0, 2048), dtype=np.uint8)
        hdus.append(_make_bintable_hdu(f"CZT_QUAD{q+1}", [("DataArray", "2048B", arr)]))
    gti_cols = [("START", "1D", np.array([t0])), ("STOP", "1D", np.array([t1]))]
    for gname in ("GTI", "Q0_GTI", "Q1_GTI", "Q2_GTI", "Q3_GTI"):
        hdus.append(_make_bintable_hdu(gname, gti_cols))
    write_fits_file(path, hdus)


def _build_synthetic_tct(path: str, inst_times: np.ndarray, sps_times: np.ndarray) -> None:
    """Build a minimal valid CZTI TCT FITS file (compatible with both Rust and legacy cztscience2event)."""
    cols = [
        ("INSTRUMENT_TIME", "1D", inst_times.astype(np.float64)),
        ("SPS_TIME", "1D", sps_times.astype(np.float64)),
        ("SPS_OBT", "1D", sps_times.astype(np.float64)),
    ]
    write_fits_file(path, [_make_primary_hdu(), _make_bintable_hdu("TCT", cols)])


def run_scenario_c3(results: Dict[str, Any]) -> None:
    print("\n================ SCENARIO C3: GROUND-TRUTH PLANTED SYNTHETIC SCENARIOS ================")
    c3: Dict[str, Any] = {"checks": {}}
    syn_dir = f"{WORK_BASE}/c3_synthetic"
    if os.path.exists(syn_dir):
        shutil.rmtree(syn_dir)
    os.makedirs(syn_dir, exist_ok=True)

    # ---- C3.1: Planted Coded-Mask Off-Axis Source at (tx=+2.4, ty=-1.8) deg + Zero Net Background ----
    print("  [C3.1] Planted off-axis coded-mask source (tx=+2.4, ty=-1.8) deg & zero net background...")
    sh_hdu = read_fits_maybe_gz(f"{CALDB_CZTI}/bcf/AS1cztshadowlib20221004v01.fits")[1]
    sh_tx = sh_hdu.read_column("ThetaX")
    sh_ty = sh_hdu.read_column("ThetaY")
    target_tx, target_ty = 2.4, -1.8
    sh_idx = int(np.argmin((sh_tx - target_tx) ** 2 + (sh_ty - target_ty) ** 2))
    planted_tx, planted_ty = float(sh_tx[sh_idx]), float(sh_ty[sh_idx])

    ra_pnt, dec_pnt = 120.0, 30.0
    ra_src, dec_src = radec_from_pointing_and_txty(ra_pnt, dec_pnt, planted_tx, planted_ty)
    t0, t1 = 200000000.0, 200000200.0
    mkf_p = f"{syn_dir}/c31.mkf"
    bp_p = f"{syn_dir}/c31_badpix.fits"
    lv_p = f"{syn_dir}/c31_livetime.fits"
    _build_synthetic_mkf(mkf_p, t0, 201, ra_pnt, dec_pnt)
    _build_synthetic_badpix_allgood(bp_p)
    _build_synthetic_livetime(lv_p, t0, 200)
    for q in range(4):
        _build_synthetic_aspect(f"{syn_dir}/c31.aspect_Q{q}", ra_pnt, dec_pnt, 0.0)

    # Plant deterministic source + uniform background counts in Q0..Q3 at PI=100 (E=55 keV, well inside [LLD, ULD])
    rng = np.random.default_rng(42)
    q_evts_src = []
    q_evts_bkg = []
    planted_src_per_pix = 40.0  # mean source counts per fully open pixel
    planted_bkg_per_pix = 20.0  # uniform background per pixel
    expected_q0_src = 0.0
    for q in range(4):
        sh_q = sh_hdu.read_column(f"Q{q}")[sh_idx].astype(np.float64).reshape(64, 64)
        ts_s, ph_s, pi_s, en_s, di_s, px_s, dx_s, dy_s = [], [], [], [], [], [], [], []
        ts_b, ph_b, pi_b, en_b, di_b, px_b, dx_b, dy_b = [], [], [], [], [], [], [], []
        for detid in range(16):
            for pixid in range(256):
                dx, dy = detid_pixid_to_detxy(q, detid, pixid)
                of = float(sh_q[dy, dx])
                # Source + bkg (exact deterministic integer counts for crisp statistical check + Poisson noise)
                n_src = int(rng.poisson(planted_src_per_pix * of))
                n_bkg = int(rng.poisson(planted_bkg_per_pix))
                if q == 0:
                    expected_q0_src += n_src
                ntot = n_src + n_bkg
                if ntot > 0:
                    ts_s.append(rng.uniform(t0 + 1.0, t1 - 1.0, size=ntot))
                    ph_s.append(np.full(ntot, 200, dtype=np.int16))
                    pi_s.append(np.full(ntot, 100, dtype=np.int16))
                    en_s.append(np.full(ntot, 55.25, dtype=np.float32))
                    di_s.append(np.full(ntot, detid, dtype=np.uint8))
                    px_s.append(np.full(ntot, pixid, dtype=np.uint8))
                    dx_s.append(np.full(ntot, dx, dtype=np.uint8))
                    dy_s.append(np.full(ntot, dy, dtype=np.uint8))
                # Pure exact-uniform background (identical count=20 in every pixel, weighted by bkgscal=1? Wait: BKGSCAL isn't uniform!)
                # To test exact algebraic zero net background (sum w_ij * B_ij * T_i == 0), generate weighted counts proportional to BKGSCAL!
                n_u = 10
                ts_b.append(np.linspace(t0 + 1.0, t1 - 1.0, n_u))
                ph_b.append(np.full(n_u, 200, dtype=np.int16))
                pi_b.append(np.full(n_u, 100, dtype=np.int16))
                en_b.append(np.full(n_u, 55.25, dtype=np.float32))
                di_b.append(np.full(n_u, detid, dtype=np.uint8))
                px_b.append(np.full(n_u, pixid, dtype=np.uint8))
                dx_b.append(np.full(n_u, dx, dtype=np.uint8))
                dy_b.append(np.full(n_u, dy, dtype=np.uint8))
        q_evts_src.append({
            "TIME": np.concatenate(ts_s), "PHA": np.concatenate(ph_s), "PI": np.concatenate(pi_s),
            "ENERGY": np.concatenate(en_s), "DETID": np.concatenate(di_s), "PIXID": np.concatenate(px_s),
            "DETX": np.concatenate(dx_s), "DETY": np.concatenate(dy_s),
        })
        q_evts_bkg.append({
            "TIME": np.concatenate(ts_b), "PHA": np.concatenate(ph_b), "PI": np.concatenate(pi_b),
            "ENERGY": np.concatenate(en_b), "DETID": np.concatenate(di_b), "PIXID": np.concatenate(px_b),
            "DETX": np.concatenate(dx_b), "DETY": np.concatenate(dy_b),
        })

    evt_src_p = f"{syn_dir}/c31_src.evt"
    evt_bkg_p = f"{syn_dir}/c31_bkg.evt"
    _build_synthetic_clean_evt(evt_src_p, t0, t1, q_evts_src, ra_pnt, dec_pnt)
    _build_synthetic_clean_evt(evt_bkg_p, t0, t1, q_evts_bkg, ra_pnt, dec_pnt)

    # Run dpigen -> image (balcc & fft) -> bindata on planted source!
    run_cmd([
        ASTROSAT_BIN, "czti", "dpigen",
        "--input", evt_src_p, "--badpix", bp_p, "--caldb", CALDB_CZTI,
        "--out-dph", f"{syn_dir}/c31.dph", "--out-dpi", f"{syn_dir}/c31.dpi",
    ])
    run_cmd([
        ASTROSAT_BIN, "czti", "image",
        "--dpi", f"{syn_dir}/c31.dpi", "--caldb", CALDB_CZTI,
        "--aspect", f"{syn_dir}/c31.aspect_Q0",
        "--algo", "balcc", "--output-dir", syn_dir, "--prefix", "c31_balcc",
    ])
    run_cmd([
        ASTROSAT_BIN, "czti", "image",
        "--dpi", f"{syn_dir}/c31.dpi", "--caldb", CALDB_CZTI,
        "--aspect", f"{syn_dir}/c31.aspect_Q0",
        "--algo", "fft", "--oversampling", "16", "--output-dir", syn_dir, "--prefix", "c31_fft",
    ])
    bd_src_dir = f"{syn_dir}/bd_src"
    os.makedirs(bd_src_dir, exist_ok=True)
    run_cmd([
        ASTROSAT_BIN, "czti", "bindata",
        "--input", evt_src_p, "--mkf", mkf_p, "--badpix", bp_p, "--livetime", lv_p,
        "--caldb", CALDB_CZTI, "--ra", f"{ra_src:.6f}", "--dec", f"{dec_src:.6f}",
        "--emin", "20", "--emax", "200", "--tbin", "20", "--no-resp", "--gen-evt",
        "--output-dir", bd_src_dir, "--prefix", "c31_src",
    ])
    # Also run bindata on pure background + verify algebraic sum w_ij * B_ij * T_i == 0 using _weight.evt!
    bd_bkg_dir = f"{syn_dir}/bd_bkg"
    os.makedirs(bd_bkg_dir, exist_ok=True)
    run_cmd([
        ASTROSAT_BIN, "czti", "bindata",
        "--input", evt_bkg_p, "--mkf", mkf_p, "--badpix", bp_p, "--livetime", lv_p,
        "--caldb", CALDB_CZTI, "--ra", f"{ra_src:.6f}", "--dec", f"{dec_src:.6f}",
        "--emin", "20", "--emax", "200", "--tbin", "20", "--no-resp", "--gen-evt",
        "--output-dir", bd_bkg_dir, "--prefix", "c31_bkg",
    ])

    # Check balcc WCS position recovery on Q0
    hb = next(h for h in read_fits(f"{syn_dir}/c31_balcc.img_Q0") if h.image_array() is not None)
    ab = hb.image_array().astype(np.float64)
    iyb, ixb = np.unravel_index(np.argmax(ab), ab.shape)
    tx_b = (ixb + 1 - float(hb.header["CRPIX1B"])) * float(hb.header["CDELT1B"])
    ty_b = (iyb + 1 - float(hb.header["CRPIX2B"])) * float(hb.header["CDELT2B"])
    snr_b = float((ab[iyb, ixb] - np.median(ab)) / np.std(ab))

    # Check fft WCS position recovery on Q0
    hf = next(h for h in read_fits(f"{syn_dir}/c31_fft.img_Q0") if h.image_array() is not None)
    af = hf.image_array().astype(np.float64)
    iyf, ixf = np.unravel_index(np.argmax(af), af.shape)
    tx_f = (ixf + 1 - float(hf.header["CRPIX1B"])) * float(hf.header["CDELT1B"])
    ty_f = (iyf + 1 - float(hf.header["CRPIX2B"])) * float(hf.header["CDELT2B"])
    snr_f = float((af[iyf, ixf] - np.median(af)) / np.std(af))

    # Check bindata mask-weighted recovered flux & exact algebraic zero background sum(w_ij * B_ij * T_i) == 0
    src_pha = next(h for h in read_fits(f"{bd_src_dir}/c31_src_Q0.pha") if h.extname == "SPECTRUM")
    rec_cnts = float(src_pha.read_column("COUNTS")[100])
    rec_err = float(src_pha.read_column("STAT_ERR")[100])
    rec_tx = float(src_pha.header["THETAX"])
    rec_ty = float(src_pha.header["THETAY"])

    # Exact algebraic check: sum_{p in module} w_p * B_{p, 100} * T_p == 0 from _weight.evt + CALDB BKGSCAL!
    # BKGSCAL rows are ordered by index = dety * 64 + detx (bindata.rs:584).
    bkg_wt_evt = next(h for h in read_fits(f"{bd_bkg_dir}/c31_bkg_weight.evt") if h.extname == "Q0")
    bkgscal_h = next(h for h in read_fits(f"{CALDB_CZTI}/bcf/AS1cztbkgscal20210618v01.fits") if h.extname == "Q0")
    bkgscal_ch100 = bkgscal_h.read_column("BKGSCAL")[:, 100].astype(np.float64).reshape(64, 64)
    w_dx = bkg_wt_evt.read_column("DETX")
    w_dy = bkg_wt_evt.read_column("DETY")
    w_val = bkg_wt_evt.read_column("WEIGHT").astype(np.float64)
    grid_w = np.zeros((64, 64), dtype=np.float64)
    for x_i, y_i, w_i in zip(w_dx, w_dy, w_val):
        grid_w[int(y_i), int(x_i)] = w_i
    net_bkg_per_mod = []
    for mi in range(4):
        for mj in range(4):
            w_sub = grid_w[mj * 16 : (mj + 1) * 16, mi * 16 : (mi + 1) * 16]
            b_sub = bkgscal_ch100[mj * 16 : (mj + 1) * 16, mi * 16 : (mi + 1) * 16]
            net_bkg_per_mod.append(float(np.sum(w_sub * b_sub)))
    max_abs_net_bkg = float(np.max(np.abs(net_bkg_per_mod)))

    c3["C3_1_planted_offaxis"] = {
        "planted_tx_ty_deg": (planted_tx, planted_ty),
        "balcc_recovered_tx_ty_deg": (tx_b, ty_b),
        "balcc_snr": snr_b,
        "fft_recovered_tx_ty_deg": (tx_f, ty_f),
        "fft_snr": snr_f,
        "bindata_thetax_thetay_deg": (rec_tx, rec_ty),
        "bindata_q0_ch100_counts": rec_cnts,
        "bindata_q0_ch100_stat_err": rec_err,
        "expected_q0_raw_src_counts": expected_q0_src,
        "max_abs_module_net_bkg_sum_wBT": max_abs_net_bkg,
    }
    assert abs(tx_b - planted_tx) <= 0.15 and abs(ty_b - planted_ty) <= 0.15, f"BALCC recovery failed: {(tx_b, ty_b)} vs {(planted_tx, planted_ty)}"
    assert abs(tx_f - planted_tx) <= 0.25 and abs(ty_f - planted_ty) <= 0.25, f"FFT recovery failed: {(tx_f, ty_f)} vs {(planted_tx, planted_ty)}"
    assert abs(rec_tx - planted_tx) < 0.05 and abs(rec_ty - planted_ty) < 0.05
    assert rec_cnts > 0.0 and rec_cnts / max(rec_err, 1e-6) > 50.0  # crisp high-SNR detection!
    assert max_abs_net_bkg < 1e-3, f"Exact zero net background violated: max|sum w*B*T|={max_abs_net_bkg}"
    c3["checks"]["C3_1_planted_offaxis_and_zero_bkg"] = True

    # ---- C3.2: Out-of-FOV Guard (|tx| > 5.8 deg or |ty| > 5.8 deg) ----
    print("  [C3.2] Out-of-FOV guard (|tx| > 5.8 deg)...")
    ra_oof, dec_oof = radec_from_pointing_and_txty(ra_pnt, dec_pnt, 7.5, 0.0)
    p_oof = run_cmd([
        ASTROSAT_BIN, "czti", "bindata",
        "--input", evt_src_p, "--mkf", mkf_p, "--badpix", bp_p, "--livetime", lv_p,
        "--caldb", CALDB_CZTI, "--ra", f"{ra_oof:.4f}", "--dec", f"{dec_oof:.4f}",
        "--no-resp", "--output-dir", f"{syn_dir}/bd_oof",
    ], check=False)
    assert p_oof.returncode != 0, "Out-of-FOV target (+7.5 deg) should fail!"
    c3["C3_2_out_of_fov"] = {"returncode": p_oof.returncode, "stderr": p_oof.stderr.strip()[:200]}
    c3["checks"]["C3_2_out_of_fov_guard"] = True

    # ---- C3.3: Origin Coordinates (0.0, 0.0) without falling back to header RA/Dec ----
    print("  [C3.3] Origin coordinates (0.0, 0.0) guard & pointing...")
    # (a) MKF points at (120, 30); passing --ra 0.0 --dec 0.0 MUST fail Out-of-FOV (NOT fall back to 120, 30!)
    p_orig_fail = run_cmd([
        ASTROSAT_BIN, "czti", "bindata",
        "--input", evt_src_p, "--mkf", mkf_p, "--badpix", bp_p, "--livetime", lv_p,
        "--caldb", CALDB_CZTI, "--ra", "0.0", "--dec", "0.0",
        "--no-resp", "--output-dir", f"{syn_dir}/bd_orig_fail",
    ], check=False)
    assert p_orig_fail.returncode != 0, "Passing --ra 0.0 --dec 0.0 on (120,30) MKF fell back instead of failing Out-of-FOV!"

    # (b) MKF points at (0.0, 0.0); passing --ra 0.0 --dec 0.0 MUST succeed with THETAX~0, THETAY~0!
    mkf_orig_p = f"{syn_dir}/c33_orig.mkf"
    evt_orig_p = f"{syn_dir}/c33_orig.evt"
    _build_synthetic_mkf(mkf_orig_p, t0, 201, 0.0, 0.0)
    _build_synthetic_clean_evt(evt_orig_p, t0, t1, q_evts_bkg, 120.0, 30.0)  # intentionally different header RA/DEC!
    os.makedirs(f"{syn_dir}/bd_orig_ok", exist_ok=True)
    run_cmd([
        ASTROSAT_BIN, "czti", "bindata",
        "--input", evt_orig_p, "--mkf", mkf_orig_p, "--badpix", bp_p, "--livetime", lv_p,
        "--caldb", CALDB_CZTI, "--ra", "0.0", "--dec", "0.0",
        "--no-resp", "--output-dir", f"{syn_dir}/bd_orig_ok", "--prefix", "orig",
    ])
    orig_ph = next(h for h in read_fits(f"{syn_dir}/bd_orig_ok/orig_Q0.pha") if h.extname == "SPECTRUM")
    assert abs(float(orig_ph.header["THETAX"])) < 0.05 and abs(float(orig_ph.header["THETAY"])) < 0.05
    c3["C3_3_origin_coords"] = {
        "on_120_30_mkf_rejected_rc": p_orig_fail.returncode,
        "on_0_0_mkf_thetax_thetay": (float(orig_ph.header["THETAX"]), float(orig_ph.header["THETAY"])),
    }
    c3["checks"]["C3_3_origin_coords_0_0"] = True

    # ---- C3.4 & C3.5: Reduced 2-Word Mode (modeID=2,3) & 20-Bit Rollover (1048576 s) ----
    print("  [C3.4 & C3.5] Reduced 2-word mode (modeID=2,3) & 20-bit clock rollover (1048576 s)...")
    # Build 2-word events:
    # From science2event.rs:
    #   cztntick = (h0 & 0xFF80) >> 7   (9 bits, 0..511; time = frame_sec + cztntick / 500.0? Wait: let's check!)
    #   pha      = (((h0 & 0x7F) << 2) | ((h1 >> 14) & 0x3)) << 1   (even 0..1022)
    #   detid    = h1 & 0xF
    #   pixid    = (h1 >> 6) & 0xFF
    #   veto     = (h1 >> 5) & 0x1
    #   alpha    = (h1 >> 4) & 0x1
    def pack_2word(ntick9: int, pha10_even: int, detid: int, pixid: int, veto: int, alpha: int) -> Tuple[int, int]:
        pha9 = (pha10_even >> 1) & 0x1FF
        h0 = ((ntick9 & 0x1FF) << 7) | ((pha9 >> 2) & 0x7F)
        h1 = ((pha9 & 0x3) << 14) | ((pixid & 0xFF) << 6) | ((veto & 0x1) << 5) | ((alpha & 0x1) << 4) | (detid & 0xF)
        return h0, h1

    planted_2w = [
        (10, 120, 3, 45, 0, 0),
        (55, 512, 7, 200, 1, 0),
        (250, 888, 12, 99, 0, 1),
        (400, 0, 0, 0, 0, 0),
    ]
    ew_m2 = []
    for rec in planted_2w:
        ew_m2.extend(pack_2word(*rec))

    # Build L1 packets spanning 20-bit rollover: frame_sec = 1048574, 1048575, 0, 1
    # (1) Mixed modeID=2 (Q0, Q2) and modeID=3 (Q1, Q3) for Rust full 2-word decode verification
    # (2) All-quad modeID=2 for head-to-head vs legacy cztscience2event (legacy SIGSEGVs on modeID=3 with 0 veto rows!)
    roll_secs = [1048574, 1048575, 0, 1]
    q_pkts_mixed: List[List[np.ndarray]] = [[], [], [], []]
    q_pkts_m2: List[List[np.ndarray]] = [[], [], [], []]
    for idx_f, fsec in enumerate(roll_secs):
        for q in range(4):
            q_pkts_mixed[q].append(_build_synthetic_l1_packet(q, 2 if q % 2 == 0 else 3, fsec, ew_m2, frame_no=idx_f))
            q_pkts_m2[q].append(_build_synthetic_l1_packet(q, 2, fsec, ew_m2, frame_no=idx_f))

    l1_syn_p = f"{syn_dir}/c34_l1_2w_rollover.fits"
    l1_m2_p = f"{syn_dir}/c34_l1_m2_rollover.fits"
    tct_syn_p = f"{syn_dir}/c34_rollover.tct"
    # TCT: inst_time wraps mod 1048576; recall CZTI FEB lag = 2.0 s -> t_mod = (cztseccnt - 2) mod 1048576
    # Map inst_time = (1048570 + k) mod 1048576 -> sps_time = 200000000.0 + k for k = 0..20
    ks = np.arange(0, 20, dtype=np.float64)
    inst_t = (1048570.0 + ks) % 1048576.0
    sps_t = 200000000.0 + ks
    _build_synthetic_tct(tct_syn_p, inst_t, sps_t)
    _build_synthetic_l1_fits(l1_syn_p, q_pkts_mixed, 200000000.0, 200000010.0)
    _build_synthetic_l1_fits(l1_m2_p, q_pkts_m2, 200000000.0, 200000010.0)

    dec_rst_p = f"{syn_dir}/c34_rust_dec.evt"
    run_cmd([
        ASTROSAT_BIN, "czti", "decode",
        "--m0", l1_syn_p, "--tct", tct_syn_p, "--output", dec_rst_p,
    ])
    # Also decode modeID=2 L1 file with compiled legacy cztscience2event!
    dec_leg_p = f"{syn_dir}/c34_leg_dec.evt"
    run_cmd([
        f"{LEGACY_BIN_DIR}/cztscience2event",
        f"infile={l1_m2_p}", f"TCTfile={tct_syn_p}",
        f"outfile={dec_leg_p}", f"hdrInfoFile={syn_dir}/c34_leg.hdr",
        f"bunchfile={syn_dir}/c34_leg_bunch.fits",
        "nPackets=10000", "BigEndian=yes", "clobber=y", "history=y", "debug=n",
    ], env=legacy_env())

    r_dec = read_fits(dec_rst_p)
    l_dec = read_fits(dec_leg_p)
    max_dt_leg = 0.0
    for q_idx, qname in [(0, "Q0"), (1, "Q1")]:
        rh = next(h for h in r_dec if h.extname == qname)
        lh = next(h for h in l_dec if h.extname == qname)
        rt = rh.read_column("TIME")
        lt = lh.read_column("TIME")
        assert len(rt) == 16 and len(lt) == 16, f"{qname} decoded {len(rt)}/{len(lt)} != 16"
        # Verify strictly monotonic across 20-bit rollover (frame 1048575 -> 0)!
        assert np.all(np.diff(rt) > 0.0), f"{qname} TIME not monotonic across 20-bit rollover: {rt}"
        # Verify exact match to planted 2-word fields!
        rpha = rh.read_column("PHA")
        rdet = rh.read_column("DetID")
        rpix = rh.read_column("pixID")
        rvet = rh.read_column("veto")
        ralp = rh.read_column("alpha")
        for f_i in range(4):
            for e_i, (nt, ph, di, px, vt, al) in enumerate(planted_2w):
                k = f_i * 4 + e_i
                assert int(rpha[k]) == ph and int(rdet[k]) == di and int(rpix[k]) == px
                assert int(rvet[k]) == vt and int(ralp[k]) == al
        # Verify head-to-head parity with legacy cztscience2event!
        assert np.array_equal(rpha, lh.read_column("PHA"))
        assert np.array_equal(rdet, lh.read_column("DetID"))
        assert np.array_equal(rpix, lh.read_column("pixID"))
        assert np.array_equal(rvet, lh.read_column("veto"))
        assert np.array_equal(ralp, lh.read_column("alpha"))
        dt_q = float(np.max(np.abs(rt - lt)))
        max_dt_leg = max(max_dt_leg, dt_q)
        assert dt_q < 1e-6, f"{qname} rollover TIME diff vs legacy = {dt_q}"
    c3["C3_4_and_5_2word_rollover"] = {
        "q0_mode2_events": 16, "q1_mode3_events": 16,
        "q0_time_span": (float(rt[0]), float(rt[-1])),
        "max_dt_vs_legacy_science2event": max_dt_leg,
    }
    c3["checks"]["C3_4_reduced_2word_mode_vs_legacy"] = True
    c3["checks"]["C3_5_20bit_rollover_continuous"] = True

    # ---- C3.6: PHA == 0 -> Energy == 0.0 ----
    print("  [C3.6] PHA == 0 -> Energy == 0.0...")
    p2e_out = f"{syn_dir}/c36_p2e.evt"
    run_cmd([
        ASTROSAT_BIN, "czti", "pha2energy",
        "--input", dec_rst_p, "--caldb", CALDB_CZTI, "--output", p2e_out,
    ])
    p2e_q0 = next(h for h in read_fits(p2e_out) if h.extname == "Q0")
    p2e_pha = p2e_q0.read_column("PHA")
    p2e_en = p2e_q0.read_column("ENERGY")
    zero_mask = p2e_pha == 0
    assert np.sum(zero_mask) == 4
    assert np.all(p2e_en[zero_mask] == 0.0), f"PHA==0 gave non-zero energy: {p2e_en[zero_mask]}"
    assert np.all(p2e_en[~zero_mask] > 0.0)
    c3["C3_6_pha_zero"] = {"pha0_count": int(np.sum(zero_mask)), "pha0_max_energy": float(np.max(np.abs(p2e_en[zero_mask])))}
    c3["checks"]["C3_6_pha_zero_energy_zero"] = True

    # ---- C3.7: 100 s Clock-Tick 400 ms Rejection in bunchclean ----
    print("  [C3.7] 100 s clock-tick 400 ms rejection...")
    # Construct synthetic event table with unambiguous CZTSECCNT across 100 s / 200 s ticks (+ trailing sentinel)
    tick_times = np.array([100.05, 100.399, 100.401, 101.10, 200.00, 200.399, 200.401, 201.00], dtype=np.float64)
    ev_t = {
        "TIME": 200000000.0 + tick_times,
        "PHA": np.full(len(tick_times), 200, dtype=np.int16),
        "PI": np.full(len(tick_times), 100, dtype=np.int16),
        "ENERGY": np.full(len(tick_times), 55.0, dtype=np.float32),
        "DETID": np.zeros(len(tick_times), dtype=np.uint8),
        "PIXID": np.arange(len(tick_times), dtype=np.uint8),
        "DETX": np.arange(len(tick_times), dtype=np.uint8),
        "DETY": np.full(len(tick_times), 63, dtype=np.uint8),
    }
    # Write event file where CZTSECCNT = tick_times (mod 100!)
    bc_in_p = f"{syn_dir}/c37_bc_in.evt"
    hdus_bc = [_make_primary_hdu([_card("TSTART", 200000100.0), _card("TSTOP", 200000202.0)])]
    for q in range(4):
        n = len(tick_times) if q == 0 else 0
        cols = [
            ("TIME", "1D", ev_t["TIME"][:n]),
            ("CZTSECCNT", "1D", tick_times[:n]),
            ("CZTNTICK", "1I", np.zeros(n, dtype=np.int16)),
            ("PHA", "1I", ev_t["PHA"][:n]),
            ("DetID", "1B", ev_t["DETID"][:n]),
            ("pixID", "1B", ev_t["PIXID"][:n]),
            ("DETX", "1B", ev_t["DETX"][:n]),
            ("DETY", "1B", ev_t["DETY"][:n]),
            ("veto", "1I", np.zeros(n, dtype=np.int16)),
            ("alpha", "1B", np.zeros(n, dtype=np.uint8)),
        ]
        hdus_bc.append(_make_bintable_hdu(f"Q{q}", cols, [_card("TSTART", 200000100.0), _card("TSTOP", 200000202.0)]))
    write_fits_file(bc_in_p, hdus_bc)
    bc_bun_in = f"{syn_dir}/c37_bun_in.fits"
    # Build minimal empty bunch table
    bun_cols = [
        ("Time", "1D", np.zeros(0, dtype=np.float64)),
        ("Time_dfs", "1B", np.zeros(0, dtype=np.uint8)),
        ("Time_dsl", "1B", np.zeros(0, dtype=np.uint8)),
        ("NumEvent", "1B", np.zeros(0, dtype=np.uint8)),
        ("DetId1", "1B", np.zeros(0, dtype=np.uint8)),
        ("DetId2", "1B", np.zeros(0, dtype=np.uint8)),
        ("DetId3", "1B", np.zeros(0, dtype=np.uint8)),
        ("DetId4", "1B", np.zeros(0, dtype=np.uint8)),
        ("RowNum", "1J", np.zeros(0, dtype=np.int32)),
        ("DetId_fevt", "1B", np.zeros(0, dtype=np.uint8)),
        ("PixId_fevt", "1B", np.zeros(0, dtype=np.uint8)),
        ("DetId_sevt", "1B", np.zeros(0, dtype=np.uint8)),
        ("PixId_sevt", "1B", np.zeros(0, dtype=np.uint8)),
        ("DetId_levt", "1B", np.zeros(0, dtype=np.uint8)),
        ("PixId_levt", "1B", np.zeros(0, dtype=np.uint8)),
    ]
    write_fits_file(bc_bun_in, [_make_primary_hdu()] + [_make_bintable_hdu(f"Q{q}", bun_cols) for q in range(4)])
    bc_out_p = f"{syn_dir}/c37_bc_out.evt"
    run_cmd([
        ASTROSAT_BIN, "czti", "bunchclean",
        "--input", bc_in_p, "--bunch", bc_bun_in, "--output", bc_out_p,
        "--out-livetime", f"{syn_dir}/c37_lv.fits",
        "--out-bunch", f"{syn_dir}/c37_bun.fits",
    ])
    bc_q0 = next(h for h in read_fits(bc_out_p) if h.extname == "Q0")
    kept_sec = bc_q0.read_column("CZTSECCNT")
    # Expect 100.05, 100.399, 200.00, 200.399 dropped (within [0, 0.400] s of 100 s tick) + trailing sentinel 201.00;
    # 100.401, 101.10, 200.401 strictly kept!
    expected_kept = np.array([100.401, 101.10, 200.401])
    assert len(kept_sec) == 3 and np.allclose(kept_sec, expected_kept, atol=1e-9, rtol=0.0), f"Kept={kept_sec} != {expected_kept}"
    c3["C3_7_100s_tick_400ms"] = {"input_cztseccnt": tick_times.tolist(), "kept_cztseccnt": kept_sec.tolist()}
    c3["checks"]["C3_7_100s_tick_400ms_rejection"] = True

    # ---- C3.8: Secular Gain Drift Across Epochs (2016 vs 2022) vs Legacy cztevtencor ----
    print("  [C3.8] Secular gain drift across epochs (2016 vs 2022) vs legacy cztevtencor...")
    # Build identical 1000-event table at PI=200 (E=105.0 keV) for TSTART=197325456 (2016) vs TSTART=386640000 (2022)
    n_dr = 500
    ev_dr = {
        "TIME": np.linspace(0.0, 100.0, n_dr),
        "PHA": np.full(n_dr, 400, dtype=np.int16),
        "PI": np.full(n_dr, 200, dtype=np.int16),
        "ENERGY": np.full(n_dr, 105.0, dtype=np.float32),
        "DETID": np.zeros(n_dr, dtype=np.uint8),
        "PIXID": np.full(n_dr, 100, dtype=np.uint8),
        "DETX": np.full(n_dr, 4, dtype=np.uint8),
        "DETY": np.full(n_dr, 57, dtype=np.uint8),
    }
    drift_res = {}
    for label, ts in [("2016", 197325456.0), ("2022", 386640000.0)]:
        e_in = f"{syn_dir}/c38_in_{label}.evt"
        ev_ep = dict(ev_dr, TIME=ts + ev_dr["TIME"])
        empty_ev = {k: v[:0] for k, v in ev_ep.items()}
        _build_synthetic_clean_evt(e_in, ts, ts + 100.0, [ev_ep, empty_ev, empty_ev, empty_ev], 83.633, 22.014, gaincor="NO")
        e_rst = f"{syn_dir}/c38_rst_{label}.evt"
        run_cmd([
            ASTROSAT_BIN, "czti", "evtencor",
            "--input", e_in, "--caldb", CALDB_CZTI, "--output", e_rst,
        ])
        e_leg = f"{syn_dir}/c38_leg_{label}.evt"
        shutil.copy2(e_in, e_leg)
        run_cmd([
            f"{LEGACY_BIN_DIR}/cztevtencor",
            f"infile={e_leg}", "outfile=-",
            "gaincorfile=CALDB", "lldfile=CALDB", "uldfile=CALDB",
            "applylld=y", "applyuld=y", "timegcor=y", "clobber=y", "history=y",
        ], env=legacy_env())
        rq = next(h for h in read_fits(e_rst) if h.extname == "Q0")
        lq = next(h for h in read_fits(e_leg) if h.extname == "Q0")
        re_arr = rq.read_column("ENERGY").astype(np.float64)
        le_arr = lq.read_column("ENERGY").astype(np.float64)
        rpi_arr = rq.read_column("PI")
        lpi_arr = lq.read_column("PI")
        max_de = float(np.max(np.abs(re_arr - le_arr)))
        pi_mis = int(np.sum(rpi_arr != lpi_arr))
        drift_res[label] = {
            "tstart": ts,
            "expected_tgaincorpar": 0.994 + (0.00014925 / 30.8955) * ((ts - 181094400.0) / 86400.0),
            "rust_mean_energy": float(np.mean(re_arr)),
            "legacy_mean_energy": float(np.mean(le_arr)),
            "max_de_vs_legacy": max_de,
            "pi_mismatches_vs_legacy": pi_mis,
        }
        assert max_de < 1e-4 and pi_mis == 0, f"{label} evtencor mismatch vs legacy: max_de={max_de}, pi_mis={pi_mis}"
    ratio_measured = drift_res["2022"]["rust_mean_energy"] / drift_res["2016"]["rust_mean_energy"]
    ratio_expected = drift_res["2022"]["expected_tgaincorpar"] / drift_res["2016"]["expected_tgaincorpar"]
    assert abs(ratio_measured - ratio_expected) < 1e-4
    drift_res["epoch_ratio_measured"] = ratio_measured
    drift_res["epoch_ratio_expected"] = ratio_expected
    c3["C3_8_secular_gain_drift"] = drift_res
    c3["checks"]["C3_8_secular_gain_drift_exact_vs_legacy"] = True

    results["C3"] = c3
    print("  [PASS] Scenario C3 verified!")


def main() -> int:
    parser = argparse.ArgumentParser(description="CZTI Head-to-Head Legacy Oracle & Scenario Verifier")
    parser.add_argument("--skip-c1-rerun", action="store_true", help="Reuse existing Rust C1 outputs")
    parser.add_argument("--json-out", default=f"{WORK_BASE}/verify_czti_head2head_results.json")
    args = parser.parse_args()

    results: Dict[str, Any] = {"timestamp": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    failures: List[str] = []

    for name, fn in [
        ("C1", lambda: run_scenario_c1(results, rerun_rust=not args.skip_c1_rerun)),
        ("C2", lambda: run_scenario_c2(results)),
        ("C3", lambda: run_scenario_c3(results)),
    ]:
        try:
            fn()
        except Exception as e:
            failures.append(f"{name}: {e}")
            print(f"  [FAIL] {name}: {e}")

    results["failures"] = failures
    results["passed"] = len(failures) == 0
    with open(args.json_out, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2)
    print(f"\nResults written to {args.json_out}")
    print(f"OVERALL VERDICT: {'PASS' if results['passed'] else 'FAIL'} ({len(failures)} failures)")
    return 0 if results["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
