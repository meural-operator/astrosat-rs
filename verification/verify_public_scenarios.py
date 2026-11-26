#!/usr/bin/env python3
"""Comprehensive Public-Release End-User Scenario, Stress & Cross-Instrument Verification Suite.

Tests the standalone installed `astrosat` binary across:
  1. End-to-end researcher workflows (chained per-stage CZTI & LAXPC CLIs vs `run`,
     declarative TOML/JSON `--config` + CLI default overrides + `--json` schemas,
     idempotency & guardrails).
  2. Ground-truth end-to-end planted physics scenarios (synthetic L1 -> L2 spectrum
     power-law + Gaussian line recovery + deadtime, and synthetic L1 -> L2 -> bary -> timing
     17.3456 Hz pulsar across 2016-2017 leap-second boundary + LEO & Earth orbit motion).
  3. Adversarial & corrupt-input stress testing (0 panics/segfaults/backtraces).
  4. Cross-instrument consistency (`inspect --full`, `validate --strict`, `fitsverify 4.22`)
     and clean-environment (`env -i`) execution with zero hidden env dependencies.

All temporary test artifacts are strictly confined to `/mnt/d/ISRO/.tmp/public_scenarios`.
"""
from __future__ import annotations

import argparse
import glob
import json
import math
import os
import re
import shutil
import struct
import subprocess
import sys
import time
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fits_reader import read_fits  # noqa: E402

ROOT_DIR = "/mnt/d/ISRO"
DEFAULT_ASTROSAT_BIN = "/mnt/d/ISRO/.tmp/public_scenarios/bin/astrosat"
FITSVERIFY_BIN = "/mnt/d/ISRO/.toolchain/bin/fitsverify"
WORK_ROOT = "/mnt/d/ISRO/.tmp/public_scenarios/work"
REPORT_JSON_DEFAULT = "/mnt/d/ISRO/verification/reports/public_scenarios.json"

CZTI_L1_DIR = "/mnt/d/ISRO/data/sample/czti/czti_sample_data/20160331_T01_112T01_9000000406_level1"
CZTI_M0 = f"{CZTI_L1_DIR}/modeM0/AS1T01_112T01_9000000406cztM0_level1.fits"
CZTI_SS = f"{CZTI_L1_DIR}/modeSS/AS1T01_112T01_9000000406cztSS_level1.fits"
CZTI_TCT = f"{CZTI_L1_DIR}/aux/AS1T01_112T01_9000000406czt_level1.tct"
CZTI_MKF = f"{CZTI_L1_DIR}/AS1T01_112T01_9000000406czt_level2.mkf"
CZTI_THR = f"{CZTI_L1_DIR}/mkfThresholds.txt"
CZTI_ORB = "/mnt/d/ISRO/data/sample/czti_orb/AS1T01_112T01_9000000406czt_level1.orb"
CZTI_CALDB = "/mnt/d/ISRO/data/caldb/czti/data/as1/czti"

LAXPC_L1_DIR = "/mnt/d/ISRO/data/sample/laxpc_l1/20160107_G02_010T01_9000000252_level1/laxpc/01509"
LAXPC_EA1 = f"{LAXPC_L1_DIR}/lxp1/modeEA/AS1G02_010T01_9000000252lxp1EA_level1.fits"
LAXPC_BB1 = f"{LAXPC_L1_DIR}/lxp1/modeBB/AS1G02_010T01_9000000252lxp1BB_level1.fits"
LAXPC_TCT = f"{LAXPC_L1_DIR}/aux/AS1G02_010T01_9000000252lxp_level1.tct"
LAXPC_MKF = f"{LAXPC_L1_DIR}/AS1G02_010T01_9000000252lxp_level1.mkf"
LAXPC_ORB = f"{LAXPC_L1_DIR}/aux/aux1/AS1G02_010T01_9000000252lxp_level1.orb"
LAXPC_CALDB_SRC = "/mnt/d/ISRO/data/caldb/laxpc"
LAXPC_REF_SRC = "/mnt/d/ISRO/data/ref_src/laxpc"
JPLEPH_405 = "/mnt/d/ISRO/data/ephem/lnxp1600p2200.405"
JPLEPH_FITS_405 = "/mnt/d/ISRO/data/ephem/JPLEPH.405"
JPLEPH_FITS_200 = "/mnt/d/ISRO/data/ephem/JPLEPH.200"

FITS_EXTS = {
    ".fits", ".fit", ".fts", ".evt", ".evn", ".gti", ".gtis", ".bti",
    ".dph", ".dpi", ".img", ".img_q0", ".img_q1", ".img_q2", ".img_q3",
    ".img_all", ".pha", ".spec", ".lc", ".lcurv", ".rsp", ".rmf", ".arf",
    ".orb", ".mkf", ".att", ".hdr", ".frt",
    ".aspect_q0", ".aspect_q1", ".aspect_q2", ".aspect_q3",
}
SUMMARY_RE = re.compile(r"Verification found (\d+) warning\(s\) and (\d+) error\(s\)")
PANIC_MARKERS = (
    "panicked at",
    "RUST_BACKTRACE",
    "stack backtrace:",
    "Segmentation fault",
    "SIGSEGV",
    "memory allocation of",
    "fatal runtime error",
)


@dataclass
class CmdResult:
    cmd: List[str]
    returncode: int
    stdout: str
    stderr: str
    elapsed_s: float
    panicked: bool


def run_cli(
    bin_path: str,
    args: Sequence[str],
    env: Optional[Dict[str, str]] = None,
    timeout: int = 240,
) -> CmdResult:
    cmd = [bin_path, *args]
    use_env = env
    if use_env is None:
        use_env = os.environ.copy()
        use_env["TMPDIR"] = "/mnt/d/ISRO/.tmp"
    t0 = time.monotonic()
    proc = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        env=use_env,
        timeout=timeout,
        check=False,
    )
    dt = time.monotonic() - t0
    combined = proc.stdout + "\n" + proc.stderr
    panicked = (proc.returncode < 0) or any(m in combined for m in PANIC_MARKERS)
    return CmdResult(
        cmd=cmd,
        returncode=proc.returncode,
        stdout=proc.stdout,
        stderr=proc.stderr,
        elapsed_s=dt,
        panicked=panicked,
    )


def require_ok(res: CmdResult, label: str) -> CmdResult:
    if res.panicked:
        raise AssertionError(f"{label} panicked/crashed (rc={res.returncode}):\n{res.stderr[:1200]}")
    if res.returncode != 0:
        raise AssertionError(
            f"{label} failed (rc={res.returncode}): {' '.join(res.cmd)}\n"
            f"STDOUT: {res.stdout[:800]}\nSTDERR: {res.stderr[:1200]}"
        )
    return res


def is_fits_file(path: str) -> bool:
    low = path.lower()
    if any(low.endswith(ext) for ext in FITS_EXTS):
        return True
    try:
        with open(path, "rb") as fh:
            return fh.read(9) == b"SIMPLE  ="
    except OSError:
        return False


def collect_fits(target: str) -> List[str]:
    if os.path.isfile(target):
        return [target]
    out: List[str] = []
    for root, _, names in os.walk(target):
        for n in sorted(names):
            p = os.path.join(root, n)
            if is_fits_file(p):
                out.append(p)
    return sorted(set(out))


def verify_fitsverify(targets: Sequence[str]) -> Dict[str, Any]:
    files: List[str] = []
    for t in targets:
        files.extend(collect_fits(t))
    files = sorted(set(files))
    if not files:
        raise AssertionError(f"No FITS files found under {targets}")
    failures = []
    warnings_total = 0
    for f in files:
        p = subprocess.run([FITSVERIFY_BIN, f], capture_output=True, text=True, check=False)
        txt = p.stdout + "\n" + p.stderr
        m = SUMMARY_RE.search(txt)
        if not m:
            failures.append({"file": f, "errors": -1, "detail": txt[-800:]})
        else:
            w, e = int(m.group(1)), int(m.group(2))
            warnings_total += w
            if e != 0:
                failures.append({"file": f, "errors": e, "warnings": w, "detail": txt[-1000:]})
    return {
        "total_files": len(files),
        "failed_files": len(failures),
        "total_warnings": warnings_total,
        "failures": failures,
        "files": files,
    }


def make_temp_laxpc_caldb(dest: str) -> str:
    if os.path.exists(dest):
        shutil.rmtree(dest)
    os.makedirs(dest, exist_ok=True)
    for src_dir in (LAXPC_CALDB_SRC,):
        for name in os.listdir(src_dir):
            src = os.path.join(src_dir, name)
            if os.path.isfile(src):
                os.symlink(src, os.path.join(dest, name))
    for name in (
        "lxp10.channel",
        "lxp20.channel",
        "lxp30.channel",
        "lxp10cshm141v2.2.channel",
        "lxp10cshm175v2.2.channel",
        "lxp10cshm197v2.2.channel",
        "lxp10cshm212v2.2.channel",
        "backlxp10jan16.pha",
        "backlxp10apr18.pha",
        "backlxp30aug17.pha",
        "backfit.lxp10jan16",
    ):
        src = os.path.join(LAXPC_REF_SRC, name)
        if os.path.isfile(src):
            os.symlink(src, os.path.join(dest, name))
    return dest


# ============================================================================
# Pure-Python minimal FITS 4.0 writer for synthetic L1 & adversarial tests
# ============================================================================

def _fits_card(key: str, val: Any, comment: str = "") -> bytes:
    k = key.ljust(8)[:8]
    if val is None:
        return k.ljust(80).encode("ascii")
    if isinstance(val, bool):
        vstr = "T" if val else "F"
        body = f"{k}= {vstr:>20}"
    elif isinstance(val, int):
        body = f"{k}= {val:>20d}"
    elif isinstance(val, float):
        body = f"{k}= {val:>20.12E}"
    elif isinstance(val, str):
        esc = val.replace("'", "''")
        lit = f"'{esc}'"
        if len(lit) < 10:
            lit = f"'{esc:<8}'"
        body = f"{k}= {lit:<20}"
    else:
        raise TypeError(f"Unsupported FITS card type: {type(val)}")
    if comment:
        body = f"{body} / {comment}"
    return body[:80].ljust(80).encode("ascii")


def _pad_2880(buf: bytes) -> bytes:
    rem = len(buf) % 2880
    if rem == 0:
        return buf
    return buf + b"\x00" * (2880 - rem)


def _pad_header_2880(cards: List[bytes]) -> bytes:
    raw = b"".join(cards) + b"END".ljust(80)
    rem = len(raw) % 2880
    if rem:
        raw += b" " * (2880 - rem)
    return raw


def _build_bintable_hdu(
    extname: str,
    columns: List[Tuple[str, str, Optional[str], Any]],
    extra_cards: Optional[List[Tuple[str, Any, str]]] = None,
) -> bytes:
    """columns: list of (name, tform, tunit, data_array_or_list).
    Supports TFORMs: 'D' (f64), 'E' (f32), 'J' (i32), 'K' (i64), 'I' (i16), '<N>B' (uint8 bytes).
    """
    nrows = len(columns[0][3]) if columns else 0
    row_bytes = 0
    parsed_cols = []
    for name, tform, tunit, data in columns:
        m = re.match(r"^(\d*)([A-Z])$", tform.strip())
        assert m, f"Unsupported TFORM {tform}"
        rep = int(m.group(1)) if m.group(1) else 1
        code = m.group(2)
        elem_sz = {"D": 8, "E": 4, "K": 8, "J": 4, "I": 2, "B": 1}[code]
        col_sz = rep * elem_sz
        row_bytes += col_sz
        parsed_cols.append((name, tform, tunit, rep, code, data))

    cards = [
        _fits_card("XTENSION", "BINTABLE", "binary table extension"),
        _fits_card("BITPIX", 8, "8-bit bytes"),
        _fits_card("NAXIS", 2, "2-dimensional binary table"),
        _fits_card("NAXIS1", row_bytes, "width of table in bytes"),
        _fits_card("NAXIS2", nrows, "number of rows in table"),
        _fits_card("PCOUNT", 0, "size of special data area"),
        _fits_card("GCOUNT", 1, "one data group"),
        _fits_card("TFIELDS", len(columns), "number of fields in each row"),
        _fits_card("EXTNAME", extname, "extension name"),
    ]
    for idx, (name, tform, tunit, _rep, _code, _data) in enumerate(parsed_cols, start=1):
        cards.append(_fits_card(f"TTYPE{idx}", name))
        cards.append(_fits_card(f"TFORM{idx}", tform))
        if tunit:
            cards.append(_fits_card(f"TUNIT{idx}", tunit))
    if extra_cards:
        for k, v, c in extra_cards:
            cards.append(_fits_card(k, v, c))

    hdr_bytes = _pad_header_2880(cards)
    if nrows == 0:
        return hdr_bytes

    chunks: List[bytes] = []
    dtype_map = {"D": ">f8", "E": ">f4", "K": ">i8", "J": ">i4", "I": ">i2", "B": "u1"}
    col_arrays = []
    for _name, _tform, _tunit, rep, code, data in parsed_cols:
        arr = np.asarray(data, dtype=dtype_map[code])
        if rep > 1:
            arr = arr.reshape(nrows, rep)
        else:
            arr = arr.reshape(nrows, 1)
        col_arrays.append(arr.view(np.uint8).reshape(nrows, -1))
    table_bytes = np.hstack(col_arrays).tobytes()
    return hdr_bytes + _pad_2880(table_bytes)


def write_synthetic_fits(
    path: str,
    primary_cards: List[Tuple[str, Any, str]],
    ext_hdus: List[bytes],
) -> None:
    cards = [
        _fits_card("SIMPLE", True, "conforms to FITS standard"),
        _fits_card("BITPIX", 8, "8-bit bytes"),
        _fits_card("NAXIS", 0, "no primary data array"),
        _fits_card("EXTEND", len(ext_hdus) > 0, "FITS dataset may contain extensions"),
    ]
    for k, v, c in primary_cards:
        cards.append(_fits_card(k, v, c))
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "wb") as fh:
        fh.write(_pad_header_2880(cards))
        for h in ext_hdus:
            fh.write(h)


# ============================================================================
# Workflow 1: Chained Per-Stage CZTI CLI Pipeline
# ============================================================================

def test_workflow1_czti(bin_path: str) -> Dict[str, Any]:
    wf_dir = os.path.join(WORK_ROOT, "wf1_czti")
    if os.path.exists(wf_dir):
        shutil.rmtree(wf_dir)
    run_def_dir = os.path.join(wf_dir, "run_default")
    run_gen_dir = os.path.join(wf_dir, "run_generalized")
    chain_def_dir = os.path.join(wf_dir, "chain_default")
    chain_gen_dir = os.path.join(wf_dir, "chain_generalized")
    for d in (run_def_dir, run_gen_dir, chain_def_dir, chain_gen_dir):
        os.makedirs(d, exist_ok=True)

    # Baseline 1: `astrosat czti run` (default + --gen-evt)
    r_def = require_ok(
        run_cli(bin_path, [
            "--json", "czti", "run",
            "--input", CZTI_L1_DIR,
            "--caldb", CZTI_CALDB,
            "--gen-evt",
            "--output", run_def_dir,
        ]),
        "czti run default",
    )
    run_def_json = json.loads(r_def.stdout)

    # Baseline 2: `astrosat czti run --algo generalized`
    r_gen = require_ok(
        run_cli(bin_path, [
            "--json", "czti", "run",
            "--input", CZTI_L1_DIR,
            "--caldb", CZTI_CALDB,
            "--algo", "generalized",
            "--output", run_gen_dir,
        ]),
        "czti run generalized",
    )
    run_gen_json = json.loads(r_gen.stdout)

    # Chained per-stage Default pipeline
    s1_evt = os.path.join(chain_def_dir, "s1_decode.evt")
    s1_bunch = os.path.join(chain_def_dir, "s1_bunch.fits")
    s1_hdr = os.path.join(chain_def_dir, "s1.hdr")
    require_ok(run_cli(bin_path, [
        "--json", "czti", "decode", "--m0", CZTI_M0, "--tct", CZTI_TCT,
        "-o", s1_evt, "--bunch", s1_bunch, "--hdr", s1_hdr,
    ]), "czti decode")

    s2_evt = os.path.join(chain_def_dir, "s2_bc.evt")
    s2_bunch = os.path.join(chain_def_dir, "s2_bc_bunch.fits")
    s2_lv = os.path.join(chain_def_dir, "s2_bc_livetime.fits")
    require_ok(run_cli(bin_path, [
        "--json", "czti", "bunchclean", "-i", s1_evt, "--bunch", s1_bunch,
        "-o", s2_evt, "--out-bunch", s2_bunch, "--out-livetime", s2_lv,
    ]), "czti bunchclean")

    s3_evt = os.path.join(chain_def_dir, "s3_p2e.evt")
    require_ok(run_cli(bin_path, [
        "--json", "czti", "pha2energy", "-i", s2_evt, "--caldb", CZTI_CALDB,
        "--ss", CZTI_SS, "--tct", CZTI_TCT, "-o", s3_evt,
    ]), "czti pha2energy")

    s4_gti = os.path.join(chain_def_dir, "s4.gti")
    require_ok(run_cli(bin_path, [
        "--json", "czti", "gtigen", "--m0", CZTI_M0, "--mkf", CZTI_MKF,
        "--mkf-thresholds", CZTI_THR, "-o", s4_gti,
    ]), "czti gtigen")

    gaas_res = require_ok(run_cli(bin_path, [
        "--json", "czti", "gaas", "--mkf", CZTI_MKF, "--caldb", CZTI_CALDB,
        "--gti", s4_gti, "--gti-type", "quad", "-o", chain_def_dir, "--prefix", "chain_def",
    ]), "czti gaas")
    gaas_json = json.loads(gaas_res.stdout)

    s6_evt = os.path.join(chain_def_dir, "s6_ds.evt")
    require_ok(run_cli(bin_path, [
        "--json", "czti", "datasel", "-i", s3_evt, "--gti", s4_gti,
        "--gti-type", "quad", "-o", s6_evt,
    ]), "czti datasel")

    s7_evt = os.path.join(chain_def_dir, "s7_clean.evt")
    s7_bp = os.path.join(chain_def_dir, "s7_badpix.fits")
    s7_lv = os.path.join(chain_def_dir, "s7_livetime.fits")
    require_ok(run_cli(bin_path, [
        "--json", "czti", "clean", "-i", s6_evt, "--caldb", CZTI_CALDB,
        "--livetime", s2_lv, "--bunch", s2_bunch, "--algo", "default",
        "--gti-type", "quad", "-o", s7_evt, "--out-badpix", s7_bp, "--out-livetime", s7_lv,
    ]), "czti clean default")

    s8_evt = os.path.join(chain_def_dir, "s8_encor.evt")
    require_ok(run_cli(bin_path, [
        "--json", "czti", "evtencor", "-i", s7_evt, "--caldb", CZTI_CALDB, "-o", s8_evt,
    ]), "czti evtencor")

    s9_evt = os.path.join(chain_def_dir, "s9_evtclean.evt")
    require_ok(run_cli(bin_path, [
        "--json", "czti", "evtclean", "-i", s8_evt, "--badpix", s7_bp, "-o", s9_evt,
    ]), "czti evtclean")

    s10_dph = os.path.join(chain_def_dir, "s10.dph")
    s10_dpi = os.path.join(chain_def_dir, "s10.dpi")
    require_ok(run_cli(bin_path, [
        "--json", "czti", "dpigen", "-i", s9_evt, "--caldb", CZTI_CALDB,
        "--badpix", s7_bp, "--emin", "0.0", "--emax", "1000.0",
        "--out-dph", s10_dph, "--out-dpi", s10_dpi,
    ]), "czti dpigen")

    img_balcc_res = require_ok(run_cli(bin_path, [
        "--json", "czti", "image", "--dpi", s10_dpi, "--aspect", chain_def_dir,
        "--caldb", CZTI_CALDB, "-o", chain_def_dir, "--prefix", "chain_balcc", "--algo", "balcc",
    ]), "czti image balcc")
    img_balcc_json = json.loads(img_balcc_res.stdout)

    img_fft_res = require_ok(run_cli(bin_path, [
        "--json", "czti", "image", "--dpi", s10_dpi, "--aspect", chain_def_dir,
        "--caldb", CZTI_CALDB, "-o", chain_def_dir, "--prefix", "chain_fft",
        "--algo", "fft", "--oversampling", "16",
    ]), "czti image fft")
    img_fft_json = json.loads(img_fft_res.stdout)

    bindata_res = require_ok(run_cli(bin_path, [
        "--json", "czti", "bindata", "-i", s9_evt, "--mkf", CZTI_MKF,
        "--caldb", CZTI_CALDB, "--badpix", s7_bp, "--livetime", s7_lv,
        "-o", chain_def_dir, "--prefix", "chain_def", "--gen-evt",
    ]), "czti bindata")
    bindata_json = json.loads(bindata_res.stdout)

    phas = [os.path.join(chain_def_dir, f"chain_def_Q{q}.pha") for q in range(4)]
    rsps = [os.path.join(chain_def_dir, f"chain_def_Q{q}.rsp") for q in range(4)]
    coadd_pha = os.path.join(chain_def_dir, "coadd_all.pha")
    coadd_rsp = os.path.join(chain_def_dir, "coadd_all.rsp")
    addspec_res = require_ok(run_cli(bin_path, [
        "--json", "czti", "addspec",
        "--pha", *phas,
        "--rsp", *rsps,
        "-o", coadd_pha,
        "--out-rsp", coadd_rsp,
    ]), "czti addspec")
    addspec_json = json.loads(addspec_res.stdout)

    bary_evt = os.path.join(chain_def_dir, "chain_bary.evt")
    run_bary_evt = os.path.join(run_def_dir, "run_bary.evt")
    require_ok(run_cli(bin_path, [
        "--json", "bary", s9_evt, "--orbit", CZTI_ORB, "-o", bary_evt,
    ]), "czti chain bary")
    require_ok(run_cli(bin_path, [
        "--json", "bary", run_def_json["clean_evt_path"], "--orbit", CZTI_ORB, "-o", run_bary_evt,
    ]), "czti run bary")

    # Chained Generalized stages (`bunchclean --generalized` -> `pha2energy` -> `datasel` -> `clean --algo generalized` -> `evtencor` -> `evtclean`)
    g_s2_evt = os.path.join(chain_gen_dir, "g_s2_bc.evt")
    g_s2_bunch = os.path.join(chain_gen_dir, "g_s2_bc_bunch.fits")
    g_s2_lv = os.path.join(chain_gen_dir, "g_s2_bc_livetime.fits")
    require_ok(run_cli(bin_path, [
        "czti", "bunchclean", "--generalized", "-i", s1_evt, "--bunch", s1_bunch,
        "-o", g_s2_evt, "--out-bunch", g_s2_bunch, "--out-livetime", g_s2_lv,
    ]), "czti bunchclean generalized")
    g_s3_evt = os.path.join(chain_gen_dir, "g_s3_p2e.evt")
    require_ok(run_cli(bin_path, [
        "czti", "pha2energy", "-i", g_s2_evt, "--caldb", CZTI_CALDB, "--ss", CZTI_SS, "--tct", CZTI_TCT, "-o", g_s3_evt,
    ]), "czti pha2energy generalized")
    g_s6_evt = os.path.join(chain_gen_dir, "g_s6_ds.evt")
    require_ok(run_cli(bin_path, [
        "czti", "datasel", "-i", g_s3_evt, "--gti", s4_gti, "--gti-type", "quad", "-o", g_s6_evt,
    ]), "czti datasel generalized")
    g_s7_evt = os.path.join(chain_gen_dir, "g_s7_clean.evt")
    g_s7_bp = os.path.join(chain_gen_dir, "g_s7_badpix.fits")
    g_s7_lv = os.path.join(chain_gen_dir, "g_s7_livetime.fits")
    require_ok(run_cli(bin_path, [
        "czti", "clean", "-i", g_s6_evt, "--caldb", CZTI_CALDB,
        "--livetime", g_s2_lv, "--bunch", g_s2_bunch, "--algo", "generalized",
        "--gti-type", "quad", "-o", g_s7_evt, "--out-badpix", g_s7_bp, "--out-livetime", g_s7_lv,
    ]), "czti clean generalized")
    g_s8_evt = os.path.join(chain_gen_dir, "g_s8_encor.evt")
    require_ok(run_cli(bin_path, [
        "czti", "evtencor", "-i", g_s7_evt, "--caldb", CZTI_CALDB, "-o", g_s8_evt,
    ]), "czti evtencor generalized")
    g_s9_evt = os.path.join(chain_gen_dir, "g_s9_evtclean.evt")
    require_ok(run_cli(bin_path, [
        "czti", "evtclean", "-i", g_s8_evt, "--badpix", g_s7_bp, "-o", g_s9_evt,
    ]), "czti evtclean generalized")

    # Compare chain_default vs run_default
    run_evt_hdus = read_fits(run_def_json["clean_evt_path"])
    chn_evt_hdus = read_fits(s9_evt)
    quad_cmp = {}
    for q in range(4):
        qn = f"Q{q}"
        rh = next(h for h in run_evt_hdus if h.extname == qn)
        ch = next(h for h in chn_evt_hdus if h.extname == qn)
        rn, cn = int(rh.header["NAXIS2"]), int(ch.header["NAXIS2"])
        if rn != cn:
            raise AssertionError(f"CZTI default {qn} event count mismatch: run={rn}, chain={cn}")
        r_time = rh.read_column("TIME")
        c_time = ch.read_column("TIME")
        r_en = rh.read_column("ENERGY")
        c_en = ch.read_column("ENERGY")
        if not np.allclose(r_time, c_time, atol=1e-9) or not np.allclose(r_en, c_en, atol=1e-5):
            raise AssertionError(f"CZTI default {qn} TIME/ENERGY mismatch between run and chain")
        quad_cmp[qn] = rn

    # Compare bary outputs
    rb_hdus = read_fits(run_bary_evt)
    cb_hdus = read_fits(bary_evt)
    for q in range(4):
        qn = f"Q{q}"
        rh = next(h for h in rb_hdus if h.extname == qn)
        ch = next(h for h in cb_hdus if h.extname == qn)
        if not np.allclose(rh.read_column("TIME"), ch.read_column("TIME"), atol=1e-9):
            raise AssertionError(f"CZTI bary {qn} TIME mismatch between run and chain")

    # Compare DPH, BalCC sky image, PHA (including QUALITY/GROUPING), LC, and RSP between run_default and chain_default
    run_dph = read_fits(run_def_json["dph_path"])
    chn_dph = read_fits(s10_dph)
    for q in range(4):
        ra = next(h for h in run_dph if h.extname == f"Q{q}").image_array()
        ca = next(h for h in chn_dph if h.extname == f"Q{q}").image_array()
        if not np.array_equal(ra, ca):
            raise AssertionError(f"CZTI DPH Q{q} mismatch between run and chain")

    for q in range(4):
        r_img_p = next(p for p in run_def_json["image_paths"] if p.endswith(f".img_Q{q}"))
        c_img_p = next(p for p in img_balcc_json["image_paths"] if p.endswith(f".img_Q{q}"))
        r_ih = next(h for h in read_fits(r_img_p) if h.extname == f"Q{q}")
        c_ih = next(h for h in read_fits(c_img_p) if h.extname == f"Q{q}")
        if not np.allclose(r_ih.image_array(), c_ih.image_array(), atol=1e-5):
            raise AssertionError(f"CZTI BalCC sky image Q{q} mismatch between run and chain")
        cdelt1 = abs(float(c_ih.header.get("CDELT1", 0.0)))
        cdelt2 = abs(float(c_ih.header.get("CDELT2", 0.0)))
        if abs(cdelt1 - 0.10) > 1e-4 or abs(cdelt2 - 0.10) > 1e-4:
            raise AssertionError(f"CZTI BalCC Q{q} WCS plate scale expected 0.10 deg/pix, got ({cdelt1}, {cdelt2})")

    for q in range(4):
        r_pha_p = next(p for p in run_def_json["spectrum_paths"] if f"_Q{q}.pha" in p)
        c_pha_p = phas[q]
        r_ph = next(h for h in read_fits(r_pha_p) if h.extname == "SPECTRUM")
        c_ph = next(h for h in read_fits(c_pha_p) if h.extname == "SPECTRUM")
        if not np.allclose(r_ph.read_column("COUNTS"), c_ph.read_column("COUNTS"), atol=1e-6):
            raise AssertionError(f"CZTI Q{q} .pha COUNTS mismatch between run and chain")
        if not np.allclose(r_ph.read_column("STAT_ERR"), c_ph.read_column("STAT_ERR"), atol=1e-6):
            raise AssertionError(f"CZTI Q{q} .pha STAT_ERR mismatch between run and chain")
        for col_name in ("QUALITY", "GROUPING"):
            if not np.array_equal(r_ph.read_column(col_name), c_ph.read_column(col_name)):
                raise AssertionError(f"CZTI Q{q} .pha {col_name} mismatch between run and chain")

        r_lc_p = next(p for p in run_def_json["lightcurve_paths"] if f"_Q{q}.lc" in p)
        c_lc_p = os.path.join(chain_def_dir, f"chain_def_Q{q}.lc")
        r_lh = next(h for h in read_fits(r_lc_p) if h.extname in ("LIGHTCURVE", "RATE"))
        c_lh = next(h for h in read_fits(c_lc_p) if h.extname in ("LIGHTCURVE", "RATE"))
        if not np.allclose(r_lh.read_column("RATE"), c_lh.read_column("RATE"), atol=1e-6):
            raise AssertionError(f"CZTI Q{q} .lc RATE mismatch between run and chain")

        r_rsp_p = next(p for p in run_def_json["response_paths"] if f"_Q{q}.rsp" in p)
        c_rsp_p = rsps[q]
        r_rsp_hdus = read_fits(r_rsp_p)
        c_rsp_hdus = read_fits(c_rsp_p)
        r_rh = next(h for h in r_rsp_hdus if h.extname == "MATRIX")
        c_rh = next(h for h in c_rsp_hdus if h.extname == "MATRIX")
        if r_rh.data_bytes.tobytes() != c_rh.data_bytes.tobytes():
            raise AssertionError(f"CZTI Q{q} .rsp MATRIX table+heap mismatch between run and chain")
        r_eb = next(h for h in r_rsp_hdus if h.extname == "EBOUNDS")
        c_eb = next(h for h in c_rsp_hdus if h.extname == "EBOUNDS")
        if not np.allclose(r_eb.read_column("E_MIN"), c_eb.read_column("E_MIN"), atol=1e-6):
            raise AssertionError(f"CZTI Q{q} .rsp EBOUNDS mismatch between run and chain")

    # Verify explicit `--ra 0.0 --dec 0.0` is honored (must fail FOV check on Crab, not silently fall back to MKF!)
    offaxis_res = run_cli(bin_path, [
        "czti", "bindata", "-i", s9_evt, "--mkf", CZTI_MKF,
        "--caldb", CZTI_CALDB, "--badpix", s7_bp, "--livetime", s7_lv,
        "-o", os.path.join(chain_def_dir, "offaxis_check"), "--prefix", "offaxis",
        "--ra", "0.0", "--dec", "0.0",
    ])
    if offaxis_res.panicked or offaxis_res.returncode == 0:
        raise AssertionError(
            f"czti bindata --ra 0.0 --dec 0.0 must reject off-axis source (rc={offaxis_res.returncode})"
        )
    # Verify omitted `--ra`/`--dec` does NOT falsely report (0.0, 0.0) in JSON report
    if bindata_json.get("ra_deg") == 0.0 and bindata_json.get("dec_deg") == 0.0:
        raise AssertionError(f"czti bindata report reported bogus (0.0, 0.0) when --ra/--dec omitted: {bindata_json}")

    # Verify addspec coadded exposure > 0 and finite spectrum
    coadd_h = next(h for h in read_fits(coadd_pha) if h.extname == "SPECTRUM")
    coadd_rates = coadd_h.read_column("COUNTS")
    if not np.all(np.isfinite(coadd_rates)) or addspec_json["exposure_sec"] <= 0:
        raise AssertionError("CZTI addspec produced invalid spectrum or zero exposure")

    # Compare generalized chain vs run_generalized (note: run_gen runs noisypixclean before datasel; chain runs datasel before clean --algo generalized, so check close agreement & non-empty)
    gen_run_hdus = read_fits(run_gen_json["clean_evt_path"])
    gen_chn_hdus = read_fits(g_s9_evt)
    gen_cmp = {}
    for q in range(4):
        qn = f"Q{q}"
        rn = int(next(h for h in gen_run_hdus if h.extname == qn).header["NAXIS2"])
        cn = int(next(h for h in gen_chn_hdus if h.extname == qn).header["NAXIS2"])
        rel = abs(rn - cn) / max(rn, 1)
        if cn == 0 or rel > 0.05:
            raise AssertionError(f"CZTI generalized {qn} count drift too high: run={rn}, chain={cn}")
        gen_cmp[qn] = {"run": rn, "chain": cn, "rel_diff": rel}

    coadd_pha_exp = os.path.join(chain_def_dir, "coadd_sumexp.pha")
    addspec_sum_json = json.loads(require_ok(run_cli(bin_path, [
        "--json", "czti", "addspec", "--pha", phas[0], phas[1], "-o", coadd_pha_exp, "--add-exposure",
    ]), "czti addspec --add-exposure").stdout)
    assert addspec_sum_json["exposure_sec"] > addspec_json["exposure_sec"] * 1.5, (
        "czti addspec --add-exposure should sum exposures instead of averaging"
    )

    # Also verify --version works cleanly
    ver_res = require_ok(run_cli(bin_path, ["--version"]), "astrosat --version")
    assert "astrosat" in ver_res.stdout.lower()

    # Verify all generated FITS files pass NASA fitsverify 4.22 (0 errors)
    fv = verify_fitsverify([run_def_dir, run_gen_dir, chain_def_dir, chain_gen_dir])
    if fv["failed_files"] != 0:
        raise AssertionError(f"CZTI fitsverify failures: {fv['failures'][:3]}")

    return {
        "passed": True,
        "default_quad_events_exact_match": quad_cmp,
        "generalized_quad_events": gen_cmp,
        "gaas_aspect_deg": gaas_json["quad_aspect_deg"],
        "image_balcc_peaks": img_balcc_json["peak_snr"],
        "image_fft_peaks": img_fft_json["peak_snr"],
        "bindata_outputs": len(bindata_json["outputs"]),
        "addspec_exposure_sec": addspec_json["exposure_sec"],
        "fitsverify_total": fv["total_files"],
        "fitsverify_warnings": fv["total_warnings"],
    }


# ============================================================================
# Workflow 2: Chained Per-Stage LAXPC CLI Pipeline
# ============================================================================

def test_workflow2_laxpc(bin_path: str) -> Dict[str, Any]:
    wf_dir = os.path.join(WORK_ROOT, "wf2_laxpc")
    if os.path.exists(wf_dir):
        shutil.rmtree(wf_dir)
    caldb_dir = make_temp_laxpc_caldb(os.path.join(wf_dir, "caldb"))
    run_dir = os.path.join(wf_dir, "run_lxp1")
    chain_dir = os.path.join(wf_dir, "chain_lxp1")
    os.makedirs(run_dir, exist_ok=True)
    os.makedirs(chain_dir, exist_ok=True)

    # Baseline `astrosat laxpc run` (lxp1, format-a, caldb)
    r_run = require_ok(run_cli(bin_path, [
        "--json", "laxpc", "run",
        "--input", LAXPC_L1_DIR,
        "--caldb", caldb_dir,
        "--units", "lxp1",
        "--format-a",
        "--output", run_dir,
    ]), "laxpc run lxp1")
    run_json = json.loads(r_run.stdout)
    u0 = run_json["unit_results"][0]

    # 1. gtigen
    gti_path = os.path.join(chain_dir, "lxp1.gti")
    gti_res = require_ok(run_cli(bin_path, [
        "--json", "laxpc", "gtigen",
        "--ea", LAXPC_EA1, "--bb", LAXPC_BB1, "--mkf", LAXPC_MKF, "--tct", LAXPC_TCT,
        "--unit", "lxp1", "-o", gti_path,
    ]), "laxpc gtigen")
    gti_json = json.loads(gti_res.stdout)

    # 2. decode
    evn_path = os.path.join(chain_dir, "lxp1.evn")
    dec_res = require_ok(run_cli(bin_path, [
        "--json", "laxpc", "decode",
        "--ea", LAXPC_EA1, "--bb", LAXPC_BB1, "--tct", LAXPC_TCT, "--gti", gti_path,
        "--unit", "lxp1", "-o", evn_path,
    ]), "laxpc decode")
    dec_json = json.loads(dec_res.stdout)

    # 3. orbit
    orb_path = os.path.join(chain_dir, "AstroSat.orb")
    orb_res = require_ok(run_cli(bin_path, [
        "--json", "laxpc", "orbit", "--mkf", LAXPC_MKF, "-o", orb_path,
    ]), "laxpc orbit")
    orb_json = json.loads(orb_res.stdout)

    # 4a. lightcurve (default settings matching `laxpc run`)
    lc_def_path = os.path.join(chain_dir, "lxp1_default.lcurv")
    lc_def_res = require_ok(run_cli(bin_path, [
        "--json", "laxpc", "lightcurve",
        "--ea", LAXPC_EA1, "--bb", LAXPC_BB1, "--mkf", LAXPC_MKF, "--tct", LAXPC_TCT,
        "--unit", "lxp1", "-o", lc_def_path,
    ]), "laxpc lightcurve default")
    lc_def_json = json.loads(lc_def_res.stdout)

    # 4b. lightcurve with `--energy-range 3:20 --tbin 0.01 --layer1-only --single-only`
    lc_filt_path = os.path.join(chain_dir, "lxp1_3_20_l1_single.lcurv")
    lc_filt_res = require_ok(run_cli(bin_path, [
        "--json", "laxpc", "lightcurve",
        "--ea", LAXPC_EA1, "--bb", LAXPC_BB1, "--mkf", LAXPC_MKF, "--tct", LAXPC_TCT,
        "--caldb", caldb_dir,
        "--unit", "lxp1",
        "--energy-range", "3:20",
        "--tbin", "0.01",
        "--layer1-only",
        "--single-only",
        "-o", lc_filt_path,
    ]), "laxpc lightcurve filtered")
    lc_filt_json = json.loads(lc_filt_res.stdout)

    # 5a. spectrum (default broadband)
    spec_path = os.path.join(chain_dir, "lxp1.spec")
    spec_res = require_ok(run_cli(bin_path, [
        "--json", "laxpc", "spectrum",
        "--ea", LAXPC_EA1, "--bb", LAXPC_BB1, "--mkf", LAXPC_MKF, "--tct", LAXPC_TCT,
        "--unit", "lxp1", "-o", spec_path,
    ]), "laxpc spectrum")
    spec_json = json.loads(spec_res.stdout)

    # 5b. spectrum with `--energy-range 3:20 --caldb <caldb_dir> --layer1-only --single-only`
    spec_filt_path = os.path.join(chain_dir, "lxp1_3_20_l1_single.spec")
    spec_filt_res = require_ok(run_cli(bin_path, [
        "--json", "laxpc", "spectrum",
        "--ea", LAXPC_EA1, "--bb", LAXPC_BB1, "--mkf", LAXPC_MKF, "--tct", LAXPC_TCT,
        "--caldb", caldb_dir, "--unit", "lxp1",
        "--energy-range", "3:20", "--layer1-only", "--single-only",
        "-o", spec_filt_path,
    ]), "laxpc spectrum filtered")
    spec_filt_json = json.loads(spec_filt_res.stdout)
    sf_rates = next(h for h in read_fits(spec_filt_path) if h.extname == "SPECTRUM").read_column("RATE")
    assert 0.0 < float(np.sum(sf_rates)) < float(spec_json["total_corrected_rate_cps"]), (
        "Filtered 3-20 keV spectrum rate should be strictly between 0 and broadband rate"
    )
    assert float(sf_rates[0]) == 0.0 and float(sf_rates[-1]) == 0.0, (
        "Channels outside 3-20 keV must be zeroed in filtered laxpc spectrum"
    )

    # 6. backshift (standard + HV-jump catalog + explicit `--bkg-peak-ch` override)
    bkg_shift_path = os.path.join(chain_dir, "lxp1_bkg_shifted.spec")
    bsh_res = require_ok(run_cli(bin_path, [
        "--json", "laxpc", "backshift",
        "--unit", "lxp1",
        "--caldb", caldb_dir,
        "--fits", spec_path,
        "--bkg-pha", os.path.join(caldb_dir, "backlxp10jan16.pha"),
        "-o", bkg_shift_path,
    ]), "laxpc backshift")
    bsh_json = json.loads(bsh_res.stdout)

    # 6b. Post-HV-jump LXP3 background (`backlxp30aug17.pha` -> catalog peak=213.0) + explicit `--bkg-peak-ch 220.0`
    bkg_hv_path = os.path.join(chain_dir, "lxp3_hv_shifted.spec")
    bsh_hv = json.loads(require_ok(run_cli(bin_path, [
        "--json", "laxpc", "backshift",
        "--unit", "lxp3",
        "--caldb", caldb_dir,
        "--mission-day", "700",
        "--bkg-pha", os.path.join(caldb_dir, "backlxp30aug17.pha"),
        "-o", bkg_hv_path,
    ]), "laxpc backshift lxp3 HV catalog").stdout)
    bsh_ovr = json.loads(require_ok(run_cli(bin_path, [
        "--json", "laxpc", "backshift",
        "--unit", "lxp3",
        "--caldb", caldb_dir,
        "--mission-day", "700",
        "--bkg-pha", os.path.join(caldb_dir, "backlxp30aug17.pha"),
        "--bkg-peak-ch", "220.0",
    ]), "laxpc backshift --bkg-peak-ch override").stdout)
    hv_dch = float(bsh_hv["gain_shift"]["delta_ch_30kev"])
    ovr_dch = float(bsh_ovr["gain_shift"]["delta_ch_30kev"])
    assert abs(hv_dch - 29.8) < 1.0, f"Expected LXP3 HV-jump catalog delta_ch ~ +29.8, got {hv_dch}"
    assert abs((hv_dch - ovr_dch) - 7.0) < 1e-4, f"Expected --bkg-peak-ch 220 vs 213 to shift delta_ch by -7.0, got {hv_dch} vs {ovr_dch}"

    # 6c. Post-day-917 LXP10 shifted `.channel` map selection (`backlxp10apr18.pha` at day 950 -> fshbk ~ 0.550028)
    bsh_cshm = json.loads(require_ok(run_cli(bin_path, [
        "--json", "laxpc", "backshift",
        "--unit", "lxp1",
        "--caldb", caldb_dir,
        "--mission-day", "950",
        "--bkg-pha", os.path.join(caldb_dir, "backlxp10apr18.pha"),
    ]), "laxpc backshift lxp10 day 950 cshm141").stdout)
    cshm_fsh = float(bsh_cshm["lcbk_shift_fraction"])
    assert abs(cshm_fsh - 0.550028) < 1e-4, (
        f"Expected LXP10 day 950 cshm141 fshbk ~ 0.550028, got {cshm_fsh}"
    )

    # 7. bary
    run_bary_evn = os.path.join(run_dir, "lxp1_bary.evn")
    chn_bary_evn = os.path.join(chain_dir, "lxp1_bary.evn")
    require_ok(run_cli(bin_path, [
        "--json", "bary", u0["event_fits_path"], "--orbit", run_json["orbit_fits_path"], "-o", run_bary_evn,
    ]), "laxpc run bary")
    require_ok(run_cli(bin_path, [
        "--json", "bary", evn_path, "--orbit", orb_path, "-o", chn_bary_evn,
    ]), "laxpc chain bary")

    # Exact parity checks between chained stages and `laxpc run`
    if abs(gti_json["net_exposure_sec"] - u0["gti"]["net_exposure_sec"]) > 1e-6:
        raise AssertionError("LAXPC gtigen net exposure mismatch vs run")
    if abs(lc_def_json["mean_rate_cps"] - u0["lightcurve_mean_rate_cps"]) > 1e-6:
        raise AssertionError("LAXPC lightcurve mean rate mismatch vs run")
    if abs(spec_json["total_corrected_rate_cps"] - u0["spectrum"]["total_corrected_rate_cps"]) > 1e-6:
        raise AssertionError("LAXPC spectrum total rate mismatch vs run")

    r_spec = next(h for h in read_fits(u0["spectrum_fits_path"]) if h.extname == "SPECTRUM").read_column("RATE")
    c_spec = next(h for h in read_fits(spec_path) if h.extname == "SPECTRUM").read_column("RATE")
    if not np.allclose(r_spec, c_spec, atol=1e-9):
        raise AssertionError("LAXPC spectrum channel rates mismatch between run and chain")

    r_evn_n = int(next(h for h in read_fits(u0["event_fits_path"]) if h.xtension == "BINTABLE" and h.extname != "GTI").header["NAXIS2"])
    c_evn_n = int(next(h for h in read_fits(evn_path) if h.xtension == "BINTABLE" and h.extname != "GTI").header["NAXIS2"])
    if r_evn_n != c_evn_n:
        raise AssertionError(f"LAXPC .evn event count mismatch: run={r_evn_n}, chain={c_evn_n}")

    r_bary_t = next(h for h in read_fits(run_bary_evn) if h.xtension == "BINTABLE" and h.extname != "GTI").read_column("TIME")
    c_bary_t = next(h for h in read_fits(chn_bary_evn) if h.xtension == "BINTABLE" and h.extname != "GTI").read_column("TIME")
    if not np.allclose(r_bary_t, c_bary_t, atol=1e-9):
        raise AssertionError("LAXPC bary TIME mismatch between run and chain")

    # Filtered lightcurve check: tbin == 0.01, rate > 0 and < broadband rate
    if abs(lc_filt_json["tbin"] - 0.01) > 1e-9:
        raise AssertionError(f"Expected tbin=0.01, got {lc_filt_json['tbin']}")
    if not (0.0 < lc_filt_json["mean_rate_cps"] < lc_def_json["mean_rate_cps"]):
        raise AssertionError(
            f"Filtered 3-20 keV L1 single rate {lc_filt_json['mean_rate_cps']} not strictly < broadband {lc_def_json['mean_rate_cps']}"
        )

    # 7b. Barycenter Format A `level2.event.fits` and verify `astrosat timing` with `--unit --layer --pi-min --pi-max`
    fmt_a_bary = os.path.join(run_dir, "level2.event.bary.fits")
    require_ok(run_cli(bin_path, [
        "--json", "bary", run_json["format_a_event_fits_path"], "--orbit", run_json["orbit_fits_path"], "-o", fmt_a_bary,
    ]), "laxpc format A bary")
    fmt_a_tim = json.loads(require_ok(run_cli(bin_path, [
        "--json", "timing", fmt_a_bary, "--period", "0.0337176", "--search-span", "0.005",
        "--fdot", "-3.693e-10", "--unit", "1", "--layer", "1", "--pi-min", "10", "--pi-max", "300",
    ]), "laxpc format A timing").stdout)
    assert 29.645 <= fmt_a_tim["best_freq_hz"] <= 29.670, f"Unexpected Crab freq on Format A: {fmt_a_tim['best_freq_hz']}"
    assert fmt_a_tim["max_z2_n"] > 100.0, f"Expected strong Crab Z2 on Format A Layer 1, got {fmt_a_tim['max_z2_n']}"

    fv = verify_fitsverify([run_dir, chain_dir])
    if fv["failed_files"] != 0:
        raise AssertionError(f"LAXPC fitsverify failures: {fv['failures'][:3]}")

    return {
        "passed": True,
        "events_exact_match": c_evn_n,
        "exposure_sec": gti_json["net_exposure_sec"],
        "lc_default_mean_cps": lc_def_json["mean_rate_cps"],
        "lc_filtered_3_20_l1_single_0p01s_cps": lc_filt_json["mean_rate_cps"],
        "spec_total_rate_cps": spec_json["total_corrected_rate_cps"],
        "spec_filtered_3_20_rate_cps": float(np.sum(sf_rates)),
        "backshift_delta_ch": bsh_json["gain_shift"]["delta_ch_30kev"],
        "backshift_rmf": bsh_json["gain_shift"]["recommended_rmf"],
        "orbit_tstart": orb_json["tstart"],
        "orbit_tstop": orb_json["tstop"],
        "fitsverify_total": fv["total_files"],
        "fitsverify_warnings": fv["total_warnings"],
    }


# ============================================================================
# Workflow 3: Declarative Config (TOML + JSON), Overrides & `--json` Schemas
# ============================================================================

def test_workflow3_configs_and_json(bin_path: str) -> Dict[str, Any]:
    wf_dir = os.path.join(WORK_ROOT, "wf3_config")
    if os.path.exists(wf_dir):
        shutil.rmtree(wf_dir)
    os.makedirs(wf_dir, exist_ok=True)
    caldb_laxpc = make_temp_laxpc_caldb(os.path.join(wf_dir, "caldb_laxpc"))

    # --- 3A. CZTI TOML & JSON + Full CLI Overrides (including default values) ---
    czti_toml = os.path.join(wf_dir, "czti_cfg.toml")
    czti_json = os.path.join(wf_dir, "czti_cfg.json")
    out_czti_toml = os.path.join(wf_dir, "out_czti_toml")
    out_czti_json = os.path.join(wf_dir, "out_czti_json")
    out_czti_ovr = os.path.join(wf_dir, "out_czti_ovr")

    # Verify `generate_eventfile = true` + `mask_weight = true` in TOML emits weighted evt, and `output_type = "spec"` + `gen_resp = false` omits .lc/.rsp
    with open(czti_toml, "w", encoding="utf-8") as fh:
        fh.write(f"""
threads = 2

[czti]
input_dir = "{CZTI_L1_DIR}"
l1_m0_file = "{CZTI_M0}"
l1_ss_file = "{CZTI_SS}"
tct_file = "{CZTI_TCT}"
mkf_file = "{CZTI_MKF}"
caldb_dir = "{CZTI_CALDB}"
output_dir = "{out_czti_toml}"
output_prefix = "toml_czti"
mkf_threshold_file = "{CZTI_THR}"
entry_stage = 1
exit_stage = 3
gti_type = "common"
algo = "generalized"
ra_deg = 83.6331
dec_deg = 22.0145
time_binsize = 5.0
emin_kev = 25.0
emax_kev = 90.0
output_type = "spec"
mask_weight = true
gen_resp = false
generate_eventfile = true
image_algo = "fft"
oversampling = 8
quads = [0, 1]
""")

    czti_cfg_dict = {
        "threads": 2,
        "czti": {
            "input_dir": CZTI_L1_DIR,
            "l1_m0_file": CZTI_M0,
            "l1_ss_file": CZTI_SS,
            "tct_file": CZTI_TCT,
            "mkf_file": CZTI_MKF,
            "caldb_dir": CZTI_CALDB,
            "output_dir": out_czti_json,
            "output_prefix": "json_czti",
            "mkf_threshold_file": CZTI_THR,
            "entry_stage": 1,
            "exit_stage": 3,
            "gti_type": "common",
            "algo": "generalized",
            "ra_deg": 83.6331,
            "dec_deg": 22.0145,
            "time_binsize": 5.0,
            "emin_kev": 25.0,
            "emax_kev": 90.0,
            "output_type": "spec",
            "mask_weight": True,
            "gen_resp": False,
            "generate_eventfile": True,
            "image_algo": "fft",
            "oversampling": 8,
            "quads": [0, 1],
        },
    }
    with open(czti_json, "w", encoding="utf-8") as fh:
        json.dump(czti_cfg_dict, fh, indent=2)

    r_ct = json.loads(require_ok(run_cli(bin_path, ["--config", czti_toml, "--json", "czti", "run"]), "czti config toml").stdout)
    r_cj = json.loads(require_ok(run_cli(bin_path, ["--config", czti_json, "--json", "czti", "run"]), "czti config json").stdout)
    assert r_ct["stage2_clean_events"] == r_cj["stage2_clean_events"], "TOML vs JSON CZTI mismatch"
    assert len(r_ct["lightcurve_paths"]) == 0, "config output_type='spec' should omit .lc"
    assert len(r_ct["response_paths"]) == 0, "config gen_resp=false should omit .rsp"
    assert r_ct["weighted_evt_path"] is not None, "config generate_eventfile=true should emit weighted evt"

    # Override EVERY CZTI config field with CLI defaults/values (entry_stage=1, exit_stage=3, output_type=both, no-gen-evt, mask-weight, resp, gti-type=quad, algo=default, image-algo=balcc)
    czti_toml_ovr_cfg = os.path.join(wf_dir, "czti_ovr.toml")
    with open(czti_toml_ovr_cfg, "w", encoding="utf-8") as fh:
        fh.write(f"""
[czti]
input_dir = "{CZTI_L1_DIR}"
caldb_dir = "{CZTI_CALDB}"
output_dir = "{out_czti_toml}"
output_prefix = "should_be_overridden"
entry_stage = 3
exit_stage = 1
gti_type = "common"
algo = "generalized"
output_type = "lc"
mask_weight = false
gen_resp = false
generate_eventfile = true
image_algo = "fft"
oversampling = 4
quads = [0]
""")
    r_co = json.loads(require_ok(run_cli(bin_path, [
        "--config", czti_toml_ovr_cfg, "--json", "czti", "run",
        "--output", out_czti_ovr,
        "--prefix", "ovr_czti",
        "--entry-stage", "1",
        "--exit-stage", "3",
        "--gti-type", "quad",
        "--algo", "default",
        "--output-type", "both",
        "--mask-weight",
        "--resp",
        "--no-gen-evt",
        "--image-algo", "balcc",
        "--oversampling", "16",
        "--quads", "0,1",
        "--tbin", "2.0",
        "--emin", "20.0",
        "--emax", "100.0",
        "--ra", "83.6331",
        "--dec", "22.0145",
    ]), "czti config CLI override").stdout)
    assert len(r_co["spectrum_paths"]) == 2, f"Expected 2 quad spectra after CLI override, got {r_co['spectrum_paths']}"
    assert len(r_co["lightcurve_paths"]) == 2, "CLI --output-type both must override config 'lc'"
    assert len(r_co["response_paths"]) == 2, "CLI --resp must override config gen_resp=false"
    assert r_co["weighted_evt_path"] is None, "CLI --no-gen-evt must override config generate_eventfile=true"

    # --- 3B. LAXPC TOML & JSON + CLI Default Overrides ---
    laxpc_toml = os.path.join(wf_dir, "laxpc_cfg.toml")
    laxpc_json = os.path.join(wf_dir, "laxpc_cfg.json")
    out_lxp_toml = os.path.join(wf_dir, "out_lxp_toml")
    out_lxp_json = os.path.join(wf_dir, "out_lxp_json")
    out_lxp_ovr = os.path.join(wf_dir, "out_lxp_ovr")

    with open(laxpc_toml, "w", encoding="utf-8") as fh:
        fh.write(f"""
threads = 2

[laxpc]
input_dir = "{LAXPC_L1_DIR}"
output_dir = "{out_lxp_toml}"
caldb_dir = "{caldb_laxpc}"
units = ["lxp1"]
tbin = 2.0
deadtime_sec = 0.000042
spectrum_channels = 256
write_format_a_combined = true
write_format_b_events = true
write_orbit_file = true
write_ascii_products = true
single_only = true
layer1_only = true
min_elv_deg = 5.0
saa_entry_lon_deg = 255.0
saa_exit_lon_deg = 5.0
merge_gap_sec = 1.0
min_duration_sec = 10.0
user_gti_only = false
""")

    laxpc_cfg_dict = {
        "threads": 2,
        "laxpc": {
            "input_dir": LAXPC_L1_DIR,
            "output_dir": out_lxp_json,
            "caldb_dir": caldb_laxpc,
            "units": ["lxp1"],
            "tbin": 2.0,
            "deadtime_sec": 0.000042,
            "spectrum_channels": 256,
            "write_format_a_combined": True,
            "write_format_b_events": True,
            "write_orbit_file": True,
            "write_ascii_products": True,
            "single_only": True,
            "layer1_only": True,
            "min_elv_deg": 5.0,
            "saa_entry_lon_deg": 255.0,
            "saa_exit_lon_deg": 5.0,
            "merge_gap_sec": 1.0,
            "min_duration_sec": 10.0,
            "user_gti_only": False,
        },
    }
    with open(laxpc_json, "w", encoding="utf-8") as fh:
        json.dump(laxpc_cfg_dict, fh, indent=2)

    r_lt = json.loads(require_ok(run_cli(bin_path, ["--config", laxpc_toml, "--json", "laxpc", "run"]), "laxpc config toml").stdout)
    r_lj = json.loads(require_ok(run_cli(bin_path, ["--config", laxpc_json, "--json", "laxpc", "run"]), "laxpc config json").stdout)
    assert r_lt["unit_results"][0]["spectrum"]["num_channels"] == 256
    assert r_lt["unit_results"][0]["gti_filtered_events"] == r_lj["unit_results"][0]["gti_filtered_events"]
    assert r_lt["format_a_event_fits_path"] is not None
    assert r_lt["orbit_fits_path"] is not None
    assert os.path.isfile(os.path.join(out_lxp_toml, "lxp1level2.pha"))

    # CLI overrides with defaults (`--no-format-a`, `--no-orbit`, `--no-ascii`, `--no-single-only`, `--no-layer1-only`, `--no-user-gti-only`, `--tbin 1.0`, `--channels 512`, `--deadtime 0.000043`)
    r_lo = json.loads(require_ok(run_cli(bin_path, [
        "--config", laxpc_toml, "--json", "laxpc", "run",
        "--output", out_lxp_ovr,
        "--no-format-a",
        "--format-b",
        "--no-orbit",
        "--no-ascii",
        "--no-single-only",
        "--no-layer1-only",
        "--no-user-gti-only",
        "--tbin", "1.0",
        "--channels", "512",
        "--deadtime", "0.000043",
        "--elv-min", "3.5",
        "--saa-entry-lon", "249.0",
        "--saa-exit-lon", "3.0",
        "--merge-gap", "0.3",
        "--min-duration", "5.0",
    ]), "laxpc config CLI override").stdout)
    assert r_lo["format_a_event_fits_path"] is None, "--no-format-a failed to override config"
    assert r_lo["orbit_fits_path"] is None, "--no-orbit failed to override config"
    assert not os.path.exists(os.path.join(out_lxp_ovr, "lxp1level2.pha")), "--no-ascii failed to override config"
    assert r_lo["unit_results"][0]["spectrum"]["num_channels"] == 512
    assert r_lo["unit_results"][0]["gti_filtered_events"] > r_lt["unit_results"][0]["gti_filtered_events"], (
        "--no-single-only --no-layer1-only should retain more events than single_only+layer1_only"
    )

    # --- 3C. Bary TOML & JSON + CLI Default Overrides (`--frame icrs`, `--ephem analytical`, `--checksum`) ---
    bary_toml = os.path.join(wf_dir, "bary_cfg.toml")
    bary_json = os.path.join(wf_dir, "bary_cfg.json")
    bary_in = r_lt["unit_results"][0]["event_fits_path"]
    orb_in = r_lt["orbit_fits_path"]
    bary_out_t = os.path.join(wf_dir, "bary_toml.evn")
    bary_out_j = os.path.join(wf_dir, "bary_json.evn")
    bary_out_o = os.path.join(wf_dir, "bary_ovr.evn")

    with open(bary_toml, "w", encoding="utf-8") as fh:
        fh.write(f"""
[bary]
orbit_file = "{orb_in}"
ra_deg = 83.6331
dec_deg = 22.0145
frame = "fk5"
ephemeris = "/nonexistent/de405.bin"
write_checksums = false
""")
    with open(bary_json, "w", encoding="utf-8") as fh:
        json.dump({
            "bary": {
                "orbit_file": orb_in,
                "ra_deg": 83.6331,
                "dec_deg": 22.0145,
                "frame": "fk5",
                "ephemeris": JPLEPH_405,
                "write_checksums": True,
            }
        }, fh, indent=2)

    # JSON config with JPL DE405 & FK5
    r_bj = json.loads(require_ok(run_cli(bin_path, [
        "--config", bary_json, "--json", "bary", bary_in, "-o", bary_out_j,
    ]), "bary config json").stdout)
    assert r_bj["frame"] == "FK5"
    assert r_bj["plephem"] == "JPL-DE405"

    # TOML config with broken ephemeris & write_checksums=false overridden by CLI defaults `--frame icrs --ephem analytical --checksum`
    r_bo = json.loads(require_ok(run_cli(bin_path, [
        "--config", bary_toml, "--json", "bary", bary_in, "-o", bary_out_o,
        "--frame", "icrs", "--ephem", "analytical", "--checksum",
    ]), "bary config CLI default override").stdout)
    assert r_bo["frame"] == "ICRS", f"Expected ICRS override, got {r_bo['frame']}"
    assert r_bo["plephem"] == "ASTROSAT-ANALYTICAL", f"Expected ASTROSAT-ANALYTICAL override, got {r_bo['plephem']}"
    bo_hdus = read_fits(bary_out_o)
    assert "CHECKSUM" in bo_hdus[1].header, "--checksum override failed to write CHECKSUM"

    # Also test TOML config with valid JPL DE405
    with open(bary_toml, "w", encoding="utf-8") as fh:
        fh.write(f"""
[bary]
orbit_file = "{orb_in}"
ra_deg = 83.6331
dec_deg = 22.0145
frame = "icrs"
ephemeris = "{JPLEPH_405}"
write_checksums = true
""")
    r_bt = json.loads(require_ok(run_cli(bin_path, [
        "--config", bary_toml, "--json", "bary", bary_in, "-o", bary_out_t,
    ]), "bary config toml").stdout)
    assert r_bt["frame"] == "ICRS"
    assert r_bt["plephem"] == "JPL-DE405"

    # Also verify HEASARC FITS `JPLEPH.405` matches binary DE405 to < 1e-9 s, HEASARC FITS `JPLEPH.200` + FK5 works,
    # and explicit `--mkf` outranks sibling `AstroSat.orb` & `cfg.orbit_file`!
    bary_out_fits405 = os.path.join(wf_dir, "bary_fits405.evn")
    bary_out_fits200 = os.path.join(wf_dir, "bary_fits200.evn")
    bary_out_mkf = os.path.join(wf_dir, "bary_mkf.evn")
    r_f405 = json.loads(require_ok(run_cli(bin_path, [
        "--json", "bary", bary_in, "--orbit", orb_in, "--ephem", JPLEPH_FITS_405, "-o", bary_out_fits405,
    ]), "bary HEASARC FITS JPLEPH.405").stdout)
    assert r_f405["plephem"] == "JPL-DE405"
    t_bin405 = next(h for h in read_fits(bary_out_t) if h.xtension == "BINTABLE" and h.extname != "GTI").read_column("TIME")
    t_fits405 = next(h for h in read_fits(bary_out_fits405) if h.xtension == "BINTABLE" and h.extname != "GTI").read_column("TIME")
    assert np.allclose(t_bin405, t_fits405, atol=1e-9), "HEASARC FITS JPLEPH.405 diverged from binary DE405!"

    r_f200 = json.loads(require_ok(run_cli(bin_path, [
        "--json", "bary", bary_in, "--orbit", orb_in, "--frame", "fk5", "--ephem", JPLEPH_FITS_200, "-o", bary_out_fits200,
    ]), "bary HEASARC FITS JPLEPH.200").stdout)
    assert r_f200["frame"] == "FK5" and r_f200["plephem"] == "JPL-DE200"

    r_mkf = json.loads(require_ok(run_cli(bin_path, [
        "--config", bary_toml, "--json", "bary", bary_in, "--mkf", LAXPC_MKF, "-o", bary_out_mkf,
    ]), "bary explicit --mkf precedence").stdout)
    assert r_mkf["orbit"] == LAXPC_MKF, f"Explicit --mkf was ignored in favor of {r_mkf['orbit']}!"
    t_mkf = next(h for h in read_fits(bary_out_mkf) if h.xtension == "BINTABLE" and h.extname != "GTI").read_column("TIME")
    assert np.allclose(t_bin405, t_mkf, atol=1e-6), "bary --mkf synthesized orbit diverged from --orbit!"

    return {
        "passed": True,
        "czti_toml_json_match": True,
        "czti_cli_overrides_ok": True,
        "laxpc_toml_json_match": True,
        "laxpc_cli_overrides_ok": True,
        "bary_toml_json_overrides_ok": True,
        "bary_fits_jpleph_405_200_and_mkf_ok": True,
    }


# ============================================================================
# Workflow 4: Idempotency & Guardrails
# ============================================================================

def test_workflow4_idempotency_and_guardrails(bin_path: str) -> Dict[str, Any]:
    wf_dir = os.path.join(WORK_ROOT, "wf4_guardrails")
    if os.path.exists(wf_dir):
        shutil.rmtree(wf_dir)
    os.makedirs(wf_dir, exist_ok=True)

    # 1. Re-running `czti evtencor` on GAINCOR='YES' `.evt` must be a clean no-op
    clean_evt = glob.glob(os.path.join(WORK_ROOT, "wf1_czti", "run_default", "*_clean.evt"))[0]
    re_encor = os.path.join(wf_dir, "re_encor.evt")
    res_ec = require_ok(run_cli(bin_path, [
        "--json", "czti", "evtencor", "-i", clean_evt, "--caldb", CZTI_CALDB, "-o", re_encor,
    ]), "czti evtencor idempotency")
    ec_json = json.loads(res_ec.stdout)
    orig_hdus = read_fits(clean_evt)
    new_hdus = read_fits(re_encor)
    for q in range(4):
        qn = f"Q{q}"
        oh = next(h for h in orig_hdus if h.extname == qn)
        nh = next(h for h in new_hdus if h.extname == qn)
        assert str(nh.header.get("GAINCOR", "")).strip().upper() == "YES"
        assert int(oh.header["NAXIS2"]) == int(nh.header["NAXIS2"]), f"{qn} count changed on re-evtencor"
        assert np.array_equal(oh.read_column("PI"), nh.read_column("PI")), f"{qn} PI shifted on re-evtencor!"
        assert np.allclose(oh.read_column("ENERGY"), nh.read_column("ENERGY"), atol=1e-6), f"{qn} ENERGY shifted!"

    # 2. Re-running `astrosat bary` on an already barycentered file (`TIMESYS='TDB'`) must fail cleanly
    bary_in = os.path.join(WORK_ROOT, "wf1_czti", "run_default", "run_bary.evt")
    bary_twice = os.path.join(wf_dir, "bary_twice.evt")
    res_bt = run_cli(bin_path, ["bary", bary_in, "--orbit", CZTI_ORB, "-o", bary_twice])
    assert not res_bt.panicked and res_bt.returncode != 0, "Double barycenter should fail cleanly"
    assert not os.path.exists(bary_twice), "Failed double barycenter must not create output"
    if "TDB" not in res_bt.stderr and "already" not in res_bt.stderr.lower():
        raise AssertionError(f"Double barycenter error unclear: {res_bt.stderr}")

    # 3. `laxpc run --energy-range 3:20` and `3-20` without `--caldb` must fail cleanly
    for er_arg in ("3-20", "3:20"):
        res_er = run_cli(bin_path, [
            "laxpc", "run", "--input", LAXPC_L1_DIR, "--units", "lxp1",
            "--energy-range", er_arg, "--output", os.path.join(wf_dir, f"er_{er_arg}"),
        ])
        assert not res_er.panicked and res_er.returncode != 0, f"laxpc run --energy-range {er_arg} without --caldb should fail!"
        err_low = res_er.stderr.lower()
        if ("caldb" not in err_low) and ("emin:emax" not in err_low) and ("energy-range" not in err_low):
            raise AssertionError(f"Unhelpful error for --energy-range {er_arg}: {res_er.stderr}")

    # 4. Also test `laxpc lightcurve --energy-range 3:20` without `--caldb`
    res_lcer = run_cli(bin_path, [
        "laxpc", "lightcurve",
        "--ea", LAXPC_EA1, "--mkf", LAXPC_MKF, "--tct", LAXPC_TCT,
        "--unit", "lxp1", "--energy-range", "3:20", "-o", os.path.join(wf_dir, "lc_nocaldb.lcurv"),
    ])
    assert not res_lcer.panicked and res_lcer.returncode != 0, "laxpc lightcurve --energy-range without --caldb should fail!"

    # 5. Check `laxpc spectrum --energy-range 3:20` (probe whether CLI silently ignores it!)
    spec_er_path = os.path.join(wf_dir, "spec_er.spec")
    res_sper = run_cli(bin_path, [
        "laxpc", "spectrum",
        "--ea", LAXPC_EA1, "--mkf", LAXPC_MKF, "--tct", LAXPC_TCT,
        "--unit", "lxp1", "--energy-range", "3:20", "-o", spec_er_path,
    ])
    spec_er_guard_ok = (not res_sper.panicked) and (res_sper.returncode != 0)

    return {
        "passed": spec_er_guard_ok,
        "evtencor_idempotent_quads": ec_json["quad_events"],
        "double_bary_rejected_cleanly": True,
        "laxpc_run_energy_range_no_caldb_rejected": True,
        "laxpc_lightcurve_energy_range_no_caldb_rejected": True,
        "laxpc_spectrum_energy_range_rejected_or_bug": {
            "guard_ok": spec_er_guard_ok,
            "returncode": res_sper.returncode,
            "note": "Bug if rc==0: `LaxpcSpectrumArgs` defines --energy-range/--pha-min/--pha-max without --caldb and `run_laxpc_spectrum` silently ignores them" if not spec_er_guard_ok else "OK",
        },
    }


# ============================================================================
# Synthetic LAXPC EA Frame Encoder (for Ground-Truth Physics Scenarios)
# ============================================================================

STBG_EPOCH_OFFSET_SEC = 1.8121823219164e8
STBG_CLOCK_HZ = 100_000.0


def encode_ea_word(fine_ticks: int, ia1: int, ia2: int, ic1: int, ic2: int = 0, k1: int = 0, k2: int = 0) -> bytes:
    b0 = (ia1 & 0x0F) | ((ia2 & 0x0F) << 4)
    b1 = (fine_ticks + 2) & 0xFF
    b2 = (ic1 >> 2) & 0xFF
    b3 = (ic2 >> 2) & 0xFF
    l1 = ic1 & 0x03
    l2 = ic2 & 0x03
    b4 = (k1 & 0x01) | (l1 << 2) | ((k2 & 0x01) << 4) | (l2 << 6)
    return bytes((b0, b1, b2, b3, b4))


def build_synthetic_ea_frames(
    t_ins_events: np.ndarray,
    phas: np.ndarray,
    anodes: np.ndarray,
    unit_id: int = 1,
    t0_offset_sec: float = STBG_EPOCH_OFFSET_SEC,
    max_words_per_frame: int = 380,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Pack strictly increasing `t_ins_events` (in seconds) into valid 2242-byte LAXPC EA frames.

    Returns (l1_times, grt_arr, seq_arr, data_arrays_2242).
    """
    ticks = np.round(t_ins_events * STBG_CLOCK_HZ).astype(np.int64)
    # Enforce strictly increasing ticks >= 1 tick apart
    for i in range(1, len(ticks)):
        if ticks[i] <= ticks[i - 1]:
            ticks[i] = ticks[i - 1] + 1

    # Split into frames so each frame holds <= max_words_per_frame 5-byte words.
    # Real LAXPC hardware emits a 0xEF coarse-time marker for EVERY 256-tick (2.56 ms) bin,
    # even when no events occur in that bin, so m2 increments strictly by +1 mod 256.
    frames_words: List[Tuple[int, List[bytes]]] = []
    cur_words: List[bytes] = []
    first_coarse = int(ticks[0]) >> 8
    last_coarse = int(ticks[-1]) >> 8
    frame_start_tick = int(ticks[0]) & ~0xFF

    ev_idx = 0
    n_ev = len(ticks)
    prev_ev_tk = int(ticks[0])

    for coarse in range(first_coarse, last_coarse + 1):
        m4 = (coarse >> 16) & 0xFF
        m3 = (coarse >> 8) & 0xFF
        m2 = coarse & 0xFF
        ef_word = bytes((0xEF, m4, m3, m2, 0))

        # Collect events belonging to this 256-tick bin
        bin_events: List[Tuple[int, int, int]] = []
        while ev_idx < n_ev and (int(ticks[ev_idx]) >> 8) == coarse:
            bin_events.append((int(ticks[ev_idx]), int(phas[ev_idx]), int(anodes[ev_idx])))
            ev_idx += 1

        need = 1 + len(bin_events)
        if cur_words and len(cur_words) + need >= max_words_per_frame - 2:
            # Append 1 dummy word right before 0xEE (ea.rs drops nad-1 on hit_eof_marker)
            cur_words.append(encode_ea_word((prev_ev_tk + 1) & 0xFF, 1, 0, 100))
            frames_words.append((frame_start_tick, cur_words))
            cur_words = []
            frame_start_tick = coarse << 8

        cur_words.append(ef_word)
        for tk_i, pha, an in bin_events:
            cur_words.append(encode_ea_word(tk_i & 0xFF, an, 0, pha))
            prev_ev_tk = tk_i

    if cur_words:
        cur_words.append(encode_ea_word((prev_ev_tk + 1) & 0xFF, 1, 0, 100))
        frames_words.append((frame_start_tick, cur_words))

    n_frames = len(frames_words)
    rows = np.zeros((n_frames, 2242), dtype=np.uint8)

    for idx, (stbg_tick, wlist) in enumerate(frames_words):
        l1_t = STBG_EPOCH_OFFSET_SEC + stbg_tick / STBG_CLOCK_HZ
        buf = bytearray(2242)
        struct.pack_into(">d", buf, 0, l1_t)
        struct.pack_into(">q", buf, 36, idx)
        buf[44] = 0xDE
        buf[45] = unit_id
        buf[46] = 0xEA
        # Store lower 32 bits in i4..i1 (bytes 50..54), i5=0 (byte 49); ea.rs reconstructs i5 via 5th-byte rollover from l1_t!
        buf[49:54] = (stbg_tick & 0xFFFFFFFF).to_bytes(5, "big")
        buf[57] = 0
        struct.pack_into(">H", buf, 58, idx & 0xFFFF)
        pos = 60
        for w in wlist:
            buf[pos:pos + 5] = w
            pos += 5
        buf[pos] = 0xEE
        struct.pack_into(">hhh", buf, 2092, 0, 0, 0)
        rows[idx] = np.frombuffer(buf, dtype=np.uint8)

    return rows


def write_synthetic_laxpc_l1_dir(
    root_dir: str,
    obs_id: str,
    ra_deg: float,
    dec_deg: float,
    t_utc_start: float,
    t_utc_stop: float,
    t_ins_events: np.ndarray,
    phas: np.ndarray,
    anodes: np.ndarray,
    orbit_pos_fn=None,
) -> Dict[str, str]:
    if os.path.exists(root_dir):
        shutil.rmtree(root_dir)
    orb_dir = os.path.join(root_dir, "00001")
    ea_dir = os.path.join(orb_dir, "lxp1", "modeEA")
    aux_dir = os.path.join(orb_dir, "aux")
    os.makedirs(ea_dir, exist_ok=True)
    os.makedirs(aux_dir, exist_ok=True)

    # Real AstroSat LAXPC TCT requires dt_sps/dt_ins in (1.000, 1.001]; use exact slope 1.00001
    tct_slope = 1.00001
    t_ins0 = t_utc_start - STBG_EPOCH_OFFSET_SEC
    # Map t_ins_events from UTC-offset space to TCT instrument-clock space so tct(t_ins) == t_utc
    t_utc_events = t_ins_events + STBG_EPOCH_OFFSET_SEC
    t_ins_clk = t_ins0 + (t_utc_events - t_utc_start) / tct_slope
    ea_rows = build_synthetic_ea_frames(
        t_ins_clk, phas, anodes, unit_id=1, t0_offset_sec=STBG_EPOCH_OFFSET_SEC
    )

    ea_path = os.path.join(ea_dir, f"AS1_{obs_id}_lxp1EA_level1.fits")
    tct_path = os.path.join(aux_dir, f"AS1_{obs_id}_lxp_level1.tct")
    mkf_path = os.path.join(orb_dir, f"AS1_{obs_id}_lxp_level1.mkf")

    p_cards = [
        ("TELESCOP", "ASTROSAT", "mission"),
        ("INSTRUME", "LAXPC1", "instrument"),
        ("OBS_ID", obs_id, "observation id"),
        ("OBJECT", "SYNTH_SRC", "target name"),
        ("RA_NOM", ra_deg, "nominal RA"),
        ("DEC_NOM", dec_deg, "nominal DEC"),
        ("RA_PNT", ra_deg, "pointing RA"),
        ("DEC_PNT", dec_deg, "pointing DEC"),
        ("TSTART", t_utc_start, "start time UTC"),
        ("TSTOP", t_utc_stop, "stop time UTC"),
        ("MJDREFI", 55197, "MJD ref int"),
        ("MJDREFF", 0.0, "MJD ref frac"),
        ("TIMESYS", "UTC", "time system"),
    ]

    # Exact AstroSat LAXPC L1 row layout (`NAXIS1 = 2242`):
    # Time(8) + GRT(28) + FrameSequenceNumber(8) + DataArray(2048) + DecodingStatusFlag(2)
    # + DecodingErrorCount(2) + DecodingErrorCorrCount(2) + ErrorLocation(72) + DecodingErrorCorrFlag(72)
    n_fr = len(ea_rows)
    ea_hdu = _build_bintable_hdu(
        "LAXPC_MODEEA",
        [
            ("Time", "D", "s", ea_rows[:, 0:8].view(">f8").reshape(n_fr)),
            ("GRT", "7J", None, ea_rows[:, 8:36].view(">i4").reshape(n_fr, 7)),
            ("FrameSequenceNumber", "K", None, ea_rows[:, 36:44].view(">i8").reshape(n_fr)),
            ("DataArray", "2048B", None, ea_rows[:, 44:2092]),
            ("DecodingStatusFlag", "I", None, ea_rows[:, 2092:2094].view(">i2").reshape(n_fr)),
            ("DecodingErrorCount", "I", None, ea_rows[:, 2094:2096].view(">i2").reshape(n_fr)),
            ("DecodingErrorCorrCount", "I", None, ea_rows[:, 2096:2098].view(">i2").reshape(n_fr)),
            ("ErrorLocation", "36I", None, ea_rows[:, 2098:2170].view(">i2").reshape(n_fr, 36)),
            ("DecodingErrorCorrFlag", "36I", None, ea_rows[:, 2170:2242].view(">i2").reshape(n_fr, 36)),
        ],
    )
    gti_hdu = _build_bintable_hdu(
        "GTI",
        [("START", "D", "s", [t_utc_start]), ("STOP", "D", "s", [t_utc_stop])],
    )
    write_synthetic_fits(ea_path, p_cards, [ea_hdu, gti_hdu])

    # TCT: exact linear mapping t_utc = t_utc_start + tct_slope * (t_ins - t_ins0)
    t_utc_grid = np.linspace(t_utc_start - 100.0, t_utc_stop + 100.0, 20)
    t_ins_grid = t_ins0 + (t_utc_grid - t_utc_start) / tct_slope
    tct_hdu = _build_bintable_hdu(
        "TCT",
        [
            ("SPS_TIME", "D", "s", t_utc_grid),
            ("SPS_OBT", "D", None, t_ins_grid),
            ("INSTRUMENT_TIME", "D", "s", t_ins_grid),
        ],
    )
    write_synthetic_fits(tct_path, p_cards, [tct_hdu])

    # MKF: 1-s cadence, clean ELV=45 deg, SAA_FLAG=0, EARTHLON=100 deg, CPM_Rate=0
    n_mkf = int(math.ceil(t_utc_stop - t_utc_start)) + 21
    mkf_t = np.linspace(t_utc_start - 10.0, t_utc_stop + 10.0, n_mkf)
    if orbit_pos_fn is None:
        posx = np.full(n_mkf, 7028.0, dtype=np.float32)
        posy = np.zeros(n_mkf, dtype=np.float32)
        posz = np.zeros(n_mkf, dtype=np.float32)
        velx = np.zeros(n_mkf, dtype=np.float32)
        vely = np.full(n_mkf, 7.53, dtype=np.float32)
        velz = np.zeros(n_mkf, dtype=np.float32)
    else:
        pos, vel = orbit_pos_fn(mkf_t)
        posx, posy, posz = pos[:, 0].astype(np.float32), pos[:, 1].astype(np.float32), pos[:, 2].astype(np.float32)
        velx, vely, velz = vel[:, 0].astype(np.float32), vel[:, 1].astype(np.float32), vel[:, 2].astype(np.float32)

    mkf_hdu = _build_bintable_hdu(
        "MKF",
        [
            ("TIME", "D", "s", mkf_t),
            ("ELV", "E", "deg", np.full(n_mkf, 45.0, dtype=np.float32)),
            ("ANG_OFFSET", "E", "deg", np.zeros(n_mkf, dtype=np.float32)),
            ("SAA_FLAG", "E", None, np.zeros(n_mkf, dtype=np.float32)),
            ("EARTHLON", "E", "deg", np.full(n_mkf, 100.0, dtype=np.float32)),
            ("EARTHLAT", "E", "deg", np.zeros(n_mkf, dtype=np.float32)),
            ("CPM_Rate", "E", None, np.zeros(n_mkf, dtype=np.float32)),
            ("POSX", "E", "km", posx),
            ("POSY", "E", "km", posy),
            ("POSZ", "E", "km", posz),
            ("VELX", "E", "km/s", velx),
            ("VELY", "E", "km/s", vely),
            ("VELZ", "E", "km/s", velz),
        ],
    )
    write_synthetic_fits(mkf_path, p_cards, [mkf_hdu])

    return {"root": root_dir, "orb_dir": orb_dir, "ea": ea_path, "tct": tct_path, "mkf": mkf_path}


# ============================================================================
# Scenario 2A: Planted Spectral Power-Law & Gaussian Line Recovery
# ============================================================================

def test_scenario2a_spectral_recovery(bin_path: str) -> Dict[str, Any]:
    sc_dir = os.path.join(WORK_ROOT, "sc2a_spectral")
    l1_dir = os.path.join(sc_dir, "l1")
    out_dir = os.path.join(sc_dir, "l2")
    if os.path.exists(sc_dir):
        shutil.rmtree(sc_dir)

    # Plant T = 50.0 s at 600 cps -> 30,000 events cleanly inside [t0+1, t0+51]
    t_utc0 = 200_000_000.0
    t_start = t_utc0
    t_stop = t_utc0 + 54.0
    gti_start = t_utc0 + 2.0
    gti_stop = t_utc0 + 52.0
    t_exp_true = gti_stop - gti_start  # 50.0 s

    # Generate evenly spaced timestamps across [gti_start-1, gti_stop+1] in t_ins = t_utc - STBG_EPOCH_OFFSET_SEC
    rng = np.random.default_rng(20261126)
    n_in_gti = 30_000
    t_ins_gti0 = gti_start - STBG_EPOCH_OFFSET_SEC
    t_ins_gti1 = gti_stop - STBG_EPOCH_OFFSET_SEC
    t_ins_gti = np.linspace(t_ins_gti0 + 0.005, t_ins_gti1 - 0.005, n_in_gti)
    t_ins_pre = np.linspace(t_ins_gti0 - 1.8, t_ins_gti0 - 0.05, 100)
    t_ins_post = np.linspace(t_ins_gti1 + 0.05, t_ins_gti1 + 1.8, 100)
    t_ins = np.concatenate([t_ins_pre, t_ins_gti, t_ins_post])

    # Planted spectral shape in 1024-ch PHA [40..800] -> 512-ch [20..400]:
    # 80% power-law (p ~ ch^-1.7) + 20% Gaussian line centered at PHA=200 (ch_512=100), sigma=8 PHA (4 ch_512)
    n_tot = len(t_ins)
    phas = np.zeros(n_tot, dtype=np.int32)
    # Build exact discrete PDF over PHA in 40..799
    ch_grid = np.arange(40, 800, dtype=np.float64)
    pdf_pl = ch_grid ** (-1.7)
    pdf_pl /= np.sum(pdf_pl)
    pdf_gauss = np.exp(-0.5 * ((ch_grid - 200.0) / 8.0) ** 2)
    pdf_gauss /= np.sum(pdf_gauss)
    pdf_total = 0.80 * pdf_pl + 0.20 * pdf_gauss
    pdf_total /= np.sum(pdf_total)
    phas = rng.choice(ch_grid.astype(np.int32), size=n_tot, p=pdf_total)
    anodes = rng.integers(1, 8, size=n_tot)

    paths = write_synthetic_laxpc_l1_dir(
        l1_dir, "SYNTH_SPEC", 83.6331, 22.0145, t_start, t_stop, t_ins, phas, anodes
    )

    user_gti_txt = os.path.join(sc_dir, "exact_gti.txt")
    with open(user_gti_txt, "w", encoding="utf-8") as fh:
        fh.write(f"{gti_start:.6f} {gti_stop:.6f}\n")

    res = json.loads(require_ok(run_cli(bin_path, [
        "--json", "laxpc", "run",
        "--input", paths["orb_dir"],
        "--units", "lxp1",
        "--channels", "512",
        "--deadtime", "0.000043",
        "--user-gti", user_gti_txt,
        "--user-gti-only",
        "--no-orbit",
        "--output", out_dir,
    ]), "synthetic spectral recovery").stdout)

    u = res["unit_results"][0]
    spec_h = next(h for h in read_fits(u["spectrum_fits_path"]) if h.extname == "SPECTRUM")
    exp_rec = float(spec_h.header["EXPOSURE"])
    rates = np.asarray(spec_h.read_column("RATE"), dtype=np.float64)
    errs = np.asarray(spec_h.read_column("STAT_ERR"), dtype=np.float64)

    # Exact analytical ground truth (exposure in instrument clock is (gti_stop - gti_start) / tct_slope):
    # Raw events in [gti_start, gti_stop] = n_in_gti = 30,000; R_obs = N / exposure_ins; tau = 43 us
    t_exp_ins_true = t_exp_true / 1.00001
    r_obs_true = n_in_gti / t_exp_ins_true
    fcd_true = 1.0 / (1.0 - 43e-6 * r_obs_true)
    r_corr_true = r_obs_true * fcd_true
    r_corr_rec = float(np.sum(rates))

    assert abs(exp_rec - t_exp_ins_true) < 1e-4, f"Exposure mismatch: rec={exp_rec}, true={t_exp_ins_true}"
    rel_rate_err = abs(r_corr_rec - r_corr_true) / r_corr_true
    assert rel_rate_err < 1e-4, f"Deadtime-corrected rate mismatch: rec={r_corr_rec:.4f}, true={r_corr_true:.4f}"

    # Verify STAT_ERR == sqrt(RATE / EXPOSURE)
    expected_errs = np.sqrt(np.maximum(rates / exp_rec, 0.0))
    assert np.allclose(errs, expected_errs, atol=1e-8), "STAT_ERR != sqrt(RATE / EXPOSURE)"

    # Verify Gaussian line peak & flux recovery in 512-ch space (PHA=200 -> ch=100, sigma=4 ch)
    # Using exact planted analytical continuum `0.80 * pdf_pl` binned to 512 channels
    ch512 = np.arange(512, dtype=np.float64)
    pdf_512_pl = np.zeros(512, dtype=np.float64)
    pdf_512_gauss = np.zeros(512, dtype=np.float64)
    for k, ch in enumerate(ch_grid.astype(int)):
        pdf_512_pl[ch // 2] += 0.80 * pdf_pl[k]
        pdf_512_gauss[ch // 2] += 0.20 * pdf_gauss[k]
    expected_rate_512 = (pdf_512_pl + pdf_512_gauss) * r_corr_true
    cont_rate_512 = pdf_512_pl * r_corr_true

    line_mask = (ch512 >= 90) & (ch512 <= 110)
    line_excess_rate = float(np.sum(rates[line_mask] - cont_rate_512[line_mask]))
    line_centroid = float(np.sum(ch512[line_mask] * (rates[line_mask] - cont_rate_512[line_mask])) / line_excess_rate)
    expected_line_rate = float(np.sum(pdf_512_gauss[line_mask]) * r_corr_true)
    line_err = float(np.sqrt(np.sum(errs[line_mask] ** 2)))
    z_line = abs(line_excess_rate - expected_line_rate) / max(line_err, 1e-6)
    assert abs(line_centroid - 100.0) < 0.8, f"Gaussian centroid off: {line_centroid:.2f} vs 100.0"
    assert z_line < 2.0, f"Line flux deviates {z_line:.2f} sigma: rec={line_excess_rate:.2f}, true={expected_line_rate:.2f}"

    # Verify chi2 goodness-of-fit across the populated band [25..380] against exact analytical PDF
    band = (ch512 >= 25) & (ch512 <= 380) & (expected_rate_512 > 0)
    red_chi2 = float(np.mean(((rates[band] - expected_rate_512[band]) / errs[band]) ** 2))
    assert red_chi2 < 1.35, f"Spectrum reduced chi2={red_chi2:.3f} vs planted analytical PDF too high"

    # Verify power-law slope on continuum region [130..350] (ch512)
    mask_pl = (ch512 >= 130) & (ch512 <= 350) & (rates > 0)
    slope, _ = np.polyfit(np.log(ch512[mask_pl]), np.log(rates[mask_pl]), 1)
    assert abs(slope - (-1.7)) < 0.12, f"Power-law slope off: {slope:.3f} vs -1.70"

    fv = verify_fitsverify([out_dir])
    assert fv["failed_files"] == 0, f"Synthetic 2A fitsverify error: {fv['failures']}"

    return {
        "passed": True,
        "exposure_true_s": t_exp_ins_true,
        "exposure_rec_s": exp_rec,
        "rate_corr_true_cps": r_corr_true,
        "rate_corr_rec_cps": r_corr_rec,
        "rate_rel_err": rel_rate_err,
        "line_centroid_ch": line_centroid,
        "line_rate_true_cps": expected_line_rate,
        "line_rate_rec_cps": line_excess_rate,
        "line_z_sigma": z_line,
        "reduced_chi2_vs_analytical": red_chi2,
        "powerlaw_slope_rec": float(slope),
        "fitsverify_total": fv["total_files"],
    }


# ============================================================================
# Scenario 2B: Planted 17.3456 Hz Pulsar Across 2016-2017 Leap Second + Orbit
# ============================================================================

def _earth_ssb_jpl_de405(jd_tdb: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Read JPL DE405 (`lnxp1600p2200.405`) directly in Python to compute exact Earth SSB pos/vel (m, m/s)."""
    with open(JPLEPH_405, "rb") as fh:
        raw = fh.read()
    # Record 0 header (size = 1018 * 8 = 8144 bytes for DE405)
    # Determine KSIZE from header: DE405 ncoeff = 1018 -> rec_size = 8144
    ncoeff = 1018
    rec_bytes = ncoeff * 8
    hdr = raw[:rec_bytes]
    # ss[0..3] at offset 2652, emrat at 2684, numde at 2692, ipt at 2696 (matches bary.rs:895-920)
    jd_start, jd_end, span = struct.unpack_from("<ddd", hdr, 2652)
    emrat = struct.unpack_from("<d", hdr, 2684)[0]
    ipt_flat = struct.unpack_from("<36i", hdr, 2696)
    ipt = [(ipt_flat[3 * i], ipt_flat[3 * i + 1], ipt_flat[3 * i + 2]) for i in range(12)]

    def eval_body(body_idx: int, rec: bytes, norm_t: float) -> Tuple[np.ndarray, np.ndarray]:
        start_1, ncf, nsub = ipt[body_idx]
        sub_span = span / nsub
        sub_idx = min(max(int(norm_t * nsub), 0), nsub - 1)
        tau = 2.0 * (norm_t * nsub - sub_idx) - 1.0
        dtau_dt = 2.0 * nsub / (span * 86400.0)
        off_base = (start_1 - 1 + sub_idx * 3 * ncf) * 8
        pos = np.zeros(3)
        vel = np.zeros(3)
        t_poly = np.zeros(ncf)
        d_poly = np.zeros(ncf)
        t_poly[0] = 1.0
        if ncf > 1:
            t_poly[1] = tau
        for k in range(2, ncf):
            t_poly[k] = 2.0 * tau * t_poly[k - 1] - t_poly[k - 2]
        d_poly[0] = 0.0
        if ncf > 1:
            d_poly[1] = 1.0
        for k in range(2, ncf):
            d_poly[k] = 2.0 * t_poly[k - 1] + 2.0 * tau * d_poly[k - 1] - d_poly[k - 2]
        for ax in range(3):
            c = np.frombuffer(rec, dtype="<f8", count=ncf, offset=off_base + ax * ncf * 8)
            pos[ax] = float(np.dot(c, t_poly))
            vel[ax] = float(np.dot(c, d_poly) * dtau_dt)
        return pos, vel

    n = len(jd_tdb)
    earth_pos_m = np.zeros((n, 3))
    earth_vel_ms = np.zeros((n, 3))
    fac = 1.0 / (1.0 + emrat)
    for i, jd in enumerate(jd_tdb):
        rec_idx = int((jd - jd_start) / span)
        rec_start = jd_start + rec_idx * span
        norm_t = (jd - rec_start) / span
        rec = raw[(rec_idx + 2) * rec_bytes : (rec_idx + 3) * rec_bytes]
        emb_p, emb_v = eval_body(2, rec, norm_t)
        moon_p, moon_v = eval_body(9, rec, norm_t)
        # Convert km -> m, km/s -> m/s
        earth_pos_m[i] = (emb_p - fac * moon_p) * 1000.0
        earth_vel_ms[i] = (emb_v - fac * moon_v) * 1000.0
    return earth_pos_m, earth_vel_ms


def test_scenario2b_leap_second_pulsar(bin_path: str) -> Dict[str, Any]:
    sc_dir = os.path.join(WORK_ROOT, "sc2b_leap_pulsar")
    l1_dir = os.path.join(sc_dir, "l1")
    out_dir = os.path.join(sc_dir, "l2")
    if os.path.exists(sc_dir):
        shutil.rmtree(sc_dir)

    # Leap second: 2017-01-01T00:00:00 UTC => MJD_UTC = 57754.0
    # AstroSat UTC = (57754.0 - 55197.0) * 86400 = 220,924,800.0 s
    t_leap_utc = (57754.0 - 55197.0) * 86400.0
    t_start_utc = t_leap_utc - 300.0  # 23:55:00 UTC on 2016-12-31
    t_stop_utc = t_leap_utc + 300.0   # 00:05:00 UTC on 2017-01-01

    ra_deg, dec_deg = 83.6331, 22.0145
    ra_r, dec_r = math.radians(ra_deg), math.radians(dec_deg)
    n_hat = np.array([math.cos(dec_r) * math.cos(ra_r), math.cos(dec_r) * math.sin(ra_r), math.sin(dec_r)])

    # LEO circular AstroSat orbit (r = 7028 km, v = 7.53 km/s, incl = 6 deg)
    r_km = 7028.0
    v_kms = 7.53
    omega = v_kms / r_km
    incl = math.radians(6.0)

    def leo_orbit(t_utc_arr: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        dt = t_utc_arr - t_start_utc
        th = omega * dt
        c, s = np.cos(th), np.sin(th)
        ci, si = math.cos(incl), math.sin(incl)
        pos = np.column_stack([r_km * c, r_km * s * ci, r_km * s * si])
        vel = np.column_stack([-v_kms * s, v_kms * c * ci, v_kms * c * si])
        return pos, vel

    # Compute Roemer + Einstein delay on SEPARATE dense pre-leap and post-leap UTC grids
    # so we NEVER interpolate across the 1-second 2017-01-01T00:00:00 UTC leap-second discontinuity!
    c_light = 299_792_458.0

    def utc_to_tdb_branch(t_utc_arr: np.ndarray, tai_offset: float) -> np.ndarray:
        t_tt = t_utc_arr + tai_offset + 32.184
        jd_tt = 55197.0 + 2400000.5 + t_tt / 86400.0
        t_mill = (jd_tt - 2451545.0) / 365250.0
        einstein = (
            1656.674564e-6 * np.sin(6283.075943033 * t_mill + 6.240054195)
            + 22.417471e-6 * np.sin(5753.384970095 * t_mill + 4.296977442)
            + 13.839792e-6 * np.sin(12566.151886066 * t_mill + 6.196904410)
        )
        jd_tdb = jd_tt + einstein / 86400.0
        e_pos_m, _ = _earth_ssb_jpl_de405(jd_tdb)
        sc_pos_km, _ = leo_orbit(t_utc_arr)
        roemer_earth = (e_pos_m @ n_hat) / c_light
        roemer_sc = (sc_pos_km * 1000.0 @ n_hat) / c_light
        return t_tt + einstein + roemer_earth + roemer_sc

    t_grid_pre = np.linspace(t_start_utc - 5.0, t_leap_utc - 1e-3, 2000)
    t_grid_post = np.linspace(t_leap_utc + 1e-3, t_stop_utc + 5.0, 2000)
    t_tdb_pre = utc_to_tdb_branch(t_grid_pre, 36.0)
    t_tdb_post = utc_to_tdb_branch(t_grid_post, 37.0)
    leap_tdb_pre = float(t_tdb_pre[-1])
    leap_tdb_post = float(t_tdb_post[0])

    # Plant f_true = 17.3456 Hz periodic pulsar in TDB across the 600 s window
    f_true = 17.3456
    # Generate ~8000 pulsed photons uniformly in cycle index (excluding the 1s leap-second gap cleanly)
    n_photons_target = 8000
    rng = np.random.default_rng(173456)
    n_cycles = int((t_tdb_post[-1] - t_tdb_pre[0] - 4.0) * f_true)
    cycle_indices = np.sort(rng.choice(n_cycles, size=n_photons_target, replace=False))
    # Sharp pulse around phase 0.30 (sigma = 0.04 cycles)
    phases_planted = 0.30 + rng.normal(0.0, 0.04, size=n_photons_target)
    t_tdb_start_phot = float(t_tdb_pre[0]) + 2.0
    t_tdb_phot = t_tdb_start_phot + (cycle_indices + phases_planted) / f_true
    # Exclude photons inside the 1-second leap-second gap [leap_tdb_pre, leap_tdb_post] so zero events clamp
    valid_mask = (t_tdb_phot <= leap_tdb_pre) | (t_tdb_phot >= leap_tdb_post)
    t_tdb_phot = np.sort(t_tdb_phot[valid_mask])
    n_photons = len(t_tdb_phot)

    # Invert t_tdb -> t_utc cleanly on exact pre-leap and post-leap branches
    mask_pre = t_tdb_phot <= leap_tdb_pre
    t_utc_phot = np.empty_like(t_tdb_phot)
    t_utc_phot[mask_pre] = np.interp(t_tdb_phot[mask_pre], t_tdb_pre, t_grid_pre)
    t_utc_phot[~mask_pre] = np.interp(t_tdb_phot[~mask_pre], t_tdb_post, t_grid_post)

    # Instrument time: t_ins = t_utc - STBG_EPOCH_OFFSET_SEC
    t_ins_phot = t_utc_phot - STBG_EPOCH_OFFSET_SEC
    phas = rng.integers(50, 400, size=n_photons)
    anodes = rng.integers(1, 8, size=n_photons)

    paths = write_synthetic_laxpc_l1_dir(
        l1_dir, "LEAP_PSR", ra_deg, dec_deg, t_start_utc, t_stop_utc,
        t_ins_phot, phas, anodes, orbit_pos_fn=leo_orbit,
    )

    # Full L1 -> L2 with user GTI covering full 600 s window -> AstroSat.orb + lxp1level2.evn
    gti_txt = os.path.join(sc_dir, "full_gti.txt")
    with open(gti_txt, "w", encoding="utf-8") as fh:
        fh.write(f"{t_start_utc + 1.0:.6f} {t_stop_utc - 1.0:.6f}\n")

    l2_res = json.loads(require_ok(run_cli(bin_path, [
        "--json", "laxpc", "run",
        "--input", paths["orb_dir"],
        "--units", "lxp1",
        "--user-gti", gti_txt,
        "--user-gti-only",
        "--output", out_dir,
    ]), "leap-second pulsar laxpc run").stdout)

    evn_path = l2_res["unit_results"][0]["event_fits_path"]
    orb_path = l2_res["orbit_fits_path"]
    bary_evn_ana = os.path.join(out_dir, "lxp1_bary_ana.evn")
    bary_evn_de = os.path.join(out_dir, "lxp1_bary_de405.evn")
    bary_evn_fits405 = os.path.join(out_dir, "lxp1_bary_jpleph405_fits.evn")

    # Run `astrosat bary` with built-in analytical, binary JPL DE405, AND HEASARC FITS JPLEPH.405
    require_ok(run_cli(bin_path, [
        "--json", "bary", evn_path, "--orbit", orb_path, "-o", bary_evn_ana,
    ]), "leap-second bary analytical")
    require_ok(run_cli(bin_path, [
        "--json", "bary", evn_path, "--orbit", orb_path, "--ephem", JPLEPH_405, "-o", bary_evn_de,
    ]), "leap-second bary DE405")
    require_ok(run_cli(bin_path, [
        "--json", "bary", evn_path, "--orbit", orb_path, "--ephem", JPLEPH_FITS_405, "-o", bary_evn_fits405,
    ]), "leap-second bary FITS JPLEPH.405")

    # Run `astrosat timing` on all three with fine frequency resolution (`df = 2e-5 Hz`)!
    tim_ana = json.loads(require_ok(run_cli(bin_path, [
        "--json", "timing", bary_evn_ana,
        "--freq-min", "17.340", "--freq-max", "17.351", "--df", "2e-5", "--harmonics", "2", "--nbins", "32",
    ]), "leap-second timing analytical").stdout)
    tim_de = json.loads(require_ok(run_cli(bin_path, [
        "--json", "timing", bary_evn_de,
        "--freq-min", "17.340", "--freq-max", "17.351", "--df", "2e-5", "--harmonics", "2", "--nbins", "32",
    ]), "leap-second timing DE405").stdout)
    tim_f405 = json.loads(require_ok(run_cli(bin_path, [
        "--json", "timing", bary_evn_fits405,
        "--freq-min", "17.340", "--freq-max", "17.351", "--df", "2e-5", "--harmonics", "2", "--nbins", "32",
    ]), "leap-second timing FITS JPLEPH.405").stdout)

    err_ana = abs(tim_ana["best_freq_hz"] - f_true)
    err_de = abs(tim_de["best_freq_hz"] - f_true)
    err_f405 = abs(tim_f405["best_freq_hz"] - f_true)
    assert err_ana < 1e-4, f"Analytical bary frequency error {err_ana:.6e} >= 1e-4 Hz"
    assert err_de < 1e-4, f"DE405 bary frequency error {err_de:.6e} >= 1e-4 Hz"
    assert err_f405 < 1e-4, f"FITS JPLEPH.405 bary frequency error {err_f405:.6e} >= 1e-4 Hz"
    assert tim_ana["max_z2_n"] > 1000.0 and tim_de["max_z2_n"] > 1000.0 and tim_f405["max_z2_n"] > 1000.0, (
        "Weak pulsar Z2 across leap second"
    )

    # Check phase coherence BEFORE vs AFTER the 2017-01-01 leap second on BOTH DE405 & analytical (plus DE405 == FITS JPLEPH.405)!
    leap_tdb_mid = 0.5 * (leap_tdb_pre + leap_tdb_post)

    def check_leap_phase(path: str, label: str) -> Tuple[float, float, float, np.ndarray]:
        ev_h = next(h for h in read_fits(path) if h.extname == "LAXPC_MODE_EA")
        t_bary = np.asarray(ev_h.read_column("TIME"), dtype=np.float64)
        pre_times = t_bary[t_bary < leap_tdb_mid]
        post_times = t_bary[t_bary >= leap_tdb_mid]
        assert len(pre_times) > 1000 and len(post_times) > 1000, f"{label}: missing pre/post leap-second events"

        def mean_phase(ts: np.ndarray) -> float:
            ph = np.mod((ts - t_bary[0]) * f_true, 1.0)
            ang = 2.0 * math.pi * ph
            return float((math.atan2(np.mean(np.sin(ang)), np.mean(np.cos(ang))) / (2.0 * math.pi)) % 1.0)

        p0 = mean_phase(pre_times)
        p1 = mean_phase(post_times)
        dp = abs(((p1 - p0 + 0.5) % 1.0) - 0.5)
        assert dp < 0.02, f"{label}: phase jump across 2017-01-01 leap second: pre={p0:.4f}, post={p1:.4f}, dphi={dp:.4f}"
        return p0, p1, dp, t_bary

    ph_pre, ph_post, dphi, t_bary_de = check_leap_phase(bary_evn_de, "DE405")
    _, _, dphi_ana, _ = check_leap_phase(bary_evn_ana, "Analytical")
    _, _, dphi_f405, t_bary_f405 = check_leap_phase(bary_evn_fits405, "FITS-JPLEPH.405")
    assert np.allclose(t_bary_de, t_bary_f405, atol=1e-9), "FITS JPLEPH.405 vs binary DE405 mismatch across leap second!"

    # Contrast: raw un-barycentered event file MUST fail frequency accuracy & have degraded chi2/Z2 across leap second!
    raw_tim = json.loads(require_ok(run_cli(bin_path, [
        "--json", "timing", evn_path,
        "--freq-min", "17.340", "--freq-max", "17.351", "--df", "2e-5", "--harmonics", "2", "--nbins", "32",
    ]), "unbarycentered timing contrast").stdout)
    assert abs(raw_tim["best_freq_hz"] - f_true) > 2e-4, "Expected un-barycentered frequency to be Doppler-shifted"
    assert raw_tim["max_chi2"] < 0.5 * tim_de["max_chi2"], "Expected barycentering to boost chi2 across leap second + orbit"

    fv = verify_fitsverify([out_dir])
    assert fv["failed_files"] == 0, f"Synthetic 2B fitsverify error: {fv['failures']}"

    return {
        "passed": True,
        "planted_freq_hz": f_true,
        "analytical_best_freq_hz": tim_ana["best_freq_hz"],
        "analytical_freq_err_hz": err_ana,
        "analytical_z2_2": tim_ana["max_z2_n"],
        "analytical_phase_jump_cycles": dphi_ana,
        "de405_best_freq_hz": tim_de["best_freq_hz"],
        "de405_freq_err_hz": err_de,
        "de405_z2_2": tim_de["max_z2_n"],
        "fits_jpleph405_freq_err_hz": err_f405,
        "fits_jpleph405_phase_jump_cycles": dphi_f405,
        "raw_unbary_z2_2": raw_tim["max_z2_n"],
        "phase_pre_leap": ph_pre,
        "phase_post_leap": ph_post,
        "phase_jump_cycles": dphi,
        "fitsverify_total": fv["total_files"],
    }


# ============================================================================
# Section 3: Adversarial & Corrupt-Input Stress Testing
# ============================================================================

def test_section3_adversarial(bin_path: str) -> Dict[str, Any]:
    adv_dir = os.path.join(WORK_ROOT, "sec3_adversarial")
    if os.path.exists(adv_dir):
        shutil.rmtree(adv_dir)
    os.makedirs(adv_dir, exist_ok=True)

    cases: List[Dict[str, Any]] = []

    def check_adversarial(name: str, args: List[str], allow_zero: bool = False) -> None:
        res = run_cli(bin_path, args)
        ok = (not res.panicked) and (allow_zero or res.returncode != 0)
        if not ok:
            raise AssertionError(
                f"Adversarial case '{name}' failed (panicked={res.panicked}, rc={res.returncode}):\n"
                f"STDERR: {res.stderr[:800]}"
            )
        cases.append({"case": name, "rc": res.returncode, "panicked": res.panicked, "msg_preview": (res.stderr or res.stdout).strip()[:120]})

    # 1. Empty 0-byte file, 100-byte garbage, header-only truncated FITS
    f_empty = os.path.join(adv_dir, "empty.fits")
    f_garbage = os.path.join(adv_dir, "garbage.fits")
    f_hdr_only = os.path.join(adv_dir, "hdr_only.fits")
    f_mid_trunc = os.path.join(adv_dir, "mid_trunc.fits")
    open(f_empty, "wb").close()
    with open(f_garbage, "wb") as fh:
        fh.write(b"NOT_FITS_" * 15)
    # Truncate real sample file mid-data and header-only
    with open(CZTI_M0, "rb") as fh:
        sample_bytes = fh.read(2880 * 4)
    with open(f_hdr_only, "wb") as fh:
        fh.write(sample_bytes[:2880])
    with open(f_mid_trunc, "wb") as fh:
        fh.write(sample_bytes[:2880 * 2 + 123])

    for label, fp in [("empty_0b", f_empty), ("garbage_135b", f_garbage), ("header_only_trunc", f_hdr_only), ("mid_block_trunc", f_mid_trunc)]:
        check_adversarial(f"inspect_{label}", ["inspect", "--full", fp])
        check_adversarial(f"validate_{label}", ["validate", "--strict", fp])
        check_adversarial(f"bary_{label}", ["bary", fp, "--orbit", CZTI_ORB, "--ra", "83.6", "--dec", "22.0", "-o", os.path.join(adv_dir, f"out_{label}.fits")])
        check_adversarial(f"timing_{label}", ["timing", fp, "--freq", "29.66"])

    # 2. Corrupt FITS headers (`NAXIS1=-5`, `BITPIX=9`, huge `PCOUNT`, huge `NAXIS2`, huge `TFIELDS`)
    def write_corrupt_fits(path: str, pri_overrides: Dict[str, Any], ext_cards: Optional[List[Tuple[str, Any, str]]] = None) -> None:
        pri = [
            _fits_card("SIMPLE", True),
            _fits_card("BITPIX", pri_overrides.get("BITPIX", 8)),
            _fits_card("NAXIS", pri_overrides.get("NAXIS", 0)),
            _fits_card("EXTEND", True),
        ]
        if "PCOUNT" in pri_overrides:
            pri.append(_fits_card("PCOUNT", pri_overrides["PCOUNT"]))
        buf = _pad_header_2880(pri)
        if ext_cards is not None:
            buf += _pad_header_2880(ext_cards) + _pad_2880(b"\x00" * 64)
        with open(path, "wb") as fh:
            fh.write(buf)

    f_bitpix9 = os.path.join(adv_dir, "bitpix9.fits")
    write_corrupt_fits(f_bitpix9, {"BITPIX": 9})

    f_naxis1_neg = os.path.join(adv_dir, "naxis1_neg5.fits")
    write_corrupt_fits(f_naxis1_neg, {}, [
        _fits_card("XTENSION", "BINTABLE"),
        _fits_card("BITPIX", 8),
        _fits_card("NAXIS", 2),
        _fits_card("NAXIS1", -5),
        _fits_card("NAXIS2", 10),
        _fits_card("PCOUNT", 0),
        _fits_card("GCOUNT", 1),
        _fits_card("TFIELDS", 1),
        _fits_card("TTYPE1", "TIME"),
        _fits_card("TFORM1", "1D"),
        _fits_card("EXTNAME", "EVENTS"),
    ])

    f_huge_pcount = os.path.join(adv_dir, "huge_pcount.fits")
    write_corrupt_fits(f_huge_pcount, {}, [
        _fits_card("XTENSION", "BINTABLE"),
        _fits_card("BITPIX", 8),
        _fits_card("NAXIS", 2),
        _fits_card("NAXIS1", 8),
        _fits_card("NAXIS2", 1),
        _fits_card("PCOUNT", 2147483647),
        _fits_card("GCOUNT", 1),
        _fits_card("TFIELDS", 1),
        _fits_card("TTYPE1", "TIME"),
        _fits_card("TFORM1", "1D"),
        _fits_card("EXTNAME", "EVENTS"),
    ])

    f_huge_naxis2 = os.path.join(adv_dir, "huge_naxis2.fits")
    write_corrupt_fits(f_huge_naxis2, {}, [
        _fits_card("XTENSION", "BINTABLE"),
        _fits_card("BITPIX", 8),
        _fits_card("NAXIS", 2),
        _fits_card("NAXIS1", 8),
        _fits_card("NAXIS2", 100000000),
        _fits_card("PCOUNT", 0),
        _fits_card("GCOUNT", 1),
        _fits_card("TFIELDS", 1),
        _fits_card("TTYPE1", "TIME"),
        _fits_card("TFORM1", "1D"),
        _fits_card("EXTNAME", "EVENTS"),
    ])

    f_huge_tfields = os.path.join(adv_dir, "huge_tfields.fits")
    write_corrupt_fits(f_huge_tfields, {}, [
        _fits_card("XTENSION", "BINTABLE"),
        _fits_card("BITPIX", 8),
        _fits_card("NAXIS", 2),
        _fits_card("NAXIS1", 8),
        _fits_card("NAXIS2", 1),
        _fits_card("PCOUNT", 0),
        _fits_card("GCOUNT", 1),
        _fits_card("TFIELDS", 9999999),
        _fits_card("EXTNAME", "EVENTS"),
    ])

    for label, fp in [
        ("bitpix9", f_bitpix9),
        ("naxis1_neg5", f_naxis1_neg),
        ("huge_pcount", f_huge_pcount),
        ("huge_naxis2", f_huge_naxis2),
        ("huge_tfields", f_huge_tfields),
    ]:
        # `inspect` reads primary metadata and skips invalid extension tables gracefully (`allow_zero=True` for huge_tfields), while `validate` MUST reject all of them
        check_adversarial(f"inspect_{label}", ["inspect", "--full", fp], allow_zero=(label == "huge_tfields"))
        check_adversarial(f"validate_{label}", ["validate", "--strict", fp])
        check_adversarial(f"bary_{label}", ["bary", fp, "--orbit", CZTI_ORB, "--ra", "83.6", "--dec", "22.0", "-o", os.path.join(adv_dir, f"bary_{label}.fits")])
        check_adversarial(f"timing_{label}", ["timing", fp, "--freq", "10.0"])
        check_adversarial(f"czti_decode_{label}", ["czti", "decode", "--m0", fp, "--tct", CZTI_TCT, "-o", os.path.join(adv_dir, f"cz_{label}.evt")])
        check_adversarial(f"laxpc_decode_{label}", ["laxpc", "decode", "--ea", fp, "--unit", "lxp1"])

    # 3. Empty GTI & inverted GTI
    f_empty_gti_fits = os.path.join(adv_dir, "empty_gti.fits")
    write_synthetic_fits(
        f_empty_gti_fits,
        [("TELESCOP", "ASTROSAT", ""), ("INSTRUME", "LAXPC", "")],
        [_build_bintable_hdu("GTI", [("START", "D", "s", []), ("STOP", "D", "s", [])])],
    )
    f_empty_gti_txt = os.path.join(adv_dir, "empty_gti.txt")
    with open(f_empty_gti_txt, "w", encoding="utf-8") as fh:
        fh.write("# empty gti\n")
    f_bad_gti_fits = os.path.join(adv_dir, "inverted_gti.fits")
    write_synthetic_fits(
        f_bad_gti_fits,
        [("TELESCOP", "ASTROSAT", "")],
        [_build_bintable_hdu("GTI", [("START", "D", "s", [200.0]), ("STOP", "D", "s", [100.0])])],
    )
    check_adversarial("laxpc_run_empty_user_gti_txt", [
        "laxpc", "run", "--input", LAXPC_L1_DIR, "--units", "lxp1",
        "--user-gti", f_empty_gti_txt, "--user-gti-only", "--output", os.path.join(adv_dir, "lxp_empty_gti"),
    ])
    # Empty user GTI FITS either fails cleanly or succeeds cleanly with 0 exposure / 0 events and 0 panics
    check_adversarial("laxpc_run_empty_user_gti_fits", [
        "laxpc", "run", "--input", LAXPC_L1_DIR, "--units", "lxp1",
        "--user-gti", f_empty_gti_fits, "--user-gti-only", "--output", os.path.join(adv_dir, "lxp_empty_gti2"),
    ], allow_zero=True)
    # Also test CZTI gtigen with empty user GTI (produces 0 exposure without panic)
    check_adversarial("czti_gtigen_empty_user_gti", [
        "czti", "gtigen", "--m0", CZTI_M0, "--mkf", CZTI_MKF,
        "--user-gti", f_empty_gti_txt, "-o", os.path.join(adv_dir, "cz_empty.gti"),
    ], allow_zero=True)
    check_adversarial("validate_inverted_gti", ["validate", "--strict", f_bad_gti_fits])

    # 4. Missing MKF/TCT/CALDB/ORBIT & non-existent paths
    missing = "/mnt/d/ISRO/.tmp/public_scenarios/does_not_exist_99999"
    check_adversarial("czti_missing_caldb", ["czti", "run", "--input", CZTI_L1_DIR, "--caldb", missing, "-o", os.path.join(adv_dir, "cz_no_caldb")])
    check_adversarial("czti_missing_m0", ["czti", "decode", "--m0", missing, "--tct", CZTI_TCT, "-o", os.path.join(adv_dir, "cz_no_m0.evt")])
    check_adversarial("czti_missing_tct", ["czti", "decode", "--m0", CZTI_M0, "--tct", missing, "-o", os.path.join(adv_dir, "cz_no_tct.evt")])
    check_adversarial("czti_missing_mkf", ["czti", "gtigen", "--m0", CZTI_M0, "--mkf", missing, "-o", os.path.join(adv_dir, "cz_no_mkf.gti")])
    check_adversarial("laxpc_missing_input", ["laxpc", "run", "--input", missing, "-o", os.path.join(adv_dir, "lxp_no_in")])
    check_adversarial("laxpc_missing_orbit_mkf", ["laxpc", "orbit", "--mkf", missing, "-o", os.path.join(adv_dir, "no_orb.orb")])
    check_adversarial("bary_missing_orbit", ["bary", LAXPC_EA1, "--orbit", missing, "--ra", "83.6", "--dec", "22.0", "-o", os.path.join(adv_dir, "bary_no_orb.fits")])
    check_adversarial("inspect_missing", ["inspect", missing])

    # 5. Invalid CLI flags & parameter bounds
    check_adversarial("invalid_flag", ["--nonexistent-flag"])
    check_adversarial("czti_bad_stage", ["czti", "run", "--input", CZTI_L1_DIR, "--caldb", CZTI_CALDB, "--entry-stage", "9", "-o", os.path.join(adv_dir, "bad_stg")])
    check_adversarial("laxpc_bad_tbin", ["laxpc", "lightcurve", "--ea", LAXPC_EA1, "--mkf", LAXPC_MKF, "--tct", LAXPC_TCT, "--unit", "lxp1", "--tbin", "-1.0", "-o", os.path.join(adv_dir, "bad_tbin.lcurv")])
    check_adversarial("laxpc_bad_channels", ["laxpc", "spectrum", "--ea", LAXPC_EA1, "--mkf", LAXPC_MKF, "--tct", LAXPC_TCT, "--unit", "lxp1", "--channels", "333", "-o", os.path.join(adv_dir, "bad_ch.spec")])
    check_adversarial("timing_bad_harmonics", ["timing", LAXPC_EA1, "--freq", "10.0", "--harmonics", "0"])
    check_adversarial("timing_bad_nbins", ["timing", LAXPC_EA1, "--freq", "10.0", "--nbins", "1"])
    check_adversarial("timing_huge_ntrials", ["timing", LAXPC_EA1, "--freq-min", "10.0", "--freq-max", "11.0", "--n-trials", "99999999"])
    check_adversarial("timing_nan_freq", ["timing", LAXPC_EA1, "--freq", "NaN"])
    check_adversarial("timing_inf_fdot", ["timing", LAXPC_EA1, "--freq", "10.0", "--fdot", "inf"])
    check_adversarial("czti_huge_oversampling", ["czti", "run", "--input", CZTI_L1_DIR, "--caldb", CZTI_CALDB, "--oversampling", "1000000", "-o", os.path.join(adv_dir, "bad_ov")])
    check_adversarial("czti_unpaired_ra", ["czti", "run", "--input", CZTI_L1_DIR, "--caldb", CZTI_CALDB, "--ra", "83.6331", "-o", os.path.join(adv_dir, "unpaired_ra")])

    return {
        "passed": True,
        "total_cases": len(cases),
        "panics": 0,
        "sample_cases": cases[:10],
    }


# ============================================================================
# Section 4: Cross-Instrument Consistency & Clean `env -i` Determinism
# ============================================================================

def test_section4_cross_instrument_and_clean_env(bin_path: str) -> Dict[str, Any]:
    # 1. Run `astrosat inspect --full` and `astrosat validate --strict` across ALL generated CZTI & LAXPC products
    product_roots = [
        os.path.join(WORK_ROOT, "wf1_czti", "run_default"),
        os.path.join(WORK_ROOT, "wf1_czti", "chain_default"),
        os.path.join(WORK_ROOT, "wf2_laxpc", "run_lxp1"),
        os.path.join(WORK_ROOT, "wf2_laxpc", "chain_lxp1"),
        os.path.join(WORK_ROOT, "sc2a_spectral", "l2"),
        os.path.join(WORK_ROOT, "sc2b_leap_pulsar", "l2"),
    ]
    all_fits: List[str] = []
    for r in product_roots:
        if os.path.exists(r):
            all_fits.extend(collect_fits(r))
    all_fits = sorted(set(all_fits))
    assert len(all_fits) >= 40, f"Expected >= 40 generated FITS files, found {len(all_fits)}"

    insp_res = require_ok(run_cli(bin_path, ["--json", "inspect", "--full", *all_fits]), "inspect --full all products")
    insp_json = json.loads(insp_res.stdout)
    bad_csum = [f["path"] for f in insp_json["files"] if not f["health"]["all_checksums_ok"]]
    # Check that NO sorted L2 event, rate, or orbit HDU has non-monotonic TIME
    # (exclude raw pre-bunchclean `s1_decode.evt`, and CZTI multi-quadrant aux tables VETOSPECTRUM/SSM/TEMP)
    monotone_critical = {"EVENTS", "LAXPC_EVENTS", "LAXPC_MODE_EA", "Q0", "Q1", "Q2", "Q3", "RATE", "ORBIT"}
    non_mono = [
        (f["path"], [h for h in f["health"]["non_monotonic_time_hdus"] if h.upper() in monotone_critical])
        for f in insp_json["files"]
        if not f["path"].endswith("s1_decode.evt")
    ]
    non_mono = [x for x in non_mono if x[1]]
    non_fin = [f["path"] for f in insp_json["files"] if f["health"]["non_finite_time_rows"] > 0]
    if bad_csum or non_mono or non_fin:
        raise AssertionError(f"Inspect health failures: bad_csum={bad_csum}, non_mono={non_mono}, non_fin={non_fin}")

    val_res = require_ok(run_cli(bin_path, ["--json", "validate", "--strict", *all_fits]), "validate --strict all products")
    val_json = json.loads(val_res.stdout)
    if not val_json["passed"] or val_json["total_errors"] != 0 or val_json["total_warnings"] != 0:
        raise AssertionError(f"validate --strict failed: errors={val_json['total_errors']}, warnings={val_json['total_warnings']}")

    # 2. Clean `env -i` execution (`PATH=/mnt/d/ISRO/.tmp/public_scenarios/bin:/usr/bin:/bin`, zero HEADAS/CALDB/HOME/TMPDIR)
    clean_env_dir = os.path.join(WORK_ROOT, "sec4_clean_env")
    if os.path.exists(clean_env_dir):
        shutil.rmtree(clean_env_dir)
    os.makedirs(os.path.join(clean_env_dir, "czti"), exist_ok=True)
    os.makedirs(os.path.join(clean_env_dir, "laxpc"), exist_ok=True)

    clean_env = {"PATH": f"{os.path.dirname(bin_path)}:/usr/bin:/bin"}
    require_ok(run_cli(bin_path, ["--help"], env=clean_env), "clean env --help")
    c_clean = json.loads(require_ok(run_cli(bin_path, [
        "--json", "czti", "run",
        "--input", CZTI_L1_DIR, "--caldb", CZTI_CALDB,
        "--output", os.path.join(clean_env_dir, "czti"),
    ], env=clean_env), "clean env czti run").stdout)
    l_clean = json.loads(require_ok(run_cli(bin_path, [
        "--json", "laxpc", "run",
        "--input", LAXPC_L1_DIR, "--units", "lxp1",
        "--output", os.path.join(clean_env_dir, "laxpc"),
    ], env=clean_env), "clean env laxpc run").stdout)

    bary_clean = os.path.join(clean_env_dir, "laxpc", "bary.evn")
    b_clean = json.loads(require_ok(run_cli(bin_path, [
        "--json", "bary", l_clean["unit_results"][0]["event_fits_path"],
        "--orbit", l_clean["orbit_fits_path"], "-o", bary_clean,
    ], env=clean_env), "clean env bary").stdout)
    t_clean = json.loads(require_ok(run_cli(bin_path, [
        "--json", "timing", bary_clean, "--freq-min", "29.65", "--freq-max", "29.665", "--fdot=-3.693e-10",
    ], env=clean_env), "clean env timing").stdout)
    require_ok(run_cli(bin_path, ["--json", "inspect", "--full", bary_clean], env=clean_env), "clean env inspect")
    require_ok(run_cli(bin_path, ["--json", "validate", "--strict", bary_clean], env=clean_env), "clean env validate")

    # Assert exact cross-environment determinism against normal-env wf1/wf2 outputs
    assert c_clean["stage2_clean_events"] == [43281, 44276, 41337, 37140], (
        f"clean env CZTI events drifted: {c_clean['stage2_clean_events']}"
    )
    assert l_clean["unit_results"][0]["gti_filtered_events"] == 8939817, (
        f"clean env LAXPC events drifted: {l_clean['unit_results'][0]['gti_filtered_events']}"
    )
    assert b_clean["events_corrected"] == 8939817, (
        f"clean env bary events drifted: {b_clean['events_corrected']}"
    )
    assert 29.650 <= t_clean["best_freq_hz"] <= 29.665, (
        f"clean env Crab frequency out of range: {t_clean['best_freq_hz']}"
    )

    return {
        "passed": True,
        "total_products_inspected_and_validated_strict": len(all_fits),
        "validate_strict_errors": val_json["total_errors"],
        "validate_strict_warnings": val_json["total_warnings"],
        "clean_env_czti_events": c_clean["stage2_clean_events"],
        "clean_env_laxpc_events": l_clean["unit_results"][0]["gti_filtered_events"],
        "clean_env_bary_events": b_clean["events_corrected"],
        "clean_env_crab_freq_hz": t_clean["best_freq_hz"],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="AstroSat Public-Release End-User Scenario Verification Suite")
    parser.add_argument("--astrosat-bin", default=DEFAULT_ASTROSAT_BIN, help="Path to standalone astrosat binary")
    parser.add_argument("--json-out", default=REPORT_JSON_DEFAULT, help="Path to write JSON summary report")
    parser.add_argument(
        "--only",
        nargs="*",
        choices=["wf1", "wf2", "wf3", "wf4", "sc2a", "sc2b", "sec3", "sec4"],
        help="Optional subset of scenarios to run",
    )
    args = parser.parse_args()

    if not os.path.isfile(args.astrosat_bin):
        print(f"ERROR: standalone binary not found at {args.astrosat_bin}", file=sys.stderr)
        return 2
    if not os.path.isfile(FITSVERIFY_BIN):
        print(f"ERROR: fitsverify not found at {FITSVERIFY_BIN}", file=sys.stderr)
        return 2

    os.makedirs(WORK_ROOT, exist_ok=True)
    selected = set(args.only) if args.only else {"wf1", "wf2", "wf3", "wf4", "sc2a", "sc2b", "sec3", "sec4"}

    tests = [
        ("wf1", "Workflow 1: Chained CZTI CLI Pipeline", test_workflow1_czti),
        ("wf2", "Workflow 2: Chained LAXPC CLI Pipeline", test_workflow2_laxpc),
        ("wf3", "Workflow 3: Declarative TOML/JSON Configs & Overrides", test_workflow3_configs_and_json),
        ("wf4", "Workflow 4: Idempotency & Guardrails", test_workflow4_idempotency_and_guardrails),
        ("sc2a", "Scenario 2A: Planted Spectral Power-Law & Gaussian Line", test_scenario2a_spectral_recovery),
        ("sc2b", "Scenario 2B: Planted 17.3456 Hz Pulsar Across 2016-2017 Leap Second + Orbit", test_scenario2b_leap_second_pulsar),
        ("sec3", "Section 3: Adversarial & Corrupt-Input Stress Testing", test_section3_adversarial),
        ("sec4", "Section 4: Cross-Instrument Strict Validation & Clean env -i", test_section4_cross_instrument_and_clean_env),
    ]

    results: Dict[str, Any] = {}
    all_passed = True
    for key, title, fn in tests:
        if key not in selected:
            continue
        t0 = time.monotonic()
        print(f"==> Running [{key}] {title} ...", flush=True)
        try:
            res = fn(args.astrosat_bin)
            res["elapsed_s"] = round(time.monotonic() - t0, 2)
            results[key] = res
            status = "PASS" if res.get("passed", False) else "FAIL"
            if not res.get("passed", False):
                all_passed = False
            print(f"[{status}] {key} ({res['elapsed_s']}s): {json.dumps(res)}", flush=True)
        except Exception as exc:
            all_passed = False
            results[key] = {
                "passed": False,
                "error": str(exc),
                "elapsed_s": round(time.monotonic() - t0, 2),
            }
            print(f"[FAIL] {key}: {exc}", file=sys.stderr, flush=True)

    summary = {"passed": all_passed, "astrosat_bin": args.astrosat_bin, "scenarios": results}
    os.makedirs(os.path.dirname(os.path.abspath(args.json_out)), exist_ok=True)
    with open(args.json_out, "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)
    print(f"\nReport written to {args.json_out} (overall passed={all_passed})")
    return 0 if all_passed else 1


if __name__ == "__main__":
    sys.exit(main())
