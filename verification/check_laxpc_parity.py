"""Compare generated LAXPC L2 outputs for Crab ObsID 9000000252 against official L2 reference."""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys
import numpy as np
from fits_reader import read_fits

REF_L2_DIR = "/mnt/d/ISRO/data/sample/laxpc_l2/20160107_G02_010T01_9000000252_level2/laxpc/01509"


def find_first(out_dir: str, patterns: list[str]) -> str | None:
    for pat in patterns:
        matches = sorted(glob.glob(os.path.join(out_dir, "**", pat), recursive=True))
        if matches:
            return matches[0]
    return None


def pearson_r(a: np.ndarray, b: np.ndarray) -> float:
    if len(a) != len(b) or len(a) == 0:
        return 0.0
    da = a - np.mean(a)
    db = b - np.mean(b)
    denom = float(np.sqrt(np.sum(da * da) * np.sum(db * db)))
    return float(np.sum(da * db) / denom) if denom > 0 else 0.0


def _check_bary_laxpc(
    bary_path: str,
    orig_path: str | None,
    label: str,
    expected_plephem: str | None,
    errors: list[str],
) -> dict:
    info: dict = {}
    if not os.path.exists(bary_path):
        errors.append(f"Missing {label} at {bary_path}")
        return info
    bh = read_fits(bary_path)
    ev_h = next((h for h in bh if h.xtension == "BINTABLE" and h.extname.upper() != "GTI"), None)
    if ev_h is None:
        errors.append(f"{label} missing event BINTABLE")
        return info
    tsys = str(ev_h.header.get("TIMESYS", "")).strip().upper()
    tref = str(ev_h.header.get("TIMEREF", "")).strip().upper()
    plephem = str(ev_h.header.get("PLEPHEM", "")).strip().upper()
    if tsys != "TDB" or tref != "SOLARSYSTEM":
        errors.append(f"{label} has TIMESYS={tsys!r}, TIMEREF={tref!r} (expected TDB/SOLARSYSTEM)")
    if expected_plephem and plephem != expected_plephem.upper():
        errors.append(f"{label} has PLEPHEM={plephem!r} (expected {expected_plephem!r})")
    b_times = np.asarray(ev_h.read_column("TIME"), dtype=np.float64)
    if len(b_times) == 0 or not np.all(np.isfinite(b_times)):
        errors.append(f"{label} has empty or non-finite TIME column")
        return info
    info["events"] = int(len(b_times))
    info["plephem"] = plephem
    if orig_path and os.path.exists(orig_path):
        oh = next((h for h in read_fits(orig_path) if h.xtension == "BINTABLE" and h.extname.upper() != "GTI"), None)
        if oh is not None:
            o_times = np.asarray(oh.read_column("TIME"), dtype=np.float64)
            if len(b_times) != len(o_times):
                errors.append(f"{label} event count {len(b_times)} != original {len(o_times)}")
            else:
                med_shift = float(np.median(b_times[:10000] - o_times[:10000]))
                info["median_shift_s"] = med_shift
                if not (10.0 <= abs(med_shift) <= 600.0):
                    errors.append(f"{label} unexpected barycentric shift {med_shift:.3f} s")
    return info


def check_laxpc(
    out_dir: str,
    bary_evn: str | None = None,
    bary_fmt_a: str | None = None,
    bary_de405_evn: str | None = None,
    parity_dir: str | None = None,
) -> tuple[bool, list[str], dict]:
    errors: list[str] = []
    metrics: dict = {"units": {}, "parity_1024ch": {}}

    unit_event_total = 0
    lxp1_evn_path: str | None = None

    # 1. Check auto-GTI Format B outputs for all 3 units
    for unit in ("lxp1", "lxp2", "lxp3"):
        unit_m: dict = {}
        evn_p = find_first(out_dir, [f"*{unit}*.evn", f"*{unit}*EA*.evn", f"*{unit}*.evt"])
        lc_p = find_first(out_dir, [f"*{unit}*.lc", f"*{unit}*.lcurv"])
        pha_p = find_first(out_dir, [f"*{unit}*.pha", f"*{unit}*.spec"])
        gti_p = find_first(out_dir, [f"*{unit}*.gti", f"*{unit}*.gtis"])
        if unit == "lxp1":
            lxp1_evn_path = evn_p

        if not evn_p:
            errors.append(f"{unit}: missing generated .evn file")
        else:
            evn_hdus = read_fits(evn_p)
            ev_hdu = next((h for h in evn_hdus if h.xtension == "BINTABLE" and h.extname.upper() not in ("GTI", "BTI")), None)
            gti_hdu = next((h for h in evn_hdus if h.extname.upper() == "GTI"), None)
            if ev_hdu is None:
                errors.append(f"{unit}: .evn missing event BINTABLE")
            else:
                n_ev = int(ev_hdu.header.get("NAXIS2", 0))
                unit_m["gen_events"] = n_ev
                unit_event_total += n_ev
                if n_ev < 7_000_000:
                    errors.append(f"{unit}: unexpectedly low event count {n_ev} on Crab 9000000252")
            if gti_hdu is None:
                errors.append(f"{unit}: .evn missing GTI HDU")
            else:
                gen_exp = float(np.sum(gti_hdu.read_column("STOP") - gti_hdu.read_column("START")))
                unit_m["gen_exposure"] = gen_exp
                if gen_exp < 2500.0 or gen_exp > 4000.0:
                    errors.append(f"{unit}: auto-GTI exposure {gen_exp:.2f}s outside expected [2500, 4000]s window")

        if not lc_p:
            errors.append(f"{unit}: missing generated .lc/.lcurv")
        else:
            lc_h = next((h for h in read_fits(lc_p) if h.xtension == "BINTABLE"), None)
            if lc_h is None:
                errors.append(f"{unit}: .lc missing BINTABLE")
            else:
                g_rate = np.asarray(lc_h.read_column("RATE"), dtype=np.float64)
                cols = {c[0].upper() for c in lc_h.bintable_columns()}
                if "FRACEXP" in cols:
                    g_frac = np.asarray(lc_h.read_column("FRACEXP"), dtype=np.float64)
                    valid = np.isfinite(g_rate) & (g_frac > 0)
                else:
                    valid = np.isfinite(g_rate)
                g_mean = float(np.mean(g_rate[valid])) if np.any(valid) else 0.0
                unit_m["gen_lc_mean_rate_cps"] = g_mean
                if g_mean < 2000.0 or g_mean > 5000.0:
                    errors.append(f"{unit}: Crab mean LC rate {g_mean:.1f} cps outside expected [2000, 5000] cps")

        if not pha_p:
            errors.append(f"{unit}: missing generated .pha/.spec")
        else:
            pha_h = next((h for h in read_fits(pha_p) if h.extname.upper() == "SPECTRUM"), None)
            if pha_h is None:
                errors.append(f"{unit}: .pha missing SPECTRUM HDU")
            else:
                cols = {c[0].upper(): c[0] for c in pha_h.bintable_columns()}
                if "COUNTS" in cols:
                    g_cnt = float(np.sum(pha_h.read_column(cols["COUNTS"])))
                    unit_m["gen_pha_sum"] = g_cnt
                    if g_cnt <= 0:
                        errors.append(f"{unit}: PHA COUNTS sum <= 0")
                elif "RATE" in cols:
                    g_sum = float(np.sum(pha_h.read_column(cols["RATE"])))
                    unit_m["gen_pha_rate_sum"] = g_sum
                    if g_sum < 2000.0:
                        errors.append(f"{unit}: PHA RATE sum {g_sum:.1f} cps unexpectedly low")

        if not gti_p:
            errors.append(f"{unit}: missing standalone .gti/.gtis file")
        else:
            gh = next((h for h in read_fits(gti_p) if h.extname.upper() == "GTI"), None)
            if gh is None or int(gh.header.get("NAXIS2", 0)) == 0:
                errors.append(f"{unit}: standalone .gtis missing non-empty GTI HDU")

        metrics["units"][unit] = unit_m

    # 2. Apples-to-apples parity run against official L2 reference (official GTI & 1024 ch, all 3 units)
    if parity_dir and os.path.isdir(parity_dir):
        for unit in ("lxp1", "lxp2", "lxp3"):
            ref_pha_p = os.path.join(REF_L2_DIR, unit, "modeEA", f"AS1G02_010T01_9000000252{unit}EA_level2.pha")
            ref_gti_p = os.path.join(REF_L2_DIR, unit, "modeEA", f"AS1G02_010T01_9000000252{unit}EA_level2.gti")
            ref_gti = read_fits(ref_gti_p)[1]
            ref_exp = float(np.sum(ref_gti.read_column("STOP") - ref_gti.read_column("START")))
            ref_pha = read_fits(ref_pha_p)[1]
            ref_counts = np.asarray(ref_pha.read_column("COUNTS"), dtype=np.float64)

            gen_pha_p = find_first(parity_dir, [f"*{unit}*.pha", f"*{unit}*.spec"])
            if not gen_pha_p:
                errors.append(f"Parity dir missing {unit} .pha/.spec")
            else:
                pha_h = next((h for h in read_fits(gen_pha_p) if h.extname.upper() == "SPECTRUM"), None)
                if pha_h is None:
                    errors.append(f"Parity {unit} .pha missing SPECTRUM HDU")
                else:
                    cols = {c[0].upper(): c[0] for c in pha_h.bintable_columns()}
                    gen_exp = float(pha_h.header.get("EXPOSURE", 0.0))
                    if "RATE" in cols:
                        gen_vec = np.asarray(pha_h.read_column(cols["RATE"]), dtype=np.float64)
                    else:
                        gen_vec = np.asarray(pha_h.read_column(cols["COUNTS"]), dtype=np.float64)
                    r = pearson_r(gen_vec, ref_counts)
                    exp_diff = abs(gen_exp - ref_exp) / max(ref_exp, 1.0)
                    metrics["parity_1024ch"][unit] = {
                        "ref_exposure": ref_exp,
                        "gen_exposure": gen_exp,
                        "exp_rel_diff": exp_diff,
                        "pearson_r": r,
                    }
                    if exp_diff > 0.01:
                        errors.append(f"Parity {unit}: exposure {gen_exp:.2f}s vs ref {ref_exp:.2f}s (>{exp_diff:.2%})")
                    if r < 0.999:
                        errors.append(f"Parity {unit}: 1024-ch spectrum Pearson r={r:.6f} < 0.999")

    # 3. Format A combined level2.event.fits
    fmt_a = find_first(out_dir, ["level2.event.fits", "*level2.event.fits"])
    if not fmt_a:
        errors.append("Missing Format A combined level2.event.fits")
    else:
        fa_hdus = read_fits(fmt_a)
        ev_h = next((h for h in fa_hdus if h.xtension == "BINTABLE" and h.extname.upper() != "GTI"), None)
        if ev_h is None:
            errors.append("Format A level2.event.fits missing event BINTABLE")
        else:
            cols = {c[0].upper() for c in ev_h.bintable_columns()}
            for req in ("TIME", "CHANNEL", "LAYER", "LAXPC_NO"):
                if req not in cols:
                    errors.append(f"Format A level2.event.fits missing column {req}")
            fa_n = int(ev_h.header.get("NAXIS2", 0))
            metrics["format_a_events"] = fa_n
            if unit_event_total > 0 and fa_n != unit_event_total:
                errors.append(f"Format A event count {fa_n} != sum of lxp1..3 ({unit_event_total})")

    # 4. AstroSat.orb
    orb_p = find_first(out_dir, ["AstroSat.orb", "*.orb"])
    if not orb_p:
        errors.append("Missing generated AstroSat.orb")
    else:
        orb_hdus = read_fits(orb_p)
        oh = next((h for h in orb_hdus if h.xtension == "BINTABLE"), None)
        if oh is None or int(oh.header.get("NAXIS2", 0)) == 0:
            errors.append("AstroSat.orb missing non-empty ORBIT BINTABLE")
        else:
            ocols = {c[0].upper() for c in oh.bintable_columns()}
            for req in ("TIME", "X", "Y", "Z", "VX", "VY", "VZ"):
                if req not in ocols:
                    errors.append(f"AstroSat.orb missing column {req}")
            metrics["orbit_rows"] = int(oh.header.get("NAXIS2", 0))

    # 5. Barycentered .evn (analytical + Format A + optional JPLEPH.405)
    if bary_evn:
        metrics["bary_evn"] = _check_bary_laxpc(bary_evn, lxp1_evn_path, "analytical bary .evn", "ASTROSAT-ANALYTICAL", errors)
    if bary_fmt_a:
        metrics["bary_fmt_a"] = _check_bary_laxpc(bary_fmt_a, fmt_a, "analytical bary Format A", "ASTROSAT-ANALYTICAL", errors)
    if bary_de405_evn:
        metrics["bary_de405_evn"] = _check_bary_laxpc(bary_de405_evn, lxp1_evn_path, "JPLEPH.405 bary .evn", "JPL-DE405", errors)

    return len(errors) == 0, errors, metrics


def main() -> int:
    parser = argparse.ArgumentParser(description="LAXPC Crab 9000000252 physics & L2 parity check")
    parser.add_argument("out_dir", help="Directory containing generated LAXPC L2 products")
    parser.add_argument("--bary-evn", default=None, help="Optional analytical barycentered LAXPC .evn path")
    parser.add_argument("--bary-fmt-a", default=None, help="Optional barycentered Format A level2.event.bary.fits path")
    parser.add_argument("--bary-de405-evn", default=None, help="Optional JPLEPH.405 barycentered LAXPC .evn path")
    parser.add_argument("--parity-dir", default=None, help="Optional 1024-ch official-GTI parity run directory")
    parser.add_argument("--json-out", default=None, help="Optional JSON summary output path")
    args = parser.parse_args()

    ok, errors, metrics = check_laxpc(args.out_dir, args.bary_evn, args.bary_fmt_a, args.bary_de405_evn, args.parity_dir)
    print(json.dumps({"passed": ok, "errors": errors, "metrics": metrics}, indent=2))
    if args.json_out:
        with open(args.json_out, "w", encoding="utf-8") as fh:
            json.dump({"passed": ok, "errors": errors, "metrics": metrics}, fh, indent=2)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
