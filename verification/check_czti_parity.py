"""Compare generated CZTI L2 outputs for ObsID 9000000406 against official L2 reference."""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys
import numpy as np
from fits_reader import read_fits

REF_DIR = "/mnt/d/ISRO/data/sample/czti/czti_sample_data/20160331_T01_112T01_9000000406_level2/modeM0"


def find_first(out_dir: str, patterns: list[str]) -> str | None:
    for pat in patterns:
        matches = sorted(glob.glob(os.path.join(out_dir, "**", pat), recursive=True))
        if matches:
            return matches[0]
    return None


def _get_asp_vals(hdus):
    for h in hdus:
        if "RA_PNT" in h.header and "DEC_PNT" in h.header:
            ra = float(h.header["RA_PNT"])
            dec = float(h.header["DEC_PNT"])
            tw = float(h.header.get("TWIST", 0.0))
            return ra, dec, tw
    asp_h = next((h for h in hdus if h.xtension == "BINTABLE" and int(h.header.get("NAXIS2", 0)) > 0), None)
    if asp_h is not None:
        cols = {c[0].upper(): c[0] for c in asp_h.bintable_columns()}
        if {"NX", "NY", "NZ"}.issubset(cols):
            nx = float(np.mean(asp_h.read_column(cols["NX"])))
            ny = float(np.mean(asp_h.read_column(cols["NY"])))
            nz = float(np.mean(asp_h.read_column(cols["NZ"])))
            ra = float(np.degrees(np.arctan2(ny, nx)) % 360.0)
            dec = float(np.degrees(np.arcsin(np.clip(nz / max(np.sqrt(nx * nx + ny * ny + nz * nz), 1e-12), -1.0, 1.0))))
            tw = 0.0
            if {"NXT", "NYT", "NZT"}.issubset(cols):
                nxt = float(np.mean(asp_h.read_column(cols["NXT"])))
                nyt = float(np.mean(asp_h.read_column(cols["NYT"])))
                tw = float(np.degrees(np.arctan2(nyt, nxt)))
            return ra, dec, tw
    return None


def _check_bary_czti(
    bary_path: str,
    orig_evt_hdus,
    label: str,
    expected_plephem: str | None,
    errors: list[str],
) -> dict:
    info: dict = {}
    if not os.path.exists(bary_path):
        errors.append(f"Missing {label} CZTI .evt at {bary_path}")
        return info
    bh = read_fits(bary_path)
    for q in ("Q0", "Q1", "Q2", "Q3"):
        qh = next((h for h in bh if h.extname == q), None)
        oh = next((h for h in orig_evt_hdus if h.extname == q), None)
        if qh is None:
            errors.append(f"{label} CZTI .evt missing {q} HDU")
            continue
        tsys = str(qh.header.get("TIMESYS", "")).strip().upper()
        tref = str(qh.header.get("TIMEREF", "")).strip().upper()
        plephem = str(qh.header.get("PLEPHEM", "")).strip().upper()
        if tsys != "TDB" or tref != "SOLARSYSTEM":
            errors.append(f"{label} CZTI {q} has TIMESYS={tsys!r}, TIMEREF={tref!r} (expected TDB/SOLARSYSTEM)")
        if expected_plephem and plephem != expected_plephem.upper():
            errors.append(f"{label} CZTI {q} has PLEPHEM={plephem!r} (expected {expected_plephem!r})")
        b_times = np.asarray(qh.read_column("Time"), dtype=np.float64)
        if not np.all(np.isfinite(b_times)):
            errors.append(f"{label} CZTI {q} contains non-finite Time values")
        if oh is not None:
            o_times = np.asarray(oh.read_column("Time"), dtype=np.float64)
            if len(b_times) != len(o_times):
                errors.append(f"{label} CZTI {q} row count {len(b_times)} != original {len(o_times)}")
            elif len(b_times) > 0:
                med_shift = float(np.median(b_times - o_times))
                info[q] = {"median_shift_s": med_shift, "plephem": plephem}
                if not (10.0 <= abs(med_shift) <= 600.0):
                    errors.append(f"{label} CZTI {q} unexpected median bary shift {med_shift:.3f} s")
    return info


def check_czti(
    out_dir: str,
    bary_evt: str | None = None,
    bary_de405_evt: str | None = None,
) -> tuple[bool, list[str], dict]:
    errors: list[str] = []
    metrics: dict = {}

    # 1a. Stage 1 bunchclean event parity (*_bc.evt)
    ref_bc_path = os.path.join(REF_DIR, "AS1T01_112T01_9000000406cztM0_level2_bc.evt")
    gen_bc_path = find_first(out_dir, ["*_bc.evt"])
    if not gen_bc_path:
        errors.append(f"Missing generated CZTI Stage 1 *_bc.evt in {out_dir}")
    else:
        ref_bc = read_fits(ref_bc_path)
        gen_bc = read_fits(gen_bc_path)
        metrics["bc_quadrants"] = {}
        for q in ("Q0", "Q1", "Q2", "Q3"):
            rh = next((h for h in ref_bc if h.extname == q), None)
            gh = next((h for h in gen_bc if h.extname == q), None)
            if gh is None or rh is None:
                errors.append(f"Stage 1 bc.evt missing HDU {q}")
                continue
            r_n = int(rh.header["NAXIS2"])
            g_n = int(gh.header["NAXIS2"])
            rel = abs(g_n - r_n) / max(r_n, 1)
            metrics["bc_quadrants"][q] = {"ref_events": r_n, "gen_events": g_n, "rel_diff": rel}
            if rel > 0.015:
                errors.append(f"Stage 1 bc.evt {q}: count {g_n} deviates >1.5% from ref {r_n} ({rel:.2%})")

    # 1b. Stage 2 clean event file
    ref_evt_path = os.path.join(REF_DIR, "AS1T01_112T01_9000000406czt_level2_quad_clean.evt")
    gen_evt_path = find_first(out_dir, ["*quad_clean.evt", "*_clean.evt", "*.evt"])
    if not gen_evt_path:
        errors.append(f"Missing generated CZTI .evt in {out_dir}")
        return False, errors, metrics

    ref_evt = read_fits(ref_evt_path)
    gen_evt = read_fits(gen_evt_path)
    metrics["evt_file"] = gen_evt_path
    metrics["quadrants"] = {}

    for q in ("Q0", "Q1", "Q2", "Q3"):
        rh = next((h for h in ref_evt if h.extname == q), None)
        gh = next((h for h in gen_evt if h.extname == q), None)
        if gh is None:
            errors.append(f"Generated CZTI .evt missing HDU {q}")
            continue
        r_n = int(rh.header["NAXIS2"]) if rh else 0
        g_n = int(gh.header["NAXIS2"])
        g_pi = np.asarray(gh.read_column("PI"), dtype=np.float64)
        g_en = np.asarray(gh.read_column("ENERGY"), dtype=np.float64)
        r_pi = np.asarray(rh.read_column("PI"), dtype=np.float64)
        r_en = np.asarray(rh.read_column("ENERGY"), dtype=np.float64)

        g_detx = np.asarray(gh.read_column("DETX"))
        g_dety = np.asarray(gh.read_column("DETY"))
        g_veto = np.asarray(gh.read_column("veto"))
        g_alpha = np.asarray(gh.read_column("alpha"))
        if g_n > 0:
            if np.min(g_detx) < 0 or np.max(g_detx) > 63 or np.min(g_dety) < 0 or np.max(g_dety) > 63:
                errors.append(f"{q}: DETX/DETY out of [0, 63] range")
            if np.any(g_veto != 0):
                errors.append(f"{q}: clean events contain non-zero veto flags ({int(np.sum(g_veto != 0))} rows)")
            if np.any(g_alpha != 0):
                errors.append(f"{q}: clean events contain non-zero alpha flags ({int(np.sum(g_alpha != 0))} rows)")
            if not np.all(np.isfinite(g_en)) or np.min(g_pi) < 0 or np.max(g_pi) >= 512:
                errors.append(f"{q}: PI/ENERGY out of valid bounds")

        rg = next((h for h in ref_evt if h.extname == f"{q}_GTI"), None)
        gg = next((h for h in gen_evt if h.extname == f"{q}_GTI"), None)
        if gg is None:
            errors.append(f"Generated CZTI .evt missing {q}_GTI HDU")
            g_exp = 0.0
        else:
            g_exp = float(np.sum(gg.read_column("STOP") - gg.read_column("START")))
        r_exp = float(np.sum(rg.read_column("STOP") - rg.read_column("START"))) if rg else 0.0

        rel_n = abs(g_n - r_n) / max(r_n, 1)
        rel_exp = abs(g_exp - r_exp) / max(r_exp, 1e-6)
        pi_diff = abs(float(np.mean(g_pi)) - float(np.mean(r_pi))) if g_n > 0 else 999.0
        en_diff = abs(float(np.mean(g_en)) - float(np.mean(r_en))) if g_n > 0 else 999.0

        metrics["quadrants"][q] = {
            "ref_events": r_n,
            "gen_events": g_n,
            "rel_events_diff": rel_n,
            "ref_exposure": r_exp,
            "gen_exposure": g_exp,
            "rel_exposure_diff": rel_exp,
            "ref_mean_pi": float(np.mean(r_pi)),
            "gen_mean_pi": float(np.mean(g_pi)) if g_n > 0 else 0.0,
            "ref_mean_energy": float(np.mean(r_en)),
            "gen_mean_energy": float(np.mean(g_en)) if g_n > 0 else 0.0,
        }
        if g_n == 0:
            errors.append(f"{q}: 0 clean events produced")
        elif rel_n > 0.18:
            errors.append(f"{q}: event count {g_n} deviates >18% from reference {r_n} ({rel_n:.1%})")
        if rel_exp > 0.02:
            errors.append(f"{q}: GTI exposure {g_exp:.2f}s deviates >2% from reference {r_exp:.2f}s")
        if pi_diff > 20.0 or en_diff > 10.0:
            errors.append(f"{q}: mean PI/ENERGY ({np.mean(g_pi):.1f}, {np.mean(g_en):.1f} keV) deviates from ref ({np.mean(r_pi):.1f}, {np.mean(r_en):.1f} keV)")

    # 2. Aspect files (.aspect_Q0..Q3)
    metrics["aspect"] = {}
    for q in ("Q0", "Q1", "Q2", "Q3"):
        ref_asp_p = os.path.join(REF_DIR, f"AS1T01_112T01_9000000406cztM0_level2.aspect_{q}")
        gen_asp_p = find_first(out_dir, [f"*.aspect_{q}", f"*aspect*{q}*"])
        if not gen_asp_p:
            errors.append(f"Missing generated aspect file for {q}")
            continue
        r_v = _get_asp_vals(read_fits(ref_asp_p))
        g_v = _get_asp_vals(read_fits(gen_asp_p))
        if not g_v:
            errors.append(f"Generated aspect_{q} missing RA/DEC/TWIST")
        else:
            metrics["aspect"][q] = {"ref": r_v, "gen": g_v}
            if r_v:
                if abs(g_v[0] - r_v[0]) > 0.25 or abs(g_v[1] - r_v[1]) > 0.25:
                    errors.append(f"aspect_{q}: RA/DEC {g_v[:2]} differs >0.25 deg from ref {r_v[:2]}")
                if abs(g_v[2] - r_v[2]) > 1.0:
                    errors.append(f"aspect_{q}: TWIST {g_v[2]:.3f} differs >1.0 deg from ref {r_v[2]:.3f}")

    # 3. DPH / DPI
    ref_dph_p = os.path.join(REF_DIR, "AS1T01_112T01_9000000406cztM0_level2_quad_clean.dph")
    gen_dph_p = find_first(out_dir, ["*.dph"])
    if not gen_dph_p:
        errors.append("Missing generated .dph")
    else:
        ref_dph = read_fits(ref_dph_p)
        gen_dph = read_fits(gen_dph_p)
        metrics["dph"] = {}
        for q in ("Q0", "Q1", "Q2", "Q3"):
            rh = next((h for h in ref_dph if h.extname == q), None)
            gh = next((h for h in gen_dph if h.extname == q), None)
            if gh is None:
                errors.append(f".dph missing HDU {q}")
                continue
            ga = gh.image_array()
            ra = rh.image_array() if rh else None
            if ga is None or ga.shape != (64, 64):
                errors.append(f".dph {q} has unexpected shape {None if ga is None else ga.shape}")
            else:
                g_sum = float(np.sum(ga))
                r_sum = float(np.sum(ra)) if ra is not None else 0.0
                metrics["dph"][q] = {"ref_sum": r_sum, "gen_sum": g_sum}
                if np.min(ga) < 0 or g_sum <= 10_000.0:
                    errors.append(f".dph {q} has invalid counts (min={float(np.min(ga))}, sum={g_sum})")

    gen_dpi_p = find_first(out_dir, ["*.dpi"])
    if not gen_dpi_p:
        errors.append("Missing generated .dpi")
    else:
        gen_dpi = read_fits(gen_dpi_p)
        full_dpi = next((h for h in gen_dpi if "FULL" in h.extname.upper()), None)
        if full_dpi is None or full_dpi.image_array() is None:
            errors.append(".dpi missing FULL_DPI IMAGE HDU")
        else:
            arr = full_dpi.image_array()
            if arr.shape != (128, 128) or not np.all(np.isfinite(arr)) or np.min(arr) < 0.0 or float(np.sum(arr)) <= 0.0:
                errors.append(f".dpi FULL_DPI invalid: shape={arr.shape}, sum={float(np.sum(arr))}")
            metrics["full_dpi_sum"] = float(np.sum(arr))

    # 4. Sky images (.img_Q0..Q3)
    img_files = sorted(glob.glob(os.path.join(out_dir, "**", "*.img_Q*"), recursive=True))
    if len(img_files) < 4:
        errors.append(f"Expected 4 quadrant .img_Q* files, found {len(img_files)}")
    else:
        metrics["img_files"] = len(img_files)
        for imf in img_files:
            hdus = read_fits(imf)
            img_h = next((h for h in hdus if h.image_array() is not None and h.image_array().size > 0), None)
            if img_h is None:
                errors.append(f"Image {imf} has no non-empty IMAGE HDU")
                continue
            arr = img_h.image_array()
            if not np.all(np.isfinite(arr)) or float(np.std(arr)) == 0.0:
                errors.append(f"Image {imf} is non-finite or constant")
            crval1 = img_h.header.get("CRVAL1")
            crval2 = img_h.header.get("CRVAL2")
            if crval1 is None or crval2 is None or abs(float(crval1) - 83.633) > 1.0 or abs(float(crval2) - 22.014) > 1.0:
                errors.append(f"Image {imf} WCS CRVAL1/2=({crval1}, {crval2}) off target Crab")

    # 5. Mask-weighted .pha, .lc, .rsp across ALL 4 quadrants Q0..Q3
    metrics["products"] = {}
    for q in ("Q0", "Q1", "Q2", "Q3"):
        q_m: dict = {}
        for ext, expected_hdus in ((".pha", ("SPECTRUM",)), (".lc", ("RATE", "LIGHTCURVE")), (".rsp", ("MATRIX", "SPECRESP MATRIX"))):
            p = find_first(out_dir, [f"*_{q}{ext}", f"*{q}*{ext}"])
            if not p:
                errors.append(f"Missing generated CZTI {q} {ext} product")
                continue
            hdus = read_fits(p)
            h = next((x for x in hdus if any(eh in x.extname.upper() for eh in expected_hdus)), None)
            if h is None:
                errors.append(f"{p} missing {expected_hdus} HDU")
                continue
            nrows = int(h.header.get("NAXIS2", 0))
            if nrows == 0:
                errors.append(f"{p} {expected_hdus} has 0 rows")
                continue
            if ext == ".lc":
                rate = np.asarray(h.read_column("RATE"), dtype=np.float64)
                frac = np.asarray(h.read_column("FRACEXP"), dtype=np.float64)
                if not np.all(np.isfinite(rate)) or np.any(frac <= 0.0) or np.any(frac > 1.0001):
                    errors.append(f"{p}: invalid RATE/FRACEXP values")
                q_m["lc_mean_rate"] = float(np.mean(rate))
            elif ext == ".pha":
                col = "RATE" if any(c[0].upper() == "RATE" for c in h.bintable_columns()) else "COUNTS"
                vals = np.asarray(h.read_column(col), dtype=np.float64)
                if not np.all(np.isfinite(vals)):
                    errors.append(f"{p}: non-finite values in {col}")
                q_m["pha_sum"] = float(np.sum(vals))
            elif ext == ".rsp":
                eb = next((x for x in hdus if x.extname.upper() == "EBOUNDS"), None)
                if eb is None or int(eb.header.get("NAXIS2", 0)) != 512:
                    errors.append(f"{p}: missing 512-row EBOUNDS HDU")
                q_m["rsp_ebins"] = nrows
        metrics["products"][q] = q_m

    # 6. Barycentered .evt (analytical + optional JPLEPH.405)
    if bary_evt:
        metrics["bary_analytical"] = _check_bary_czti(bary_evt, gen_evt, "analytical bary", "ASTROSAT-ANALYTICAL", errors)
    if bary_de405_evt:
        metrics["bary_de405"] = _check_bary_czti(bary_de405_evt, gen_evt, "JPLEPH.405 bary", "JPL-DE405", errors)

    return len(errors) == 0, errors, metrics


def main() -> int:
    parser = argparse.ArgumentParser(description="CZTI 9000000406 physics & L2 parity check")
    parser.add_argument("out_dir", help="Directory containing generated CZTI L2 products")
    parser.add_argument("--bary-evt", default=None, help="Optional analytical barycentered CZTI .evt path")
    parser.add_argument("--bary-de405-evt", default=None, help="Optional JPLEPH.405 barycentered CZTI .evt path")
    parser.add_argument("--json-out", default=None, help="Optional JSON summary output path")
    args = parser.parse_args()

    ok, errors, metrics = check_czti(args.out_dir, args.bary_evt, args.bary_de405_evt)
    print(json.dumps({"passed": ok, "errors": errors, "metrics": metrics}, indent=2))
    if args.json_out:
        with open(args.json_out, "w", encoding="utf-8") as fh:
            json.dump({"passed": ok, "errors": errors, "metrics": metrics}, fh, indent=2)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
