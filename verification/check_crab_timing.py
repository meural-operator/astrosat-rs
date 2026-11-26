"""Independent Crab pulsar (~29.66 Hz) timing verification on barycentered LAXPC 9000000252 events."""
from __future__ import annotations

import argparse
import json
import math
import sys
import numpy as np
from fits_reader import read_fits


def load_bary_events(path: str, hdu_name: str | None = None) -> tuple[np.ndarray, dict]:
    hdus = read_fits(path)
    target_hdu = None
    for h in hdus:
        if h.xtension != "BINTABLE":
            continue
        if hdu_name and h.extname.upper() == hdu_name.upper():
            target_hdu = h
            break
        cols = {c[0].upper() for c in h.bintable_columns()}
        if "TIME" in cols and h.extname.upper() not in ("GTI", "BTI"):
            target_hdu = h
            break
    if target_hdu is None:
        raise ValueError(f"No event BINTABLE with TIME column found in {path}")
    times = np.asarray(target_hdu.read_column("TIME"), dtype=np.float64)
    # Check for layer/channel filtering if present (optional, keeps SNR high)
    cols = {c[0].upper(): c[0] for c in target_hdu.bintable_columns()}
    if "LAYER" in cols:
        layers = target_hdu.read_column(cols["LAYER"])
        # Keep valid non-negative times
        mask = np.isfinite(times) & (layers >= 0)
        times = times[mask]
    else:
        times = times[np.isfinite(times)]
    times.sort()
    return times, target_hdu.header


def z2_n_stat(phases: np.ndarray, n_harmonics: int = 2) -> float:
    n = float(len(phases))
    if n == 0:
        return 0.0
    total = 0.0
    two_pi_phi = 2.0 * math.pi * phases
    for k in range(1, n_harmonics + 1):
        c = float(np.sum(np.cos(k * two_pi_phi)))
        s = float(np.sum(np.sin(k * two_pi_phi)))
        total += c * c + s * s
    return (2.0 / n) * total


def epoch_fold_chi2(phases: np.ndarray, nbins: int = 32) -> tuple[float, np.ndarray]:
    counts, _ = np.histogram(phases, bins=nbins, range=(0.0, 1.0))
    mean = float(np.mean(counts))
    if mean <= 0:
        return 0.0, counts
    chi2 = float(np.sum((counts - mean) ** 2 / mean))
    return chi2, counts


def search_crab(
    times: np.ndarray,
    f_min: float = 29.650,
    f_max: float = 29.665,
    fdot: float = -3.693e-10,
    n_freq: int = 600,
    nbins: int = 32,
) -> dict:
    if len(times) < 100:
        raise ValueError(f"Too few events ({len(times)}) for Crab timing check")
    t0 = float(times[0])
    dt = times - t0
    dt2 = 0.5 * fdot * (dt * dt)

    freqs = np.linspace(f_min, f_max, n_freq)
    best_f = float(freqs[0])
    best_chi2 = -1.0
    mean = float(len(times)) / float(nbins)

    # Fast C-level bincount epoch-folding scan over all events (~2s for 8.9M events)
    for f in freqs:
        phases = np.mod(f * dt + dt2, 1.0)
        bins = (phases * nbins).astype(np.int64)
        bins = np.clip(bins, 0, nbins - 1)
        counts = np.bincount(bins, minlength=nbins)
        chi2 = float(np.sum((counts - mean) ** 2 / mean))
        if chi2 > best_chi2:
            best_chi2 = chi2
            best_f = float(f)

    # Fine exact Z_2^2 refinement around peak on all events
    step = (f_max - f_min) / max(n_freq - 1, 1)
    fine_freqs = np.linspace(best_f - 2 * step, best_f + 2 * step, 41)
    best_z2 = -1.0
    for f in fine_freqs:
        phases = np.mod(f * dt + dt2, 1.0)
        z2 = z2_n_stat(phases, n_harmonics=2)
        if z2 > best_z2:
            best_z2 = z2
            best_f = float(f)

    best_phases = np.mod(best_f * dt + dt2, 1.0)
    chi2, profile = epoch_fold_chi2(best_phases, nbins=nbins)
    p_max = int(np.max(profile))
    p_min = int(np.min(profile))
    p_mean = float(np.mean(profile))
    pulsed_frac = (p_max - p_min) / (p_max + p_min) if (p_max + p_min) > 0 else 0.0

    return {
        "n_events": int(len(times)),
        "t0": t0,
        "duration_s": float(times[-1] - times[0]),
        "best_freq_hz": best_f,
        "fdot_hz_s": fdot,
        "z2_2": best_z2,
        "chi2_32bins": chi2,
        "profile_min": p_min,
        "profile_max": p_max,
        "profile_mean": p_mean,
        "pulsed_amplitude_ratio": pulsed_frac,
        "profile": profile.tolist(),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Independent Crab pulsar timing check")
    parser.add_argument("bary_evt", help="Barycentered LAXPC event FITS file")
    parser.add_argument("--hdu", default=None, help="Optional event HDU name")
    parser.add_argument("--min-z2", type=float, default=1000.0, help="Minimum Z_2^2 threshold for detection")
    parser.add_argument("--json-out", default=None, help="Optional JSON output path")
    args = parser.parse_args()

    times, hdr = load_bary_events(args.bary_evt, args.hdu)
    timesys = str(hdr.get("TIMESYS", "")).strip().upper()
    timeref = str(hdr.get("TIMEREF", "")).strip().upper()
    print(f"Loaded {len(times)} events from {args.bary_evt} (TIMESYS={timesys!r}, TIMEREF={timeref!r})")

    errors = []
    if timesys != "TDB":
        errors.append(f"Expected TIMESYS='TDB' in barycentered file, got {timesys!r}")
    if timeref != "SOLARSYSTEM":
        errors.append(f"Expected TIMEREF='SOLARSYSTEM' in barycentered file, got {timeref!r}")

    res = search_crab(times)
    res["timesys"] = timesys
    res["timeref"] = timeref

    print(
        f"Crab search result: f={res['best_freq_hz']:.6f} Hz, "
        f"Z_2^2={res['z2_2']:.2f}, chi2(32 bins)={res['chi2_32bins']:.2f}, "
        f"profile min/mean/max={res['profile_min']}/{res['profile_mean']:.1f}/{res['profile_max']}"
    )

    # Render ASCII pulse profile (2 cycles)
    prof = np.array(res["profile"], dtype=float)
    lo, hi = prof.min(), prof.max()
    if hi > lo:
        bars = " .:-=+*#%@"
        line = "".join(bars[min(9, int(9.99 * (v - lo) / (hi - lo)))] for v in prof)
        print(f"Pulse profile (2 cycles): |{line}{line}|")

    if not (29.655 <= res["best_freq_hz"] <= 29.662):
        errors.append(f"Peak frequency {res['best_freq_hz']:.6f} Hz outside expected Jan-2016 Crab [29.655, 29.662] Hz window")
    if res["z2_2"] < args.min_z2:
        errors.append(f"Z_2^2={res['z2_2']:.2f} below detection threshold {args.min_z2}")
    if res["chi2_32bins"] < 1000.0:
        errors.append(f"chi2(32 bins)={res['chi2_32bins']:.2f} below detection threshold 1000.0")
    if res["pulsed_amplitude_ratio"] < 0.05:
        errors.append(f"Pulsed amplitude ratio {res['pulsed_amplitude_ratio']:.4f} < 0.05")

    if args.json_out:
        with open(args.json_out, "w", encoding="utf-8") as fh:
            json.dump({"passed": len(errors) == 0, "errors": errors, **res}, fh, indent=2)

    if errors:
        for e in errors:
            print(f"FAIL: {e}", file=sys.stderr)
        return 1
    print("PASS: Clear Crab pulsar detection at ~29.66 Hz on barycentered events!")
    return 0


if __name__ == "__main__":
    sys.exit(main())
