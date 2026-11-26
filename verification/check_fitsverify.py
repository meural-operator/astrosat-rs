"""Run fitsverify on all FITS files under given paths and require 0 errors."""
from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
from typing import Dict, List

FITSVERIFY_BIN = "/mnt/d/ISRO/.toolchain/bin/fitsverify"
SUMMARY_RE = re.compile(r"Verification found (\d+) warning\(s\) and (\d+) error\(s\)")
FITS_EXTS = {
    ".fits", ".fit", ".fts",
    ".evt", ".evn", ".gti", ".gtis", ".bti", ".dph", ".dpi",
    ".img", ".img_q0", ".img_q1", ".img_q2", ".img_q3",
    ".pha", ".spec", ".lc", ".lcurv", ".rsp", ".rmf", ".arf",
    ".orb", ".mkf", ".att", ".hdr", ".frt", ".aspect_q0", ".aspect_q1", ".aspect_q2", ".aspect_q3",
}


def is_fits_candidate(path: str) -> bool:
    low = path.lower()
    if any(low.endswith(ext) for ext in FITS_EXTS):
        return True
    try:
        with open(path, "rb") as f:
            return f.read(9) == b"SIMPLE  ="
    except OSError:
        return False


def collect_files(targets: List[str]) -> List[str]:
    files: List[str] = []
    for t in targets:
        if os.path.isfile(t):
            files.append(t)
        elif os.path.isdir(t):
            for root, _, names in os.walk(t):
                for name in sorted(names):
                    full = os.path.join(root, name)
                    if is_fits_candidate(full):
                        files.append(full)
    return sorted(set(files))


def verify_one(path: str) -> Dict[str, object]:
    proc = subprocess.run(
        [FITSVERIFY_BIN, path],
        capture_output=True,
        text=True,
        check=False,
    )
    out = proc.stdout + "\n" + proc.stderr
    m = SUMMARY_RE.search(out)
    if not m:
        return {
            "file": path,
            "ok": False,
            "warnings": -1,
            "errors": -1,
            "returncode": proc.returncode,
            "detail": out.strip()[-1000:],
        }
    warnings = int(m.group(1))
    errors = int(m.group(2))
    return {
        "file": path,
        "ok": (errors == 0),
        "warnings": warnings,
        "errors": errors,
        "returncode": proc.returncode,
        "detail": out.strip() if errors > 0 else "",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate FITS files with HEASARC fitsverify")
    parser.add_argument("targets", nargs="+", help="FITS files or directories to scan")
    parser.add_argument("--json-out", help="Optional path to write JSON results")
    parser.add_argument("--strict-warnings", action="store_true", help="Fail on warnings too")
    args = parser.parse_args()

    if not os.path.exists(FITSVERIFY_BIN):
        print(f"ERROR: fitsverify binary not found at {FITSVERIFY_BIN}", file=sys.stderr)
        return 2

    files = collect_files(args.targets)
    if not files:
        print("ERROR: No FITS files found to verify", file=sys.stderr)
        return 2

    results = []
    failed = 0
    for f in files:
        res = verify_one(f)
        if args.strict_warnings and int(res["warnings"]) > 0:
            res["ok"] = False
        results.append(res)
        status = "PASS" if res["ok"] else "FAIL"
        print(f"[{status}] {f} (warnings={res['warnings']}, errors={res['errors']})")
        if not res["ok"]:
            failed += 1
            if res["detail"]:
                print(res["detail"])

    if args.json_out:
        os.makedirs(os.path.dirname(os.path.abspath(args.json_out)), exist_ok=True)
        with open(args.json_out, "w", encoding="utf-8") as fh:
            json.dump({"total": len(results), "failed": failed, "results": results}, fh, indent=2)

    print(f"\nfitsverify summary: {len(results) - failed}/{len(results)} passed (0 errors required)")
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
