"""Verify reproducible benchmark execution and cross-check README.md reported metrics."""
from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys


ROW_RE = re.compile(
    r"^\|\s*`([^`]+)`\s*\|\s*(\d+)\s*\|\s*([\d.]+)\s*\|\s*([\d.]+)\s*\|\s*([\d.]+)\s*\|\s*([\d.]+)\s*\|\s*([\d.]+)x\s*\|\s*([\d.]+)\s*\|"
)


def parse_readme_table(readme_text: str) -> dict[tuple[str, int], dict[str, float]]:
    rows: dict[tuple[str, int], dict[str, float]] = {}
    for line in readme_text.splitlines():
        m = ROW_RE.match(line.strip())
        if not m:
            continue
        bench_id = m.group(1)
        threads = int(m.group(2))
        rows[(bench_id, threads)] = {
            "wall_ms": float(m.group(3)),
            "median_ms": float(m.group(4)),
            "events_per_sec": float(m.group(5)),
            "mb_per_sec": float(m.group(6)),
            "speedup_vs_1t": float(m.group(7)),
            "peak_rss_mib": float(m.group(8)),
        }
    return rows


def run_and_check(bench_script: str, readme_path: str, json_out: str, md_out: str) -> int:
    if not os.path.exists(bench_script):
        print(f"ERROR: {bench_script} not found", file=sys.stderr)
        return 2
    if not os.path.exists(readme_path):
        print(f"ERROR: {readme_path} not found", file=sys.stderr)
        return 2

    os.makedirs(os.path.dirname(os.path.abspath(json_out)), exist_ok=True)
    proc = subprocess.run(
        ["bash", bench_script, "--json-out", json_out, "--md-out", md_out],
        capture_output=True,
        text=True,
        check=False,
    )
    print(proc.stdout)
    if proc.stderr:
        print(proc.stderr, file=sys.stderr)
    if proc.returncode != 0:
        print(f"FAIL: {bench_script} exited with code {proc.returncode}", file=sys.stderr)
        return 1

    with open(json_out, "r", encoding="utf-8") as fh:
        bench_data = json.load(fh)
    with open(readme_path, "r", encoding="utf-8") as fh:
        readme_text = fh.read()

    placeholders = re.findall(r"\b(TBD|TODO|PLACEHOLDER|XXX|Pending)\b", readme_text, flags=re.IGNORECASE)
    if placeholders:
        print(f"FAIL: README.md contains placeholder markers: {placeholders}", file=sys.stderr)
        return 1

    readme_rows = parse_readme_table(readme_text)
    if not readme_rows:
        print("FAIL: Could not parse Detailed Throughput table in README.md", file=sys.stderr)
        return 1

    measured = {(r["id"], int(r["threads"])): r for r in bench_data.get("results", [])}
    errors: list[str] = []

    for key, r_row in sorted(readme_rows.items()):
        if key not in measured:
            errors.append(f"Missing measured benchmark for {key} listed in README.md")
            continue
        m_row = measured[key]
        m_wall = float(m_row["wall_ms"])
        r_wall = r_row["wall_ms"]
        wall_ratio = max(m_wall, r_wall) / max(min(m_wall, r_wall), 1e-6)

        m_rss = float(m_row["peak_rss_mib"])
        r_rss = r_row["peak_rss_mib"]
        rss_ratio = max(m_rss, r_rss) / max(min(m_rss, r_rss), 1e-6)

        m_eps = float(m_row["events_per_sec"])
        r_eps = r_row["events_per_sec"]
        eps_ratio = max(m_eps, r_eps) / max(min(m_eps, r_eps), 1e-6)

        print(
            f"  {key[0]:24s} T={key[1]:2d}: wall={m_wall:8.2f} ms (README {r_wall:8.2f}, x{wall_ratio:.2f}) | "
            f"ev/s={m_eps:10.0f} (README {r_eps:10.0f}) | RSS={m_rss:6.1f} MiB (README {r_rss:6.1f})"
        )
        # Allow 2.5x wall/throughput variance (WSL2 I/O & cache jitter) and 1.5x RSS variance (with 64 MiB floor for small baselines)
        if wall_ratio > 2.5 or eps_ratio > 2.5:
            errors.append(f"{key}: wall/throughput diverges >2.5x from README ({m_wall:.1f} ms vs {r_wall:.1f} ms)")
        if rss_ratio > 1.5 and abs(m_rss - r_rss) > 64.0:
            errors.append(f"{key}: peak RSS diverges >1.5x and >64 MiB from README ({m_rss:.1f} MiB vs {r_rss:.1f} MiB)")

    # Verify scaling exists on compute-heavy bary_synth
    b1 = measured.get(("bary_synth", 1))
    bmax = max((r for (bid, _), r in measured.items() if bid == "bary_synth"), key=lambda x: x["threads"], default=None)
    if b1 and bmax and bmax["threads"] > 1:
        if float(bmax["speedup_vs_1t"]) < 2.0:
            errors.append(f"bary_synth scaling at max threads unexpectedly low: {bmax['speedup_vs_1t']:.2f}x")

    if errors:
        for e in errors:
            print(f"FAIL: {e}", file=sys.stderr)
        return 1

    print(f"PASS: All {len(readme_rows)} README.md benchmark entries reproduced cleanly!")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description="Verify astrosat-rs benchmark suite and README.md")
    parser.add_argument("--bench-script", default="/mnt/d/ISRO/scripts/run_benchmarks.sh")
    parser.add_argument("--readme", default="/mnt/d/ISRO/README.md")
    parser.add_argument("--json-out", default="/mnt/d/ISRO/.tmp/verify_bench.json")
    parser.add_argument("--md-out", default="/mnt/d/ISRO/.tmp/verify_bench.md")
    args = parser.parse_args()
    return run_and_check(args.bench_script, args.readme, args.json_out, args.md_out)


if __name__ == "__main__":
    sys.exit(main())
