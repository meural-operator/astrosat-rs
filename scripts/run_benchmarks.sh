#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

if [[ -f "$SCRIPT_DIR/env.sh" ]]; then
  # shellcheck disable=SC1091
  source "$SCRIPT_DIR/env.sh"
fi

JSON_OUT="${BENCH_JSON_OUT:-$ROOT_DIR/.tmp/bench_results.json}"
MD_OUT="${BENCH_MD_OUT:-$ROOT_DIR/.tmp/bench_results.md}"
THREADS="${BENCH_THREADS:-1,4,16,max}"
SUITE="${BENCH_SUITE:-all}"
ITERATIONS="${BENCH_ITERATIONS:-3}"
REAL_ITERATIONS="${BENCH_REAL_ITERATIONS:-1}"

usage() {
  cat <<EOF
Usage: $(basename "$0") [OPTIONS]

Run the reproducible astrosat-rs benchmark suite in release mode.

Options:
  --suite <name>            Suite target (default: all)
  --threads <list>          Thread sweep (default: 1,4,16,max)
  --iterations <N>          Iterations for micro/synth benchmarks (default: 3)
  --real-iterations <N>     Iterations for real-dataset E2E runs (default: 1)
  --json-out <path>         JSON report path (default: $ROOT_DIR/.tmp/bench_results.json)
  --md-out <path>           Markdown report path (default: $ROOT_DIR/.tmp/bench_results.md)
  -h, --help                Show this help
EOF
}

require_arg() {
  if [[ "$1" -lt 2 ]]; then
    echo "ERROR: option $2 requires a value" >&2
    usage; exit 1
  fi
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --suite) require_arg $# "$1"; SUITE="$2"; shift 2 ;;
    --threads) require_arg $# "$1"; THREADS="$2"; shift 2 ;;
    --iterations) require_arg $# "$1"; ITERATIONS="$2"; shift 2 ;;
    --real-iterations) require_arg $# "$1"; REAL_ITERATIONS="$2"; shift 2 ;;
    --json-out) require_arg $# "$1"; JSON_OUT="$2"; shift 2 ;;
    --md-out) require_arg $# "$1"; MD_OUT="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

echo "==> Building astrosat-bench (--release)..."
cargo build --release -p astrosat-bench --manifest-path "$ROOT_DIR/Cargo.toml"

BIN="${CARGO_TARGET_DIR:-$ROOT_DIR/target}/release/astrosat-bench"
if [[ ! -x "$BIN" ]]; then
  echo "ERROR: release binary not found at $BIN" >&2
  exit 1
fi

echo "==> Running suite=$SUITE threads=$THREADS..."
"$BIN" \
  --suite "$SUITE" \
  --threads "$THREADS" \
  --iterations "$ITERATIONS" \
  --real-iterations "$REAL_ITERATIONS" \
  --json-out "$JSON_OUT" \
  --md-out "$MD_OUT"

echo "==> Benchmark artifacts written:"
echo "    Markdown: $MD_OUT"
echo "    JSON:     $JSON_OUT"
