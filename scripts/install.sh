#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

if [[ -f "$SCRIPT_DIR/env.sh" ]]; then
  # shellcheck disable=SC1091
  source "$SCRIPT_DIR/env.sh"
fi

PREFIX="${PREFIX:-${HOME}/.local}"
SKIP_SMOKE=0

usage() {
  cat <<EOF
Zero-patch standalone installer for astrosat-rs.

Usage: $(basename "$0") [--prefix <dir>] [--skip-smoke]

Options:
  --prefix <dir>   Install directory prefix (default: \$PREFIX or \$HOME/.local)
  --skip-smoke     Skip post-install binary verification checks
  -h, --help       Show this help

Environment variables:
  PREFIX           Alternative way to specify installation prefix (no sudo required)
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --prefix)
      if [[ $# -lt 2 ]]; then
        echo "ERROR: --prefix requires a directory argument" >&2
        usage; exit 1
      fi
      PREFIX="$2"; shift 2 ;;
    --prefix=*) PREFIX="${1#*=}"; shift ;;
    --skip-smoke) SKIP_SMOKE=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

if ! command -v cargo >/dev/null 2>&1; then
  echo "ERROR: cargo not found on PATH. Install Rust (https://rustup.rs) and retry." >&2
  exit 1
fi

echo "==> Building standalone astrosat binary (--release, zero C/CFITSIO/HEASoft deps)..."
cargo build --release -p astrosat-cli --manifest-path "$ROOT_DIR/Cargo.toml"

TARGET_DIR="${CARGO_TARGET_DIR:-$ROOT_DIR/target}"
SRC_BIN="$TARGET_DIR/release/astrosat"
if [[ ! -x "$SRC_BIN" ]]; then
  echo "ERROR: built binary not found at $SRC_BIN" >&2
  exit 1
fi

BIN_DIR="$PREFIX/bin"
DST_BIN="$BIN_DIR/astrosat"
mkdir -p "$BIN_DIR"
install -m 0755 "$SRC_BIN" "$DST_BIN"
echo "==> Installed: $DST_BIN"

if [[ "$SKIP_SMOKE" -eq 0 ]]; then
  echo "==> Verifying standalone binary..."
  "$DST_BIN" --help >/dev/null
  "$DST_BIN" czti --help >/dev/null
  "$DST_BIN" czti addspec --help >/dev/null
  "$DST_BIN" laxpc --help >/dev/null
  "$DST_BIN" laxpc backshift --help >/dev/null
  "$DST_BIN" bary --help >/dev/null
  "$DST_BIN" inspect --help >/dev/null
  "$DST_BIN" validate --help >/dev/null
  "$DST_BIN" timing --help >/dev/null

  if command -v ldd >/dev/null 2>&1; then
    echo "==> Shared library audit (ldd):"
    LDD_OUT="$(ldd "$DST_BIN" 2>&1 || true)"
    echo "$LDD_OUT" | sed 's/^/    /'
    if echo "$LDD_OUT" | grep -qiE 'cfitsio|heasoft|glog|pil|gfortran'; then
      echo "ERROR: unexpected external scientific C/Fortran runtime dependency detected!" >&2
      exit 1
    fi
  fi
  echo "==> Verification passed (zero external CFITSIO/HEASoft/glog dependencies)."
fi

case ":${PATH}:" in
  *":${BIN_DIR}:"*) ;;
  *)
    echo ""
    echo "NOTE: Add $BIN_DIR to your PATH to invoke 'astrosat' directly:"
    echo "  export PATH=\"$BIN_DIR:\$PATH\""
    ;;
esac
