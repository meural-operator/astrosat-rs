#!/usr/bin/env bash
set -euo pipefail

# Independent end-to-end verification suite for astrosat-rs.
# Always keeps all caches, builds, and temporary outputs on D:.

SKIP_BENCHMARKS=0
SKIP_CARGO_TESTS=0
for arg in "$@"; do
    case "$arg" in
        --skip-benchmarks) SKIP_BENCHMARKS=1 ;;
        --skip-cargo-tests) SKIP_CARGO_TESTS=1 ;;
        *) echo "Unknown option: $arg" >&2; exit 2 ;;
    esac
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
if [[ -f "${ROOT_DIR}/scripts/env.sh" ]]; then
    # shellcheck disable=SC1091
    source "${ROOT_DIR}/scripts/env.sh"
fi

DATA_DIR="${ASTROSAT_TEST_DATA:-/mnt/d/ISRO/data}"
PY="${PYTHON:-python3}"
FITSVERIFY="${FITSVERIFY_BIN:-$(command -v fitsverify || echo /mnt/d/ISRO/.toolchain/bin/fitsverify)}"
VERIFY_TMP="${VERIFY_TMP:-${TMPDIR:-/tmp}/astrosat_verify_run}"
FRESH_PREFIX="${FRESH_PREFIX:-${TMPDIR:-/tmp}/astrosat_fresh_install}"
REPORT_DIR="${ROOT_DIR}/verification/reports"
JPLEPH405="${DATA_DIR}/ephem/JPLEPH.405"

rm -rf "${VERIFY_TMP}" "${FRESH_PREFIX}"
mkdir -p "${VERIFY_TMP}/czti" "${VERIFY_TMP}/laxpc" "${VERIFY_TMP}/laxpc_parity" "${REPORT_DIR}"

echo "=== [1/6] C: Storage & Environment Audit ==="
df -h / /mnt/c /mnt/d
for var_name in CARGO_HOME RUSTUP_HOME CARGO_TARGET_DIR TMPDIR; do
    val="${!var_name:-}"
    if [[ "${val}" != /mnt/d/ISRO* ]]; then
        echo "FAIL: ${var_name}=${val} is not under /mnt/d/ISRO!" >&2
        exit 1
    fi
done
for bad_dir in "/home/diat/.cargo" "/home/diat/.rustup"; do
    if [[ -d "${bad_dir}" ]]; then
        echo "WARNING: Unexpected directory on C:/WSL root: ${bad_dir}"
    fi
done
echo "CARGO_HOME=${CARGO_HOME} RUSTUP_HOME=${RUSTUP_HOME} CARGO_TARGET_DIR=${CARGO_TARGET_DIR} TMPDIR=${TMPDIR}"

echo "=== [2/6] Workspace Tests, Clippy & Release Build ==="
cd "${ROOT_DIR}"
if [[ "${SKIP_CARGO_TESTS}" -eq 0 ]]; then
    cargo test --workspace
else
    echo "SKIPPED cargo test --workspace (--skip-cargo-tests)"
fi
cargo clippy --workspace --all-targets -- -D warnings
cargo build --release

ASTROSAT_BIN="${CARGO_TARGET_DIR}/release/astrosat"
if [[ ! -x "${ASTROSAT_BIN}" ]]; then
    echo "FAIL: Missing standalone binary at ${ASTROSAT_BIN}" >&2
    exit 1
fi

echo "=== [3/6] Standalone Binary ldd & Zero-Patch Install Smoke Test ==="
LDD_OUT="$(ldd "${ASTROSAT_BIN}")"
echo "${LDD_OUT}"
if echo "${LDD_OUT}" | grep -iE 'cfitsio|heasoft|glog|pil|lapack|gfortran'; then
    echo "FAIL: astrosat binary links forbidden external C/Fortran libraries!" >&2
    exit 1
fi

if [[ -x "${ROOT_DIR}/scripts/install.sh" ]]; then
    PREFIX="${FRESH_PREFIX}" bash "${ROOT_DIR}/scripts/install.sh" --prefix "${FRESH_PREFIX}"
    INSTALLED_BIN="${FRESH_PREFIX}/bin/astrosat"
else
    echo "FAIL: scripts/install.sh not found or not executable" >&2
    exit 1
fi

if [[ ! -x "${INSTALLED_BIN}" ]]; then
    echo "FAIL: Installed binary missing at ${INSTALLED_BIN}" >&2
    exit 1
fi

"${INSTALLED_BIN}" --help >/dev/null
for sub in czti laxpc bary inspect validate timing; do
    "${INSTALLED_BIN}" "${sub}" --help >/dev/null
done
for csub in run decode bunchclean pha2energy gtigen gaas datasel clean evtencor evtclean dpigen image bindata addspec; do
    "${INSTALLED_BIN}" czti "${csub}" --help >/dev/null
done
for lsub in run decode gtigen lightcurve spectrum backshift orbit; do
    "${INSTALLED_BIN}" laxpc "${lsub}" --help >/dev/null
done
echo "PASS: Zero-patch install & all 28 subcommand --help checks succeeded."

echo "=== [4/6] End-to-End CZTI 9000000406 Pipeline + addspec + bary + fitsverify + Parity ==="
CZTI_L1_DIR="/mnt/d/ISRO/data/sample/czti/czti_sample_data/20160331_T01_112T01_9000000406_level1"
CZTI_CALDB="/mnt/d/ISRO/data/caldb/czti/data/as1/czti"
CZTI_ORB="${CZTI_L1_DIR}/aux/aux1/AS1T01_112T01_9000000406czt_level1.orb"
if [[ ! -f "${CZTI_ORB}" ]]; then
    CZTI_ORB="/mnt/d/ISRO/data/sample/czti_orb/AS1T01_112T01_9000000406czt_level1.orb"
fi

"${INSTALLED_BIN}" czti run \
    --input "${CZTI_L1_DIR}" \
    --caldb "${CZTI_CALDB}" \
    --gen-evt \
    --output "${VERIFY_TMP}/czti"

CZTI_CLEAN_EVT="$(find "${VERIFY_TMP}/czti" -name '*quad_clean.evt' | head -n 1)"
if [[ -z "${CZTI_CLEAN_EVT}" ]]; then
    echo "FAIL: Missing generated CZTI *quad_clean.evt" >&2
    exit 1
fi

"${INSTALLED_BIN}" czti addspec \
    --pha "${VERIFY_TMP}/czti/"*_Q0.pha "${VERIFY_TMP}/czti/"*_Q1.pha "${VERIFY_TMP}/czti/"*_Q2.pha "${VERIFY_TMP}/czti/"*_Q3.pha \
    --rsp "${VERIFY_TMP}/czti/"*_Q0.rsp "${VERIFY_TMP}/czti/"*_Q1.rsp "${VERIFY_TMP}/czti/"*_Q2.rsp "${VERIFY_TMP}/czti/"*_Q3.rsp \
    --out-pha "${VERIFY_TMP}/czti/combined_Q0123.pha" \
    --out-rsp "${VERIFY_TMP}/czti/combined_Q0123.rsp"

"${INSTALLED_BIN}" bary \
    --input "${CZTI_CLEAN_EVT}" \
    --orbit "${CZTI_ORB}" \
    --output "${VERIFY_TMP}/czti/czti_bary.evt"

CZTI_DE405_ARG=()
if [[ -f "${JPLEPH405}" ]]; then
    "${INSTALLED_BIN}" bary \
        --input "${CZTI_CLEAN_EVT}" \
        --orbit "${CZTI_ORB}" \
        --ephem "${JPLEPH405}" \
        --output "${VERIFY_TMP}/czti/czti_bary_de405.evt"
    "${INSTALLED_BIN}" validate --strict "${VERIFY_TMP}/czti/czti_bary_de405.evt" >/dev/null
    CZTI_DE405_ARG=(--bary-de405-evt "${VERIFY_TMP}/czti/czti_bary_de405.evt")
fi

"${INSTALLED_BIN}" validate --strict "${VERIFY_TMP}/czti/czti_bary.evt" "${VERIFY_TMP}/czti/combined_Q0123.pha" "${VERIFY_TMP}/czti/combined_Q0123.rsp" >/dev/null
"${PY}" "${ROOT_DIR}/verification/check_fitsverify.py" --strict-warnings "${VERIFY_TMP}/czti" --json-out "${REPORT_DIR}/czti_fitsverify.json"
"${PY}" "${ROOT_DIR}/verification/check_czti_parity.py" "${VERIFY_TMP}/czti" --bary-evt "${VERIFY_TMP}/czti/czti_bary.evt" "${CZTI_DE405_ARG[@]}" --json-out "${REPORT_DIR}/czti_parity.json"

echo "=== [5/6] End-to-End LAXPC Crab 9000000252 Pipeline + fitsverify + Parity + Crab Timing ==="
LAXPC_L1_DIR="/mnt/d/ISRO/data/sample/laxpc_l1/20160107_G02_010T01_9000000252_level1/laxpc/01509"
LAXPC_L2_REF="/mnt/d/ISRO/data/sample/laxpc_l2/20160107_G02_010T01_9000000252_level2/laxpc/01509"
LAXPC_CALDB="/mnt/d/ISRO/data/caldb/laxpc"

"${INSTALLED_BIN}" laxpc run \
    --input "${LAXPC_L1_DIR}" \
    --caldb "${LAXPC_CALDB}" \
    --format-a \
    --output "${VERIFY_TMP}/laxpc"

for unit in lxp1 lxp2 lxp3; do
    "${INSTALLED_BIN}" laxpc run \
        --input "${LAXPC_L1_DIR}" \
        --caldb "${LAXPC_CALDB}" \
        --units "${unit}" \
        --channels 1024 \
        --user-gti "${LAXPC_L2_REF}/${unit}/modeEA/AS1G02_010T01_9000000252${unit}EA_level2.gti" \
        --user-gti-only \
        --no-orbit \
        --no-format-b \
        --output "${VERIFY_TMP}/laxpc_parity"
done

LAXPC_EVN="$(find "${VERIFY_TMP}/laxpc" -name '*lxp1*.evn' | head -n 1)"
LAXPC_FMT_A="$(find "${VERIFY_TMP}/laxpc" -name 'level2.event.fits' | head -n 1)"
LAXPC_ORB="$(find "${VERIFY_TMP}/laxpc" -name 'AstroSat.orb' | head -n 1)"
if [[ -z "${LAXPC_EVN}" || -z "${LAXPC_FMT_A}" || -z "${LAXPC_ORB}" ]]; then
    echo "FAIL: Missing LAXPC outputs (evn=${LAXPC_EVN}, fmt_a=${LAXPC_FMT_A}, orb=${LAXPC_ORB})" >&2
    exit 1
fi

"${INSTALLED_BIN}" bary \
    --input "${LAXPC_EVN}" \
    --orbit "${LAXPC_ORB}" \
    --output "${VERIFY_TMP}/laxpc/laxpc_bary.evn"
"${INSTALLED_BIN}" bary \
    --input "${LAXPC_FMT_A}" \
    --orbit "${LAXPC_ORB}" \
    --output "${VERIFY_TMP}/laxpc/level2.event.bary.fits"

LAXPC_DE405_ARG=()
CRAB_TIMING_FILE="${VERIFY_TMP}/laxpc/laxpc_bary.evn"
if [[ -f "${JPLEPH405}" ]]; then
    "${INSTALLED_BIN}" bary \
        --input "${LAXPC_EVN}" \
        --orbit "${LAXPC_ORB}" \
        --ephem "${JPLEPH405}" \
        --output "${VERIFY_TMP}/laxpc/laxpc_bary_de405.evn"
    "${INSTALLED_BIN}" validate --strict "${VERIFY_TMP}/laxpc/laxpc_bary_de405.evn" >/dev/null
    LAXPC_DE405_ARG=(--bary-de405-evn "${VERIFY_TMP}/laxpc/laxpc_bary_de405.evn")
    CRAB_TIMING_FILE="${VERIFY_TMP}/laxpc/laxpc_bary_de405.evn"
fi

"${INSTALLED_BIN}" inspect --full "${VERIFY_TMP}/laxpc/laxpc_bary.evn" >/dev/null
"${INSTALLED_BIN}" validate --strict "${VERIFY_TMP}/laxpc/laxpc_bary.evn" "${VERIFY_TMP}/laxpc/level2.event.bary.fits" >/dev/null
"${PY}" "${ROOT_DIR}/verification/check_fitsverify.py" --strict-warnings "${VERIFY_TMP}/laxpc" --json-out "${REPORT_DIR}/laxpc_fitsverify.json"
"${PY}" "${ROOT_DIR}/verification/check_laxpc_parity.py" "${VERIFY_TMP}/laxpc" \
    --bary-evn "${VERIFY_TMP}/laxpc/laxpc_bary.evn" \
    --bary-fmt-a "${VERIFY_TMP}/laxpc/level2.event.bary.fits" \
    "${LAXPC_DE405_ARG[@]}" \
    --parity-dir "${VERIFY_TMP}/laxpc_parity" \
    --json-out "${REPORT_DIR}/laxpc_parity.json"
"${PY}" "${ROOT_DIR}/verification/check_crab_timing.py" "${CRAB_TIMING_FILE}" --json-out "${REPORT_DIR}/crab_timing.json"
"${INSTALLED_BIN}" timing "${CRAB_TIMING_FILE}" --freq-min 29.65 --freq-max 29.665 --fdot=-3.693e-10

echo "=== [6/6] Benchmark Suite Verification ==="
if [[ "${SKIP_BENCHMARKS}" -eq 1 ]]; then
    echo "SKIPPED (--skip-benchmarks)"
else
    "${PY}" "${ROOT_DIR}/verification/check_benchmarks.py" --bench-script "${ROOT_DIR}/scripts/run_benchmarks.sh" --readme "${ROOT_DIR}/README.md"
fi

echo "ALL VERIFICATION CHECKS PASSED!"
