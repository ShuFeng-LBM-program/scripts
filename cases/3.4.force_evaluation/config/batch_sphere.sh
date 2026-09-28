#!/bin/bash
# Batch runner for sphere flow cases across resolutions (D/dx) and Reynolds numbers.
# Each case's stdout/stderr is captured to log.txt inside its own OutputDir.
#
# Usage:
#   ./batch_sphere.sh                       # default: U=0.03 batch, all dx x Re
#   ./batch_sphere.sh u003                  # U=0.03 batch, all dx x Re
#   ./batch_sphere.sh u001                  # U=0.01 batch (legacy, see note below)
#   ./batch_sphere.sh u003 25               # U=0.03, only dx=25
#   ./batch_sphere.sh u003 25 10 50         # U=0.03, dx=25 with Re=10,50
#   ./batch_sphere.sh 20                    # default U=0.03, dx=20 only
#
# Note on U=0.01 (u001): those JSONs live in sphere_dx{N}/ at the top level
# (not inside a sphere_bgk_u001/ folder) and their OutputDir is also
# ./sphere_dx{N}/... Re-running u001 works as expected.
set -u
set -o pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BIN_DIR="/home/liangmj/repo/LBM-Program/amrlbm/bin"
RUN_EXE="${BIN_DIR}/run"

# Working directory contains ico_sphere_d1_sdv4.stl (so JSON's STL path resolves).
WORK_DIR="${SCRIPT_DIR}"

DX_VALUES=(20 25 32)
RE_VALUES=(5 10 50 100)
UCASE="u003"   # default to the active U=0.03 batch

# Parse optional U selector (first arg starting with "u"), then dx/Re filters.
while [ "$#" -gt 0 ]; do
    case "$1" in
        u001|u003)
            UCASE="$1"; shift
            ;;
        *)
            break
            ;;
    esac
done

if [ "$#" -ge 1 ]; then
    DX_VALUES=("$1"); shift
fi
if [ "$#" -ge 1 ]; then
    RE_VALUES=("$@")
fi

# u001 JSONs live in sphere_dx{N}/ (top level); u003 JSONs in sphere_bgk_u003/sphere_dx{N}/
if [ "${UCASE}" = "u001" ]; then
    CASE_PREFIX=""
else
    CASE_PREFIX="sphere_bgk_${UCASE}/"
fi

cd "${WORK_DIR}"

if [ ! -x "${RUN_EXE}" ]; then
    echo "ERROR: executable not found or not executable: ${RUN_EXE}" >&2
    exit 1
fi
if [ ! -f "./ico_sphere_d1_sdv4.stl" ]; then
    echo "ERROR: ico_sphere_d1_sdv4.stl missing in ${WORK_DIR}" >&2
    exit 1
fi

echo "=== Batch sphere flow (U case: ${UCASE}) ==="
echo "WorkDir  : ${WORK_DIR}"
echo "RunExe   : ${RUN_EXE}"
echo "CaseRoot : ${CASE_PREFIX:-.}"
echo "D/dx     : ${DX_VALUES[*]}"
echo "Re       : ${RE_VALUES[*]}"
echo

total_fail=0
for dx in "${DX_VALUES[@]}"; do
    for Re in "${RE_VALUES[@]}"; do
        JSON="${CASE_PREFIX}sphere_dx${dx}/sphereRe${Re}_dx${dx}.json"
        if [ ! -f "${JSON}" ]; then
            echo "WARN: JSON not found, skipping: ${JSON}" >&2
            continue
        fi

        OUTDIR=$(python3 -c "import json,sys;print(json.load(open(sys.argv[1]))['Setup']['Output']['OutputDir'])" "${JSON}")
        if [ -z "${OUTDIR}" ]; then
            echo "ERROR: could not read OutputDir from ${JSON}" >&2
            exit 1
        fi

        mkdir -p "${OUTDIR}"
        echo ">>> ${UCASE} dx=${dx} Re=${Re}  ->  ${OUTDIR}/log.txt"
        START_TIME=$(date +%s)

        if "${RUN_EXE}" "${JSON}" 2>&1 | tee "${OUTDIR}/log.txt"; then
            status="OK"
        else
            status="FAIL(code=${PIPESTATUS[0]})"
            total_fail=$((total_fail + 1))
        fi

        END_TIME=$(date +%s)
        ELAPSED=$(awk -v s="${START_TIME}" -v e="${END_TIME}" 'BEGIN{printf "%.1f min",(e-s)/60.0}')
        echo "<<< ${UCASE} dx=${dx} Re=${Re}  ${status}  (${ELAPSED})"
        echo
    done
done

echo "=== Batch finished. Failures: ${total_fail} ==="
exit ${total_fail}
