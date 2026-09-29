#!/bin/bash
# Batch runner for cylinder flow cases across resolutions (D/dx) and Reynolds numbers.
# Each case's stdout/stderr is captured to log.txt inside its own OutputDir.
#
# Usage:
#   ./batch_cylinder.sh                       # default: U=0.03 batch, all dx x Re
#   ./batch_cylinder.sh u003                  # U=0.03 batch, all dx x Re
#   ./batch_cylinder.sh u003_y20D             # U=0.03, expanded Y domain [-10,10]
#   ./batch_cylinder.sh u003_y20D_spg0        # U=0.03, Y=[-10,10], NO sponge layer
#   ./batch_cylinder.sh u001                  # U=0.01 batch (legacy, see note below)
#   ./batch_cylinder.sh u003 25               # U=0.03, only dx=25
#   ./batch_cylinder.sh u003_y20D_spg0 25 40 50 100  # u003_y20D_spg0, dx=25 Re=40,50,100
#   ./batch_cylinder.sh 25                    # default U=0.03, dx=25 only
#
# Note on U=0.01 (u001): those JSONs were moved into cylinder_bgk_u001/ for
# archival, but their internal OutputDir still points to ./cylinder_dx{N}/...
# (the original pre-move location). Re-running u001 will write outputs there,
# not next to the JSONs. Edit the JSONs' OutputDir if you want them co-located.
set -u
set -o pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BIN_DIR="/home/liangmj/repo/LBM-Program/amrlbm/bin"
RUN_EXE="${BIN_DIR}/run"

# Working directory contains 04_Cylinder.stl (so JSON's "./04_Cylinder.stl" resolves).
WORK_DIR="${SCRIPT_DIR}"

DX_VALUES=(50 40 32 25)
RE_VALUES=(20 40 50 100 150 200)
UCASE="u003"   # default to the active U=0.03 batch

# Parse optional U selector (first arg starting with "u"), then dx/Re filters.
while [ "$#" -gt 0 ]; do
    case "$1" in
        u*)
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

CASE_ROOT="cylinder_bgk_${UCASE}"

cd "${WORK_DIR}"

if [ ! -x "${RUN_EXE}" ]; then
    echo "ERROR: executable not found or not executable: ${RUN_EXE}" >&2
    exit 1
fi
if [ ! -f "./04_Cylinder.stl" ]; then
    echo "ERROR: 04_Cylinder.stl missing in ${WORK_DIR}" >&2
    exit 1
fi
if [ ! -d "${CASE_ROOT}" ]; then
    echo "ERROR: case root not found: ${CASE_ROOT}" >&2
    exit 1
fi

echo "=== Batch cylinder flow (U case: ${UCASE}) ==="
echo "WorkDir  : ${WORK_DIR}"
echo "RunExe   : ${RUN_EXE}"
echo "CaseRoot : ${CASE_ROOT}"
echo "D/dx     : ${DX_VALUES[*]}"
echo "Re       : ${RE_VALUES[*]}"
echo

total_fail=0
for dx in "${DX_VALUES[@]}"; do
    for Re in "${RE_VALUES[@]}"; do
        JSON="${CASE_ROOT}/cylinder_dx${dx}/cylinderRe${Re}_${dx}dx.json"
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

        # Run; capture combined output to both terminal and per-case log.txt
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
