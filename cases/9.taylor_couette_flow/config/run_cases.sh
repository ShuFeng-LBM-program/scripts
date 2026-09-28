#!/bin/bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BIN_DIR="$(dirname "$SCRIPT_DIR")"

CASES=("50dx_noGeoUpdate" "40dx_noGeoUpdate" "32dx_noGeoUpdate" "25dx_noGeoUpdate")

cd "$BIN_DIR"

for case in "${CASES[@]}"; do
    JSON="$SCRIPT_DIR/$case/taylor_couette.json"
    LOG="$SCRIPT_DIR/$case/log.txt"

    echo "=========================================="
    echo "Running $case ..."
    echo "=========================================="

    ./run "$JSON" >& "$LOG"

    echo "$case finished."
    echo ""
done

echo "All cases completed."
