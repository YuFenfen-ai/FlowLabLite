#!/usr/bin/env bash
# run_local.sh – Run FlowLabLite locally and extract JSON results.
#
# Usage:
#   bash run_local.sh                  # run and save results.json
#   bash run_local.sh results_500.json # run and save to custom filename
#
# The JSON file can then be loaded in cmd/main/local_viewer.html for
# browser-based visualization of speed and pressure fields.

set -euo pipefail

OUTPUT_FILE="${1:-results.json}"

echo "[run_local] Building and running FlowLabLite locally..."
echo "[run_local] Target: wasm (local Wasmtime runtime)"

# Run the simulation and capture full output
FULL_OUTPUT=$(moon run cmd/main --target wasm 2>&1)

echo "[run_local] Simulation completed. Extracting JSON..."

# Extract JSON block between markers
echo "$FULL_OUTPUT" \
  | sed -n '/===JSON_DATA_START===/,/===JSON_DATA_END===/p' \
  | grep -v '===JSON_DATA' \
  > "$OUTPUT_FILE"

# Validate JSON
if node -e "JSON.parse(require('fs').readFileSync('$OUTPUT_FILE','utf8'))" 2>/dev/null; then
  POINTS=$(node -e "const d=JSON.parse(require('fs').readFileSync('$OUTPUT_FILE','utf8'));console.log(d.grid.length)")
  echo "[run_local] Valid JSON saved to: $OUTPUT_FILE ($POINTS grid points)"
else
  echo "[run_local] WARNING: JSON validation failed, but file was written to: $OUTPUT_FILE"
fi

echo "[run_local] Open cmd/main/local_viewer.html in a browser and load $OUTPUT_FILE"
