#!/usr/bin/env bash
# run_local.sh – Build and run FlowLabLite locally, extract JSON results.
#
# Usage:
#   bash run_local.sh                        # both solvers, save to results.json
#   bash run_local.sh my_output.json         # custom output filename
#
# The JSON file can then be loaded in cmd/main/local_viewer.html:
#   - Drag & drop results.json onto the viewer page.
#   - Use the "Chorin / SIMPLE" tabs to switch between solver views.
#
# Solver configuration (steps / iterations / timing on-off):
#   Edit the constants at the top of cmd/main/main.mbt:
#     let run_chorin     : Bool = true    # include Chorin solver
#     let local_nt       : Int  = nt      # Chorin steps (default: 500)
#     let run_simple     : Bool = true    # include SIMPLE solver
#     let local_simple_n : Int  = 100     # SIMPLE iterations
#     let timing_enabled : Bool = true    # print [Timing] lines

set -euo pipefail

OUTPUT_FILE="${1:-results.json}"

echo "[run_local] Building and running FlowLabLite locally..."
echo "[run_local] Target: wasm (Wasmtime runtime)"

# Run the simulation and capture full output (warnings go to stderr → captured via 2>&1)
FULL_OUTPUT=$(moon run cmd/main --target wasm 2>&1)

echo "[run_local] Simulation completed. Extracting JSON..."

# Extract JSON block between markers (strip marker lines)
echo "$FULL_OUTPUT" \
  | sed -n '/===JSON_DATA_START===/,/===JSON_DATA_END===/p' \
  | grep -v '===JSON_DATA' \
  > "$OUTPUT_FILE"

# Validate JSON and report grid size
if node -e "JSON.parse(require('fs').readFileSync('$OUTPUT_FILE','utf8'))" 2>/dev/null; then
  CHORIN_PTS=$(node -e "
    const d=JSON.parse(require('fs').readFileSync('$OUTPUT_FILE','utf8'));
    console.log((d.grid||[]).length);
  " 2>/dev/null || echo "?")
  SIMPLE_PTS=$(node -e "
    const d=JSON.parse(require('fs').readFileSync('$OUTPUT_FILE','utf8'));
    console.log((d.simple_grid||[]).length);
  " 2>/dev/null || echo "0")
  echo "[run_local] Valid JSON saved to: $OUTPUT_FILE"
  echo "[run_local]   Chorin grid  : $CHORIN_PTS points"
  echo "[run_local]   SIMPLE grid  : $SIMPLE_PTS points"
else
  echo "[run_local] WARNING: JSON validation failed. File written: $OUTPUT_FILE"
fi

echo "[run_local] Open cmd/main/local_viewer.html in a browser and load $OUTPUT_FILE"
