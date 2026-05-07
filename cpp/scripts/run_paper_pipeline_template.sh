#!/usr/bin/env bash
set -euo pipefail

# Run the seven paper approaches on one graph and summarize the paper-style metrics.
# Suggested thresholds from the manuscript:
#   MIT:    DX_THRESHOLD=1e-3 ETA_THRESHOLD=1
#   FR079:  DX_THRESHOLD=1e-4 ETA_THRESHOLD=0.6
#   CSAIL:  DX_THRESHOLD=1e-5 ETA_THRESHOLD=0.95
#   Intel:  DX_THRESHOLD=1e-6 ETA_THRESHOLD=0.72
#   FRH:    DX_THRESHOLD=1e-7 ETA_THRESHOLD=0.45
# MIT-P additionally requires generating the MIT graph with intermittent position priors.

BUILD_DIR="${BUILD_DIR:-build}"
EXE="${EXE:-$BUILD_DIR/islam_run_benchmark}"
GRAPH="${1:-path/to/graph.g2o}"
OUT_DIR="${2:-paper_outputs}"
DATASET_LABEL="${DATASET_LABEL:-$(basename "$GRAPH")}" 
MAX_BATCHES="${MAX_BATCHES:-100000000}"
MAX_GN_ITER="${MAX_GN_ITER:-10}"
SELECTIVE_ALPHA="${SELECTIVE_ALPHA:-0.3}"
LC_GAP="${LC_GAP:-5}"
ETA_THRESHOLD="${ETA_THRESHOLD:-1.0}"
DX_THRESHOLD="${DX_THRESHOLD:-1e-3}"
BACKEND="${BACKEND:-custom-sparse-expanding}"
RUN_SUMMARY="${RUN_SUMMARY:-1}"

mkdir -p "$OUT_DIR"
if [[ ! -x "$EXE" ]]; then
  echo "missing executable: $EXE" >&2
  echo "Build first, or set BUILD_DIR/EXE." >&2
  exit 1
fi

run_one() {
  local label="$1"
  local gn_iter="$2"
  local gating="$3"
  local use_spo="$4"
  local csv="$OUT_DIR/${label}.csv"
  echo "running $label: gn_iter=$gn_iter gating=$gating use_spo=$use_spo backend=$BACKEND -> $csv"
  "$EXE" \
    "$GRAPH" "$csv" \
    "$MAX_BATCHES" "$gn_iter" "$SELECTIVE_ALPHA" "$LC_GAP" \
    "$ETA_THRESHOLD" "$gating" "$use_spo" "$BACKEND" "$DX_THRESHOLD"
}

run_one "GNi-SPO-IGG" "$MAX_GN_ITER" "IGG"  "1"
run_one "GNi-SPO-LCG" "$MAX_GN_ITER" "LCG"  "1"
run_one "GNi-SPO"     "$MAX_GN_ITER" "none" "1"
run_one "GNi-IGG"     "$MAX_GN_ITER" "IGG"  "0"
run_one "GNi-LCG"     "$MAX_GN_ITER" "LCG"  "0"
run_one "GNi"         "$MAX_GN_ITER" "none" "0"
run_one "GN1"         "1"            "none" "0"

if [[ "$RUN_SUMMARY" == "1" ]]; then
  env PYTHONNOUSERSITE=1 PYTHONPATH= python3 -S "$(dirname "$0")/summarize_paper_results.py" "$OUT_DIR" --dataset "$DATASET_LABEL"
fi

echo "paper pipeline completed; inspect $OUT_DIR"
