#!/usr/bin/env bash
set -euo pipefail
BUILD_DIR="${BUILD_DIR:-build}"
DATASET="${1:-path/to/graph.g2o}"
OUT_DIR="${2:-benchmark_outputs}"
MAX_BATCHES="${MAX_BATCHES:-200}"
MAX_GN_ITER="${MAX_GN_ITER:-10}"
SELECTIVE_ALPHA="${SELECTIVE_ALPHA:-0.3}"
LC_GAP="${LC_GAP:-5}"
ETA_THRESHOLD="${ETA_THRESHOLD:-1.0}"
DX_THRESHOLD="${DX_THRESHOLD:-1e-3}"
GATING_MODE="${GATING_MODE:-IGG}"
USE_SPO="${USE_SPO:-1}"
mkdir -p "$OUT_DIR"
if [[ ! -x "$BUILD_DIR/islam_run_benchmark" ]]; then
  echo "missing executable: $BUILD_DIR/islam_run_benchmark" >&2
  exit 1
fi
run_one() {
  local label="$1"
  local backend="$2"
  local csv="$OUT_DIR/${label}.csv"
  echo "running $label backend=$backend -> $csv"
  "$BUILD_DIR/islam_run_benchmark" \
    "$DATASET" "$csv" \
    "$MAX_BATCHES" "$MAX_GN_ITER" "$SELECTIVE_ALPHA" "$LC_GAP" \
    "$ETA_THRESHOLD" "$GATING_MODE" "$USE_SPO" "$backend" "$DX_THRESHOLD"
  if command -v awk >/dev/null 2>&1; then
    awk -F, 'NR==1 {for (i=1;i<=NF;i++) idx[$i]=i; next} NR==2 {print "  csv backend_requested=" $idx["backend_requested"] ", backend_actual=" $idx["backend_actual"]}' "$csv" || true
  fi
}
run_one dense_eigen dense-eigen
run_one cholmod_when_available cholmod
run_one custom_sparse_expanding custom-sparse-expanding
echo "backend comparison completed; inspect $OUT_DIR"
