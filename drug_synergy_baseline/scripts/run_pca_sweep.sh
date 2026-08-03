#!/usr/bin/env bash
# PCA-dimension sweep: 10-fold CV per (compressed dim x split), plus a no-genes point per split.
# Run from the drug_synergy_baseline/ directory. Uses the non-iCloud venv to avoid dataless hangs.
#
#   bash scripts/run_pca_sweep.sh
#
# Progress is logged per run under /tmp/sweep_<run>.log; this script also tees a summary.
set -u

BASE="/Users/arnoldcheskis/Documents/Projects/drug_discovery/drug_synergy_baseline"
PY="$HOME/.venvs/drug_synergy/bin/python"
DATA="../data/data_compression/pca_sweep"
DIMS=(1 2 4 8 16 32 58)
SPLITS=(random cell_line drug_and_cell_line)
EPOCHS=35
cd "$BASE" || exit 1

run_one () {  # $1=output-dir  shift; rest = extra train args
  local out="$1"; shift
  if [ -f "$out/cv_metrics.json" ]; then
    echo "[skip] $out already has cv_metrics.json"
    return 0
  fi
  echo "[run ] $out"
  "$PY" -u -m src.train \
    --output-dir "$out" \
    --cv-folds 10 --cv-seeds 42 --max-samples 10000 --epochs "$EPOCHS" \
    "$@" > "/tmp/sweep_$(basename "$out").log" 2>&1
  if grep -q "Saved CV summary" "/tmp/sweep_$(basename "$out").log"; then
    echo "[ ok ] $out"
  else
    echo "[FAIL] $out  (see /tmp/sweep_$(basename "$out").log)"
  fi
}

for split in "${SPLITS[@]}"; do
  for K in "${DIMS[@]}"; do
    run_one "outputs/sweep/pca${K}_${split}" \
      --split-strategy "$split" \
      --synergy-path "$DATA/drugcomb_pca${K}.csv" \
      --cell-expression-path "$DATA/drugcomb_pca${K}.csv"
  done
  # no-genes reference point (dim 0): same rows, gene features off.
  run_one "outputs/sweep/nogenes_${split}" \
    --split-strategy "$split" \
    --synergy-path "$DATA/drugcomb_pca8.csv" \
    --no-use-gene-expression
done

echo "=== SWEEP COMPLETE ==="
