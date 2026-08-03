#!/usr/bin/env bash
# Single-split (NOT 10-fold) version of the PCA-dimension sweep, for the same splits.
# Uses the same per-dim feature CSVs; --cv-folds 1 -> one train/val/test split per run.
# Run from drug_synergy_baseline/ with the non-iCloud venv.
set -u

BASE="/Users/arnoldcheskis/Documents/Projects/drug_discovery/drug_synergy_baseline"
PY="$HOME/.venvs/drug_synergy/bin/python"
DATA="../data/data_compression/pca_sweep"
DIMS=(1 2 4 8 16 32 58)
SPLITS=(random cell_line drug_and_cell_line)
EPOCHS=35
cd "$BASE" || exit 1

run_one () {
  local out="$1"; shift
  if [ -f "$out/metrics.json" ]; then echo "[skip] $out"; return 0; fi
  echo "[run ] $out"
  "$PY" -u -m src.train \
    --output-dir "$out" \
    --cv-folds 1 --seed 42 --max-samples 10000 --epochs "$EPOCHS" \
    "$@" > "/tmp/single_$(basename "$out").log" 2>&1
  if grep -q "Saved metrics to" "/tmp/single_$(basename "$out").log"; then echo "[ ok ] $out"; else echo "[FAIL] $out"; fi
}

for split in "${SPLITS[@]}"; do
  for K in "${DIMS[@]}"; do
    run_one "outputs/sweep_single/pca${K}_${split}" \
      --split-strategy "$split" \
      --synergy-path "$DATA/drugcomb_pca${K}.csv" \
      --cell-expression-path "$DATA/drugcomb_pca${K}.csv"
  done
  run_one "outputs/sweep_single/nogenes_${split}" \
    --split-strategy "$split" \
    --synergy-path "$DATA/drugcomb_pca8.csv" \
    --no-use-gene-expression
done
echo "=== SINGLE-SPLIT SWEEP COMPLETE ==="
