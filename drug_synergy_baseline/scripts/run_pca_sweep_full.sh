#!/usr/bin/env bash
# Full-dataset (no --max-samples cap) version of the PCA-dimension 10-fold CV sweep,
# using codex's EXACT original model hyperparameters (epochs=10, lr=0.001,
# hidden_dims=[512,256], dropout=0.2, batch_size=64) instead of the practical_research
# macro preset. IMPORTANT CAVEAT: codex's original baseline was a single train/val/test
# split, never 10-fold CV -- there is no "codex CV run" to reproduce. This script only
# matches the *model* hyperparameters; the CV *protocol* itself is necessarily new.
# Writes to outputs/sweep_full/ so the existing 10k sweep in outputs/sweep/ (which used
# practical_research) is left untouched for comparison. Uses the non-iCloud venv.
#
#   bash /Users/arnoldcheskis/Documents/Projects/drug_discovery/drug_synergy_baseline/scripts/run_pca_sweep_full.sh
#
# Progress is logged per run under /tmp/sweep_full_<run>.log.
set -u

BASE="/Users/arnoldcheskis/Documents/Projects/drug_discovery/drug_synergy_baseline"
PY="$HOME/.venvs/drug_synergy/bin/python"
DATA="../data/data_compression/pca_sweep"
DIMS=(1 2 4 8 16 32 58)
SPLITS=(random cell_line drug_and_cell_line)
cd "$BASE" || exit 1

run_one () {
  local out="$1"; shift
  if [ -f "$out/cv_metrics.json" ]; then
    echo "[skip] $out already has cv_metrics.json"
    return 0
  fi
  echo "[run ] $out"
  "$PY" -u -m src.train \
    --output-dir "$out" \
    --cv-folds 10 --cv-seeds 42 \
    --epochs 10 --lr 0.001 --hidden-dims 512 256 --dropout 0.2 --batch-size 64 \
    "$@" > "/tmp/sweep_full_$(basename "$out").log" 2>&1
  if grep -q "Saved CV summary" "/tmp/sweep_full_$(basename "$out").log"; then
    echo "[ ok ] $out"
  else
    echo "[FAIL] $out  (see /tmp/sweep_full_$(basename "$out").log)"
  fi
}

for split in "${SPLITS[@]}"; do
  for K in "${DIMS[@]}"; do
    run_one "outputs/sweep_full/pca${K}_${split}" \
      --split-strategy "$split" \
      --synergy-path "$DATA/drugcomb_pca${K}.csv" \
      --cell-expression-path "$DATA/drugcomb_pca${K}.csv"
  done
  run_one "outputs/sweep_full/nogenes_${split}" \
    --split-strategy "$split" \
    --synergy-path "$DATA/drugcomb_pca8.csv" \
    --no-use-gene-expression
done

echo "=== FULL-DATASET CV SWEEP COMPLETE ==="
