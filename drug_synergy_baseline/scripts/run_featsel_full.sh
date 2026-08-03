#!/usr/bin/env bash
# Full-dataset (no --max-samples cap) 10-fold CV for the raw feature-selection views
# (627/3171/23808), using codex's EXACT original hyperparameters (epochs=10, lr=0.001,
# hidden_dims=[512,256], dropout=0.2, batch_size=64) -- matching run_pca_sweep_full.sh.
# Writes into outputs/sweep_full/ (same dir as the PCA CV sweep) so
# compression_rmse_plot.py picks these up as the dashed feature-selection points
# alongside the PCA curve, for the CV (not single-split) plots.
#
#   bash /Users/arnoldcheskis/Documents/Projects/drug_discovery/drug_synergy_baseline/scripts/run_featsel_full.sh
#
# Progress is logged per run under /tmp/featsel_full_cv_<run>.log.
set -u

BASE="/Users/arnoldcheskis/Documents/Projects/drug_discovery/drug_synergy_baseline"
PY="$HOME/.venvs/drug_synergy/bin/python"
PKL="../data/data_compression/source_data/drugcomb.pkl"
VIEWS=(2 1 0); declare -A DIM=( [2]=627 [1]=3171 [0]=23808 )
SPLITS=(random cell_line drug_and_cell_line)
cd "$BASE" || exit 1

run_one () {
  local out="$1"; shift
  if [ -f "$out/cv_metrics.json" ]; then echo "[skip] $out"; return 0; fi
  echo "[run ] $out"
  "$PY" -u -m src.train \
    --output-dir "$out" \
    --cv-folds 10 --cv-seeds 42 \
    --epochs 10 --lr 0.001 --hidden-dims 512 256 --dropout 0.2 --batch-size 64 \
    "$@" > "/tmp/featsel_full_cv_$(basename "$out").log" 2>&1
  if grep -q "Saved CV summary" "/tmp/featsel_full_cv_$(basename "$out").log"; then echo "[ ok ] $out"; else echo "[FAIL] $out"; fi
}

for split in "${SPLITS[@]}"; do
  for v in "${VIEWS[@]}"; do
    d=${DIM[$v]}
    run_one "outputs/sweep_full/featsel${d}_${split}" \
      --split-strategy "$split" \
      --synergy-path "$PKL" \
      --fallback-pickle-path "$PKL" \
      --cell-feature-view "$v"
  done
done
echo "=== FULL-DATASET FEATSEL CV COMPLETE ==="
