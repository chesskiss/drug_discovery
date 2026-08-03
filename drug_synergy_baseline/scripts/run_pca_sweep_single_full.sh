#!/usr/bin/env bash
# Full-dataset (no --max-samples cap) single-split version of the PCA sweep, using
# the EXACT hyperparameters from codex's original baseline runs (outputs/genes_on_random
# etc.): epochs=10, lr=0.001, hidden_dims=[512,256], dropout=0.2, batch_size=64,
# train/val fractions 0.8/0.1. These do NOT match the practical_research macro preset
# (35 epochs, lr=0.0001, [1024,512], dropout=0.3) -- explicit flags below override it.
# Writes to outputs/sweep_single_full/ so outputs/sweep_single/ (10k, practical_research
# settings) is untouched.
#
#   bash /Users/arnoldcheskis/Documents/Projects/drug_discovery/drug_synergy_baseline/scripts/run_pca_sweep_single_full.sh
set -u

BASE="/Users/arnoldcheskis/Documents/Projects/drug_discovery/drug_synergy_baseline"
PY="$HOME/.venvs/drug_synergy/bin/python"
DATA="../data/data_compression/pca_sweep"
DIMS=(1 2 4 8 16 32 58)
SPLITS=(random cell_line drug_and_cell_line)
cd "$BASE" || exit 1

run_one () {
  local out="$1"; shift
  if [ -f "$out/metrics.json" ]; then echo "[skip] $out"; return 0; fi
  echo "[run ] $out"
  "$PY" -u -m src.train \
    --output-dir "$out" \
    --cv-folds 1 --seed 42 \
    --epochs 10 --lr 0.001 --hidden-dims 512 256 --dropout 0.2 --batch-size 64 \
    --train-fraction 0.8 --val-fraction 0.1 \
    "$@" > "/tmp/single_full_$(basename "$out").log" 2>&1
  if grep -q "Saved metrics to" "/tmp/single_full_$(basename "$out").log"; then echo "[ ok ] $out"; else echo "[FAIL] $out"; fi
}

for split in "${SPLITS[@]}"; do
  for K in "${DIMS[@]}"; do
    run_one "outputs/sweep_single_full/pca${K}_${split}" \
      --split-strategy "$split" \
      --synergy-path "$DATA/drugcomb_pca${K}.csv" \
      --cell-expression-path "$DATA/drugcomb_pca${K}.csv"
  done
  run_one "outputs/sweep_single_full/nogenes_${split}" \
    --split-strategy "$split" \
    --synergy-path "$DATA/drugcomb_pca8.csv" \
    --no-use-gene-expression
done
echo "=== FULL-DATASET SINGLE-SPLIT SWEEP COMPLETE ==="
