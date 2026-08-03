#!/usr/bin/env bash
# Feature-selection (pickle views) single-split runs at MATCHED sweep settings
# (10k samples, 35 epochs), extending the compression axis to 627 / 3171 / 23808.
# These are the "baseline comparison" high-dim points, re-run so they are directly
# comparable to the PCA sweep. Run from drug_synergy_baseline/.
set -u

BASE="/Users/arnoldcheskis/Documents/Projects/drug_discovery/drug_synergy_baseline"
PY="$HOME/.venvs/drug_synergy/bin/python"
PKL="../data/data_compression/source_data/drugcomb.pkl"
# view index -> effective dim
VIEWS=(2 1 0); declare -A DIM=( [2]=627 [1]=3171 [0]=23808 )
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
    "$@" > "/tmp/$(basename "$out").log" 2>&1
  if grep -q "Saved metrics to" "/tmp/$(basename "$out").log"; then echo "[ ok ] $out"; else echo "[FAIL] $out"; fi
}

for split in "${SPLITS[@]}"; do
  for v in "${VIEWS[@]}"; do
    d=${DIM[$v]}
    run_one "outputs/sweep_single/featsel${d}_${split}" \
      --split-strategy "$split" \
      --synergy-path "$PKL" \
      --fallback-pickle-path "$PKL" \
      --cell-feature-view "$v"
  done
done
echo "=== FEATSEL SINGLE-SPLIT COMPLETE ==="
