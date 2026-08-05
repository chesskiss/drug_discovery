#!/usr/bin/env bash
# Fast wiring check for the KEGG bio-context pipeline before committing days of compute.
# Dims 0 / 8 / 336 -- 336 is included deliberately so the WIDEST MLP0 is exercised.
# Tiny sample cap, 2 epochs, single split. Expect a couple of minutes.
#
#   bash /Users/arnoldcheskis/Documents/Projects/drug_discovery/drug_synergy_biocontext/scripts/smoke_kegg.sh
set -u

BASE="/Users/arnoldcheskis/Documents/Projects/drug_discovery/drug_synergy_biocontext"
PY="$HOME/.venvs/drug_synergy/bin/python"
DIMS=(0 8 336)
cd "$BASE" || exit 1

for K in "${DIMS[@]}"; do
  out="outputs/smoke_kegg/mlp0${K}_random"
  echo "[run ] $out"
  extra=()
  [ "$K" = "0" ] && extra+=(--no-use-gene-expression)
  "$PY" -u -m src.train \
    --output-dir "$out" \
    --bio-context kegg \
    --mlp0-out-dim "$K" \
    --split-strategy random \
    --cv-folds 1 --seed 42 --max-samples 2000 --epochs 2 \
    "${extra[@]+"${extra[@]}"}" > "/tmp/smoke_kegg_mlp0${K}.log" 2>&1
  if grep -q "Saved metrics to" "/tmp/smoke_kegg_mlp0${K}.log"; then
    echo "[ ok ] $out"
  else
    echo "[FAIL] $out  (see /tmp/smoke_kegg_mlp0${K}.log)"; tail -5 "/tmp/smoke_kegg_mlp0${K}.log"
  fi
done

echo "=== generating plot ==="
"$PY" -m src.biocontext_rmse_plot \
  --sweep-dir outputs/smoke_kegg \
  --eval-name "KEGG smoke (single split)" \
  --metric rmse \
  --output outputs/training_curves/smoke_kegg.png

echo "=== KEGG SMOKE COMPLETE ==="
echo "  metrics: $BASE/outputs/smoke_kegg/mlp0{0,8,336}_random/metrics.json"
echo "  plot   : $BASE/outputs/training_curves/smoke_kegg.png"
