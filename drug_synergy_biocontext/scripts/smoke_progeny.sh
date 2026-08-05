#!/usr/bin/env bash
# Fast end-to-end smoke test for the PROGENy bio-context pipeline.
# Tiny sample cap, 2 epochs, single split, 3 MLP0 widths -> then a plot.
# Expect a couple of minutes total.
#
#   bash /Users/arnoldcheskis/Documents/Projects/drug_discovery/drug_synergy_biocontext/scripts/smoke_progeny.sh
set -u

BASE="/Users/arnoldcheskis/Documents/Projects/drug_discovery/drug_synergy_biocontext"
PY="$HOME/.venvs/drug_synergy/bin/python"
DIMS=(0 2 14)
cd "$BASE" || exit 1

for K in "${DIMS[@]}"; do
  out="outputs/smoke/mlp0${K}_random"
  echo "[run ] $out"
  extra=()
  [ "$K" = "0" ] && extra+=(--no-use-gene-expression)
  "$PY" -u -m src.train \
    --output-dir "$out" \
    --bio-context progeny \
    --mlp0-out-dim "$K" \
    --split-strategy random \
    --cv-folds 1 --seed 42 --max-samples 2000 --epochs 2 \
    "${extra[@]+"${extra[@]}"}" > "/tmp/smoke_mlp0${K}.log" 2>&1
  if grep -q "Saved metrics to" "/tmp/smoke_mlp0${K}.log"; then
    echo "[ ok ] $out"
  else
    echo "[FAIL] $out  (see /tmp/smoke_mlp0${K}.log)"; tail -5 "/tmp/smoke_mlp0${K}.log"
  fi
done

echo "=== generating plot ==="
"$PY" -m src.biocontext_rmse_plot \
  --sweep-dir outputs/smoke \
  --eval-name "smoke (single split)" \
  --metric rmse \
  --output outputs/training_curves/smoke_progeny.png

echo "=== SMOKE COMPLETE ==="
echo "  metrics: $BASE/outputs/smoke/mlp0{0,2,14}_random/metrics.json"
echo "  plot   : $BASE/outputs/training_curves/smoke_progeny.png"
