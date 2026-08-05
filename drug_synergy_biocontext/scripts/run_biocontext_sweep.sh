#!/usr/bin/env bash
# Core MLP0-width sweep driver for the bio-context model.
# Everything is env-parameterised so the same script serves 10k/full data, CV/single-split,
# and progeny/kegg/progeny_kegg. The thin wrappers in this directory just set the env.
#
#   DIMS         space-separated MLP0 widths            (default "0 1 2 4 8 14")
#   SPLITS       space-separated split strategies       (default all three)
#   BIO_CONTEXT  progeny | kegg | progeny_kegg          (default progeny)
#   CV_FOLDS     10 = cross-validation, 1 = single split(default 10)
#   MAX_SAMPLES  row cap; EMPTY STRING = full dataset   (default 10000)
#   EPOCHS LR HIDDEN_DIMS DROPOUT BATCH_SIZE            (default: train.py macro preset)
#   OUT_ROOT     output directory root                  (default outputs/<BIO_CONTEXT>)
set -u

BASE="/Users/arnoldcheskis/Documents/Projects/drug_discovery/drug_synergy_biocontext"
PY="$HOME/.venvs/drug_synergy/bin/python"

DIMS=(${DIMS:-0 1 2 4 8 14})
SPLITS=(${SPLITS:-random cell_line drug_and_cell_line})
BIO_CONTEXT="${BIO_CONTEXT:-progeny}"
CV_FOLDS="${CV_FOLDS:-10}"
MAX_SAMPLES="${MAX_SAMPLES-10000}"   # note: ${VAR-default}, so MAX_SAMPLES="" means FULL dataset
OUT_ROOT="${OUT_ROOT:-outputs/$BIO_CONTEXT}"
LOG_PREFIX="${LOG_PREFIX:-biocontext}"
cd "$BASE" || exit 1

# Optional hyperparameter overrides; unset -> train.py's macro preset decides.
COMMON_ARGS=()
[ -n "${EPOCHS:-}" ]      && COMMON_ARGS+=(--epochs "$EPOCHS")
[ -n "${LR:-}" ]          && COMMON_ARGS+=(--lr "$LR")
[ -n "${DROPOUT:-}" ]     && COMMON_ARGS+=(--dropout "$DROPOUT")
[ -n "${BATCH_SIZE:-}" ]  && COMMON_ARGS+=(--batch-size "$BATCH_SIZE")
[ -n "${HIDDEN_DIMS:-}" ] && COMMON_ARGS+=(--hidden-dims ${HIDDEN_DIMS})
# Only pass --max-samples when non-empty; omitting it trains on the FULL dataset.
[ -n "$MAX_SAMPLES" ]     && COMMON_ARGS+=(--max-samples "$MAX_SAMPLES")

if [ "$CV_FOLDS" -gt 1 ]; then
  EVAL_ARGS=(--cv-folds "$CV_FOLDS" --cv-seeds 42)
  DONE_FILE="cv_metrics.json"
  DONE_MARKER="Saved CV summary"
else
  EVAL_ARGS=(--cv-folds 1 --seed 42)
  DONE_FILE="metrics.json"
  DONE_MARKER="Saved metrics to"
fi

echo "=== sweep: context=$BIO_CONTEXT cv_folds=$CV_FOLDS max_samples=${MAX_SAMPLES:-FULL} dims=(${DIMS[*]}) ==="

run_one () {
  local out="$1"; shift
  if [ -f "$out/$DONE_FILE" ]; then echo "[skip] $out"; return 0; fi
  local log="/tmp/${LOG_PREFIX}_$(basename "$out").log"
  echo "[run ] $out"
  "$PY" -u -m src.train \
    --output-dir "$out" \
    "${EVAL_ARGS[@]}" ${COMMON_ARGS[@]+"${COMMON_ARGS[@]}"} \
    "$@" > "$log" 2>&1
  if grep -q "$DONE_MARKER" "$log"; then
    echo "[ ok ] $out"
  else
    echo "[FAIL] $out  (see $log)"; tail -3 "$log"
  fi
}

for split in "${SPLITS[@]}"; do
  for K in "${DIMS[@]}"; do
    args=(--split-strategy "$split" --bio-context "$BIO_CONTEXT" --mlp0-out-dim "$K")
    # K=0 is the no-cell-line-branch control.
    [ "$K" = "0" ] && args+=(--no-use-gene-expression)
    run_one "$OUT_ROOT/mlp0${K}_${split}" "${args[@]}"
  done
done

echo "=== SWEEP COMPLETE: $OUT_ROOT ==="
