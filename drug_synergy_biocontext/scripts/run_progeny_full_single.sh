#!/usr/bin/env bash
# PROGENy, FULL dataset (no row cap), SINGLE split, all three splits.
# Same codex-matched hyperparameters as the CV variant; much faster (one train per point).
#
#   bash /Users/arnoldcheskis/Documents/Projects/drug_discovery/drug_synergy_biocontext/scripts/run_progeny_full_single.sh
set -u
export BIO_CONTEXT=progeny
export DIMS="${DIMS:-0 1 2 4 8 14}"
export SPLITS="${SPLITS:-random cell_line drug_and_cell_line}"
export CV_FOLDS=1
export MAX_SAMPLES=""          # empty => full dataset
export EPOCHS=10 LR=0.001 HIDDEN_DIMS="512 256" DROPOUT=0.2 BATCH_SIZE=64
export OUT_ROOT="outputs/progeny_full_single"
export LOG_PREFIX="progeny_full_single"
exec bash "$(dirname "$0")/run_biocontext_sweep.sh"
