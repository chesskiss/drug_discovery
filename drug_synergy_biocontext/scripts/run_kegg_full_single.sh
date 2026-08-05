#!/usr/bin/env bash
# KEGG, FULL dataset (no row cap), SINGLE split, all three splits.
# Same hyperparameters as run_kegg_full_cv.sh -- one train per point instead of ten,
# so this is the cheap early read (hours, not days) on whether the wide end degrades.
#
#   bash /Users/arnoldcheskis/Documents/Projects/drug_discovery/drug_synergy_biocontext/scripts/run_kegg_full_single.sh
set -u
export BIO_CONTEXT=kegg
export DIMS="${DIMS:-0 1 2 4 8 14 32 58 128 256 336}"
export SPLITS="${SPLITS:-random cell_line drug_and_cell_line}"
export CV_FOLDS=1
export MAX_SAMPLES=""          # empty => full dataset
export EPOCHS=10 LR=0.001 HIDDEN_DIMS="512 256" DROPOUT=0.2 BATCH_SIZE=64
export OUT_ROOT="outputs/kegg_full_single"
export LOG_PREFIX="kegg_full_single"
exec bash "$(dirname "$0")/run_biocontext_sweep.sh"
