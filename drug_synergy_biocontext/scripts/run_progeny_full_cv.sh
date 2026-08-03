#!/usr/bin/env bash
# PROGENy, FULL dataset (no row cap), 10-fold CV, all three splits.
#
# Hyperparameters match the codex baseline (epochs=10, lr=0.001, hidden=[512,256],
# dropout=0.2, batch=64) so these results are directly comparable to the full-dataset
# PCA plots in drug_synergy_baseline/outputs/training_curves/compression_*_full.png.
#
#   bash /Users/arnoldcheskis/Documents/Projects/drug_discovery/drug_synergy_biocontext/scripts/run_progeny_full_cv.sh
set -u
export BIO_CONTEXT=progeny
export DIMS="${DIMS:-0 1 2 4 8 14}"
export SPLITS="${SPLITS:-random cell_line drug_and_cell_line}"
export CV_FOLDS=10
export MAX_SAMPLES=""          # empty => full dataset
export EPOCHS=10 LR=0.001 HIDDEN_DIMS="512 256" DROPOUT=0.2 BATCH_SIZE=64
export OUT_ROOT="outputs/progeny_full_cv"
export LOG_PREFIX="progeny_full_cv"
exec bash "$(dirname "$0")/run_biocontext_sweep.sh"
