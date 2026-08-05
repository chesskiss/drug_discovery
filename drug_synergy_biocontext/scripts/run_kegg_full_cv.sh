#!/usr/bin/env bash
# KEGG, FULL dataset (no row cap), 10-fold CV, all three splits.
#
# Hyperparameters are IDENTICAL to the completed PROGENy full-CV run
# (outputs/progeny_full_cv), so the two are directly comparable:
#   epochs=10, lr=0.001, hidden_dims=[512,256], dropout=0.2, batch_size=64
#
# Dim list keeps 14 (PROGENy's ceiling) so every PROGENy point has an exact KEGG
# counterpart, and adds 336 = all usable KEGG pathways (no compression), the KEGG
# analogue of PCA-58.
#
# 11 dims x 3 splits = 33 CV runs. At ~85 min/run (measured on PROGENy) expect ~47h.
#
#   bash /Users/arnoldcheskis/Documents/Projects/drug_discovery/drug_synergy_biocontext/scripts/run_kegg_full_cv.sh
set -u
export BIO_CONTEXT=kegg
export DIMS="${DIMS:-0 1 2 4 8 14 32 58 128 256 336}"
export SPLITS="${SPLITS:-random cell_line drug_and_cell_line}"
export CV_FOLDS=10
export MAX_SAMPLES=""          # empty => full dataset
export EPOCHS=10 LR=0.001 HIDDEN_DIMS="512 256" DROPOUT=0.2 BATCH_SIZE=64
export OUT_ROOT="outputs/kegg_full_cv"
export LOG_PREFIX="kegg_full_cv"
exec bash "$(dirname "$0")/run_biocontext_sweep.sh"
