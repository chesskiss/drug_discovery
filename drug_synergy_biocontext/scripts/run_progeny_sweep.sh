#!/usr/bin/env bash
# PROGENy, 10k-row subset, 10-fold CV, all three splits.
# Uses train.py's macro preset (practical_research) -- comparable to the earlier 10k
# PCA sweeps in drug_synergy_baseline/outputs/sweep/.
# For the full dataset use run_progeny_full_cv.sh / run_progeny_full_single.sh instead.
#
#   bash /Users/arnoldcheskis/Documents/Projects/drug_discovery/drug_synergy_biocontext/scripts/run_progeny_sweep.sh
set -u
export BIO_CONTEXT=progeny
export DIMS="${DIMS:-0 1 2 4 8 14}"
export SPLITS="${SPLITS:-random cell_line drug_and_cell_line}"
export CV_FOLDS=10
export MAX_SAMPLES="${MAX_SAMPLES-10000}"
export OUT_ROOT="${OUT_ROOT:-outputs/progeny}"
export LOG_PREFIX="progeny"
exec bash "$(dirname "$0")/run_biocontext_sweep.sh"
