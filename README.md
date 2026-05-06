# Drug Discovery Project Map

This repo studies drug synergy prediction, with the current focus on improving the cell-line / gene-expression branch while keeping a controlled baseline.

## Architecture

- [drug_synergy_baseline](/Users/arnoldcheskis/Documents/Projects/drug_discovery/drug_synergy_baseline): main training pipeline. `src/train.py` runs the baseline MLP, `dataset.py` builds splits/datasets, `data_loading.py` reads DrugComb-derived data, `model.py` defines the model, `results/` tracks canonical experiments.
- [data_compression](/Users/arnoldcheskis/Documents/Projects/drug_discovery/data_compression): cell-line feature preprocessing and compression. Includes the `zscore_var_pca_128` pipeline and data inspection utilities.
- [cell_line_difficulty](/Users/arnoldcheskis/Documents/Projects/drug_discovery/cell_line_difficulty): side analysis that scores how easy/hard each cell line appears for synergy discovery.
- [reverse_engineering](/Users/arnoldcheskis/Documents/Projects/drug_discovery/reverse_engineering): placeholder area for model reverse-engineering / attribution work.
- [sota_reference_model](/Users/arnoldcheskis/Documents/Projects/drug_discovery/sota_reference_model): separate reference implementation kept apart from the controlled baseline.
- [RESEARCH_JOURNAL.md](/Users/arnoldcheskis/Documents/Projects/drug_discovery/RESEARCH_JOURNAL.md): detailed research log. [RESEARCH_PLAN_AND_MILESTONES.md](/Users/arnoldcheskis/Documents/Projects/drug_discovery/RESEARCH_PLAN_AND_MILESTONES.md): planning document.

## Implemented So Far

- End-to-end DeepSynergy-style baseline for `Synergy_ZIP` is working.
- Baseline supports `random`, `drug`, `cell_line`, `drug_pair`, and `drug_and_cell_line` splits.
- Baseline supports no-gene, raw/filtered built-in cell-line features, and PCA128-compressed features.
- Canonical experiment tracking is in [drug_synergy_baseline/results/baseline_experiments.csv](/Users/arnoldcheskis/Documents/Projects/drug_discovery/drug_synergy_baseline/results/baseline_experiments.csv) with a compact summary in [baseline_summary.csv](/Users/arnoldcheskis/Documents/Projects/drug_discovery/drug_synergy_baseline/results/baseline_summary.csv).
- Compression pipeline for `z-score -> variance top-k -> PCA128` exists under [data_compression/zscore_var_pca_128](/Users/arnoldcheskis/Documents/Projects/drug_discovery/data_compression/zscore_var_pca_128).
- Cell-line-only descriptive and predictive diagnostics are implemented in [cell_line_difficulty](/Users/arnoldcheskis/Documents/Projects/drug_discovery/cell_line_difficulty).
- OCA / attribution scripts exist in the baseline package for saved runs and CV folds.
- 10-fold CV OCA for the combined PCA128 baseline is under [drug_synergy_baseline/outputs/long_cv10_10k_pca128_drug_and_cell_line_practical/oca_cv](/Users/arnoldcheskis/Documents/Projects/drug_discovery/drug_synergy_baseline/outputs/long_cv10_10k_pca128_drug_and_cell_line_practical/oca_cv). It aggregates fold-wise component importance and stability; per-fold OCA artifacts live under each `fold_runs/.../oca/` directory.

## Baseline Output Name Meanings

Outputs under `drug_synergy_baseline/outputs/` are intentionally descriptive:

- `short_*`: single train/val/test run from the canonical short comparison matrix.
- `long_cv10_10k_*`: 10-fold evaluation on a 10k subset.
- `random`: ordinary random row split.
- `drug`: unseen-drug split.
- `cell_line`: unseen-cell-line split.
- `drug_and_cell_line`: combined group-aware split; test rows contain held-out drugs or held-out cell lines.
- `no_genes`: drug-only baseline.
- `pca128`: PCA-compressed cell-line features.
- `practical`: run used the `practical_research` macro preset.

Examples:

- `short_no_genes_random_practical`: short single-split, drug-only, random split, practical preset.
- `long_cv10_10k_pca128_practical`: 10-fold CV on 10k rows, PCA128 features, random split, practical preset.
- `long_cv10_10k_pca128_drug_and_cell_line_practical`: 10-fold group-aware evaluation on 10k rows with PCA128 features and the combined unseen-drug + unseen-cell-line split.
