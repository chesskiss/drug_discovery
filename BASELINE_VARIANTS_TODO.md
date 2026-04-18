# Baseline Variant TODO

Stage ID: `1.1`

Goal: establish how strong and stable the current baseline is before changing gene compression.

This file tracks the current research stage:

- `Stage 1.1` = baseline variant diagnostics before any new gene-compression model is introduced
- `Stage 1.3A` = first explicit compression export: `z-score -> variance top-k -> PCA128`

## Completed Setup

- [x] Extend the baseline runner so it can compare:
  - gene expression on vs off
  - CellLine view `0/1/2` from the pickle (`23808 / 3171 / 627`)
  - split strategies: `random`, `cell_line`, `drug`, `drug_pair`
  - single split and `k`-fold CV runs

## Experiments To Run

- [ ] `Genes vs No Genes`
  - Run baseline with raw genes:
    - `cd /Users/arnoldcheskis/Documents/Projects/drug_discovery/drug_synergy_baseline`
    - `uv run python -m src.train --output-dir outputs/genes_on_random --split-strategy random --cell-feature-view 0 --epochs 10 --seed 42`
  - Run the same baseline without genes:
    - `cd /Users/arnoldcheskis/Documents/Projects/drug_discovery/drug_synergy_baseline`
    - `uv run python -m src.train --output-dir outputs/genes_off_random --split-strategy random --no-use-gene-expression --epochs 10 --seed 42`
  - Compare:
    - `outputs/genes_on_random/metrics.json`
    - `outputs/genes_off_random/metrics.json`
  - Decision:
    - if genes do not beat no-genes consistently, do not move to compression claims yet

- [ ] `Current Compression Baseline`
  - Compare the built-in cell-line views directly:
    - raw: `--cell-feature-view 0`
    - filtered: `--cell-feature-view 1`
    - compact/pathway-like: `--cell-feature-view 2`
  - Suggested commands:
    - `uv run python -m src.train --output-dir outputs/view0_raw --split-strategy random --cell-feature-view 0 --epochs 10 --seed 42`
    - `uv run python -m src.train --output-dir outputs/view1_filtered --split-strategy random --cell-feature-view 1 --epochs 10 --seed 42`
    - `uv run python -m src.train --output-dir outputs/view2_compact --split-strategy random --cell-feature-view 2 --epochs 10 --seed 42`
  - Compare:
    - `test_rmse`
    - `test_pearson`
    - gap to `test_mean_baseline_mse`

- [ ] `Cross-Validation`
  - Run 10-fold CV on a smaller subset first for speed:
    - `uv run python -m src.train --output-dir outputs/cv10_view1_10k --split-strategy random --cell-feature-view 1 --cv-folds 10 --stratified-cv --cv-seeds 42 --max-samples 10000 --epochs 5`
  - If runtime is acceptable, repeat on the full aligned dataset:
    - `uv run python -m src.train --output-dir outputs/cv10_view1_full --split-strategy random --cell-feature-view 1 --cv-folds 10 --stratified-cv --cv-seeds 42 --epochs 10`
  - Review:
    - `outputs/cv10_view1_10k/cv_metrics.json`
    - `outputs/cv10_view1_full/cv_metrics.json`

- [ ] `Test-Set Variants`
  - Random row split:
    - `uv run python -m src.train --output-dir outputs/split_random --split-strategy random --cell-feature-view 1 --epochs 10 --seed 42`
  - Unseen cell-line split:
    - `uv run python -m src.train --output-dir outputs/split_cell_line --split-strategy cell_line --cell-feature-view 1 --epochs 10 --seed 42`
  - Unseen drug split:
    - `uv run python -m src.train --output-dir outputs/split_drug --split-strategy drug --cell-feature-view 1 --epochs 10 --seed 42`
  - Optional stricter pair split:
    - `uv run python -m src.train --output-dir outputs/split_drug_pair --split-strategy drug_pair --cell-feature-view 1 --epochs 10 --seed 42`
  - Deliverable:
    - one table comparing performance drop from `random` to `cell_line` and `drug`

- [ ] `Training Stability / Hyperparameter Sanity`
  - Check whether current settings are too aggressive:
    - `uv run python -m src.train --output-dir outputs/sanity_lr_1e3 --split-strategy random --cell-feature-view 1 --epochs 10 --lr 1e-3 --seed 42`
    - `uv run python -m src.train --output-dir outputs/sanity_lr_3e4 --split-strategy random --cell-feature-view 1 --epochs 10 --lr 3e-4 --seed 42`
    - `uv run python -m src.train --output-dir outputs/sanity_lr_1e4 --split-strategy random --cell-feature-view 1 --epochs 20 --lr 1e-4 --seed 42`
  - Inspect:
    - train/val curves in each `metrics.json` history block
    - whether validation loss is still decreasing at the last epoch

## Reporting Format

- [ ] Build `baseline_variant_summary.csv` with columns:
  - `experiment_name`
  - `split_strategy`
  - `use_gene_expression`
  - `cell_feature_view`
  - `epochs`
  - `lr`
  - `seed`
  - `test_rmse`
  - `test_pearson`
  - `test_spearman`
  - `test_mean_baseline_mse`
  - `notes`

- [ ] Add a short interpretation block to `RESEARCH_JOURNAL.md` after the first full comparison round:
  - whether genes help at all
  - which cell-line view is the best current baseline
  - whether the model fails mainly on unseen drugs or unseen cell lines
  - whether current hyperparameters are obviously unstable

## Next Compression Step

- [ ] `Stage 1.3A: z-score -> variance top-k -> PCA128`
  - Build compressed export from the repository root:
    - `/Users/arnoldcheskis/Documents/Projects/drug_discovery/drug_synergy_baseline/.venv/bin/python -m data_compression.zscore_var_pca_128.build`
  - Copy the generated file back into baseline data when ready:
    - `cp /Users/arnoldcheskis/Documents/Projects/drug_discovery/data_compression/zscore_var_pca_128/data/drugcomb.csv /Users/arnoldcheskis/Documents/Projects/drug_discovery/drug_synergy_baseline/data/drugcomb.csv`
  - Run baseline on the compressed file:
    - `cd /Users/arnoldcheskis/Documents/Projects/drug_discovery/drug_synergy_baseline`
    - `uv run python -m src.train --synergy-path data/drugcomb.csv --cell-expression-path data/drugcomb.csv --cell-feature-view 0 --output-dir outputs/pca128_random --split-strategy random --epochs 10 --seed 42`
  - Record:
    - output gene dimension
    - test RMSE / Pearson / Spearman
    - whether the compressed branch still collapses to constant predictions
