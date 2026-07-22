# Drug Synergy Baseline

Minimal starting point for dataset utilities and baseline integration.

## Research Stages

- `Stage 1.0`: data grounding and initial baseline pipeline
- `Stage 1.1`: baseline variant diagnostics before gene-compression changes
- `Stage 1.2`: modular gene encoder replacement
- `Stage 1.3`: compression baseline comparisons
- `Stage 1.3A`: `z-score -> variance top-k -> PCA128` compression export

## Setup

Install `uv`, then from `drug_synergy_baseline/` run:

```bash
uv sync
```

Run scripts with:

```bash
uv run python -m <module>
```

## Hyperparameter Macros

The baseline runner now reads its default training hyperparameters from [macros.toml](/Users/arnoldcheskis/Documents/Projects/drug_discovery/drug_synergy_baseline/macros.toml).

Available presets:

- `practical_research`: stable working default for current experiments
- `deep_synergy_reference`: literature-inspired preset based on DeepSynergy
- `fast_debug`: cheap smoke-test preset

Use a preset:

```bash
uv run python -m src.train --macro-preset practical_research
```

Override individual values on top of the preset:

```bash
uv run python -m src.train --macro-preset practical_research --lr 3e-4 --epochs 20
```

The `deep_synergy_reference` preset is adapted from the DeepSynergy paper, which reported that a small learning rate (`1e-5`), conic hidden layers, and dropout were important for performance on this task family. Source: [DeepSynergy (PMC)](https://pmc.ncbi.nlm.nih.gov/articles/PMC5925774/) and [PubMed](https://pubmed.ncbi.nlm.nih.gov/29253077/).

## Results Tracking

Use the canonical experiment tracker in [results/README.md](/Users/arnoldcheskis/Documents/Projects/drug_discovery/drug_synergy_baseline/results/README.md).

[`results/baseline_summary.csv`](/Users/arnoldcheskis/Documents/Projects/drug_discovery/drug_synergy_baseline/results/baseline_summary.csv) — one row per model (`no_genes`, `pca128`) comparing test MSE/Pearson across split strategies and the 10-fold CV runs, including the shared-holdout columns when a run was trained with `--holdout-mode additional`.

Refresh the tables after any run completes:

```bash
uv run python -m src.results_tracking
```

Backfill regression plots from existing saved predictions without retraining:

```bash
uv run python -m src.backfill_regression_plots --run long_cv10_10k_pca128_drug_and_cell_line_practical --kind cv
```

Backfill all available regression plots:

```bash
uv run python -m src.backfill_regression_plots --kind all
```

## OCA / Attribution

Run component-level occlusion attribution on a saved single-run checkpoint:

```bash
uv run python -m src.oca \
  --model-path outputs/pca128_random/baseline_mlp.pt \
  --config-path outputs/pca128_random/config.json
```

Artifacts are written under `outputs/<run>/oca/`:

- `component_importance.csv`
- `local_explanations.csv`
- `component_importance_topk.png`
- `component_importance_head_tail_summary.png`
- `local_explanations_heatmap.png`
- `oca_summary.json`

Meaning:

- `component_importance.csv`: component ranking for that run
- `local_explanations.csv`: per-sample occlusion deltas for selected rows
- `component_importance_topk.png`: top helpful components
- `component_importance_head_tail_summary.png`: helpful head plus harmful tail
- `local_explanations_heatmap.png`: selected rows versus important components
- `oca_summary.json`: OCA run metadata

Run fold-wise OCA for a CV directory:

```bash
uv run python -m src.oca_cv \
  --cv-output-dir outputs/long_cv10_10k_pca128_drug_and_cell_line_practical
```

This expects each fold run to contain `baseline_mlp.pt`. Older CV outputs created before this checkpoint-saving patch must be rerun first.

Each fold directory also contains:

- `baseline_mlp.pt`: trained fold model weights
- `config.json`: exact fold settings
- `metrics.json`: fold-level evaluation metrics
- `history.csv`: train/val loss per epoch
- `test_predictions.csv`: fold test predictions
- `val_predictions.csv`: fold validation predictions
- `oca/component_importance.csv`: component ranking for that fold
- `oca/local_explanations.csv`: per-sample occlusion deltas for selected rows
- `oca/component_importance_topk.png`: top helpful components for that fold
- `oca/component_importance_head_tail_summary.png`: helpful head plus harmful tail for that fold
- `oca/local_explanations_heatmap.png`: selected rows versus important components for that fold
- `oca/oca_summary.json`: fold OCA metadata

CV OCA aggregate artifacts are written under `outputs/<cv_run>/oca_cv/`:

- `oca_cv_component_importance_per_fold.csv`
- `oca_cv_component_summary.csv`
- `oca_cv_topk_stability.png`
- `oca_cv_fold_component_heatmap.png`
- `oca_cv_topk_frequency.png`
- `oca_cv_summary.json`

Meaning:

- `oca_cv_component_importance_per_fold.csv`: all fold OCA results stacked together
- `oca_cv_component_summary.csv`: mean and spread of component importance across folds
- `oca_cv_topk_stability.png`: helpful head, compressed zero block, and harmful tail with mean dots and min/max intervals
- `oca_cv_fold_component_heatmap.png`: fold-by-component heatmap for the top aggregate components
- `oca_cv_topk_frequency.png`: how often each component appears in the fold top-k
- `oca_cv_summary.json`: aggregate OCA metadata

Current `oca_cv_topk_stability.png` semantics:

- dot: mean `mean_delta_squared_error` across folds
- vertical black interval: min to max fold value for that component
- zero-valued components are compressed into one block

## Run pipeline

```bash
uv run python -m src.train --synergy-path data/drugcomb.csv --fallback-pickle-path ../data_compression/source_data/drugcomb.pkl --output-dir outputs
```

## Top-Level Program Flow

```text
Raw DrugComb / pickle data
        |
        v
data_compression/
download / inspect / compress cell-line features
        |
        v
drug_synergy_baseline/data/drugcomb.csv
prepared training table
        |
        v
drug_synergy_baseline/src/train.py
load data -> split -> build datasets -> train/eval model
        |
        v
outputs/
metrics, predictions, checkpoints, plots
        |
        v
results_tracking.py
canonical experiment tables
```

## Main Baseline Flow

```text
CLI command
uv run python -m src.train
        |
        v
[train.py]
parse args + macro preset
        |
        v
[dataset.py]
load aligned rows + build split
        |
        v
[data_loading.py]
read CSV / pickle and cell-line features
        |
        v
[dataset.py]
convert rows into PyTorch datasets
        |
        v
[model.py]
build DeepSynergyMLP
        |
        v
[train.py]
train loop / validation / test / CV
        |
        v
[training_artifacts.py]
save metrics, history, curves, regression plots
        |
        v
outputs/<run_name>/
```

## Stage 1.1: Baseline Variant Diagnostics

This is the current step.

Goal:

- measure whether gene expression helps at all
- compare the existing `CellLine` feature views
- compare easier vs harder test splits
- check whether the baseline is stable enough before changing compression

Run all commands from `/Users/arnoldcheskis/Documents/Projects/drug_discovery/drug_synergy_baseline`.

Default data behavior after you place a prepared compressed file back into baseline `data/`:

- synergy rows: `data/drugcomb.csv`
- gene/cell-line source: `../data_compression/source_data/drugcomb.pkl`
- use the pickle for gene-on experiments because CSV `CellLine` payloads may be truncated

### Current runner features

The training entrypoint now supports:

- `--use-gene-expression` or `--no-use-gene-expression`
- `--gene-feature-set raw|filtered|compact`
- `--cell-feature-view 0|1|2`
- `--split-strategy random|cell_line|drug|drug_pair|drug_and_cell_line`
- `--cv-folds <k>`
- `--cv-seeds <seed ...>`
- `--stratified-cv`
- `--holdout-test-fraction <f>`, `--holdout-seed <n>`, `--holdout-mode instead|additional`
- `--train-fraction`, `--val-fraction`
- `--max-samples`

Cross-validation supports:

- `--split-strategy random`
- `--split-strategy drug_and_cell_line`

For `drug_and_cell_line`, CV is group-aware repeated holdout rather than classic row-disjoint CV.

#### Shared global holdout test set

By default every fold builds its own cold-start test set, so each fold is graded on a *different* set of
rows. To grade all folds on the *same* reserved set, use a global holdout that is carved out once, before
any fold splitting, using the same cold drug + cold cell-line exclusion logic:

- `--holdout-test-fraction <f>` — group fraction reserved as the shared test set. `0` disables (default).
- `--holdout-seed <n>` — seed for the carve. Deliberately separate from `--cv-seeds` so the holdout is
  identical across every CV seed.
- `--holdout-mode instead|additional`
  - `instead` (default): the shared holdout **is** each fold's test set; the fold's own group bucket is used
    for validation only.
  - `additional`: each fold keeps its own cold test set **and** is also scored on the shared holdout as an
    extra metric.

> **Group fraction ≠ row fraction.** Because a row is reserved when *either* drug *or* its cell line is
> held out, `--holdout-test-fraction 0.1` removes roughly `1 − 0.9³ ≈ 27%` of *rows*, not 10%. The runner
> prints the realized row fraction and records it in `cv_metrics.json` (`global_holdout`). This is the same
> `or`-based behavior the per-fold splits already have.

Extra outputs when a holdout is enabled: `global_holdout_groups.json` (the exact reserved drugs/cell lines
and counts), `global_holdout_test.csv` (the reserved rows), and `cv_holdout_ensemble.csv` plus an ensemble
regression plot (per-fold predictions on the shared holdout, averaged into a mean/std ensemble). In
`additional` mode each fold dir also gets `holdout_predictions.csv`.

```bash
uv run python -m src.train \
  --output-dir outputs/cv10_drug_and_cell_line_10k_holdout \
  --split-strategy drug_and_cell_line \
  --cell-feature-view 1 \
  --cv-folds 10 \
  --cv-seeds 42 \
  --max-samples 10000 \
  --epochs 5 \
  --holdout-test-fraction 0.1 \
  --holdout-mode instead
```

### Fastest first runs

Compare genes on vs off:

```bash
uv run python -m src.train \
  --output-dir outputs/genes_on_random \
  --split-strategy random \
  --cell-feature-view 0 \
  --epochs 10 \
  --seed 42

uv run python -m src.train \
  --output-dir outputs/genes_off_random \
  --split-strategy random \
  --no-use-gene-expression \
  --epochs 10 \
  --seed 42
```

Compare current built-in cell-line views:

```bash
uv run python -m src.train --output-dir outputs/view0_raw --split-strategy random --gene-feature-set raw --epochs 10 --seed 42
uv run python -m src.train --output-dir outputs/view1_filtered --split-strategy random --gene-feature-set filtered --epochs 10 --seed 42
uv run python -m src.train --output-dir outputs/view2_compact --split-strategy random --gene-feature-set compact --epochs 10 --seed 42
```

Try the filtered-gene baseline explicitly:

```bash
uv run python -m src.train \
  --output-dir outputs/genes_filtered_random \
  --split-strategy random \
  --gene-feature-set filtered \
  --epochs 10 \
  --seed 42
```

Compare held-out split strategies:

```bash
uv run python -m src.train --output-dir outputs/split_random --split-strategy random --cell-feature-view 1 --epochs 10 --seed 42
uv run python -m src.train --output-dir outputs/split_cell_line --split-strategy cell_line --cell-feature-view 1 --epochs 10 --seed 42
uv run python -m src.train --output-dir outputs/split_drug --split-strategy drug --cell-feature-view 1 --epochs 10 --seed 42
uv run python -m src.train --output-dir outputs/split_drug_and_cell_line --split-strategy drug_and_cell_line --cell-feature-view 1 --epochs 10 --seed 42
uv run python -m src.train --output-dir outputs/split_drug_pair --split-strategy drug_pair --cell-feature-view 1 --epochs 10 --seed 42
```

For `drug_and_cell_line`, rows are routed with `test > val > train` priority:

- if either drug or the cell line is in the test bucket, the row goes to test
- else if either drug or the cell line is in the val bucket, the row goes to val
- else the row goes to train

Run 10-fold CV on a smaller subset first:

```bash
uv run python -m src.train \
  --output-dir outputs/cv10_view1_10k \
  --split-strategy random \
  --cell-feature-view 1 \
  --cv-folds 10 \
  --stratified-cv \
  --cv-seeds 42 \
  --max-samples 10000 \
  --epochs 5
```

Run group-aware drug + cell-line CV on a smaller subset:

```bash
uv run python -m src.train \
  --output-dir outputs/cv10_drug_and_cell_line_10k \
  --split-strategy drug_and_cell_line \
  --cell-feature-view 1 \
  --cv-folds 10 \
  --cv-seeds 42 \
  --max-samples 10000 \
  --epochs 5
```

### Outputs to inspect

Single-run outputs:

- `outputs/<run_name>/metrics.json`
- `outputs/<run_name>/config.json`
- `outputs/<run_name>/history.csv`
- `outputs/<run_name>/val_predictions.csv`
- `outputs/<run_name>/test_predictions.csv`
- `outputs/training_curves/<run_name>/loss_curve.png`
- `outputs/training_curves/<run_name>/test_regression.png`

Cross-validation outputs:

- `outputs/<run_name>/cv_metrics.json`
- `outputs/<run_name>/cv_runs.csv`
- `outputs/<run_name>/cv_test_predictions.csv`
- `outputs/<run_name>/fold_runs/<fold_name>/history.csv`
- `outputs/training_curves/<run_name>/fold_runs/<fold_name>/loss_curve.png`
- `outputs/training_curves/<run_name>/fold_runs/<fold_name>/test_regression.png`
- `outputs/training_curves/<run_name>/cv_test_regression.png`

### What to compare

For Stage `1.1`, compare these fields first:

- `test_rmse`
- `test_pearson`
- `test_spearman`
- `test_mean_baseline_mse`
- training history in `metrics.json`

Regression plot meaning:

- x-axis = true test `Synergy_ZIP`
- y-axis = model prediction for the same point
- each dot = one held-out test example
- dashed diagonal = perfect prediction
- tighter concentration around the diagonal means better calibration and lower error
- a horizontal band means the model is collapsing toward near-constant predictions

Run single-sample prediction with:

```bash
uv run python -m src.predict --row-idx 0
```

or

```bash
uv run python -m src.predict --smiles-a "..." --smiles-b "..." --cell-line "786-0"
```

## Stage 1.3A: Data Compression Export

This stage moves dataset download and data-compression work out of the baseline package into the root-level [data_compression](/Users/arnoldcheskis/Documents/Projects/drug_discovery/data_compression) module.

Current method:

- method id: `1.3A`
- pipeline: `raw CellLine[0] -> variance top-k (3000) -> z-score -> PCA -> pad to 128`

Important:

- the raw CSV now lives at [drugcomb_raw.csv](/Users/arnoldcheskis/Documents/Projects/drug_discovery/data_compression/zscore_var_pca_128/data/drugcomb_raw.csv)
- because CSV `CellLine` payloads are truncated, the compression pipeline reads raw vectors from [drugcomb.pkl](/Users/arnoldcheskis/Documents/Projects/drug_discovery/data_compression/source_data/drugcomb.pkl)
- the pipeline output is [drugcomb.csv](/Users/arnoldcheskis/Documents/Projects/drug_discovery/data_compression/zscore_var_pca_128/data/drugcomb.csv)
- once validated, that output can be copied or moved back into `drug_synergy_baseline/data/drugcomb.csv` for baseline training

Run the compression export from the repository root:

```bash
/Users/arnoldcheskis/Documents/Projects/drug_discovery/drug_synergy_baseline/.venv/bin/python \
  -m data_compression.zscore_var_pca_128.build
```

Then place the result into baseline data:

```bash
cp /Users/arnoldcheskis/Documents/Projects/drug_discovery/data_compression/zscore_var_pca_128/data/drugcomb.csv \
   /Users/arnoldcheskis/Documents/Projects/drug_discovery/drug_synergy_baseline/data/drugcomb.csv
```

Train against the compressed export:

```bash
cd /Users/arnoldcheskis/Documents/Projects/drug_discovery/drug_synergy_baseline
uv run python -m src.train \
  --synergy-path data/drugcomb.csv \
  --cell-expression-path data/drugcomb.csv \
  --cell-feature-view 0 \
  --output-dir outputs/pca128_random \
  --split-strategy random \
  --epochs 10 \
  --seed 42
```

## Existing utilities

### Visualize data

```bash
cd /Users/arnoldcheskis/Documents/Projects/drug_discovery
/Users/arnoldcheskis/Documents/Projects/drug_discovery/drug_synergy_baseline/.venv/bin/python -m data_compression.visualize \
  --input data_compression/zscore_var_pca_128/data/drugcomb_raw.csv \
  --output-dir outputs/visualization \
  --preprocess
```

### Download dataset

```bash
cd /Users/arnoldcheskis/Documents/Projects/drug_discovery
/Users/arnoldcheskis/Documents/Projects/drug_discovery/drug_synergy_baseline/.venv/bin/python -m data_compression.download_drugcomb_tdc
```

This downloader saves:

- `data_compression/zscore_var_pca_128/data/drugcomb_raw.csv`
- `data_compression/zscore_var_pca_128/data/reduced_drugcomb.csv`
- `data_compression/zscore_var_pca_128/data/tdc_cell_features.csv`

Important:

- `tdc_cell_features.csv` is a TDC cell-feature table
- it is not the same thing as the embedded 3-view `CellLine` payload export
- if you want to compare `CellLine[0]`, `CellLine[1]`, and `CellLine[2]`, you still need a non-truncated source for those payloads

## Baseline model

We keep a separate SOTA reference implementation at [sota_reference_model](/Users/arnoldcheskis/Documents/Projects/drug_discovery/sota_reference_model), based on MatchMaker ([tastanlab/matchmaker](https://github.com/tastanlab/matchmaker)).

For now:
- it is separate from the lightweight `src/` baseline
- it preserves its own data expectations and TensorFlow stack
- it should be compared against the current baseline, not merged into the same runner
