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

Refresh the tables after any run completes:

```bash
uv run python -m src.results_tracking
```

## Run pipeline

```bash
uv run python -m src.train --synergy-path data/drugcomb.csv --fallback-pickle-path data/drugcomb.pkl --output-dir outputs
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
- gene/cell-line source: `data/drugcomb.pkl`
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
- `--train-fraction`, `--val-fraction`
- `--max-samples`

Cross-validation supports:

- `--split-strategy random`
- `--split-strategy drug_and_cell_line`

For `drug_and_cell_line`, CV is group-aware repeated holdout rather than classic row-disjoint CV.

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

Cross-validation outputs:

- `outputs/<run_name>/cv_metrics.json`
- `outputs/<run_name>/cv_runs.csv`
- `outputs/<run_name>/cv_test_predictions.csv`
- `outputs/<run_name>/fold_runs/<fold_name>/history.csv`
- `outputs/training_curves/<run_name>/fold_runs/<fold_name>/loss_curve.png`

### What to compare

For Stage `1.1`, compare these fields first:

- `test_rmse`
- `test_pearson`
- `test_spearman`
- `test_mean_baseline_mse`
- training history in `metrics.json`

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
- because CSV `CellLine` payloads are truncated, the compression pipeline reads raw vectors from [drugcomb.pkl](/Users/arnoldcheskis/Documents/Projects/drug_discovery/drug_synergy_baseline/data/drugcomb.pkl)
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
