# Research Journal

## Objective

Main objective: improve the gene expression encoder for drug synergy prediction, with emphasis on compression quality and biological relevance.

Current strategy:

1. Get a working end-to-end baseline
2. Replace only the gene expression branch
3. Compare encoder variants while keeping the downstream pipeline fixed

Current stage map:

- `Stage 1.0`: data grounding and initial pipeline validation
- `Stage 1.1`: baseline variant diagnostics
- `Stage 1.2`: modular gene encoder replacement
- `Stage 1.3`: compression baseline comparisons
- `Stage 1.3A`: `z-score -> variance top-k -> PCA128` export pipeline

## Dataset Notes

Primary dataset:

- `data_compression/zscore_var_pca_128/data/drugcomb_raw.csv`
- `drug_synergy_baseline/data/drugcomb.pkl`

Primary target:

- `Synergy_ZIP`

Why ZIP:

- standard synergy target
- reasonably balanced compared with some alternatives
- appropriate first regression label

Important columns in the synergy table:

- `Drug1_ID`: drug A identifier
- `Drug2_ID`: drug B identifier
- `Cell_Line_ID`: cell line name / ID
- `CSS`: combination sensitivity score
- `Synergy_ZIP`: main regression label
- `Synergy_Bliss`: alternative synergy metric
- `Synergy_Loewe`: alternative synergy metric
- `Synergy_HSA`: alternative synergy metric
- `Drug1`: SMILES string for drug A
- `Drug2`: SMILES string for drug B
- `CellLine`: embedded cell-line feature payload

Current layout note:

- `drugcomb_raw.csv` is the raw export file kept under `data_compression/`
- `drugcomb.pkl` is the non-truncated source for embedded `CellLine` arrays
- future compressed exports are written as `data_compression/.../data/drugcomb.csv` and can then be moved into baseline `data/`

## CellLine Interpretation

The `CellLine` payload contains 3 aligned views of the same cell line representation.

Current confirmed lengths:

- `23808`: raw/high-dimensional view
- `3171`: filtered/compressed view
- `627`: smaller pathway-like / summarized view

Interpretation for now:

- same index should correspond to the same gene/feature across samples
- `2.3` appears often and may represent a floor / low-expression value
- this suggests a future sparse or threshold-based encoding path may be useful

## What We Implemented

### Stage 1.0A. Minimal baseline pipeline

Implemented a DeepSynergy-style baseline in `src/`:

- drug A features
- drug B features
- raw gene expression vector
- concatenation
- MLP predictor

Relevant files:

- `src/data_loading.py`
- `src/dataset.py`
- `src/model.py`
- `src/train.py`
- `src/predict.py`

What this achieved:

- training works end-to-end
- prediction works for a single sample
- dataset loading is based on TDC DrugComb data, not MatchMaker preprocessing

Stored outputs:

- `outputs/baseline_mlp.pt`
- `outputs/config.json`
- `outputs/metrics.json`
- `outputs/val_predictions.csv`

Current baseline metrics from `outputs/metrics.json`:

- `train_samples = 237678`
- `val_samples = 29709`
- `test_samples = 29711`
- `drug_dim = 256`
- `gene_dim = 23808`
- `test_mse = 29.0515`
- `test_rmse = 5.3899`

Important note:

- current baseline is functional, but performance is not yet clearly better than a trivial baseline
- this is acceptable for `Stage 1.0` because the main goal was pipeline validation

### Stage 1.0B. Visualization / data inspection

Relevant outputs:

- `outputs/visualization/summary.json`
- `outputs/visualization/target_distribution.png`
- `outputs/visualization/top_cells.png`
- `outputs/visualization/top_pair_frequency_hist.png`

Figures:

![Target Distribution](outputs/visualization/target_distribution.png)

![Top Cells](outputs/visualization/top_cells.png)

![Top Pair Frequency](outputs/visualization/top_pair_frequency_hist.png)

### Stage 1.0C. Cell-line-only difficulty analysis

Implemented a separate module:

- `../cell_line_difficulty/`

Purpose:

- ignore the two drugs
- aggregate observed ZIP values per cell line
- estimate which cell lines appear easier or harder to find synergistic combinations for

Relevant outputs:

- `../cell_line_difficulty/outputs/ranked_by_ease.csv`
- `../cell_line_difficulty/outputs/ranked_by_difficulty.csv`
- `../cell_line_difficulty/outputs/summary.json`

Current summary:

- total rows analyzed: `297098`
- total cell lines: `59`
- easiest examples: `HL-60(TB)`, `MOLT-4`, `NCI-H460`
- hardest examples: `UO-31`, `HOP-92`, `MDA-MB-231`

Interpretation:

- this is not true biological “treatability”
- it is an observational screen-level heuristic
- it mixes cell-line susceptibility with which drugs and combinations were tested

### Stage 1.0D. Cell-line-only predictive model

Built a predictive model on top of the cell-line-only difficulty target.

Purpose:

- use cell-line features only
- predict how easy it may be to find synergistic drug combinations for that cell line

Inputs:

- target source: `drug_synergy_baseline/data/drugcomb.csv`
- feature source: non-truncated `CellLine` payload source
- feature view used: `CellLine[1]`
- feature dimension: `3171`

Modeling setup:

- one supervised sample per cell line
- total samples: `59`
- target: composite `ease_score`
- evaluation: LOOCV
- models compared:
  - ridge regression
  - tiny MLP
  - fold-mean baseline

Run command:

```bash
python3 -m cell_line_difficulty.src.cell_line_difficulty.predict_cli --synergy-path drug_synergy_baseline/data/drugcomb.csv --pickle-path drug_synergy_baseline/data/drugcomb.pkl --output-dir cell_line_difficulty/outputs
```

Relevant outputs:

- `../cell_line_difficulty/outputs/predictive_dataset.csv`
- `../cell_line_difficulty/outputs/loocv_predictions.csv`
- `../cell_line_difficulty/outputs/predictive_metrics.json`

Current predictive metrics from `../cell_line_difficulty/outputs/predictive_metrics.json`:

- `sample_count = 59`
- `feature_dimension = 3171`
- `ridge rmse = 0.7631`
- `ridge mae = 0.6001`
- `ridge pearson = 0.5094`
- `ridge spearman = 0.3778`
- `mlp rmse = 0.8370`
- `fold-mean baseline rmse = 0.9006`

How to interpret this:

- the target `ease_score` has `std ≈ 0.893`
- so `rmse = 0.7631` is not small in an absolute sense
- however, ridge improves over the baseline:
  - `0.9006 -> 0.7631`
- this suggests the filtered cell-line features contain real predictive signal
- ridge outperforming the tiny MLP is expected here:
  - only `59` samples
  - high-dimensional input
  - linear regularization is more stable

Important caution:

- the fold-mean baseline correlation is not meaningful in LOOCV
- because the leave-one-out training mean is mechanically anti-correlated with the held-out target
- RMSE/MAE are the more useful comparison here

## Current Stage

### Stage 1.1. Baseline variant diagnostics

This is the current step before any new gene-compression model is introduced.

Goal:

- determine whether gene expression helps over a drug-only baseline
- determine which current `CellLine` view is the strongest baseline
- determine whether the model mainly fails on unseen cell lines or unseen drugs
- determine whether the present training setup is stable enough to use as a reference

Implemented support in `drug_synergy_baseline/src/train.py`:

- gene expression on/off:
  - `--use-gene-expression`
  - `--no-use-gene-expression`
- cell-line feature-view selection:
  - `--cell-feature-view 0|1|2`
- split strategy selection:
  - `--split-strategy random|cell_line|drug|drug_pair`
- cross-validation:
  - `--cv-folds`
  - `--cv-seeds`
  - `--stratified-cv`

Primary outputs for this stage:

- `drug_synergy_baseline/outputs/<run_name>/metrics.json`

Current findings:

- raw genes (`CellLine[0]`, `23808`) produced weak regression quality:
  - `test_mse ~= 29.05`
  - near-zero correlation
  - predictions collapsed toward a constant
- no-gene baseline was materially stronger:
  - `test_mse ~= 20.25`
  - `test_pearson ~= 0.55`
- filtered built-in view (`CellLine[1]`, `3171`) still did not fix the collapse

Interpretation:

- the current issue is not just dimensionality in the abstract
- the present gene branch is likely unstable or poorly conditioned relative to the amount of unique cell-line information
- that justifies moving to an explicit compression pipeline before testing more complex encoders

### Stage 1.3A. z-score -> variance top-k -> PCA128

Implemented a separate root-level module:

- `data_compression/`

Current method directory:

- `data_compression/zscore_var_pca_128/`

Pipeline:

1. start from raw `CellLine[0]` (`23808`)
2. rank features by raw weighted variance
3. keep the top `3000`
4. z-score the selected features
5. fit weighted PCA
6. export a baseline-compatible `drugcomb.csv` with one compressed `CellLine` view

Important constraint:

- this dataset currently has only `59` unique cell lines
- so PCA cannot provide `128` informative dimensions
- the pipeline writes the informative PCA dimensions first and zero-pads the rest to reach `128`

Files:

- raw CSV: `data_compression/zscore_var_pca_128/data/drugcomb_raw.csv`
- compressed output: `data_compression/zscore_var_pca_128/data/drugcomb.csv`
- metadata: `data_compression/zscore_var_pca_128/data/metadata.json`

Status:

- compression export pipeline implemented
- compressed baseline run still pending
- `drug_synergy_baseline/outputs/<run_name>/test_predictions.csv`
- `drug_synergy_baseline/outputs/<run_name>/cv_metrics.json`
- `drug_synergy_baseline/outputs/<run_name>/cv_runs.csv`

Immediate runs for this stage:

```bash
cd /Users/arnoldcheskis/Documents/Projects/drug_discovery/drug_synergy_baseline

uv run python -m src.train --output-dir outputs/genes_on_random --split-strategy random --cell-feature-view 0 --epochs 10 --seed 42
uv run python -m src.train --output-dir outputs/genes_off_random --split-strategy random --no-use-gene-expression --epochs 10 --seed 42

uv run python -m src.train --output-dir outputs/view0_raw --split-strategy random --cell-feature-view 0 --epochs 10 --seed 42
uv run python -m src.train --output-dir outputs/view1_filtered --split-strategy random --cell-feature-view 1 --epochs 10 --seed 42
uv run python -m src.train --output-dir outputs/view2_compact --split-strategy random --cell-feature-view 2 --epochs 10 --seed 42
```

## Immediate Research Questions

### 1. Baseline quality

- Why is the baseline near trivial-performance territory?
- Is the raw gene representation too noisy?
- Is the simple drug representation too weak?
- Are we normalizing gene expression appropriately?

### 2. Gene expression meaning

- What exactly does each index correspond to?
- Can we recover gene names / metadata?
- Which preprocessing steps were applied before these views were produced?

### 3. Compression strategy

The main research direction remains:

- better compression of gene expression
- without losing biological meaning
- while improving downstream synergy prediction

### 4. Cell-line-only predictor value

- Can the filtered `3171`-feature view predict screening difficulty better than a trivial baseline?
- Does this remain true if the target is changed from composite `ease_score` to:
  - mean ZIP
  - high-hit rate (`ZIP > 10`)?
- Is the signal biological, or mostly a consequence of how the screen was designed?

## Planned Next Steps

### Stage 1.2. Modular encoder

Replace raw gene input with:

- `z_cell = encoder(gene_expr)`

Keep fixed:

- drug representation branch
- fusion style
- predictor head

First version:

- MLP encoder

Later versions:

- autoencoder bottleneck
- graph-informed encoder

### Sparse-expression idea

One concrete idea to test:

- keep only values `> 2.3`
- also keep their original indices

Possible structures:

- sparse index-value list
- dictionary / hash table
- graph-ready node-feature representation

Why this may help:

- removes obvious low-expression floor values
- preserves which genes remain active
- may reduce noise before compression

Open issues:

- need to confirm biological meaning of `2.3`
- need mapping from index to gene identity
- need to compare against standard variance filtering

## Milestones

### Completed

- TDC DrugComb data loaded
- baseline pipeline implemented
- prediction CLI implemented
- cell-line-only difficulty module implemented
- cell-line-only predictive model implemented and tested

### In progress

- understand why baseline performance is weak
- understand gene-expression feature semantics

### Next milestone

- modular gene encoder with fixed downstream model
- Result:
  - `________________`

### After that

- test expression filtering / sparse encoding
- Result:
  - `________________`

### Later

- autoencoder baseline
- graph-informed encoder
- biological prior integration
- Result:
  - `________________`

## Useful Paths

- Baseline code: `src/`
- SOTA reference model: `sota_reference_model/`
- Main synergy dataset: `drug_synergy_baseline/data/drugcomb.csv`
- Baseline metrics: `outputs/metrics.json`
- Baseline predictions: `outputs/val_predictions.csv`
- Visual summaries: `outputs/visualization/`
- Cell-line analysis module: `../cell_line_difficulty/`
- Cell-line analysis outputs: `../cell_line_difficulty/outputs/`
- Predictive metrics: `../cell_line_difficulty/outputs/predictive_metrics.json`
- Predictive dataset: `../cell_line_difficulty/outputs/predictive_dataset.csv`
- Predictive LOOCV predictions: `../cell_line_difficulty/outputs/loocv_predictions.csv`
