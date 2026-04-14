# Research Journal

## Objective

Main objective: improve the gene expression encoder for drug synergy prediction, with emphasis on compression quality and biological relevance.

Current strategy:

1. Get a working end-to-end baseline
2. Replace only the gene expression branch
3. Compare encoder variants while keeping the downstream pipeline fixed

## Dataset Notes

Primary dataset:

- `data/drugcomb_synergy.csv`
- `data/drugcomb_with_smiles.csv`
- `data/drugcomb.pkl`

Primary target:

- `Synergy_ZIP`

Why ZIP:

- standard synergy target
- reasonably balanced compared with some alternatives
- appropriate first regression label

Important columns in `data/drugcomb_synergy.csv`:

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

### 1. Minimal baseline pipeline

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
- this is acceptable for Step 1A because the main goal was pipeline validation

### 2. Visualization / data inspection

Relevant outputs:

- `outputs/visualization/summary.json`
- `outputs/visualization/target_distribution.png`
- `outputs/visualization/top_cells.png`
- `outputs/visualization/top_pair_frequency_hist.png`

Figures:

![Target Distribution](outputs/visualization/target_distribution.png)

![Top Cells](outputs/visualization/top_cells.png)

![Top Pair Frequency](outputs/visualization/top_pair_frequency_hist.png)

### 3. Cell-line-only difficulty analysis

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

### 4. Cell-line-only predictive model

Built a predictive model on top of the cell-line-only difficulty target.

Purpose:

- use cell-line features only
- predict how easy it may be to find synergistic drug combinations for that cell line

Inputs:

- target source: `data/drugcomb_synergy.csv`
- feature source: `data/drugcomb.pkl`
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
python3 -m cell_line_difficulty.src.cell_line_difficulty.predict_cli --synergy-path drug_synergy_baseline/data/drugcomb_synergy.csv --pickle-path drug_synergy_baseline/data/drugcomb.pkl --output-dir cell_line_difficulty/outputs
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

### Step 1B: modular encoder

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
- MatchMaker reference: `matchmaker/`
- Main dataset: `data/drugcomb_synergy.csv`
- Baseline metrics: `outputs/metrics.json`
- Baseline predictions: `outputs/val_predictions.csv`
- Visual summaries: `outputs/visualization/`
- Cell-line analysis module: `../cell_line_difficulty/`
- Cell-line analysis outputs: `../cell_line_difficulty/outputs/`
- Predictive metrics: `../cell_line_difficulty/outputs/predictive_metrics.json`
- Predictive dataset: `../cell_line_difficulty/outputs/predictive_dataset.csv`
- Predictive LOOCV predictions: `../cell_line_difficulty/outputs/loocv_predictions.csv`
