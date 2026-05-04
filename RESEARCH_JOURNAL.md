# Research Journal

## Objective

Main objective: improve the gene-expression / cell-line encoder for drug synergy prediction while keeping the downstream baseline controlled.

Current strategy:

1. Build a reproducible weak baseline.
2. Quantify whether genes help before changing the gene encoder.
3. Compare random, unseen-drug, unseen-cell-line, and combined unseen-drug + unseen-cell-line splits.
4. Compress or replace the gene branch only after the baseline behavior is understood.
5. Reverse-engineer the trained MLP to identify which gene inputs or compressed components affect predictions.

## Stage Map

- `Stage 1.0`: data grounding and initial baseline pipeline
- `Stage 1.1`: baseline variant diagnostics
- `Stage 1.1A`: canonical results tracking and short split matrix
- `Stage 1.1B`: group-aware drug + cell-line CV
- `Stage 1.2`: modular gene encoder replacement
- `Stage 1.3`: compression baseline comparisons
- `Stage 1.3A`: `z-score -> variance top-k -> PCA128` export pipeline
- `Stage 1.4`: MLP reverse-engineering and attribution

## Dataset Notes

Primary data sources:

- `data_compression/zscore_var_pca_128/data/drugcomb_raw.csv`
- `data_compression/zscore_var_pca_128/data/drugcomb.csv`
- `data_compression/source_data/drugcomb.pkl`

Main target:

- `Synergy_ZIP`

Why ZIP:

- standard synergy target
- appropriate first regression label
- keeps the task comparable across baseline variants

Important columns:

- `Drug1_ID`: drug A identifier
- `Drug2_ID`: drug B identifier
- `Cell_Line_ID`: cell-line ID
- `Drug1`: SMILES for drug A
- `Drug2`: SMILES for drug B
- `CellLine`: embedded cell-line feature payload
- `Synergy_ZIP`: main regression target
- `Synergy_Bliss`, `Synergy_Loewe`, `Synergy_HSA`: alternative synergy targets

Confirmed `CellLine` views:

- `23808`: raw/high-dimensional view
- `3171`: filtered view
- `627`: compact/pathway-like view

Important constraint:

- the dataset has many synergy rows but only `59` unique cell lines
- this makes raw high-dimensional gene input easy to overfit
- this also means PCA128 cannot have 128 fully informative independent components

## Completed Work

### Stage 1.0A. Minimal Baseline Pipeline

Implemented a DeepSynergy-style baseline in `drug_synergy_baseline/src/`.

Inputs:

- drug A features
- drug B features
- optional gene-expression/cell-line vector

Model:

- concatenate inputs
- MLP regression head
- target = `Synergy_ZIP`

Relevant files:

- `drug_synergy_baseline/src/data_loading.py`
- `drug_synergy_baseline/src/dataset.py`
- `drug_synergy_baseline/src/model.py`
- `drug_synergy_baseline/src/train.py`
- `drug_synergy_baseline/src/predict.py`

What this achieved:

- training works end-to-end
- prediction works for a single sample
- dataset loading is based on TDC DrugComb
- MatchMaker/SOTA code is separated from the baseline direction

Initial raw-gene baseline:

- train samples: `237678`
- validation samples: `29709`
- test samples: `29711`
- drug feature dim: `256`
- gene dim: `23808`
- test MSE: `29.0515`
- test RMSE: `5.3899`

Interpretation:

- the raw-gene baseline is functional but weak
- predictions collapsed toward an almost constant output
- this motivated explicit gene compression and stronger diagnostics

### Stage 1.0B. Visualization / Data Inspection

Generated visual inspection artifacts for target distribution and dataset structure.

Earlier relevant outputs:

- `outputs/visualization/summary.json`
- `outputs/visualization/target_distribution.png`
- `outputs/visualization/top_cells.png`
- `outputs/visualization/top_pair_frequency_hist.png`

Current direction:

- data visualization and manipulation scripts belong under `data_compression/`, not under `drug_synergy_baseline/src/`
- model outputs stay under `drug_synergy_baseline/outputs/`

### Stage 1.0C. Cell-Line-Only Difficulty Analysis

Implemented a separate module:

- `cell_line_difficulty/`

Purpose:

- ignore the two drugs
- aggregate observed ZIP values per cell line
- estimate which cell lines appear easier or harder to find synergistic combinations for

Current summary:

- total rows analyzed: `297098`
- total cell lines: `59`
- easier examples: `HL-60(TB)`, `MOLT-4`, `NCI-H460`
- harder examples: `UO-31`, `HOP-92`, `MDA-MB-231`

Interpretation:

- this is not true biological treatability
- it mixes cell-line susceptibility with which drugs/combinations were tested
- it remains useful as a diagnostic side analysis

### Stage 1.0D. Cell-Line-Only Predictive Model

Built a small predictive model for cell-line screening difficulty.

Setup:

- one supervised sample per cell line
- total samples: `59`
- feature view: `CellLine[1]`
- feature dimension: `3171`
- target: composite `ease_score`
- evaluation: LOOCV

Current metrics:

- ridge RMSE: `0.7631`
- ridge MAE: `0.6001`
- ridge Pearson: `0.5094`
- ridge Spearman: `0.3778`
- tiny MLP RMSE: `0.8370`
- fold-mean baseline RMSE: `0.9006`

Interpretation:

- filtered cell-line features contain some signal
- ridge outperforming tiny MLP is expected with only `59` samples
- RMSE/MAE are more meaningful than fold-mean baseline correlation in LOOCV

## Stage 1.1. Baseline Variant Diagnostics

Goal:

- determine whether gene expression helps over a drug-only baseline
- determine whether existing `CellLine` views help
- determine whether PCA compression helps
- determine whether failures are mainly on unseen drugs, unseen cell lines, or both

Implemented support:

- `--use-gene-expression`
- `--no-use-gene-expression`
- `--cell-feature-view 0|1|2`
- `--gene-feature-set raw|filtered|compact`
- `--split-strategy random|cell_line|drug|drug_pair|drug_and_cell_line`
- `--cv-folds`
- `--cv-seeds`
- `--stratified-cv`
- `--macro-preset`
- `--max-samples`

Artifacts saved per single run:

- `metrics.json`
- `config.json`
- `history.csv`
- `val_predictions.csv`
- `test_predictions.csv`
- `baseline_mlp.pt`
- loss curve under `outputs/training_curves/`

Artifacts saved per CV run:

- `cv_metrics.json`
- `cv_runs.csv`
- `cv_test_predictions.csv`
- fold-level `history.csv`
- fold-level loss curves

## Stage 1.1A. Results Tracking

Canonical results file:

- `drug_synergy_baseline/results/baseline_experiments.csv`

Derived summary:

- `drug_synergy_baseline/results/baseline_summary.csv`

Current comparison policy:

- source of truth is the long CSV
- each row is one attempted run
- `raw_genes` and `filtered_genes` are archived negatives
- formal comparison set is `no_genes` vs `pca128`

Current summary values:

| Model | Random MSE | Drug MSE | Cell-Line MSE | Drug + Cell-Line MSE | Random CV10 10k MSE |
|---|---:|---:|---:|---:|---:|
| no genes | `20.2502` | `27.9074` | `21.8703` | `24.7152` | `16.9083` |
| PCA128 | `19.0419` | `25.0928` | `36.0284` | `31.9951` | `14.9237` |

Current summary values with Pearson:

| Model | Random Pearson | Drug Pearson | Cell-Line Pearson | Drug + Cell-Line Pearson | Random CV10 10k Pearson |
|---|---:|---:|---:|---:|---:|
| no genes | `0.5505` | `0.0316` | `0.5614` | `0.3800` | `0.4778` |
| PCA128 | `0.5892` | `0.1239` | `0.4094` | `0.3032` | `0.5650` |

Interpretation so far:

- raw genes do not help the current MLP
- built-in filtered genes do not help the current MLP
- PCA128 improves random split and unseen-drug split
- PCA128 hurts unseen-cell-line split and combined unseen-drug + unseen-cell-line split
- no-gene baseline remains strong enough that gene claims must be OOD-specific, not just random-split-specific

Rows still planned:

- `short_no_genes_random_practical`
- `short_pca128_random_practical`
- `long_cv10_10k_no_genes_drug_and_cell_line_practical`
- `long_cv10_10k_pca128_drug_and_cell_line_practical`

## Stage 1.1B. Split Semantics

Implemented split strategies:

- `random`: random row split
- `drug`: hold out unique drugs
- `cell_line`: hold out unique cell lines
- `drug_pair`: hold out unique drug pairs
- `drug_and_cell_line`: hold out unique drugs and cell lines together

Important detail for drug splits:

- both drug columns are checked
- if a drug is excluded from train, it cannot appear as either `Drug1` or `Drug2` in train

Important detail for `drug_and_cell_line`:

- rows are not dropped
- each row is routed by `test > val > train`
- if either drug or the cell line belongs to the test bucket, row goes to test
- else if either drug or the cell line belongs to the validation bucket, row goes to validation
- else row goes to train

Why this matters:

- train is clean from validation/test drugs and cell lines
- validation and test are harder than random row splits
- mixed rows are preserved rather than thrown away

CV status:

- random stratified CV is implemented
- `drug_and_cell_line` group-aware CV is implemented
- group-aware CV is repeated holdout, not row-disjoint classic CV
- validation fold is `(i + 1) % cv_folds`

## Stage 1.3A. z-score -> Variance Top-k -> PCA128

Implemented root-level module:

- `data_compression/`

Current method directory:

- `data_compression/zscore_var_pca_128/`

Pipeline:

1. start from raw `CellLine[0]` (`23808`)
2. rank features by weighted variance
3. keep top `3000`
4. z-score the selected features
5. fit weighted PCA
6. export a baseline-compatible `drugcomb.csv`
7. pad to `128` dimensions where rank is limited by the number of unique cell lines

Files:

- raw CSV: `data_compression/zscore_var_pca_128/data/drugcomb_raw.csv`
- compressed output: `data_compression/zscore_var_pca_128/data/drugcomb.csv`
- metadata: `data_compression/zscore_var_pca_128/data/metadata.json`

Result:

- random split improved vs no genes:
  - no genes MSE: `20.2502`
  - PCA128 MSE: `19.0419`
- random CV10 10k improved vs no genes:
  - no genes MSE: `16.9083`
  - PCA128 MSE: `14.9237`
- combined OOD split worsened:
  - no genes MSE: `24.7152`
  - PCA128 MSE: `31.9951`

Current interpretation:

- PCA128 may capture useful cell-line signal when test rows share the same distribution as train
- PCA128 is not yet robust for unseen cell lines
- compression alone is not enough; the model likely needs better regularization, grouping, or biological structure

## Stage 1.4. Reverse-Engineering / Attribution Plan

Motivation:

- the next phase should not only ask whether genes improve MSE
- it should also ask what the MLP is using from the gene representation

Primary method: occlusion / perturbation attribution.

This matches the meeting note about masking the input and measuring the change in output error.

Terminology note:

- the meeting note said `OCA`
- this journal currently interprets that as occlusion-based contribution analysis / occlusion attribution
- if the intended acronym was a different formal method, the implementation still starts from the same perturb-and-measure principle

Basic procedure:

1. train and freeze a baseline model
2. choose a gene feature, PCA component, or feature group
3. mask it, zero it, replace it with the train mean, or permute it across samples
4. rerun predictions
5. measure output change:
   - delta prediction
   - delta absolute error
   - delta MSE
6. rank features by the measured effect

Why this is better than first looking at individual neurons:

- hidden neurons are not guaranteed to correspond to meaningful biological concepts
- MLP hidden representations are not identifiable; equivalent functions can use rotated/rescaled hidden states
- a neuron can mix many unrelated signals
- activation magnitude does not necessarily equal causal importance
- correlated gene features can make neuron-level interpretation misleading

Additional attribution methods to consider:

- permutation importance
- integrated gradients
- gradient times input
- SHAP-style analysis on PCA/compressed features
- PCA loading analysis to map important components back to raw gene dimensions
- group/pathway occlusion after gene metadata is available

Attention/gating idea from the meeting:

- add an attention or gating mechanism over genes, gene groups, or compressed components
- inspect the learned gates/attention weights as a weak interpretability signal
- modify attention temperature/scaling in the `Q * K` interaction to control sharpness
- in standard scaled dot-product attention, larger temperature/scale divisor makes attention softer; smaller divisor makes attention sharper
- if the goal is noisy or robust selection, explicit logit noise, dropout, or entropy regularization may be cleaner than only changing the scale factor

Open caution:

- attention is not automatically explanation
- learned weights must be validated with perturbation tests

## Immediate Research Questions

### 1. Do genes help?

Current answer:

- raw and filtered genes do not help
- PCA128 helps on random and unseen-drug splits
- PCA128 hurts unseen-cell-line and combined split

Remaining answer needed:

- rerun random short rows under `practical_research`
- complete group-aware CV10 10k

### 2. Is the model learning cell-line biology or split artifacts?

Current evidence:

- PCA128 helps random/CV random
- PCA128 hurts unseen-cell-line

Interpretation:

- the model may be exploiting cell-line identities or train/test overlap patterns
- unseen-cell-line evaluation is the critical test for biological generalization

### 3. Is PCA128 the right compression?

Current answer:

- PCA128 is a useful first compression baseline
- it is not sufficient as a final gene encoder

Next compression candidates:

- `z-score -> Var3k -> PCA512 -> AE128`
- `z-score -> BioFilter1k-3k -> AE128`
- `raw -> AE128`, lower priority because it is hardest to stabilize

### 4. What should be interpreted first?

Priority order:

1. no-gene vs PCA128 prediction differences on the same rows
2. PCA component occlusion
3. map important PCA components back to raw top-variance gene dimensions
4. raw-gene or pathway-level occlusion once gene metadata is available
5. neuron-level inspection only as a secondary diagnostic

## Milestones

### Completed

- TDC DrugComb data loaded
- baseline pipeline implemented
- prediction CLI implemented
- cell-line-only difficulty module implemented
- cell-line-only predictive model implemented
- raw/filtered/no-gene/PCA128 comparisons started
- canonical results tracking implemented
- loss curves and per-epoch histories implemented
- `drug_and_cell_line` split implemented
- group-aware `drug_and_cell_line` CV implemented
- PCA128 compression pipeline implemented

### In Progress

- finish canonical short random reruns under `practical_research`
- finish group-aware CV10 10k runs
- interpret whether PCA128 is useful beyond random splits

### Next

- Stage `1.4`: implement occlusion/perturbation attribution for PCA128
- Stage `1.2`: implement modular gene encoder branch
- Stage `1.3`: test AE-based compression after PCA128 attribution

## Useful Paths

- Baseline code: `drug_synergy_baseline/src/`
- Baseline results: `drug_synergy_baseline/results/`
- Baseline outputs: `drug_synergy_baseline/outputs/`
- Training curves: `drug_synergy_baseline/outputs/training_curves/`
- Hyperparameter presets: `drug_synergy_baseline/macros.toml`
- Data compression module: `data_compression/`
- PCA128 method: `data_compression/zscore_var_pca_128/`
- SOTA/reference model: `sota_reference_model/`
- Cell-line analysis module: `cell_line_difficulty/`
