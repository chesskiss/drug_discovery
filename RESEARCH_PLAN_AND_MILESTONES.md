# Research Plan And Milestones

## Goal

Improve the cell-line / gene-expression representation used for drug synergy prediction on TDC DrugComb.

The current research strategy is deliberately staged:

1. Establish a weak but reliable baseline.
2. Measure whether genes help under multiple train/test split assumptions.
3. Replace or compress only the gene branch while keeping the rest of the model fixed.
4. Use attribution / reverse-engineering methods to understand what the model uses from the gene input.

## Stage Map

- `Stage 1.0`: data grounding and initial baseline pipeline
- `Stage 1.1`: baseline variant diagnostics
- `Stage 1.1A`: canonical results tracking and split-matrix evaluation
- `Stage 1.1B`: group-aware CV for unseen drugs and unseen cell lines
- `Stage 1.2`: modular gene encoder replacement
- `Stage 1.3`: compression baseline comparisons
- `Stage 1.3A`: first explicit compression export with `z-score -> variance top-k -> PCA128`
- `Stage 1.4`: MLP reverse-engineering and gene-attribution analysis

## Current Baseline Direction

The baseline is intentionally simple:

- featurize `Drug1` and `Drug2`
- optionally append a cell-line gene-expression vector
- train one MLP regressor on `Synergy_ZIP`

This is not meant to beat every published method. It is meant to provide a controlled reference where only the cell-line/gene representation changes.

Important comparison principle:

- compare models under the same macro preset
- compare genes vs no genes under the same split
- compare random, unseen-drug, unseen-cell-line, and combined unseen-drug + unseen-cell-line splits
- reserve expensive CV for the most important candidates

## Dataset Notes

Primary DrugComb fields:

- `Drug1`: SMILES for drug A
- `Drug2`: SMILES for drug B
- `Drug1_ID`: drug A identifier
- `Drug2_ID`: drug B identifier
- `Cell_Line_ID`: cell-line identifier
- `CellLine`: embedded cell-line feature payload
- `Synergy_ZIP`: main regression target
- `Synergy_Bliss`, `Synergy_Loewe`, `Synergy_HSA`: alternative synergy metrics

Confirmed `CellLine` views:

- `CellLine[0]`: `23808` raw/high-dimensional features
- `CellLine[1]`: `3171` filtered features
- `CellLine[2]`: `627` compact/pathway-like features

Important limitation:

- this DrugComb extract has only `59` unique cell lines
- many rows share the same cell-line vectors
- this makes high-dimensional gene encoders easy to overfit and makes naive MLP interpretation fragile

## Completed

### Stage 1.0A. Data Grounding And Visualization

Status: done.

Completed:

- loaded TDC DrugComb data
- confirmed the important columns and target
- generated visualization outputs for target distribution and dataset structure
- moved data visualization/manipulation scripts out of the baseline package and into `data_compression/`

Relevant outputs:

- `data_compression/visualize.py`
- `data_compression/outputs/` or prior visualization outputs, depending on the run

### Stage 1.0B. Baseline Pipeline

Status: done.

Completed:

- implemented DeepSynergy-style baseline training
- added prediction CLI
- added deterministic drug featurization
- added optional gene branch
- added train/validation/test metrics
- added artifact saving:
  - `metrics.json`
  - `config.json`
  - `history.csv`
  - prediction CSVs
  - model checkpoint
  - train-vs-validation loss curve

Key files:

- `drug_synergy_baseline/src/data_loading.py`
- `drug_synergy_baseline/src/dataset.py`
- `drug_synergy_baseline/src/model.py`
- `drug_synergy_baseline/src/train.py`
- `drug_synergy_baseline/src/predict.py`
- `drug_synergy_baseline/src/training_artifacts.py`

Initial raw-gene baseline:

- gene setting: raw `CellLine[0]`
- gene dim: `23808`
- split: random
- test MSE: `29.0515`
- test RMSE: `5.3899`
- behavior: near-constant predictions / near-zero correlation

### Stage 1.0C. Cell-Line-Only Difficulty Analysis

Status: done.

Completed:

- created `cell_line_difficulty/`
- ignored molecules and ranked cell lines by observed ZIP-derived ease/difficulty
- built a small cell-line-only prediction task

Current LOOCV result:

- samples: `59` cell lines
- features: `3171`
- ridge RMSE: `0.7631`
- tiny MLP RMSE: `0.8370`
- fold-mean baseline RMSE: `0.9006`

Interpretation:

- the filtered cell-line representation contains some predictive signal
- ridge beating MLP is expected because there are only `59` unique cell lines
- this task is diagnostic, not a direct measure of biological treatability

### Stage 1.1. Baseline Variant Diagnostics

Status: mostly done; remaining work is interpretation and a few planned long-CV rows.

Implemented runner support:

- genes on/off
- `CellLine` view selection
- gene feature set aliases: raw, filtered, compact
- split strategies:
  - `random`
  - `drug`
  - `cell_line`
  - `drug_pair`
  - `drug_and_cell_line`
- cross-validation:
  - random stratified CV
  - group-aware `drug_and_cell_line` CV
- macro presets through `drug_synergy_baseline/macros.toml`
- canonical results tracking

Important results so far:

| Model | Split | Test MSE | Pearson |
|---|---:|---:|---:|
| raw genes | random | `29.0515` | near zero / undefined |
| filtered genes | random | `29.0467` | near zero / undefined |
| no genes | random | `20.2502` | `0.5505` |
| PCA128 genes | random | `19.0419` | `0.5892` |
| no genes | drug | `27.9074` | `0.0316` |
| PCA128 genes | drug | `25.0928` | `0.1239` |
| no genes | cell line | `21.8703` | `0.5614` |
| PCA128 genes | cell line | `36.0284` | `0.4094` |
| no genes | drug + cell line | `24.7152` | `0.3800` |
| PCA128 genes | drug + cell line | `31.9951` | `0.3032` |
| no genes | random CV10 10k | `16.9083` | `0.4778` |
| PCA128 genes | random CV10 10k | `14.9237` | `0.5650` |

Current interpretation:

- raw `23808` genes do not help the current MLP baseline
- built-in filtered `3171` genes also do not help
- PCA128 improves random split and unseen-drug split
- PCA128 hurts unseen-cell-line and combined unseen-drug + unseen-cell-line splits
- no-gene baseline is still very competitive, so any gene encoder must prove it helps out-of-distribution, not only on random rows

### Stage 1.1A. Results Tracking

Status: done.

Canonical source-of-truth file:

- `drug_synergy_baseline/results/baseline_experiments.csv`

Derived human-readable summary:

- `drug_synergy_baseline/results/baseline_summary.csv`

Tracker command:

```bash
cd /Users/arnoldcheskis/Documents/Projects/drug_discovery/drug_synergy_baseline
uv run python -m src.results_tracking
```

Current canonical comparison set:

- `no_genes`
- `pca128`

Archived negatives:

- `raw_genes`
- `filtered_genes`

### Stage 1.1B. Split Strategy Work

Status: implemented; long group-aware CV rows are planned but not completed.

Implemented split strategies:

- `random`: normal row split
- `drug`: held-out drugs; if a drug is assigned to validation/test, rows containing it are routed away from train
- `cell_line`: held-out cell lines
- `drug_pair`: held-out drug pairs
- `drug_and_cell_line`: held-out drugs and held-out cell lines together

`drug_and_cell_line` semantics:

- split unique drugs into train/validation/test buckets
- split unique cell lines into train/validation/test buckets
- route each row by priority:
  - `test` if either drug or the cell line is in a test bucket
  - `val` if no test membership exists but either drug or the cell line is in a validation bucket
  - `train` otherwise

Reason for priority routing:

- train remains clean from validation/test drugs and cell lines
- mixed rows are not dropped
- this preserves more data while still testing generalization to unseen groups

Group-aware CV:

- supported for `--split-strategy drug_and_cell_line`
- fold `i` uses drug fold `i` and cell-line fold `i` as test
- validation uses fold `(i + 1) % cv_folds`
- this is group-aware repeated holdout CV, not classic row-disjoint CV
- rows can appear in test in more than one fold because each row has multiple group memberships

### Stage 1.3A. PCA128 Compression Export

Status: done and tested.

Module:

- `data_compression/`

Method:

- `data_compression/zscore_var_pca_128/`

Pipeline:

1. read raw `CellLine[0]` vectors from `drugcomb.pkl`
2. rank genes/features by weighted variance
3. keep top `3000`
4. z-score selected features
5. fit weighted PCA
6. export one compressed `CellLine` view compatible with the baseline

Result:

- output dim: `128`
- effective informative PCA rank is limited by `59` unique cell lines
- PCA dimensions beyond the effective rank are zero-padded
- random split MSE improved from no-genes `20.2502` to PCA128 `19.0419`
- random CV10 10k improved from no-genes `16.9083` to PCA128 `14.9237`
- OOD cell-line and combined OOD split results are worse for PCA128 than no-genes

## Remaining Work

### Immediate

1. Complete canonical short random reruns under `practical_research`

Rows still marked planned:

- `short_no_genes_random_practical`
- `short_pca128_random_practical`

Why:

- the current random values are legacy 10-epoch runs
- the other short-phase rows use the current `practical_research` preset
- rerunning makes the matrix internally consistent

2. Complete group-aware CV10 10k rows

Rows still marked planned:

- `long_cv10_10k_no_genes_drug_and_cell_line_practical`
- `long_cv10_10k_pca128_drug_and_cell_line_practical`

Why:

- random CV answers stability under random row sampling
- group-aware CV answers generalization to unseen drugs and unseen cell lines

3. Write interpretation after the remaining rows finish

Questions to answer:

- Do genes help on random split?
- Do genes help on unseen drugs?
- Do genes help on unseen cell lines?
- Do genes help when both drugs and cell lines are unseen?
- Is PCA128 a real improvement or mainly a random-split improvement?

### Stage 1.2. Modular Gene Encoder

Status: planned.

Replace direct gene input with:

```text
z_cell = encoder(gene_expr)
prediction = MLP(drug_a, drug_b, z_cell)
```

Keep fixed:

- drug featurization
- train/test splits
- macro preset
- downstream fusion/predictor as much as possible

Candidate encoders:

- MLP encoder
- autoencoder bottleneck
- denoising autoencoder
- sparse / thresholded gene encoder
- graph-informed encoder using biological priors

### Stage 1.3. Additional Compression Baselines

Status: planned.

Candidate methods from the compression TODO:

- `z-score -> Var3k -> PCA512 -> AE128`
- `z-score -> BioFilter1k-3k -> AE128`
- `raw -> AE128`

Priority:

1. `z-score -> Var3k -> PCA512 -> AE128`: best practical upgrade from current PCA128.
2. `z-score -> BioFilter1k-3k -> AE128`: more interpretable but depends on reliable gene mapping.
3. `raw -> AE128`: hardest to stabilize; avoid until simpler methods are exhausted.

### Stage 1.4. MLP Reverse-Engineering And Attribution

Status: planned for the upcoming weeks.

Goal:

- understand which gene dimensions or compressed components influence the MLP prediction
- distinguish real cell-line signal from artifacts of the split, scaling, or repeated cell-line identities

Primary analysis: occlusion / perturbation attribution.

Terminology note:

- the meeting note said `OCA`
- in this plan, that is treated as occlusion-based contribution analysis / occlusion attribution
- if the instructor meant a more specific named method, this section should be renamed later, but the implementation idea is the same: perturb inputs and measure prediction/error change

Procedure:

1. train a frozen baseline model
2. mask, zero, replace, or permute one gene/group/component at a time
3. measure the change in output error or prediction:
   - delta MSE
   - delta absolute error
   - delta prediction for selected examples
4. rank features by how much the perturbation changes performance

Why this is useful:

- directly tests whether a feature matters to the trained predictor
- works for raw genes, PCA components, or learned embeddings
- can be extended from single features to gene groups/pathways

Important caveats:

- correlated genes make single-feature occlusion hard to interpret
- masking can create out-of-distribution inputs
- PCA components are not directly biological unless mapped back through loadings
- repeated cell-line vectors can make attribution look stronger than it really is

Neuron-level analysis is lower priority.

Reason:

- hidden neurons are not stable semantic units
- equivalent networks can represent the same function with rotated or rescaled hidden activations
- individual neurons can mix many unrelated signals
- dead/saturated activations can look important or unimportant for misleading reasons

Complementary methods:

- permutation importance
- integrated gradients
- gradient times input
- SHAP-style analysis on compressed features, if runtime allows
- component loadings for PCA128
- group/pathway occlusion after gene-index mapping is available

Attention/gating idea:

- add a small attention or gating module over gene groups/components
- tune attention temperature / scaling in the `Q * K` interaction to control how sharp or soft the attention is
- note: in standard scaled dot-product attention, dividing by a larger temperature makes attention softer, while a smaller divisor makes it sharper
- if the goal is to force robust selection under noise, explicit attention dropout, logit noise, or entropy regularization may be cleaner than only changing the scale factor

## Can `BASELINE_VARIANTS_TODO.md` Be Deleted?

Yes, after this file is committed, the useful information from `BASELINE_VARIANTS_TODO.md` has been consolidated here.

What was preserved:

- Stage `1.1` goal
- genes vs no-genes comparison
- raw/filtered/compact `CellLine` view comparison
- random, drug, cell-line, drug-pair split idea
- CV-on-10k-first principle
- training-stability/hyperparameter sanity principle
- reporting-table idea, now replaced by `baseline_experiments.csv`
- PCA128 compression next step, now completed as Stage `1.3A`

What is obsolete:

- old 10-epoch command list
- old `baseline_variant_summary.csv` plan
- `view1` as the main competitive model
- old assumption that only random CV exists

If you want to keep an audit trail, rename it to:

```text
archive/BASELINE_VARIANTS_TODO_2026-04-18.md
```

If you want less clutter, it is safe to delete after confirming this document is sufficient.

## Working Commands

Refresh results:

```bash
cd /Users/arnoldcheskis/Documents/Projects/drug_discovery/drug_synergy_baseline
uv run python -m src.results_tracking
```

Run no-genes single split:

```bash
uv run python -m src.train \
  --macro-preset practical_research \
  --output-dir outputs/short_no_genes_random_practical \
  --split-strategy random \
  --no-use-gene-expression \
  --seed 42
```

Run PCA128 single split:

If `drug_synergy_baseline/data/drugcomb.csv` is missing, first place the PCA128 export into the baseline data directory:

```bash
cp /Users/arnoldcheskis/Documents/Projects/drug_discovery/data_compression/zscore_var_pca_128/data/drugcomb.csv \
   /Users/arnoldcheskis/Documents/Projects/drug_discovery/drug_synergy_baseline/data/drugcomb.csv
```

```bash
uv run python -m src.train \
  --macro-preset practical_research \
  --output-dir outputs/short_pca128_random_practical \
  --split-strategy random \
  --cell-expression-path data/drugcomb.csv \
  --cell-feature-view 0 \
  --seed 42
```

Run no-genes group-aware CV:

```bash
uv run python -m src.train \
  --macro-preset practical_research \
  --output-dir outputs/long_cv10_10k_no_genes_drug_and_cell_line_practical \
  --split-strategy drug_and_cell_line \
  --no-use-gene-expression \
  --cv-folds 10 \
  --cv-seeds 42 \
  --max-samples 10000 \
  --seed 42
```

Run PCA128 group-aware CV:

```bash
uv run python -m src.train \
  --macro-preset practical_research \
  --output-dir outputs/long_cv10_10k_pca128_drug_and_cell_line_practical \
  --split-strategy drug_and_cell_line \
  --cell-expression-path data/drugcomb.csv \
  --cell-feature-view 0 \
  --cv-folds 10 \
  --cv-seeds 42 \
  --max-samples 10000 \
  --seed 42
```

## Longer-Term Direction

Research directions after the current baseline matrix:

- map gene indices to gene names
- connect features to pathways or PPI/transcription networks
- test autoencoder and denoising-autoencoder compression
- test biological filtering before learned compression
- evaluate whether gene encoders improve OOD splits, not only random splits
- reverse-engineer the trained MLP with occlusion/perturbation attribution
- compare against the separated SOTA/reference model only after the baseline is stable
