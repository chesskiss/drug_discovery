# Results Tracking

This directory holds the experiment tracking tables for the baseline-comparison funnel.

## Canonical Files

- [baseline_experiments.csv](/Users/arnoldcheskis/Documents/Projects/drug_discovery/drug_synergy_baseline/results/baseline_experiments.csv)
  - source of truth
  - one row per run
  - includes completed historical runs and planned rows
- [baseline_summary.csv](/Users/arnoldcheskis/Documents/Projects/drug_discovery/drug_synergy_baseline/results/baseline_summary.csv)
  - small derived table for quick reading
  - prefers the canonical short-phase rows when they exist

## Refresh The Tables

From `drug_synergy_baseline/`:

```bash
uv run python -m src.results_tracking
```

Run this after each experiment completes.

## Canonical Short-Phase Runs

Use one frozen preset for the 8-run matrix:

- `--macro-preset practical_research`

### No genes

```bash
uv run python -m src.train \
  --macro-preset practical_research \
  --output-dir outputs/short_no_genes_random_practical \
  --split-strategy random \
  --no-use-gene-expression \
  --seed 42

uv run python -m src.train \
  --macro-preset practical_research \
  --output-dir outputs/short_no_genes_drug_practical \
  --split-strategy drug \
  --no-use-gene-expression \
  --seed 42

uv run python -m src.train \
  --macro-preset practical_research \
  --output-dir outputs/short_no_genes_cell_line_practical \
  --split-strategy cell_line \
  --no-use-gene-expression \
  --seed 42

uv run python -m src.train \
  --macro-preset practical_research \
  --output-dir outputs/short_no_genes_drug_and_cell_line_practical \
  --split-strategy drug_and_cell_line \
  --no-use-gene-expression \
  --seed 42
```

### PCA128 genes

```bash
uv run python -m src.train \
  --macro-preset practical_research \
  --output-dir outputs/short_pca128_random_practical \
  --split-strategy random \
  --cell-expression-path data/drugcomb.csv \
  --cell-feature-view 0 \
  --seed 42

uv run python -m src.train \
  --macro-preset practical_research \
  --output-dir outputs/short_pca128_drug_practical \
  --split-strategy drug \
  --cell-expression-path data/drugcomb.csv \
  --cell-feature-view 0 \
  --seed 42

uv run python -m src.train \
  --macro-preset practical_research \
  --output-dir outputs/short_pca128_cell_line_practical \
  --split-strategy cell_line \
  --cell-expression-path data/drugcomb.csv \
  --cell-feature-view 0 \
  --seed 42

uv run python -m src.train \
  --macro-preset practical_research \
  --output-dir outputs/short_pca128_drug_and_cell_line_practical \
  --split-strategy drug_and_cell_line \
  --cell-expression-path data/drugcomb.csv \
  --cell-feature-view 0 \
  --seed 42
```

## Combined Split Semantics

The `drug_and_cell_line` split uses `test > val > train` priority routing.

- split drugs into `train/val/test`
- split cell lines into `train/val/test`
- row goes to `test` if either drug or the cell line is in the test bucket
- else row goes to `val` if either drug or the cell line is in the val bucket
- else row goes to `train`

This keeps train free of validation/test drugs and cell lines while avoiding row drops.

## Long-Phase Run

Random CV and group-aware CV are tracked separately.

### Random CV

```bash
uv run python -m src.train \
  --macro-preset practical_research \
  --output-dir outputs/long_cv10_10k_no_genes_practical \
  --split-strategy random \
  --no-use-gene-expression \
  --cv-folds 10 \
  --stratified-cv \
  --cv-seeds 42 \
  --max-samples 10000 \
  --seed 42
```

```bash
uv run python -m src.train \
  --macro-preset practical_research \
  --output-dir outputs/long_cv10_10k_pca128_practical \
  --split-strategy random \
  --cell-expression-path data/drugcomb.csv \
  --cell-feature-view 0 \
  --cv-folds 10 \
  --stratified-cv \
  --cv-seeds 42 \
  --max-samples 10000 \
  --seed 42
```

### Drug + Cell-Line Group-Aware CV

This is group-aware repeated holdout CV, not classic row-disjoint CV. A row can be test in more than one fold if its drug and cell-line memberships land in different held-out folds.

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
