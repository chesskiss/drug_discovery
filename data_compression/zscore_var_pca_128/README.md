# Stage 1.3A: z-score -> VarTopK -> PCA128

This pipeline prepares a compressed `drugcomb.csv` for the baseline model.

## Pipeline

1. Load raw `CellLine[0]` vectors from `drug_synergy_baseline/data/drugcomb.pkl`.
2. Rank raw features by weighted variance across cell lines.
3. Keep the top `k` features.
4. Z-score the selected features.
5. Fit weighted PCA.
6. Export a new `drugcomb.csv` where `CellLine` contains a single compressed view.

## Important Constraint

DrugComb in this workspace has only `59` unique cell lines. That means PCA cannot provide `128` informative components.

This implementation therefore:

- computes the maximum effective PCA rank from the available cell lines
- writes those informative PCA components first
- zero-pads the remaining dimensions so the exported `CellLine[0]` still has length `128`

The exact values are written to `metadata.json`.

## Why Variance Ranking Happens Before Z-score

If you z-score every feature first, almost every non-constant feature ends up with variance `1`, so variance ranking becomes degenerate.

For that reason, this implementation ranks features by raw weighted variance first, then z-scores the selected subset before PCA.

## Run

From the repository root:

```bash
/Users/arnoldcheskis/Documents/Projects/drug_discovery/drug_synergy_baseline/.venv/bin/python \
  -m data_compression.zscore_var_pca_128.build
```

Custom example:

```bash
/Users/arnoldcheskis/Documents/Projects/drug_discovery/drug_synergy_baseline/.venv/bin/python \
  -m data_compression.zscore_var_pca_128.build \
  --input-pickle drug_synergy_baseline/data/drugcomb.pkl \
  --output-csv data_compression/zscore_var_pca_128/data/drugcomb.csv \
  --variance-top-k 3000 \
  --pca-components 128
```

## Output

- `data_compression/zscore_var_pca_128/data/drugcomb.csv`
- `data_compression/zscore_var_pca_128/data/metadata.json`
