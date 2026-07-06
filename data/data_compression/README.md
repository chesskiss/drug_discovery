# Data Compression

This module now owns dataset download, dataset inspection, and cell-line compression pipelines.

## Scope

- download raw DrugComb exports
- keep raw and reduced dataset artifacts out of the baseline runner
- build compressed `CellLine` representations for baseline experiments
- store dataset-visualization artifacts that describe the data itself, not model outputs

## Current Method

Current stage id:

- `Stage 1.3A`

Current pipeline id:

- `zscore_var_pca_128`

### Implemented pipeline:

1. start from raw `CellLine[0]` with `23808` dimensions
2. rank raw features by variance
3. keep top `3000`
4. z-score the selected features
5. run PCA
6. export a baseline-compatible `drugcomb.csv`

### Important constraint:

- this dataset currently has only `59` unique cell lines
- so PCA cannot yield `128` informative components
- the current pipeline therefore produces `58` effective PCA components and zero-pads the remaining dimensions to reach `128`

See the pipeline metadata in [metadata.json](/Users/arnoldcheskis/Documents/Projects/drug_discovery/data_compression/zscore_var_pca_128/data/metadata.json).

## Directory Layout

- raw dataset:
  - [drugcomb_raw.csv](/Users/arnoldcheskis/Documents/Projects/drug_discovery/data_compression/zscore_var_pca_128/data/drugcomb_raw.csv)
- compressed dataset:
  - [drugcomb.csv](/Users/arnoldcheskis/Documents/Projects/drug_discovery/data_compression/zscore_var_pca_128/data/drugcomb.csv)
- compression metadata:
  - [metadata.json](/Users/arnoldcheskis/Documents/Projects/drug_discovery/data_compression/zscore_var_pca_128/data/metadata.json)
- data-visualization artifacts:
  - [summary.json](/Users/arnoldcheskis/Documents/Projects/drug_discovery/data_compression/visualization/summary.json)
  - [target_distribution.png](/Users/arnoldcheskis/Documents/Projects/drug_discovery/data_compression/visualization/target_distribution.png)
  - [top_cells.png](/Users/arnoldcheskis/Documents/Projects/drug_discovery/data_compression/visualization/top_cells.png)
  - [top_pair_frequency_hist.png](/Users/arnoldcheskis/Documents/Projects/drug_discovery/data_compression/visualization/top_pair_frequency_hist.png)

## Current Results

Reference values discussed so far:

- raw genes baseline:
  - `Test MSE = 29.051481`
- compressed genes baseline (`zscore_var_pca_128`):
  - `Test MSE = 19.041908`
- no cell-line genes baseline:
  - `Test MSE = 20.250234`
- trivial mean baseline on the compressed-gene run:
  - `Val baseline MSE = 28.866857`
  - `Test baseline MSE = 29.050099`

Interpretation:

- compression helped substantially versus raw genes
- compressed genes currently beat the no-gene baseline by about `1.21` MSE
- this is an improvement, but not yet a large enough margin to claim the compression problem is solved

## Commands

Initialize the module environment:

```bash
cd /Users/arnoldcheskis/Documents/Projects/drug_discovery/data_compression
uv sync
```

Download DrugComb into this module:

```bash
cd /Users/arnoldcheskis/Documents/Projects/drug_discovery/data_compression
uv run python -m download_drugcomb_tdc
```

This writes by default to:

- `zscore_var_pca_128/data/drugcomb_raw.csv`
- `zscore_var_pca_128/data/reduced_drugcomb.csv`
- optionally `zscore_var_pca_128/data/tdc_cell_features.csv`

Build the compressed export from inside this directory:

```bash
cd /Users/arnoldcheskis/Documents/Projects/drug_discovery/data_compression
uv run python -m zscore_var_pca_128.build
```

Important:

- the compression step requires the full embedded `CellLine` payload source
- by default it expects:
  - `/Users/arnoldcheskis/Documents/Projects/drug_discovery/data_compression/source_data/drugcomb.pkl`
- `drugcomb_raw.csv` alone is not enough for the PCA pipeline because CSV `CellLine` payloads may be truncated

Run the data visualization utility:

```bash
cd /Users/arnoldcheskis/Documents/Projects/drug_discovery/data_compression
uv run python -m visualize \
  --input zscore_var_pca_128/data/drugcomb_raw.csv \
  --output-dir visualization \
  --preprocess
```

Verify whether a CSV contains truncated `CellLine` payloads:

```bash
cd /Users/arnoldcheskis/Documents/Projects/drug_discovery/data_compression
uv run python -m verify_cellline_payload \
  --input zscore_var_pca_128/data/drugcomb_raw.csv
```
