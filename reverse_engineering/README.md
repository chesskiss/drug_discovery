# Reverse Engineering

Standalone module for reverse-engineering the baseline MLP and future gene/cell-line encoders.

This module is intended to hold:

- attribution CLIs
- occlusion / perturbation analysis
- component-to-gene backprojection utilities
- interpretation plots and summary exports

## Current Focus

The first target is the current baseline MLP used in `drug_synergy_baseline/`.

Near-term analyses:

- PCA-component occlusion for the `pca128` model
- prediction/error delta under masking
- global feature importance summaries
- later: backprojection from important PCA components to original selected features

## Run

From `drug_discovery/`:

```bash
python3 -m reverse_engineering.src.reverse_engineering.cli --help
```

or with the baseline environment:

```bash
cd drug_synergy_baseline
uv run python -m reverse_engineering.src.reverse_engineering.cli --help
```

## Planned Outputs

- `reverse_engineering/outputs/<run_name>/attribution_scores.csv`
- `reverse_engineering/outputs/<run_name>/summary.json`
- `reverse_engineering/outputs/<run_name>/plots/`
