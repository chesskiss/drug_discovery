# Cell Line Difficulty

Standalone analysis module for ranking DrugComb cell lines by observed synergy outcomes while ignoring molecule identity.

## Run

From `drug_discovery/`:

```bash
python3 -m cell_line_difficulty.src.cell_line_difficulty.cli
```

or with the baseline environment:

```bash
cd drug_synergy_baseline
uv run python -m cell_line_difficulty.src.cell_line_difficulty.cli --input data/drugcomb_synergy.csv --output-dir ../cell_line_difficulty/outputs
```

## Outputs

- `outputs/ranked_by_ease.csv`
- `outputs/ranked_by_difficulty.csv`
- `outputs/summary.json`

## Predictive Model

Build a cell-line-only predictive dataset and run leave-one-cell-line-out evaluation with ridge and a tiny MLP:

```bash
python3 -m cell_line_difficulty.src.cell_line_difficulty.predict_cli
```

Compare the three `CellLine` views on the same task:

```bash
python3 -m cell_line_difficulty.src.cell_line_difficulty.predict_cli --compare-views --models ridge
```

Predictive outputs:

- `outputs/predictive_dataset.csv`
- `outputs/loocv_predictions.csv`
- `outputs/predictive_metrics.json`
- `outputs/view_comparison.csv`
