# Drug Synergy Baseline

Minimal starting point for dataset utilities and baseline integration.

## Setup

Install `uv`, then from `drug_synergy_baseline/` run:

```bash
uv sync
```

Run scripts with:

```bash
uv run python -m <module>
```

## Run pipeline

```bash
uv run python -m src.train --synergy-path data/drugcomb_synergy.csv --output-dir outputs
```

Run single-sample prediction with:

```bash
uv run python -m src.predict --row-idx 0
```

or

```bash
uv run python -m src.predict --smiles-a "..." --smiles-b "..." --cell-line "786-0"
```

## Existing utilities

### Visualize data

```bash
uv run python -m src.visualize \
  --input data/drugcomb.csv \
  --output-dir outputs/visualization \
  --preprocess
```

### Download dataset

```bash
uv add PyTDC
uv run python -m src.download_drugcomb_tdc --save-dir data
```

## Baseline model

We use MatchMaker ([tastanlab/matchmaker](https://github.com/tastanlab/matchmaker)) as the baseline implementation for drug synergy prediction.

For now:
- it is included as a reference
- it will later be integrated into this pipeline
