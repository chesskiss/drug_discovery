# bio_context

Small biological prior-knowledge resources used as side-context (pathway signatures,
TF regulons, gene sets, etc.), as opposed to the drug/cell-line data under
`data/data_compression`.

**Goal:** collectively, these datasets should cover as much as possible of the ~24k genes that
appear in a single cell-line expression vector in the DrugComb NCI-60 dataset (TDC). Each
dataset's README notes how many unique genes it contributes.

## Layout

Each dataset lives in its own subfolder:

```
bio_context/
  pyproject.toml        # shared deps (requests, pandas) for all download scripts
  <dataset>/
    download_<dataset>.py
    README.md            # what it is, source, how to (re)download
    data/                 # downloaded output (committed if small; noted in the dataset README if not)
```

## Conventions

- Download scripts are plain, dependency-light Python (`requests` + `pandas`), runnable via
  `uv run --project data/bio_context <dataset>/download_<dataset>.py`.
- Scripts hit official REST endpoints directly rather than pulling in heavyweight client
  libraries, unless the dataset genuinely needs one.
- Output is written to `<dataset>/data/` and is safe to delete and re-download at any time.

## Datasets

- [`progeny/`](progeny/README.md) — PROGENy pathway-responsive gene weights: 14 signaling
  pathways, 1,295 unique genes, 72 KB. Has statistically-derived per-gene `weight`.
- [`kegg/`](kegg/README.md) — KEGG gene-to-pathway membership: 372 human pathways (metabolic,
  signaling, disease, etc.), 9,416 unique genes, ~5 MB. No statistical weights — instead
  per-(gene,pathway) structural metrics: `category`, `is_enzyme`, `degree`, `betweenness`,
  `is_articulation`, and a combined `importance` (see dataset README for caveats).

