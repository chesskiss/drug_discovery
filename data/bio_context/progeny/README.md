# PROGENy

**PROGENy** (Pathway RespOnsive GENes) provides gene weight signatures for **14 signaling
pathways** (Androgen, EGFR, Estrogen, Hypoxia, JAK-STAT, MAPK, NFkB, p53, PI3K, TGFb, TNFa,
Trail, VEGF, WNT), derived from a large compendium of publicly available perturbation
experiments. Each pathway is a set of target genes with an associated weight (sign/direction
of response) and p-value, and is commonly used to infer pathway activity from bulk or
single-cell expression data.

- Homepage: https://saezlab.github.io/progeny/
- Source repo: https://github.com/saezlab/progeny
- Data source used here: [OmniPath](https://omnipathdb.org/) REST API (the same source the
  official `decoupler`/`decoupleR` packages wrap) — no heavyweight client library required.

## Download

```bash
uv run --project data/bio_context data/bio_context/progeny/download_progeny.py
```

Options:

- `--top` (default 100): number of top genes per pathway to keep (by p-value).
- `--thr-padj` (default 0.05): significance threshold for keeping an interaction.
- `--license` (default `academic`): OmniPath license tier (`academic`, `commercial`, `nonprofit`).

Output: `data/progeny_<organism>_top<top>.csv` with columns `source` (pathway), `target` (gene),
`weight`, `padj`. Only `organism=human` is supported (no ortholog translation); for mouse, use
the `decoupler` package's `dc.op.progeny(organism='mouse')`.

With the defaults (top 100), the file is ~70 KB / 1,400 rows, so it's committed to the repo.
Re-run the script any time to refresh it.

## Contents (default download)

- **File:** `data/progeny_human_top100.csv` — 72 KB, 1,400 rows.
- **14 pathways** × **100 genes each** (all p-values < 0.05): Androgen, EGFR, Estrogen, Hypoxia,
  JAK-STAT, MAPK, NFkB, p53, PI3K, TGFb, TNFa, Trail, VEGF, WNT.
- **1,295 unique genes** across the 1,400 rows (100 per pathway, with 105 genes shared by more
  than one pathway).
- **Columns:** `source` (pathway), `target` (gene symbol), `weight` (signed pathway responsiveness),
  `padj` (adjusted p-value).
- Note: 2 of the 1,400 rows are exact duplicates (`PRNP` in MAPK, `PIDD1` in p53, same weight),
  so MAPK and p53 have 99 distinct genes rather than 100.

## Alignment to the TDC gene axis

`align_progeny_to_tdc.py` places these weights onto the 23,808-dim TDC NCI-60 gene axis, so a
cell line's expression vector `x` (`CellLine[0]`) compresses to a 14-dim pathway-activity vector
via `z = W @ x`.

```bash
# prerequisites (build once, shared by all datasets)
uv run --project data/bio_context tdc_gene_index/build_tdc_gene_index.py
uv run --project data/bio_context _alias/build_alias_map.py

uv run --project data/bio_context progeny/align_progeny_to_tdc.py
```

Matching strategy (symbols are unstable, Entrez ids are not):
1. direct symbol match against [`tdc_gene_index`](../tdc_gene_index) → `matched_via=symbol`
2. else symbol → Entrez via the NCBI [`_alias`](../_alias) map → Entrez match → `matched_via=alias`
3. else dropped (logged)

**Result: 1,295/1,295 genes matched (100%)** — 1,246 by symbol, **49 recovered by alias** (HGNC
renames such as `ACP3`→`ACPP`, `ATP5F1E`→`ATP5E`, `CEP43`→`FGFR1OP`). 14/14 pathways retained.

Outputs:
- `data/progeny_tdc_weights.npz` — `W` `[14 x 23808]` float32 (1,398 nonzero), `pathways`,
  plus `gene_symbol`/`entrez` per column.
- `data/progeny_tdc_alignment.csv` — audit trail: `pathway, gene_symbol, entrez, tdc_idx, weight,
  matched_via`.

The script warns loudly if two genes collide on the same `(pathway, tdc_idx)` with *conflicting*
weights (would silently lose data). For PROGENy the only collisions are the harmless exact
duplicates noted above.

Verified end-to-end: `X @ W.T` over all 59 cell lines yields a finite `[59 x 14]` activity matrix
with real cross-cell-line variation (JAK-STAT varies most, Trail least).
