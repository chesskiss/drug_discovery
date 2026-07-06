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
