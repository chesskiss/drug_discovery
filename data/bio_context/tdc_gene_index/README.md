# TDC gene index

The gene axis for TDC DrugComb's NCI-60 cell-line vectors — i.e. **which gene sits at each
position** of the 23,808-long `CellLine[0]` expression vector.

This is the join key that makes every other dataset in `bio_context` usable: TDC's vectors are
anonymous float arrays (no labels anywhere in the pickle, the repo, or PyTDC), so without this
index a pathway weight matrix and an expression vector don't share an axis.

## Source (verified, not assumed)

`CellLine[0]` is **bit-for-bit identical** to NCI CellMiner's **RNA-seq composite expression**
dataset, row-for-row:

- File: `nci60_RNA__RNA_seq_composite_expression.zip` (~13 MB) from
  [CellMiner downloads](https://discover.nci.nih.gov/cellminer/loadDownload.do)
- Inside: `output/RNA__RNA_seq_composite_expression.xls`, header on row 11 (0-indexed 10),
  columns `Gene name d` (HGNC symbol) and `Entrez gene id e` — both populated for all 23,808 rows.
- Verified by matching cell lines by name and comparing values: **maxabsdiff = 0.0** across 54
  cell lines. Row *i* of this file == position *i* of the TDC vector.

> **The TDC paper/docs are wrong here.** They describe the gene view as "5-platform microarray,
> 25,723 genes". It is actually RNA-seq composite with 23,808 genes. The 5-platform file exists
> but has 24,190 rows and does not match. Do not trust the paper on this point.
>
> The other two views do match the docs: `CellLine[1]` = 3,171 **proteomic** features,
> `CellLine[2]` = 627 **microRNA** features (i.e. they are *not* "filtered/compact" gene views,
> as previously assumed elsewhere in the repo).

## Build

```bash
uv run --project data/bio_context tdc_gene_index/build_tdc_gene_index.py
```

Output: `data/tdc_gene_index.csv` — columns `idx` (0..23807, = vector position), `gene_symbol`,
`entrez`.

The script **verifies value-identity against `data/data_compression/source_data/drugcomb.pkl`**
on every run and fails loudly if CellMiner ever changes the file such that the row order stops
matching — the index is only trustworthy while that assertion holds.

## Usage

Join a gene-keyed dataset onto `idx` to place its weights on the TDC axis. Prefer joining on
`entrez` (stable) over `gene_symbol` (renamed over time by HGNC) — see [`_alias/`](../_alias)
for resolving arbitrary/legacy symbols to Entrez ids.
