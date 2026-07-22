"""Build the ordered gene index for TDC DrugComb's NCI-60 cell-line vectors.

TDC's `CellLine[0]` vectors (23,808 floats, in
`data/data_compression/source_data/drugcomb.pkl`) are anonymous -- no gene labels.
They are bit-for-bit identical to NCI CellMiner's "RNA-seq composite expression"
processed dataset, row-for-row. This script downloads that CellMiner file and emits a
table mapping each vector position (0..23807) to its gene symbol and Entrez id, so any
gene-keyed dataset (PROGENy, KEGG, ...) can be aligned onto the TDC expression axis.

IMPORTANT: the TDC paper/docs say the gene view is "5-platform microarray, 25,723
genes". That is WRONG -- it is RNA-seq composite, 23,808 genes. This script *verifies*
the true source by asserting value-identity against the pickle, and fails loudly if the
CellMiner file ever changes so the row order no longer matches.

Usage:
    uv run --project data/bio_context tdc_gene_index/build_tdc_gene_index.py
"""

from __future__ import annotations

import io
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
import requests

CELLMINER_URL = (
    "https://discover.nci.nih.gov/cellminer/download/processeddataset/"
    "nci60_RNA__RNA_seq_composite_expression.zip"
)
XLS_HEADER_ROW = 10  # gene symbols start on spreadsheet row 11 (0-indexed 10)
N_GENES = 23808

REPO_ROOT = Path(__file__).resolve().parents[3]
DRUGCOMB_PKL = REPO_ROOT / "data/data_compression/source_data/drugcomb.pkl"
OUT_DIR = Path(__file__).parent / "data"


def download_cellminer_rnaseq() -> pd.DataFrame:
    resp = requests.get(CELLMINER_URL, timeout=300)
    resp.raise_for_status()
    z = zipfile.ZipFile(io.BytesIO(resp.content))
    member = next(n for n in z.namelist() if n.endswith(".xls") and "METADATA" not in n)
    return pd.read_excel(io.BytesIO(z.read(member)), header=XLS_HEADER_ROW)


def verify_against_tdc(cellminer: pd.DataFrame) -> int:
    """Assert CellMiner rows match TDC CellLine[0] value-for-value. Returns #cells checked."""
    if not DRUGCOMB_PKL.exists():
        print(f"WARNING: {DRUGCOMB_PKL} not found -- skipping value-identity verification.")
        return 0

    cell_cols = [c for c in cellminer.columns if ":" in str(c)]
    short_to_col = {c.split(":", 1)[1].strip(): c for c in cell_cols}

    tdc = pd.read_pickle(DRUGCOMB_PKL)
    cl = tdc.drop_duplicates("Cell_Line_ID")[["Cell_Line_ID", "CellLine"]].set_index("Cell_Line_ID")

    checked = 0
    for name in cl.index:
        if name not in short_to_col:
            continue
        tdc_vec = np.asarray(cl.loc[name, "CellLine"][0], dtype=float)
        cm_vec = cellminer[short_to_col[name]].to_numpy(dtype=float)
        if not np.allclose(tdc_vec, cm_vec, atol=1e-6, equal_nan=True):
            raise AssertionError(
                f"Value mismatch for cell line {name!r}: CellMiner RNA-seq file no longer "
                f"matches TDC CellLine[0] row order. The gene index would be WRONG."
            )
        checked += 1

    if checked == 0:
        raise AssertionError("No cell lines could be matched by name -- cannot verify gene order.")
    print(f"Verified value-identity for {checked} cell lines (maxabsdiff = 0). Gene order confirmed.")
    return checked


def main() -> None:
    print("Downloading CellMiner RNA-seq composite expression (~13 MB)...")
    cellminer = download_cellminer_rnaseq()

    if len(cellminer) != N_GENES:
        raise AssertionError(
            f"Expected {N_GENES} genes, got {len(cellminer)} -- CellMiner file changed; "
            f"the TDC alignment assumption must be re-verified before trusting this index."
        )

    verify_against_tdc(cellminer)

    index = pd.DataFrame(
        {
            "idx": np.arange(len(cellminer)),
            "gene_symbol": cellminer["Gene name d"].astype(str).str.strip(),
            "entrez": cellminer["Entrez gene id e"].astype("Int64"),
        }
    )

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out_path = OUT_DIR / "tdc_gene_index.csv"
    index.to_csv(out_path, index=False)
    print(
        f"Saved {len(index)} rows to {out_path} "
        f"({index['gene_symbol'].nunique()} unique symbols, {index['entrez'].notna().sum()} with Entrez id)"
    )


if __name__ == "__main__":
    main()
