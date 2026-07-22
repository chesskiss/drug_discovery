"""Build a gene-symbol-alias -> Entrez id map from NCBI gene_info.

Gene symbols get renamed over time by HGNC (e.g. ACP3 <- ACPP, ATP5F1E <- ATP5E), but
Entrez ids are stable. A dataset keyed by current symbols (PROGENy, DoRothEA, ...) will
miss genes whose symbol was renamed relative to whatever the target uses. This map lets
us resolve ANY symbol -- current or historical synonym -- to its stable Entrez id, so
alignment can join on Entrez instead of fragile symbols.

Source: NCBI `Homo_sapiens.gene_info.gz`, columns `GeneID`, `Symbol`, `Synonyms`
(pipe-separated). We explode both the primary symbol and every synonym into rows.

Usage:
    uv run --project data/bio_context _alias/build_alias_map.py
"""

from __future__ import annotations

import gzip
import io
from pathlib import Path

import pandas as pd
import requests

GENE_INFO_URL = "https://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz"
OUT_DIR = Path(__file__).parent / "data"


def build_alias_map() -> pd.DataFrame:
    resp = requests.get(GENE_INFO_URL, timeout=300)
    resp.raise_for_status()
    with gzip.open(io.BytesIO(resp.content), "rt") as fh:
        gene_info = pd.read_csv(fh, sep="\t", usecols=["GeneID", "Symbol", "Synonyms"], dtype=str)

    rows = []
    for _, r in gene_info.iterrows():
        entrez = r["GeneID"]
        aliases = {r["Symbol"]}
        if r["Synonyms"] and r["Synonyms"] != "-":
            aliases.update(r["Synonyms"].split("|"))
        for a in aliases:
            a = a.strip()
            if a and a != "-":
                rows.append((a, entrez))

    alias = pd.DataFrame(rows, columns=["alias_symbol", "entrez"]).drop_duplicates()
    return alias


def main() -> None:
    print("Downloading NCBI Homo_sapiens.gene_info.gz (~20 MB)...")
    alias = build_alias_map()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out_path = OUT_DIR / "alias_to_entrez.csv"
    alias.to_csv(out_path, index=False)
    ambiguous = alias["alias_symbol"].duplicated(keep=False).sum()
    print(
        f"Saved {len(alias)} alias->entrez rows to {out_path} "
        f"({alias['alias_symbol'].nunique()} unique symbols, {alias['entrez'].nunique()} genes; "
        f"{ambiguous} rows have a symbol mapping to >1 gene)"
    )


if __name__ == "__main__":
    main()
