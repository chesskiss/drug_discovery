"""Download KEGG gene-to-pathway membership (+ per-pathway graph degree) for human.

KEGG PATHWAY (https://www.kegg.jp/pathway/map01100) publishes its pathway maps as
diagrams, but the underlying gene/pathway associations are available as plain
tables through the KEGG REST API (https://rest.kegg.jp), with no need to parse
the graph/image. This script pulls:

  1. link/pathway/hsa  -> every human gene linked to every KEGG pathway it's in
  2. list/pathway/hsa   -> pathway id -> pathway name
  3. list/hsa           -> gene id -> primary gene symbol

and joins them into one long-format table.

It then downloads each pathway's own KGML graph (the machine-readable form of the
diagram, one file per pathway -- NOT the combined map01100 overview) and adds a
`degree` column: how many reaction/relation edges that gene's node touches within
that specific pathway's graph. This is a topological signal (is this gene directly
wired into a reaction/relation step in this pathway, and how many), not a
statistically-derived effect size like PROGENy's `weight` -- KEGG has no equivalent
of that. Caveats: a KGML node can bundle several isozyme genes together (they all
get the same degree), and degree doesn't distinguish edge direction or sign
(activation vs inhibition).

Usage:
    uv run --project data/bio_context kegg/download_kegg.py
"""

from __future__ import annotations

import time
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd
import requests

BASE_URL = "https://rest.kegg.jp"
OUT_DIR = Path(__file__).parent / "data"
REQUEST_DELAY_S = 0.1  # be polite to the free public API across ~372 requests


def _get_tsv(path: str, columns: list[str]) -> pd.DataFrame:
    resp = requests.get(f"{BASE_URL}/{path}", timeout=120)
    resp.raise_for_status()
    return pd.read_csv(pd.io.common.BytesIO(resp.content), sep="\t", names=columns)


def fetch_kegg_human_pathways() -> pd.DataFrame:
    gene_pathway = _get_tsv("link/pathway/hsa", ["gene_id", "pathway_id"])
    pathway_names = _get_tsv("list/pathway/hsa", ["pathway_id", "pathway_name"])
    # link/pathway/hsa ids are prefixed ("path:hsa00010"), list/pathway/hsa ids are not
    # ("hsa00010") -- normalize so the merge below actually matches.
    pathway_names["pathway_id"] = "path:" + pathway_names["pathway_id"]
    genes = _get_tsv("list/hsa", ["gene_id", "gene_type", "position", "gene_desc"])
    genes["gene_symbol"] = genes["gene_desc"].str.split(",", n=1).str[0].str.strip()

    df = gene_pathway.merge(pathway_names, on="pathway_id", how="left")
    df = df.merge(genes[["gene_id", "gene_symbol"]], on="gene_id", how="left")
    return df[["gene_id", "gene_symbol", "pathway_id", "pathway_name"]]


def fetch_pathway_gene_degrees(pathway_id: str) -> dict[str, int]:
    """Per-gene edge count within one pathway's own KGML graph.

    `pathway_id` like "hsa00010" (no "path:" prefix). Returns {gene_id ("hsa:123"): degree}.
    Only counts `gene`-type entries (real hsa: ids), not `ortholog`/`compound`/`map` entries.
    """
    resp = requests.get(f"{BASE_URL}/get/{pathway_id}/kgml", timeout=60)
    if resp.status_code != 200 or not resp.content:
        return {}
    root = ET.fromstring(resp.content)

    gene_entries = {e.get("id"): e.get("name").split() for e in root.findall("entry") if e.get("type") == "gene"}
    edge_count: dict[str, int] = dict.fromkeys(gene_entries, 0)

    for reaction in root.findall("reaction"):
        eid = reaction.get("id")
        if eid in edge_count:
            edge_count[eid] += 1

    for relation in root.findall("relation"):
        for eid in (relation.get("entry1"), relation.get("entry2")):
            if eid in edge_count:
                edge_count[eid] += 1

    gene_degree: dict[str, int] = {}
    for eid, gene_ids in gene_entries.items():
        for gid in gene_ids:
            gene_degree[gid] = gene_degree.get(gid, 0) + edge_count[eid]
    return gene_degree


def add_pathway_degrees(df: pd.DataFrame) -> pd.DataFrame:
    pathway_ids = sorted(df["pathway_id"].unique())
    degree_rows = []
    for i, pid in enumerate(pathway_ids, 1):
        bare_id = pid.removeprefix("path:")
        gene_degree = fetch_pathway_gene_degrees(bare_id)
        for gene_id, degree in gene_degree.items():
            degree_rows.append({"gene_id": gene_id, "pathway_id": pid, "degree": degree})
        if i % 50 == 0 or i == len(pathway_ids):
            print(f"  fetched degrees for {i}/{len(pathway_ids)} pathways")
        time.sleep(REQUEST_DELAY_S)

    degrees = pd.DataFrame(degree_rows)
    merged = df.merge(degrees, on=["gene_id", "pathway_id"], how="left")
    merged["degree"] = merged["degree"].fillna(0).astype(int)
    return merged


def main() -> None:
    df = fetch_kegg_human_pathways()
    print(f"Fetched {len(df)} gene-pathway rows, computing per-pathway graph degree...")
    df = add_pathway_degrees(df)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out_path = OUT_DIR / "kegg_human_gene_pathways.csv"
    df.to_csv(out_path, index=False)
    print(
        f"Saved {len(df)} gene-pathway rows "
        f"({df['gene_id'].nunique()} unique genes, {df['pathway_id'].nunique()} unique pathways) "
        f"to {out_path}"
    )
    print(f"Rows with degree=0 (no reaction/relation edge found, e.g. bundled ortholog-only or lookup miss): "
          f"{(df['degree'] == 0).sum()} / {len(df)}")


if __name__ == "__main__":
    main()
