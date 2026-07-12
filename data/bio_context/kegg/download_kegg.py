"""Download KEGG gene-to-pathway membership + per-pathway graph importance metrics (human).

KEGG PATHWAY (https://www.kegg.jp/pathway/map01100) publishes its pathway maps as
diagrams, but the underlying data is available as plain tables through the KEGG REST
API (https://rest.kegg.jp), with no need to parse the graph/image. This script pulls:

  1. link/pathway/hsa   -> every human gene linked to every KEGG pathway it's in
  2. list/pathway/hsa    -> pathway id -> pathway name
  3. list/hsa            -> gene id -> primary gene symbol
  4. link/enzyme/hsa     -> which human genes have an EC/enzyme link (is_enzyme flag)
  5. get/br:br08901/json -> BRITE top-level category per pathway (Metabolism, etc.)

and joins them into one long-format table (one row per gene-pathway pair).

It then downloads each pathway's own KGML graph (372 individual files, NOT the combined
map01100 overview) and, per pathway, builds a networkx graph to compute per-gene
structural metrics:

  - degree          : number of reaction/relation edges touching that gene's node
  - betweenness     : networkx betweenness_centrality (normalized 0-1) of the gene's node
  - is_articulation : True if removing the gene's node disconnects the pathway graph

and a combined `importance` score. These are topological signals (how wired-in / how much
of a bottleneck a gene is within a pathway), NOT statistically-derived effect sizes like
PROGENy's `weight` -- KEGG has no equivalent of that. Caveats: a KGML node can bundle
several isozyme genes together (they all get the same metrics), and the metrics don't
distinguish edge direction or sign (activation vs inhibition).

Graph construction note: metabolic pathways connect enzymes only THROUGH shared compound
nodes (enzyme -> compound -> enzyme), so the graph includes compound/ortholog/group nodes
(not just gene nodes) for correct connectivity; metrics are only reported for gene nodes.

Usage:
    uv run --project data/bio_context kegg/download_kegg.py
"""

from __future__ import annotations

import re
import time
import xml.etree.ElementTree as ET
from pathlib import Path

import networkx as nx
import pandas as pd
import requests

BASE_URL = "https://rest.kegg.jp"
OUT_DIR = Path(__file__).parent / "data"
REQUEST_DELAY_S = 0.1  # be polite to the free public API across ~372 requests

# Entry types kept as graph nodes (for correct connectivity). "map" links to other
# pathways and would create artificial shortcuts/hubs, so it is excluded.
GRAPH_NODE_TYPES = {"gene", "compound", "ortholog", "group"}


def _get_tsv(path: str, columns: list[str]) -> pd.DataFrame:
    resp = requests.get(f"{BASE_URL}/{path}", timeout=120)
    resp.raise_for_status()
    return pd.read_csv(pd.io.common.BytesIO(resp.content), sep="\t", names=columns)


def fetch_membership_table() -> pd.DataFrame:
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


def fetch_enzyme_gene_ids() -> set[str]:
    """Human genes that have an EC/enzyme link (for the is_enzyme flag)."""
    enzyme = _get_tsv("link/enzyme/hsa", ["gene_id", "ec_id"])
    return set(enzyme["gene_id"])


def fetch_pathway_categories() -> dict[str, str]:
    """map-number (5 digits, e.g. '00010') -> BRITE top-level category name."""
    resp = requests.get(f"{BASE_URL}/get/br:br08901/json", timeout=120)
    resp.raise_for_status()
    data = resp.json()

    lookup: dict[str, str] = {}

    def walk(node: dict, depth: int = 0, top: str | None = None) -> None:
        name = node.get("name", "")
        if depth == 1:
            top = name
        children = node.get("children")
        if children is None:
            m = re.match(r"(\d{5})\s", name)
            if m and top is not None:
                lookup[m.group(1)] = top
        else:
            for child in children:
                walk(child, depth + 1, top)

    walk(data)
    return lookup


def fetch_pathway_gene_metrics(pathway_id: str) -> dict[str, dict]:
    """Per-gene structural metrics within one pathway's own KGML graph.

    `pathway_id` like "hsa00010" (no "path:" prefix). Returns
    {gene_id ("hsa:123"): {"degree", "betweenness", "is_articulation"}}.
    Metrics are only reported for `gene`-type entries (real hsa: ids), but the graph
    itself includes compound/ortholog/group nodes for correct connectivity.
    """
    resp = requests.get(f"{BASE_URL}/get/{pathway_id}/kgml", timeout=60)
    if resp.status_code != 200 or not resp.content:
        return {}
    root = ET.fromstring(resp.content)

    # entry id -> (type, [gene ids])   for every graph-node entry
    entries: dict[str, tuple[str, list[str]]] = {}
    gene_entry_ids: set[str] = set()
    for e in root.findall("entry"):
        etype = e.get("type")
        if etype not in GRAPH_NODE_TYPES:
            continue
        names = e.get("name", "").split()
        entries[e.get("id")] = (etype, names)
        if etype == "gene":
            gene_entry_ids.add(e.get("id"))

    graph = nx.Graph()
    graph.add_nodes_from(entries)

    # reaction edges: the reaction id equals the catalyzing entry id; connect it to each
    # substrate/product compound entry (enzyme -> compound), so enzymes sharing a compound
    # become connected through it.
    for reaction in root.findall("reaction"):
        rid = reaction.get("id")
        if rid not in entries:
            continue
        for side in ("substrate", "product"):
            for node in reaction.findall(side):
                cid = node.get("id")
                if cid in entries:
                    graph.add_edge(rid, cid)

    # relation edges: entry1 -- entry2 (mostly gene-to-gene in signaling pathways)
    for relation in root.findall("relation"):
        a, b = relation.get("entry1"), relation.get("entry2")
        if a in entries and b in entries:
            graph.add_edge(a, b)

    betweenness = nx.betweenness_centrality(graph) if graph.number_of_edges() else {}
    articulation = set()
    for component in nx.connected_components(graph):
        if len(component) > 2:
            sub = graph.subgraph(component)
            articulation.update(nx.articulation_points(sub))

    # aggregate per entry, then fan out to each gene id inside that entry (isozyme bundle)
    gene_metrics: dict[str, dict] = {}
    for eid in gene_entry_ids:
        _etype, gene_ids = entries[eid]
        deg = graph.degree(eid)
        btw = betweenness.get(eid, 0.0)
        art = eid in articulation
        for gid in gene_ids:
            prev = gene_metrics.get(gid)
            # a gene may appear in multiple boxes within one map; keep the strongest node
            if prev is None or deg > prev["degree"]:
                gene_metrics[gid] = {"degree": deg, "betweenness": btw, "is_articulation": art}
    return gene_metrics


def add_graph_metrics(df: pd.DataFrame) -> pd.DataFrame:
    pathway_ids = sorted(df["pathway_id"].unique())
    rows = []
    for i, pid in enumerate(pathway_ids, 1):
        bare_id = pid.removeprefix("path:")
        for gene_id, m in fetch_pathway_gene_metrics(bare_id).items():
            rows.append({"gene_id": gene_id, "pathway_id": pid, **m})
        if i % 50 == 0 or i == len(pathway_ids):
            print(f"  computed graph metrics for {i}/{len(pathway_ids)} pathways")
        time.sleep(REQUEST_DELAY_S)

    metrics = pd.DataFrame(rows)
    merged = df.merge(metrics, on=["gene_id", "pathway_id"], how="left")
    merged["degree"] = merged["degree"].fillna(0).astype(int)
    merged["betweenness"] = merged["betweenness"].fillna(0.0)
    merged["is_articulation"] = merged["is_articulation"].fillna(False).astype(bool)
    return merged


def add_importance(df: pd.DataFrame) -> pd.DataFrame:
    # local saturating importance from degree, blended with global bottleneck signal;
    # a provable articulation point is a hard bottleneck and overrides both.
    local = 1.0 - 1.0 / (1.0 + df["degree"])
    df["importance"] = df[["betweenness"]].assign(local=local).max(axis=1)
    df.loc[df["is_articulation"], "importance"] = 1.0
    return df


def main() -> None:
    df = fetch_membership_table()
    print(f"Fetched {len(df)} gene-pathway rows.")

    enzyme_ids = fetch_enzyme_gene_ids()
    df["is_enzyme"] = df["gene_id"].isin(enzyme_ids)

    categories = fetch_pathway_categories()
    df["category"] = df["pathway_id"].str.extract(r"hsa(\d{5})")[0].map(categories)

    print("Computing per-pathway graph metrics (degree, betweenness, articulation)...")
    df = add_graph_metrics(df)
    df = add_importance(df)

    df = df[[
        "gene_id", "gene_symbol", "pathway_id", "pathway_name", "category",
        "is_enzyme", "degree", "betweenness", "is_articulation", "importance",
    ]]

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out_path = OUT_DIR / "kegg_human_gene_pathways.csv"
    df.to_csv(out_path, index=False)
    print(
        f"Saved {len(df)} rows "
        f"({df['gene_id'].nunique()} unique genes, {df['pathway_id'].nunique()} pathways, "
        f"{df['is_enzyme'].sum()} enzyme rows, {df['is_articulation'].sum()} articulation rows) "
        f"to {out_path}"
    )


if __name__ == "__main__":
    main()
