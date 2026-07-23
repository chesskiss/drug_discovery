"""Align KEGG pathway membership onto the TDC NCI-60 gene axis.

Produces a weight matrix `W` of shape [n_pathways x 23808] such that, for a cell line's
TDC expression vector `x` (CellLine[0], length 23808), `z = W @ x` is a per-pathway
activity vector -- the interpretable low-dim "bio-program" representation we want. Each
cell holds KEGG's per-(gene, pathway) `importance` score (a structural/topological
signal: `max(1 - 1/(1+degree), betweenness)`, forced to 1.0 at articulation points --
see download_kegg.py / kegg/README.md), NOT a statistically-derived effect size like
PROGENy's `weight`.

Gene matching is the mirror image of PROGENy's (see tdc_gene_index/ and _alias/): KEGG's
`gene_id` is `hsa:<entrez>`, so Entrez is embedded and the primary match is a direct
number-to-number join -- no symbols, no alias map. Only the residual misses fall back to
symbol/alias resolution.
  1. direct: KEGG `gene_id` entrez == an Entrez in the TDC gene index -> matched_via=direct
  2. symbol: else KEGG gene_symbol == a symbol in the TDC gene index    -> matched_via=symbol
  3. alias : else KEGG symbol -> Entrez (NCBI alias map) -> Entrez in TDC -> matched_via=alias
  4. else  : unmatched (dropped, logged)

Prereqs (run first):
    uv run --project data/bio_context tdc_gene_index/build_tdc_gene_index.py
    uv run --project data/bio_context _alias/build_alias_map.py
    uv run --project data/bio_context kegg/download_kegg.py
Usage:
    uv run --project data/bio_context kegg/align_kegg_to_tdc.py
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

BIO_CONTEXT = Path(__file__).resolve().parents[1]
KEGG_CSV = Path(__file__).parent / "data" / "kegg_human_gene_pathways.csv"
TDC_INDEX_CSV = BIO_CONTEXT / "tdc_gene_index" / "data" / "tdc_gene_index.csv"
ALIAS_CSV = BIO_CONTEXT / "_alias" / "data" / "alias_to_entrez.csv"
OUT_DIR = Path(__file__).parent / "data"
N_GENES = 23808
VALUE_COL = "importance"  # which KEGG per-(gene, pathway) metric populates W


def resolve_symbol_to_tdc_idx(symbols: pd.Series, tdc_index: pd.DataFrame, alias: pd.DataFrame) -> pd.DataFrame:
    """Map each symbol to a TDC row index (symbol-direct, then NCBI-alias fallback).

    Returns cols: gene_symbol, tdc_idx, entrez, matched_via (symbol | alias | unmatched).
    Used here only for the residual genes that miss the primary direct-Entrez join.
    """
    sym_to_idx = tdc_index.dropna(subset=["gene_symbol"]).drop_duplicates("gene_symbol").set_index("gene_symbol")["idx"]
    # build entrez -> idx from the FILTERED frame's own entrez column (not the full-length
    # tdc_index["entrez"], which would misalign / raise if entrez ever had a null or dup).
    tdc_ent = tdc_index.dropna(subset=["entrez"]).drop_duplicates("entrez").copy()
    tdc_ent["entrez"] = tdc_ent["entrez"].astype("Int64")
    entrez_to_idx = tdc_ent.set_index("entrez")["idx"]
    # symbol -> entrez, keep only symbols that map unambiguously to a single gene
    alias_unique = alias.drop_duplicates("alias_symbol", keep=False)
    sym_to_entrez = alias_unique.set_index("alias_symbol")["entrez"].astype("Int64")

    idx_to_entrez = tdc_index.set_index("idx")["entrez"]
    out = []
    for s in symbols:
        if s in sym_to_idx.index:
            idx = int(sym_to_idx[s])
            out.append((s, idx, idx_to_entrez.get(idx, pd.NA), "symbol"))
            continue
        ent = sym_to_entrez.get(s, pd.NA)
        if pd.notna(ent) and ent in entrez_to_idx.index:
            out.append((s, int(entrez_to_idx[ent]), ent, "alias"))
        else:
            out.append((s, -1, pd.NA, "unmatched"))
    return pd.DataFrame(out, columns=["gene_symbol", "tdc_idx", "entrez", "matched_via"])


def resolve_kegg_genes(genes: pd.DataFrame, tdc_index: pd.DataFrame, alias: pd.DataFrame) -> pd.DataFrame:
    """Resolve each unique KEGG gene to a TDC row index, entrez-first.

    `genes`: one row per gene, cols gene_id ("hsa:<entrez>") and gene_symbol (KEGG's raw
    label). Returns one row per gene_id: gene_id, tdc_idx (Int64, NA if unmatched),
    matched_via (direct | symbol | alias | unmatched).
    """
    tdc_ent = tdc_index.dropna(subset=["entrez"]).drop_duplicates("entrez").copy()
    tdc_ent["entrez"] = tdc_ent["entrez"].astype("Int64")
    entrez_to_idx = tdc_ent.set_index("entrez")["idx"]

    genes = genes.copy()
    kegg_entrez = genes["gene_id"].str.removeprefix("hsa:").astype("Int64")
    genes["tdc_idx"] = kegg_entrez.map(entrez_to_idx).astype("Int64")
    genes["matched_via"] = np.where(genes["tdc_idx"].notna(), "direct", "unmatched")

    miss = genes["matched_via"] == "unmatched"
    if miss.any():
        fb = resolve_symbol_to_tdc_idx(genes.loc[miss, "gene_symbol"], tdc_index, alias)
        rescued = (fb["matched_via"] != "unmatched").to_numpy()
        miss_pos = genes.index[miss]
        genes.loc[miss_pos[rescued], "tdc_idx"] = fb.loc[rescued, "tdc_idx"].astype("Int64").to_numpy()
        genes.loc[miss_pos, "matched_via"] = fb["matched_via"].to_numpy()

    return genes[["gene_id", "tdc_idx", "matched_via"]]


def build_weight_matrix(matched_dedup: pd.DataFrame, pathway_ids: list[str], value_col: str) -> np.ndarray:
    """Scatter one value per (pathway, gene) cell into W [n_pathways x 23808], vectorized.

    `matched_dedup` must already be unique on (pathway_id, tdc_idx) so the fancy-index
    write is collision-free (deterministic, not relying on NumPy duplicate-write order).
    """
    pw_to_row = {p: i for i, p in enumerate(pathway_ids)}
    W = np.zeros((len(pathway_ids), N_GENES), dtype=np.float32)
    rows = matched_dedup["pathway_id"].map(pw_to_row).to_numpy()
    cols = matched_dedup["tdc_idx"].astype(int).to_numpy()
    W[rows, cols] = matched_dedup[value_col].to_numpy(dtype=np.float32)
    return W


def main() -> None:
    for p in (KEGG_CSV, TDC_INDEX_CSV, ALIAS_CSV):
        if not p.exists():
            raise FileNotFoundError(f"Missing prerequisite: {p}. See module docstring for build steps.")

    kegg = pd.read_csv(KEGG_CSV)
    tdc_index = pd.read_csv(TDC_INDEX_CSV)
    tdc_index["entrez"] = tdc_index["entrez"].astype("Int64")
    alias = pd.read_csv(ALIAS_CSV)

    genes = kegg.drop_duplicates("gene_id")[["gene_id", "gene_symbol"]]
    resolved = resolve_kegg_genes(genes, tdc_index, alias)
    if not resolved["gene_id"].is_unique:
        raise AssertionError("resolve_kegg_genes returned duplicate gene_id rows -- merge would multiply rows.")

    aligned = kegg.merge(resolved, on="gene_id", how="left")
    # clean, unambiguous TDC labels looked up via tdc_idx (KEGG's own gene_symbol is
    # corrupted for ~625 genes; see kegg/README.md -- keep it only for traceability).
    idx_to_sym = tdc_index.set_index("idx")["gene_symbol"]
    idx_to_ent = tdc_index.set_index("idx")["entrez"]
    aligned["tdc_gene_symbol"] = aligned["tdc_idx"].map(idx_to_sym)
    aligned["tdc_entrez"] = aligned["tdc_idx"].map(idx_to_ent)

    matched = aligned[aligned["matched_via"] != "unmatched"].copy()
    pathways = sorted(kegg["pathway_id"].unique())
    name_map = kegg.drop_duplicates("pathway_id").set_index("pathway_id")["pathway_name"]

    # Two genes can land on the same (pathway, tdc_idx) -- e.g. a paralog rescued by symbol
    # onto a gene already claimed directly. Identical importance collapses harmlessly, but
    # conflicting values would be silently lost, so surface those loudly.
    collisions = matched[matched.duplicated(["pathway_id", "tdc_idx"], keep=False)]
    conflicting = [
        (pw, idx, sorted(g[VALUE_COL].unique()))
        for (pw, idx), g in collisions.groupby(["pathway_id", "tdc_idx"])
        if g[VALUE_COL].nunique() > 1
    ]
    if conflicting:
        print(f"WARNING: {len(conflicting)} (pathway, gene) collisions with CONFLICTING {VALUE_COL} -- last wins:")
        for pw, idx, vals in conflicting[:10]:
            print(f"    {pw} @ tdc_idx {idx}: {vals}")
    elif len(collisions):
        print(f"Note: {len(collisions)} rows collide on (pathway, gene) but agree on {VALUE_COL} (harmless duplicates).")

    matched_dedup = matched.drop_duplicates(["pathway_id", "tdc_idx"], keep="last")
    W = build_weight_matrix(matched_dedup, pathways, VALUE_COL)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        OUT_DIR / "kegg_tdc_weights.npz",
        W=W,
        pathways=np.array(pathways),
        pathway_names=np.array([name_map[p] for p in pathways]),
        gene_symbol=tdc_index["gene_symbol"].to_numpy(),
        entrez=tdc_index["entrez"].astype(float).to_numpy(),
    )
    audit = aligned.rename(
        columns={"gene_symbol": "kegg_gene_symbol", "tdc_gene_symbol": "gene_symbol", "tdc_entrez": "entrez"}
    )[[
        "pathway_id", "pathway_name", "gene_id", "kegg_gene_symbol", "gene_symbol", "entrez",
        "tdc_idx", "matched_via", "category", "is_enzyme", "degree", "betweenness",
        "is_articulation", "importance",
    ]]
    audit.to_csv(OUT_DIR / "kegg_tdc_alignment.csv", index=False)

    # coverage report (per unique gene, then per source row)
    via = resolved["matched_via"].value_counts()
    n_genes = len(resolved)
    n_matched = n_genes - via.get("unmatched", 0)
    unmatched_sym = genes.merge(resolved, on="gene_id")
    unmatched_sym = sorted(unmatched_sym.loc[unmatched_sym.matched_via == "unmatched", "gene_symbol"].astype(str))
    print(f"KEGG genes: {n_genes} unique")
    print(f"  matched by direct entrez: {via.get('direct', 0)}")
    print(f"  matched by symbol       : {via.get('symbol', 0)}")
    print(f"  matched by alias        : {via.get('alias', 0)}")
    print(f"  unmatched                : {via.get('unmatched', 0)}  {unmatched_sym[:10]}")
    print(f"  -> matched {n_matched}/{n_genes} = {100*n_matched/n_genes:.1f}%")
    print(f"     (unmatched = KEGG Entrez id absent from CellMiner's panel: non-coding RNAs,")
    print(f"      readthrough/fusion transcripts, mito genes, newer/uncharacterized ids)")
    print(f"Source rows: {len(matched)}/{len(kegg)} = {100*len(matched)/len(kegg):.1f}% matched at (gene, pathway) grain")
    populated = len(matched_dedup)
    print(
        f"W shape {W.shape}; pathways retained {matched['pathway_id'].nunique()}/{len(pathways)}; "
        f"populated cells {populated}, nonzero {int((W != 0).sum())} "
        f"(gap = genes with {VALUE_COL}=0, i.e. degree-0 / edgeless-pathway members)"
    )


if __name__ == "__main__":
    main()
