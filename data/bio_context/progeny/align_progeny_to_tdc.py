"""Align PROGENy pathway weights onto the TDC NCI-60 gene axis.

Produces a weight matrix `W` of shape [n_pathways x 23808] such that, for a cell line's
TDC expression vector `x` (CellLine[0], length 23808), `z = W @ x` is a per-pathway
activity vector -- the interpretable low-dim "bio-program" representation we want.

Gene matching (see tdc_gene_index/ and _alias/):
  1. direct: PROGENy `target` symbol == a symbol in the TDC gene index  -> matched_via=symbol
  2. alias : else PROGENy symbol -> Entrez (NCBI alias map) -> Entrez in TDC index -> matched_via=alias
  3. else  : unmatched (dropped, logged)

Prereqs (run first):
    uv run --project data/bio_context tdc_gene_index/build_tdc_gene_index.py
    uv run --project data/bio_context _alias/build_alias_map.py
Usage:
    uv run --project data/bio_context progeny/align_progeny_to_tdc.py
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

BIO_CONTEXT = Path(__file__).resolve().parents[1]
PROGENY_CSV = Path(__file__).parent / "data" / "progeny_human_top100.csv"
TDC_INDEX_CSV = BIO_CONTEXT / "tdc_gene_index" / "data" / "tdc_gene_index.csv"
ALIAS_CSV = BIO_CONTEXT / "_alias" / "data" / "alias_to_entrez.csv"
OUT_DIR = Path(__file__).parent / "data"
N_GENES = 23808


def resolve_to_tdc_idx(symbols: pd.Series, tdc_index: pd.DataFrame, alias: pd.DataFrame) -> pd.DataFrame:
    """Map each symbol to a TDC row index. Returns cols: gene_symbol, tdc_idx, entrez, matched_via."""
    sym_to_idx = tdc_index.dropna(subset=["gene_symbol"]).drop_duplicates("gene_symbol").set_index("gene_symbol")["idx"]
    entrez_to_idx = (
        tdc_index.dropna(subset=["entrez"]).drop_duplicates("entrez").set_index(tdc_index["entrez"].astype("Int64"))["idx"]
    )
    # symbol -> entrez, keep only symbols that map unambiguously to a single gene
    alias_unique = alias.drop_duplicates("alias_symbol", keep=False)
    sym_to_entrez = alias_unique.set_index("alias_symbol")["entrez"].astype("Int64")

    out = []
    for s in symbols:
        if s in sym_to_idx.index:
            out.append((s, int(sym_to_idx[s]), tdc_index.loc[sym_to_idx[s] == tdc_index["idx"], "entrez"].iloc[0], "symbol"))
            continue
        ent = sym_to_entrez.get(s, pd.NA)
        if pd.notna(ent) and ent in entrez_to_idx.index:
            idx = int(entrez_to_idx[ent])
            out.append((s, idx, ent, "alias"))
        else:
            out.append((s, -1, pd.NA, "unmatched"))
    return pd.DataFrame(out, columns=["gene_symbol", "tdc_idx", "entrez", "matched_via"])


def main() -> None:
    for p in (PROGENY_CSV, TDC_INDEX_CSV, ALIAS_CSV):
        if not p.exists():
            raise FileNotFoundError(f"Missing prerequisite: {p}. See module docstring for build steps.")

    progeny = pd.read_csv(PROGENY_CSV)  # source, target, weight, padj
    tdc_index = pd.read_csv(TDC_INDEX_CSV)
    tdc_index["entrez"] = tdc_index["entrez"].astype("Int64")
    alias = pd.read_csv(ALIAS_CSV)

    resolved = resolve_to_tdc_idx(pd.Series(progeny["target"].unique()), tdc_index, alias)
    aligned = progeny.merge(resolved, left_on="target", right_on="gene_symbol", how="left")

    matched = aligned[aligned["matched_via"] != "unmatched"].copy()
    pathways = sorted(progeny["source"].unique())
    pw_to_row = {p: i for i, p in enumerate(pathways)}

    # Two source rows can land on the same (pathway, tdc_idx) -- either duplicate rows in the
    # source data, or two symbols resolving to the same gene. Identical weights collapse
    # harmlessly, but conflicting weights would be silently lost, so surface those loudly.
    collisions = matched[matched.duplicated(["source", "tdc_idx"], keep=False)]
    conflicting = [
        (pw, idx, sorted(g["weight"].unique()))
        for (pw, idx), g in collisions.groupby(["source", "tdc_idx"])
        if g["weight"].nunique() > 1
    ]
    if conflicting:
        print(f"WARNING: {len(conflicting)} (pathway, gene) collisions with CONFLICTING weights -- last wins:")
        for pw, idx, weights in conflicting[:10]:
            print(f"    {pw} @ tdc_idx {idx}: {weights}")
    elif len(collisions):
        print(f"Note: {len(collisions)} rows collide on (pathway, gene) but agree on weight (harmless duplicates).")

    W = np.zeros((len(pathways), N_GENES), dtype=np.float32)
    for _, r in matched.iterrows():
        W[pw_to_row[r["source"]], int(r["tdc_idx"])] = r["weight"]

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        OUT_DIR / "progeny_tdc_weights.npz",
        W=W,
        pathways=np.array(pathways),
        gene_symbol=tdc_index["gene_symbol"].to_numpy(),
        entrez=tdc_index["entrez"].astype(float).to_numpy(),
    )
    audit = aligned[["source", "target", "entrez", "tdc_idx", "weight", "matched_via"]].rename(
        columns={"source": "pathway", "target": "gene_symbol"}
    )
    audit.to_csv(OUT_DIR / "progeny_tdc_alignment.csv", index=False)

    # coverage report (per unique gene)
    via = resolved["matched_via"].value_counts()
    n_genes = len(resolved)
    print(f"PROGENy genes: {n_genes} unique")
    print(f"  matched by symbol: {via.get('symbol', 0)}")
    print(f"  matched by alias : {via.get('alias', 0)}")
    print(f"  unmatched        : {via.get('unmatched', 0)}  {sorted(resolved.loc[resolved.matched_via=='unmatched','gene_symbol'])[:10]}")
    print(f"  -> matched {n_genes - via.get('unmatched', 0)}/{n_genes} = {100*(n_genes-via.get('unmatched',0))/n_genes:.1f}%")
    print(f"W shape {W.shape}; pathways retained {matched['source'].nunique()}/{len(pathways)}; nonzero entries {int((W!=0).sum())}")


if __name__ == "__main__":
    main()
