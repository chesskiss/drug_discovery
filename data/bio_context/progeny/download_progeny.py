"""Download PROGENy pathway-responsive gene weights.

PROGENy (Pathway RespOnsive GENes) provides gene weight signatures for 14
signaling pathways, derived from a large compendium of perturbation
experiments. See https://saezlab.github.io/progeny/ for background.

This script fetches the resource directly from the OmniPath REST API
(the same source the official `decoupler`/`decoupleR` packages wrap), so no
heavyweight client library is required.

Usage:
    uv run --project data/bio_context progeny/download_progeny.py
    uv run --project data/bio_context progeny/download_progeny.py --organism mouse --top 100
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import requests

OMNIPATH_URL = "https://omnipathdb.org/annotations?databases=PROGENy"
OUT_DIR = Path(__file__).parent / "data"


def fetch_progeny(organism: str = "human", license: str = "academic") -> pd.DataFrame:
    """Fetch raw PROGENy annotations from OmniPath and pivot into long format.

    Returns a dataframe with columns: source (pathway), target (gene), weight, padj.
    """
    if organism != "human":
        raise NotImplementedError(
            "Only 'human' is supported here (no ortholog translation). "
            "For mouse, use the decoupler package's dc.op.progeny(organism='mouse')."
        )

    resp = requests.get(f"{OMNIPATH_URL}&license={license}", timeout=120)
    resp.raise_for_status()

    raw = pd.read_csv(pd.io.common.BytesIO(resp.content), sep="\t")
    wide = raw.pivot(index=["genesymbol", "record_id"], columns="label", values="value").reset_index()
    wide = wide.rename(columns={"genesymbol": "target", "pathway": "source"})
    wide["weight"] = wide["weight"].astype(float)
    wide["p_value"] = wide["p_value"].astype(float)
    wide = wide.rename(columns={"p_value": "padj"})
    return wide[["source", "target", "weight", "padj"]]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--organism", default="human", choices=["human"])
    parser.add_argument("--license", default="academic", choices=["academic", "commercial", "nonprofit"])
    parser.add_argument("--thr-padj", type=float, default=0.05, help="Significance threshold to keep an interaction.")
    parser.add_argument("--top", type=int, default=100, help="Top N genes per pathway (by p-value) to keep.")
    parser.add_argument("--out", type=Path, default=None, help="Output CSV path (default: data/progeny_<organism>_top<top>.csv)")
    args = parser.parse_args()

    df = fetch_progeny(organism=args.organism, license=args.license)
    df = df[df["padj"] < args.thr_padj]
    df = (
        df.sort_values("padj")
        .groupby("source", group_keys=False)
        .head(args.top)
        .sort_values(["source", "padj"])
        .reset_index(drop=True)
    )

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out_path = args.out or OUT_DIR / f"progeny_{args.organism}_top{args.top}.csv"
    df.to_csv(out_path, index=False)
    print(f"Saved {len(df)} pathway-gene interactions ({df['source'].nunique()} pathways) to {out_path}")


if __name__ == "__main__":
    main()
