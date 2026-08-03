from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from .oca_plots import plot_component_importance_head_tail_summary


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Backfill OCA head/tail summary plots from component_importance.csv")
    parser.add_argument("--oca-dir", type=str, required=True, help="Directory containing component_importance.csv")
    parser.add_argument("--head-k", type=int, default=10)
    parser.add_argument("--tail-k", type=int, default=5)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    oca_dir = Path(args.oca_dir)
    csv_path = oca_dir / "component_importance.csv"
    if not csv_path.exists():
        raise FileNotFoundError(f"Missing component importance CSV: {csv_path}")

    output_path = oca_dir / "component_importance_head_tail_summary.png"
    component_importance = pd.read_csv(csv_path)
    plot_component_importance_head_tail_summary(
        component_importance,
        output_path,
        head_k=args.head_k,
        tail_k=args.tail_k,
    )
    print(f"[oca-backfill] Saved summary plot to {output_path}")


if __name__ == "__main__":
    main()
