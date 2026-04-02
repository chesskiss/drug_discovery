from __future__ import annotations

import argparse
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Download DrugComb synergy and cell features via TDC")
    parser.add_argument(
        "--save-dir",
        type=str,
        default="data",
        help="Directory where downloaded files will be written",
    )
    return parser.parse_args()


def download_drugcomb_with_genes(save_dir: str = "data") -> None:
    out_dir = Path(save_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    try:
        from tdc.multi_pred import DrugSyn
    except ImportError as exc:
        raise ImportError("TDC is required for this script. Install it with: uv add PyTDC") from exc

    print("Downloading DrugComb dataset from TDC...")
    data = DrugSyn(name="DrugComb")

    df = data.get_data()
    df_path = out_dir / "drugcomb_synergy.csv"
    df.to_csv(df_path, index=False)
    print(f"Saved synergy data to {df_path}")

    print("Downloading gene expression (cell line features)...")
    cell_expr = data.get_cell_features()
    expr_path = out_dir / "cell_line_gene_expression.csv"
    cell_expr.to_csv(expr_path)
    print(f"Saved gene expression data to {expr_path}")

    print("\nSanity checks:")
    print(f"Synergy shape: {df.shape}")
    print(f"Gene expression shape: {cell_expr.shape}")

    cell_col = None
    for candidate in ("cell_line", "Cell_Line_ID", "cell_line_id"):
        if candidate in df.columns:
            cell_col = candidate
            break

    if cell_col is None:
        print("WARNING: Could not detect a cell-line column in the synergy table.")
        return

    example_cell = str(df[cell_col].iloc[0])
    expr_index = cell_expr.index.astype(str)
    expr_columns = cell_expr.columns.astype(str)
    if example_cell in expr_index or example_cell in expr_columns:
        print(f"Example cell line '{example_cell}' found in gene expression table.")
    else:
        print(f"WARNING: '{example_cell}' not found in gene expression table.")


def main() -> None:
    args = parse_args()
    download_drugcomb_with_genes(save_dir=args.save_dir)


if __name__ == "__main__":
    main()
