from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Download DrugComb synergy rows and TDC cell-feature tables."
    )
    parser.add_argument(
        "--save-dir",
        type=str,
        default="data_compression/zscore_var_pca_128/data",
        help="Directory where downloaded files will be written",
    )
    parser.add_argument(
        "--raw-filename",
        type=str,
        default="drugcomb_raw.csv",
        help="Filename for the raw DrugComb table returned by TDC",
    )
    parser.add_argument(
        "--reduced-filename",
        type=str,
        default="reduced_drugcomb.csv",
        help="Filename for the reduced DrugComb table with normalized column names",
    )
    parser.add_argument(
        "--cell-features-filename",
        type=str,
        default="tdc_cell_features.csv",
        help="Filename for the TDC cell-feature table. This is not the embedded 3-view `CellLine` payload export.",
    )
    return parser.parse_args()


def _pick_first_existing(columns: list[str], candidates: list[str], field_name: str) -> str:
    normalized = {column.lower(): column for column in columns}
    for candidate in candidates:
        match = normalized.get(candidate.lower())
        if match is not None:
            return match
    raise ValueError(f"Could not detect {field_name}. Tried: {candidates}")


def _build_reduced_frame(df: pd.DataFrame) -> pd.DataFrame:
    columns = list(df.columns)
    smiles_a_col = _pick_first_existing(columns, ["smiles_a", "Drug1", "drug1", "drug_a"], "smiles_a")
    smiles_b_col = _pick_first_existing(columns, ["smiles_b", "Drug2", "drug2", "drug_b"], "smiles_b")
    cell_col = _pick_first_existing(columns, ["cell_line", "Cell_Line_ID", "CellLine", "cell"], "cell line")
    target_col = _pick_first_existing(columns, ["Synergy_ZIP", "target", "zip", "ZIP"], "target")

    reduced = pd.DataFrame(
        {
            "smiles_a": df[smiles_a_col].astype(str),
            "smiles_b": df[smiles_b_col].astype(str),
            "cell_line": df[cell_col].astype(str),
            "target": pd.to_numeric(df[target_col], errors="coerce"),
        }
    )
    return reduced.dropna(subset=["smiles_a", "smiles_b", "cell_line", "target"]).reset_index(drop=True)


def download_drugcomb_with_genes(
    save_dir: str = "data_compression/zscore_var_pca_128/data",
    *,
    raw_filename: str = "drugcomb_raw.csv",
    reduced_filename: str = "reduced_drugcomb.csv",
    cell_features_filename: str = "tdc_cell_features.csv",
) -> None:
    out_dir = Path(save_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    try:
        from tdc.multi_pred import DrugSyn
    except ImportError as exc:
        raise ImportError("TDC is required for this script. Install it with: uv add PyTDC") from exc

    print("Downloading DrugComb dataset from TDC...")
    data = DrugSyn(name="DrugComb")

    df = data.get_data()
    df_path = out_dir / raw_filename
    df.to_csv(df_path, index=False)
    print(f"Saved raw DrugComb data to {df_path}")

    reduced_df = _build_reduced_frame(df)
    reduced_path = out_dir / reduced_filename
    reduced_df.to_csv(reduced_path, index=False)
    print(f"Saved reduced DrugComb data to {reduced_path}")

    print("Downloading TDC cell features...")
    get_cell_features = getattr(data, "get_cell_features", None)
    if callable(get_cell_features):
        cell_expr = get_cell_features()
        expr_path = out_dir / cell_features_filename
        cell_expr.to_csv(expr_path)
        print(f"Saved gene expression data to {expr_path}")
        print("Note: this is the TDC cell-feature table, not the embedded 3-view `CellLine` payload export.")
    else:
        cell_expr = None
        print("TDC object does not expose `get_cell_features()` in this installed version.")
        print("Skipped cell-feature download. Raw and reduced DrugComb tables were still saved.")

    print("\nSanity checks:")
    print(f"Synergy shape: {df.shape}")
    print(f"Reduced shape: {reduced_df.shape}")
    if cell_expr is not None:
        print(f"Gene expression shape: {cell_expr.shape}")

    cell_col = None
    for candidate in ("cell_line", "Cell_Line_ID", "cell_line_id"):
        if candidate in df.columns:
            cell_col = candidate
            break

    if cell_col is None:
        print("WARNING: Could not detect a cell-line column in the synergy table.")
        return

    if cell_expr is not None:
        example_cell = str(df[cell_col].iloc[0])
        expr_index = cell_expr.index.astype(str)
        expr_columns = cell_expr.columns.astype(str)
        if example_cell in expr_index or example_cell in expr_columns:
            print(f"Example cell line '{example_cell}' found in gene expression table.")
        else:
            print(f"WARNING: '{example_cell}' not found in gene expression table.")


def main() -> None:
    args = parse_args()
    download_drugcomb_with_genes(
        save_dir=args.save_dir,
        raw_filename=args.raw_filename,
        reduced_filename=args.reduced_filename,
        cell_features_filename=args.cell_features_filename,
    )


if __name__ == "__main__":
    main()
