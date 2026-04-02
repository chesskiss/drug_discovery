from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from .data_loading import detect_schema


DRUG_A_FALLBACK = ["Drug1", "Drug1_ID", "drug1", "drug_a", "drug_a_id"]
DRUG_B_FALLBACK = ["Drug2", "Drug2_ID", "drug2", "drug_b", "drug_b_id"]
CELL_FALLBACK = ["CellLine", "Cell_Line_ID", "cell_line", "cell", "cellline"]
TARGET_FALLBACK = ["Synergy_ZIP", "ZIP", "zip", "Synergy_Bliss", "Synergy_HSA", "Synergy_Loewe", "CSS"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Visualize DrugComb-style tabular data")
    parser.add_argument("--input", type=str, default="outputs/drugcomb.csv", help="Input file (csv/tsv/txt/parquet)")
    parser.add_argument("--output-dir", type=str, default="outputs/visualization", help="Directory for plots and summary")
    parser.add_argument(
        "--preprocess",
        action="store_true",
        help="Apply light preprocessing before visualization",
    )
    return parser.parse_args()


def load_table(path: Path) -> pd.DataFrame:
    suffix = path.suffix.lower()
    if suffix == ".csv":
        return pd.read_csv(path)
    if suffix in {".tsv", ".txt"}:
        return pd.read_csv(path, sep="\t")
    if suffix == ".parquet":
        return pd.read_parquet(path)
    raise ValueError(f"Unsupported input format: {path}")


def save_json(path: Path, payload: dict) -> None:
    with open(path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)


def pick_first_existing(columns: list[str], candidates: list[str], field_name: str) -> str:
    normalized = {c.lower(): c for c in columns}
    for candidate in candidates:
        key = candidate.lower()
        if key in normalized:
            return normalized[key]
    raise ValueError(f"Could not detect {field_name}. Candidates: {candidates}")


def detect_columns_for_visualization(df: pd.DataFrame) -> tuple[str, str, str, str]:
    columns = list(df.columns)
    try:
        schema = detect_schema(df)
        return schema.smiles_a_col, schema.smiles_b_col, schema.cell_col, schema.target_col
    except Exception:
        drug_a_col = pick_first_existing(columns, DRUG_A_FALLBACK, "drug_a")
        drug_b_col = pick_first_existing(columns, DRUG_B_FALLBACK, "drug_b")
        cell_col = pick_first_existing(columns, CELL_FALLBACK, "cell")
        target_col = pick_first_existing(columns, TARGET_FALLBACK, "target")
        return drug_a_col, drug_b_col, cell_col, target_col


def main() -> None:
    args = parse_args()

    input_path = Path(args.input)
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    df_raw = load_table(input_path)
    rows_before = len(df_raw)

    print(f"Loaded: {input_path}")
    print(f"Shape: {df_raw.shape}")
    print("Columns:", list(df_raw.columns))

    drug_a_col, drug_b_col, cell_col, target_col = detect_columns_for_visualization(df_raw)

    df = df_raw.copy()

    if args.preprocess:
        df = df.dropna(subset=[drug_a_col, drug_b_col, cell_col, target_col]).copy()
        df[target_col] = pd.to_numeric(df[target_col], errors="coerce")
        df = df.dropna(subset=[target_col]).copy()

    rows_after = len(df)

    summary = {
        "input": str(input_path),
        "rows_before": rows_before,
        "rows_after": rows_after,
        "cols": len(df.columns),
        "detected_columns": {
            "drug_a": drug_a_col,
            "drug_b": drug_b_col,
            "cell_line": cell_col,
            "target": target_col,
        },
        "missing_rate_top10": df_raw.isna().mean().sort_values(ascending=False).head(10).to_dict(),
        "target_describe": pd.to_numeric(df[target_col], errors="coerce").describe().to_dict(),
        "unique_cells": int(df[cell_col].nunique()),
        "unique_drug_a": int(df[drug_a_col].nunique()),
        "unique_drug_b": int(df[drug_b_col].nunique()),
    }

    save_json(out_dir / "summary.json", summary)

    plt.figure(figsize=(7, 4))
    plt.hist(pd.to_numeric(df[target_col], errors="coerce").dropna(), bins=50)
    plt.title(f"Target Distribution: {target_col}")
    plt.xlabel(target_col)
    plt.ylabel("Count")
    plt.tight_layout()
    plt.savefig(out_dir / "target_distribution.png", dpi=150)
    plt.close()

    top_cells = df[cell_col].value_counts().head(20)
    plt.figure(figsize=(10, 5))
    top_cells.sort_values().plot(kind="barh")
    plt.title("Top 20 Cell Lines by Sample Count")
    plt.xlabel("Count")
    plt.ylabel("Cell line")
    plt.tight_layout()
    plt.savefig(out_dir / "top_cells.png", dpi=150)
    plt.close()

    pair_count = df.groupby([drug_a_col, drug_b_col], dropna=False).size().sort_values(ascending=False)
    top_pair_counts = pair_count.head(30)
    plt.figure(figsize=(8, 4))
    plt.hist(top_pair_counts.values, bins=min(20, len(top_pair_counts)))
    plt.title("Histogram of Top 30 Pair Frequencies")
    plt.xlabel("Count per pair")
    plt.ylabel("Frequency")
    plt.tight_layout()
    plt.savefig(out_dir / "top_pair_frequency_hist.png", dpi=150)
    plt.close()

    print(f"Detected columns: drug_a={drug_a_col}, drug_b={drug_b_col}, cell={cell_col}, target={target_col}")
    print(f"Saved summary: {out_dir / 'summary.json'}")
    print(f"Saved plots in: {out_dir}")


if __name__ == "__main__":
    main()
