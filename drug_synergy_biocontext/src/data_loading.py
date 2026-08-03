from __future__ import annotations

import ast
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd


SMILES_A_CANDIDATES = ["smiles_a", "Drug1", "drug1", "drug_a"]
SMILES_B_CANDIDATES = ["smiles_b", "Drug2", "drug2", "drug_b"]
CELL_CANDIDATES = ["cell_line", "Cell_Line_ID", "CellLine", "cell"]
TARGET_CANDIDATES = ["Synergy_ZIP", "target", "zip", "ZIP"]


@dataclass(frozen=True)
class DrugCombSchema:
    smiles_a_col: str
    smiles_b_col: str
    cell_col: str
    target_col: str


def _pick_first_existing(columns: list[str], candidates: list[str], field_name: str) -> str:
    normalized = {column.lower(): column for column in columns}
    for candidate in candidates:
        match = normalized.get(candidate.lower())
        if match is not None:
            return match
    raise ValueError(f"Could not detect {field_name}. Tried: {candidates}")


def detect_schema(df: pd.DataFrame) -> DrugCombSchema:
    columns = list(df.columns)
    return DrugCombSchema(
        smiles_a_col=_pick_first_existing(columns, SMILES_A_CANDIDATES, "smiles_a"),
        smiles_b_col=_pick_first_existing(columns, SMILES_B_CANDIDATES, "smiles_b"),
        cell_col=_pick_first_existing(columns, CELL_CANDIDATES, "cell line"),
        target_col=_pick_first_existing(columns, TARGET_CANDIDATES, "target"),
    )


def load_table(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    suffix = path.suffix.lower()
    if suffix == ".csv":
        return pd.read_csv(path)
    if suffix in {".tsv", ".txt"}:
        return pd.read_csv(path, sep="\t")
    if suffix == ".pkl":
        return pd.read_pickle(path)
    raise ValueError(f"Unsupported file format: {path}")


def load_synergy_table(path: str | Path) -> pd.DataFrame:
    df = load_table(path)
    schema = detect_schema(df)

    frame = pd.DataFrame(
        {
            "smiles_a": df[schema.smiles_a_col].astype(str),
            "smiles_b": df[schema.smiles_b_col].astype(str),
            "cell_line": df[schema.cell_col].astype(str),
            "target": pd.to_numeric(df[schema.target_col], errors="coerce"),
        }
    )
    frame = frame.dropna(subset=["smiles_a", "smiles_b", "cell_line", "target"]).reset_index(drop=True)
    return frame


def load_cell_expression_table(path: str | Path) -> dict[str, np.ndarray]:
    df = pd.read_csv(path)
    cell_col = _pick_first_existing(list(df.columns), CELL_CANDIDATES, "cell line")

    feature_df = df.drop(columns=[cell_col]).copy()
    feature_df = feature_df.apply(pd.to_numeric, errors="coerce")
    feature_df = feature_df.dropna(axis=1, how="all")

    expr = {}
    for cell_line, values in zip(df[cell_col].astype(str), feature_df.to_numpy(dtype=np.float32), strict=False):
        expr[cell_line] = values
    return expr


def _parse_cellline_payload(value: object) -> list[np.ndarray]:
    if isinstance(value, (list, tuple)):
        return [np.asarray(item, dtype=np.float32) for item in value]

    if isinstance(value, str):
        if "..." in value:
            raise ValueError(
                "Embedded CellLine CSV text is truncated with ellipsis. "
                "A non-truncated source (for example the original pickle/export) is required for views 0 and 1."
            )
        try:
            parsed = eval(value, {"__builtins__": {}}, {"array": np.array})  # noqa: S307
        except Exception as exc:  # pragma: no cover - defensive parse guard
            raise ValueError("Could not parse embedded CellLine payload from CSV text.") from exc
        if isinstance(parsed, (list, tuple)):
            return [np.asarray(item, dtype=np.float32) for item in parsed]

    raise TypeError("Unsupported CellLine payload; expected list/tuple or a parseable string representation.")


def load_embedded_cellline_table(
    path: str | Path,
    *,
    feature_view_index: int = 0,
) -> dict[str, np.ndarray]:
    df = load_table(path)
    columns = list(df.columns)
    cell_col = _pick_first_existing(columns, ["Cell_Line_ID", "cell_line", "CellLine"], "cell line")
    embedded_col = _pick_first_existing(columns, ["CellLine"], "embedded cell features")

    lookup: dict[str, np.ndarray] = {}
    for _, row in df[[cell_col, embedded_col]].dropna().iterrows():
        cell_line = str(row[cell_col])
        features = _parse_cellline_payload(row[embedded_col])
        if feature_view_index >= len(features):
            raise IndexError(
                f"Requested CellLine view index {feature_view_index}, but only {len(features)} views are available."
            )
        lookup.setdefault(cell_line, np.asarray(features[feature_view_index], dtype=np.float32))
    return lookup


def load_expression_lookup(
    cell_expression_path: str | Path | None,
    fallback_pickle_path: str | Path | None = None,
    feature_view_index: int = 0,
) -> dict[str, np.ndarray]:
    if cell_expression_path is not None and Path(cell_expression_path).exists():
        preview = load_table(cell_expression_path).head(1)
        if "CellLine" in preview.columns:
            lookup = load_embedded_cellline_table(
                cell_expression_path,
                feature_view_index=feature_view_index,
            )
            if lookup:
                return lookup
        return load_cell_expression_table(cell_expression_path)

    if fallback_pickle_path is not None and Path(fallback_pickle_path).exists():
        lookup = load_embedded_cellline_table(
            fallback_pickle_path,
            feature_view_index=feature_view_index,
        )
        if lookup:
            return lookup

    raise FileNotFoundError(
        "Could not load cell expression features. Expected a numeric cell-expression table "
        "or a file with embedded `CellLine` payloads."
    )


def get_cellline_array_lengths(value: object) -> list[int]:
    if isinstance(value, (list, tuple)):
        return [len(np.asarray(item)) for item in value]

    if isinstance(value, str):
        if "..." in value:
            raise ValueError("CellLine CSV text is truncated with ellipsis; use the pickle-backed row instead.")
        parsed = ast.literal_eval(value)
        if isinstance(parsed, (list, tuple)):
            return [len(np.asarray(item)) for item in parsed]

    raise TypeError("Unsupported CellLine value; expected a list/tuple of arrays or a non-truncated string.")


def get_first_cellline_array_lengths(
    csv_path: str | Path,
    pickle_path: str | Path | None = None,
) -> list[int]:
    csv_df = pd.read_csv(csv_path, nrows=1)
    first_value = csv_df.iloc[0]["CellLine"]

    try:
        return get_cellline_array_lengths(first_value)
    except ValueError:
        if pickle_path is None:
            raise
        pickle_df = load_table(pickle_path)
        return get_cellline_array_lengths(pickle_df.iloc[0]["CellLine"])
