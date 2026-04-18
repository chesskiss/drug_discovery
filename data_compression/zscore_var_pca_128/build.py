from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


DEFAULT_INPUT_PICKLE = "drug_synergy_baseline/data/drugcomb.pkl"
DEFAULT_INPUT_RAW_CSV = "data_compression/zscore_var_pca_128/data/drugcomb_raw.csv"
DEFAULT_OUTPUT_CSV = "data_compression/zscore_var_pca_128/data/drugcomb.csv"
DEFAULT_METADATA_JSON = "data_compression/zscore_var_pca_128/data/metadata.json"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a compressed DrugComb export with z-score + variance selection + PCA."
    )
    parser.add_argument(
        "--input-pickle",
        type=str,
        default=DEFAULT_INPUT_PICKLE,
        help="Non-truncated source with embedded CellLine payloads.",
    )
    parser.add_argument(
        "--input-raw-csv",
        type=str,
        default=DEFAULT_INPUT_RAW_CSV,
        help="Raw CSV path kept for provenance. It is not used for raw CellLine vectors because CSV payloads are truncated.",
    )
    parser.add_argument(
        "--output-csv",
        type=str,
        default=DEFAULT_OUTPUT_CSV,
        help="Where to write the compressed DrugComb CSV.",
    )
    parser.add_argument(
        "--metadata-json",
        type=str,
        default=DEFAULT_METADATA_JSON,
        help="Where to write pipeline metadata.",
    )
    parser.add_argument(
        "--raw-view-index",
        type=int,
        default=0,
        help="Which CellLine view to treat as the raw input.",
    )
    parser.add_argument(
        "--variance-top-k",
        type=int,
        default=3000,
        help="How many raw features to keep before PCA.",
    )
    parser.add_argument(
        "--pca-components",
        type=int,
        default=128,
        help="Requested output dimensionality. Extra dimensions are zero-padded if PCA rank is smaller.",
    )
    return parser.parse_args()


def _extract_view_vector(value: object, view_index: int) -> np.ndarray:
    if not isinstance(value, (list, tuple)):
        raise TypeError("Expected CellLine payload to already be a list/tuple from the pickle source.")
    if view_index >= len(value):
        raise IndexError(f"Requested CellLine view {view_index}, but only {len(value)} views are available.")
    return np.asarray(value[view_index], dtype=np.float32)


def _build_unique_cell_matrix(df: pd.DataFrame, raw_view_index: int) -> tuple[list[str], np.ndarray, np.ndarray]:
    grouped = df.groupby("Cell_Line_ID", sort=False, dropna=True)

    cell_lines: list[str] = []
    vectors: list[np.ndarray] = []
    counts: list[int] = []

    for cell_line, group in grouped:
        first_payload = group.iloc[0]["CellLine"]
        vector = _extract_view_vector(first_payload, raw_view_index)
        cell_lines.append(str(cell_line))
        vectors.append(vector)
        counts.append(int(len(group)))

    matrix = np.vstack(vectors).astype(np.float32)
    weights = np.asarray(counts, dtype=np.float64)
    return cell_lines, matrix, weights


def _weighted_mean_and_var(matrix: np.ndarray, weights: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    normalized = weights / weights.sum()
    mean = np.sum(matrix * normalized[:, None], axis=0)
    centered = matrix - mean
    var = np.sum((centered**2) * normalized[:, None], axis=0)
    return mean.astype(np.float32), var.astype(np.float32)


def _weighted_standardize(matrix: np.ndarray, weights: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    mean, var = _weighted_mean_and_var(matrix, weights)
    std = np.sqrt(var).astype(np.float32)
    std[std == 0] = 1.0
    scaled = (matrix - mean) / std
    return scaled.astype(np.float32), mean, std


def _weighted_pca_scores(
    scaled_matrix: np.ndarray,
    weights: np.ndarray,
    requested_components: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, int]:
    normalized = weights / weights.sum()
    centered = scaled_matrix - np.sum(scaled_matrix * normalized[:, None], axis=0)
    weighted_centered = centered * np.sqrt(normalized)[:, None]

    _, singular_values, vt = np.linalg.svd(weighted_centered, full_matrices=False)
    max_rank = int(min(requested_components, vt.shape[0], max(0, scaled_matrix.shape[0] - 1)))

    if max_rank == 0:
        raise ValueError("PCA rank is zero. There are not enough distinct cell lines to fit the pipeline.")

    components = vt[:max_rank].astype(np.float32)
    scores = centered @ components.T
    explained_variance = (singular_values[:max_rank] ** 2).astype(np.float32)
    return scores.astype(np.float32), components, explained_variance, max_rank


def _pad_scores(scores: np.ndarray, output_dim: int) -> np.ndarray:
    if scores.shape[1] >= output_dim:
        return scores[:, :output_dim]
    padded = np.zeros((scores.shape[0], output_dim), dtype=np.float32)
    padded[:, : scores.shape[1]] = scores
    return padded


def _serialize_single_view(vector: np.ndarray) -> str:
    values = ",".join(f"{float(value):.8g}" for value in vector.astype(np.float32))
    return f"[array([{values}])]"


def build_compressed_export(
    input_pickle: str | Path,
    output_csv: str | Path,
    metadata_json: str | Path,
    *,
    input_raw_csv: str | Path,
    raw_view_index: int,
    variance_top_k: int,
    pca_components: int,
) -> dict[str, object]:
    input_pickle = Path(input_pickle)
    input_raw_csv = Path(input_raw_csv)
    output_csv = Path(output_csv)
    metadata_json = Path(metadata_json)

    if not input_pickle.exists():
        raise FileNotFoundError(f"Pickle source not found: {input_pickle}")

    df = pd.read_pickle(input_pickle)
    required_columns = {"Cell_Line_ID", "CellLine"}
    missing = required_columns.difference(df.columns)
    if missing:
        raise ValueError(f"Pickle source is missing required columns: {sorted(missing)}")

    cell_lines, raw_matrix, weights = _build_unique_cell_matrix(df, raw_view_index)
    _, raw_variance = _weighted_mean_and_var(raw_matrix, weights)

    non_constant_indices = np.flatnonzero(raw_variance > 0)
    if len(non_constant_indices) == 0:
        raise ValueError("All raw CellLine features are constant.")

    top_k = min(int(variance_top_k), len(non_constant_indices))
    ranked = non_constant_indices[np.argsort(raw_variance[non_constant_indices])[::-1]]
    selected_indices = np.sort(ranked[:top_k])

    selected_matrix = raw_matrix[:, selected_indices]
    scaled_matrix, selected_mean, selected_std = _weighted_standardize(selected_matrix, weights)
    scores, components, explained_variance, effective_components = _weighted_pca_scores(
        scaled_matrix,
        weights,
        pca_components,
    )
    padded_scores = _pad_scores(scores, pca_components)

    compressed_lookup = {
        cell_line: padded_scores[idx]
        for idx, cell_line in enumerate(cell_lines)
    }

    output_df = df.copy()
    output_df["CellLine"] = output_df["Cell_Line_ID"].astype(str).map(
        lambda cell_line: _serialize_single_view(compressed_lookup[str(cell_line)])
    )

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    metadata_json.parent.mkdir(parents=True, exist_ok=True)
    output_df.to_csv(output_csv, index=False)

    metadata = {
        "method_id": "1.3A",
        "pipeline_name": "zscore_var_pca_128",
        "input_pickle": str(input_pickle),
        "input_raw_csv": str(input_raw_csv),
        "output_csv": str(output_csv),
        "raw_view_index": int(raw_view_index),
        "raw_feature_dim": int(raw_matrix.shape[1]),
        "unique_cell_lines": int(len(cell_lines)),
        "row_count": int(len(df)),
        "variance_top_k_requested": int(variance_top_k),
        "variance_top_k_used": int(len(selected_indices)),
        "pca_components_requested": int(pca_components),
        "effective_pca_components": int(effective_components),
        "output_feature_dim": int(padded_scores.shape[1]),
        "zero_padded_dimensions": int(padded_scores.shape[1] - effective_components),
        "notes": [
            "Variance ranking is computed on the raw features before z-score; otherwise variance selection is degenerate.",
            "PCA rank is limited by the number of distinct cell lines in DrugComb. Remaining output dimensions are zero-padded.",
            "The exported CSV keeps a single compressed CellLine view in baseline-compatible `[array([...])]` format.",
        ],
        "selected_feature_indices_preview": selected_indices[:20].tolist(),
        "explained_variance_preview": explained_variance[:20].astype(float).tolist(),
        "selected_feature_mean_preview": selected_mean[:10].astype(float).tolist(),
        "selected_feature_std_preview": selected_std[:10].astype(float).tolist(),
        "component_matrix_shape": [int(x) for x in components.shape],
    }

    with metadata_json.open("w", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2)

    return metadata


def main() -> None:
    args = parse_args()
    metadata = build_compressed_export(
        input_pickle=args.input_pickle,
        input_raw_csv=args.input_raw_csv,
        output_csv=args.output_csv,
        metadata_json=args.metadata_json,
        raw_view_index=args.raw_view_index,
        variance_top_k=args.variance_top_k,
        pca_components=args.pca_components,
    )
    print(f"[compression] Wrote compressed DrugComb CSV to {metadata['output_csv']}")
    print(f"[compression] Unique cell lines: {metadata['unique_cell_lines']}")
    print(f"[compression] Raw feature dim: {metadata['raw_feature_dim']}")
    print(f"[compression] Top-k selected: {metadata['variance_top_k_used']}")
    print(f"[compression] Effective PCA components: {metadata['effective_pca_components']}")
    print(f"[compression] Output feature dim: {metadata['output_feature_dim']}")


if __name__ == "__main__":
    main()
