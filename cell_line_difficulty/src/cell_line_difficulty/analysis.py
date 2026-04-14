from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd


REQUIRED_COLUMNS = ("Cell_Line_ID", "Synergy_ZIP")


def _validate_input_columns(df: pd.DataFrame) -> None:
    missing = [column for column in REQUIRED_COLUMNS if column not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")


def _normalize_series(series: pd.Series) -> pd.Series:
    std = float(series.std(ddof=0))
    if std == 0 or np.isnan(std):
        return pd.Series(np.zeros(len(series), dtype=float), index=series.index)
    return (series - float(series.mean())) / std


def _compute_reliability_weight(counts: pd.Series) -> pd.Series:
    max_count = float(counts.max())
    if max_count <= 0:
        return pd.Series(np.zeros(len(counts), dtype=float), index=counts.index)
    return np.sqrt(counts / max_count)


def _build_summary(input_path: str, metrics: pd.DataFrame) -> dict:
    return {
        "input": input_path,
        "n_cell_lines": int(len(metrics)),
        "n_total_rows": int(metrics["n_rows"].sum()),
        "zip_gt_0_threshold": 0.0,
        "zip_gt_10_threshold": 10.0,
        "score_definition": {
            "ease_score": "reliability_weight * mean(z_mean_zip, z_high_synergy_rate_zip_gt_10, z_positive_rate_zip_gt_0)",
            "difficulty_score": "-ease_score",
        },
        "guardrails": [
            "This is an observational screening heuristic, not a biological causality estimate.",
            "Drug identity is intentionally ignored, so rankings mix cell susceptibility with screen coverage.",
            "Low-sample cell lines are down-weighted via the reliability term.",
        ],
        "top_5_easiest": metrics.sort_values("ease_rank").head(5)["cell_line"].tolist(),
        "top_5_hardest": metrics.sort_values("difficulty_rank").head(5)["cell_line"].tolist(),
    }


def analyze_cell_line_difficulty(input_path: str | Path) -> pd.DataFrame:
    df = pd.read_csv(input_path)
    _validate_input_columns(df)

    analysis_df = df.loc[:, ["Cell_Line_ID", "Synergy_ZIP"]].copy()
    analysis_df["Synergy_ZIP"] = pd.to_numeric(analysis_df["Synergy_ZIP"], errors="coerce")
    analysis_df = analysis_df.dropna(subset=["Cell_Line_ID", "Synergy_ZIP"]).copy()

    grouped = analysis_df.groupby("Cell_Line_ID")["Synergy_ZIP"]
    metrics = grouped.agg(["count", "mean", "median", "std", "max"]).rename(
        columns={
            "count": "n_rows",
            "mean": "mean_zip",
            "median": "median_zip",
            "std": "std_zip",
            "max": "max_zip",
        }
    )
    metrics["std_zip"] = metrics["std_zip"].fillna(0.0)
    metrics["q75_zip"] = grouped.quantile(0.75)
    metrics["q90_zip"] = grouped.quantile(0.90)
    metrics["positive_rate_zip_gt_0"] = grouped.apply(lambda values: float((values > 0).mean()))
    metrics["high_synergy_rate_zip_gt_10"] = grouped.apply(lambda values: float((values > 10).mean()))
    metrics = metrics.reset_index(names="cell_line")

    metrics["reliability_weight"] = _compute_reliability_weight(metrics["n_rows"])
    metrics["z_mean_zip"] = _normalize_series(metrics["mean_zip"])
    metrics["z_positive_rate_zip_gt_0"] = _normalize_series(metrics["positive_rate_zip_gt_0"])
    metrics["z_high_synergy_rate_zip_gt_10"] = _normalize_series(metrics["high_synergy_rate_zip_gt_10"])

    metrics["ease_signal"] = (
        metrics["z_mean_zip"]
        + metrics["z_positive_rate_zip_gt_0"]
        + metrics["z_high_synergy_rate_zip_gt_10"]
    ) / 3.0
    metrics["ease_score"] = metrics["ease_signal"] * metrics["reliability_weight"]
    metrics["difficulty_score"] = -metrics["ease_score"]
    metrics["ease_rank"] = metrics["ease_score"].rank(method="dense", ascending=False).astype(int)
    metrics["difficulty_rank"] = metrics["difficulty_score"].rank(method="dense", ascending=False).astype(int)

    return metrics.sort_values(["ease_rank", "cell_line"]).reset_index(drop=True)


def save_analysis_outputs(metrics: pd.DataFrame, input_path: str | Path, output_dir: str | Path) -> None:
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    ranked_by_ease = metrics.sort_values(["ease_rank", "cell_line"]).reset_index(drop=True)
    ranked_by_difficulty = metrics.sort_values(["difficulty_rank", "cell_line"], ascending=[True, True]).reset_index(drop=True)

    ranked_by_ease.to_csv(output_path / "ranked_by_ease.csv", index=False)
    ranked_by_difficulty.to_csv(output_path / "ranked_by_difficulty.csv", index=False)

    summary = _build_summary(str(input_path), metrics)
    with open(output_path / "summary.json", "w", encoding="utf-8") as file:
        json.dump(summary, file, indent=2)

