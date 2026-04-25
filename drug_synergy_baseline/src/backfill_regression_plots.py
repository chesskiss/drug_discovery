from __future__ import annotations

import argparse
import math
from pathlib import Path

import pandas as pd

from .training_artifacts import build_regression_figure_path, save_regression_plot


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Backfill test-set regression plots from existing prediction CSV artifacts."
    )
    parser.add_argument(
        "--outputs-dir",
        type=str,
        default="outputs",
        help="Root outputs directory to scan.",
    )
    parser.add_argument(
        "--run",
        dest="runs",
        action="append",
        default=None,
        help="Specific run directory name under outputs/ to process. Repeatable. If omitted, scan all runs.",
    )
    parser.add_argument(
        "--kind",
        choices=["single", "cv", "all"],
        default="all",
        help="Which prediction artifact types to process.",
    )
    return parser.parse_args()


def compute_mse(frame: pd.DataFrame) -> float:
    if frame.empty:
        return float("nan")
    return float(((frame["y_pred"] - frame["y_true"]) ** 2).mean())


def compute_pearson(frame: pd.DataFrame) -> float:
    if frame.empty:
        return float("nan")
    return float(frame["y_pred"].corr(frame["y_true"], method="pearson"))


def normalize_prediction_frame(frame: pd.DataFrame) -> pd.DataFrame:
    if "y_true" not in frame.columns or "y_pred" not in frame.columns:
        raise ValueError("Prediction CSV must contain `y_true` and `y_pred` columns.")
    return frame[["y_true", "y_pred"]].apply(pd.to_numeric, errors="coerce").dropna()


def backfill_single_run(run_dir: Path) -> bool:
    predictions_path = run_dir / "test_predictions.csv"
    if not predictions_path.exists():
        return False

    frame = normalize_prediction_frame(pd.read_csv(predictions_path))
    if frame.empty:
        return False

    plot_path = build_regression_figure_path(run_dir, run_dir.name)
    save_regression_plot(
        plot_path,
        targets=frame["y_true"].tolist(),
        predictions=frame["y_pred"].tolist(),
        run_label=run_dir.name,
        mse=compute_mse(frame),
        pearson=compute_pearson(frame),
    )
    print(f"Saved single-run regression plot to {plot_path}")
    return True


def backfill_cv_run(run_dir: Path) -> bool:
    predictions_path = run_dir / "cv_test_predictions.csv"
    if not predictions_path.exists():
        return False

    frame = normalize_prediction_frame(pd.read_csv(predictions_path))
    if frame.empty:
        return False

    plot_path = build_regression_figure_path(run_dir, "cv_aggregate", aggregate=True)
    save_regression_plot(
        plot_path,
        targets=frame["y_true"].tolist(),
        predictions=frame["y_pred"].tolist(),
        run_label=run_dir.name,
        mse=compute_mse(frame),
        pearson=compute_pearson(frame),
    )
    print(f"Saved CV regression plot to {plot_path}")
    return True


def resolve_run_dirs(outputs_dir: Path, runs: list[str] | None) -> list[Path]:
    if runs:
        return [outputs_dir / run for run in runs]
    return sorted(path for path in outputs_dir.iterdir() if path.is_dir() and path.name != "training_curves")


def main() -> None:
    args = parse_args()
    outputs_dir = Path(args.outputs_dir)
    run_dirs = resolve_run_dirs(outputs_dir, args.runs)

    processed = 0
    skipped = 0

    for run_dir in run_dirs:
        matched = False

        if args.kind in {"single", "all"}:
            matched = backfill_single_run(run_dir) or matched

        if args.kind in {"cv", "all"}:
            matched = backfill_cv_run(run_dir) or matched

        if matched:
            processed += 1
        else:
            skipped += 1

    print(f"Processed {processed} run directories.")
    if skipped:
        print(f"Skipped {skipped} run directories with no matching prediction CSVs.")


if __name__ == "__main__":
    main()
