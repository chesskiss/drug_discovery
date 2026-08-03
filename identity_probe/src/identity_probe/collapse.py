"""P5 — prediction collapse.

The most direct test of "did the model stop using drugs and collapse to a
per-cell-line constant?". Reads only the already-saved predictions from trained
runs: no representation matrices, no pickle, no retraining.

Key metric
----------
``identity_r2``: the between-cell-line variance share of the model's OWN output.
Formally the R^2 of the "predict every row with its cell line's mean prediction"
predictor, evaluated against ``y_pred``. 1.0 means the model emits a single
constant per cell line (fully collapsed); 0.0 means cell line explains nothing.

The same quantity computed on ``y_true`` (``identity_r2_true``) is the calibration
ceiling — how much of the *real* signal is genuinely between-cell-line. A model
whose ``identity_r2`` far exceeds ``identity_r2_true`` has over-committed to
cell-line identity.

Interpretation subtlety
-----------------------
On a cold (cell_line / drug_and_cell_line) split the test cell lines are UNSEEN.
A high ``identity_r2`` there does not mean the model successfully looked up a
memorised value; it means the model learned a gene-vector -> offset mapping and is
applying it confidently (and wrongly) to cell lines it has never seen.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd


def _variance_shares(values: np.ndarray, groups: np.ndarray) -> tuple[float, float, float]:
    """Return (identity_r2, within_group_std, overall_std) for ``values`` grouped by ``groups``.

    ``identity_r2`` is 1 - SS_within / SS_total, i.e. the fraction of variance in
    ``values`` explained by group membership alone.
    """
    overall_std = float(np.std(values))
    total_ss = float(np.sum((values - values.mean()) ** 2))
    if total_ss <= 0:
        # Degenerate: the model emitted a single global constant.
        return float("nan"), 0.0, overall_std

    frame = pd.DataFrame({"v": values, "g": groups})
    group_mean = frame.groupby("g", sort=False)["v"].transform("mean")
    within_ss = float(np.sum((frame["v"] - group_mean) ** 2))
    identity_r2 = 1.0 - within_ss / total_ss

    # Mean over groups of the within-group std (how much prediction moves when
    # only the drug pair changes).
    within_std = float(frame.groupby("g", sort=False)["v"].std(ddof=0).mean())
    return identity_r2, within_std, overall_std


def _fold_metrics(fold: pd.DataFrame) -> dict[str, float]:
    cells = fold["cell_line"].to_numpy()
    pred = fold["y_pred"].to_numpy(dtype=float)
    true = fold["y_true"].to_numpy(dtype=float)

    pred_id_r2, pred_within_std, pred_std = _variance_shares(pred, cells)
    true_id_r2, true_within_std, true_std = _variance_shares(true, cells)

    return {
        "identity_r2": pred_id_r2,
        "identity_r2_true": true_id_r2,
        "within_cell_pred_std": pred_within_std,
        "within_cell_true_std": true_within_std,
        "overall_pred_std": pred_std,
        "overall_true_std": true_std,
        "drug_sensitivity_ratio": pred_within_std / pred_std if pred_std > 0 else float("nan"),
        "pred_std_ratio": pred_std / true_std if true_std > 0 else float("nan"),
        "n_test_cell_lines": float(len(np.unique(cells))),
        "n_rows": float(len(fold)),
    }


def load_run_predictions(run_dir: Path) -> pd.DataFrame | None:
    """Load a run's predictions, preferring the aggregated CSV over per-fold files."""
    aggregate = run_dir / "cv_test_predictions.csv"
    if aggregate.exists():
        return pd.read_csv(aggregate)

    fold_files = sorted(run_dir.glob("fold_runs/*/test_predictions.csv"))
    if not fold_files:
        return None
    frames = []
    for path in fold_files:
        frame = pd.read_csv(path)
        # Recover the fold index from the directory name (cv_seed_42_fold_7 -> 7).
        frame["cv_fold"] = int(path.parent.name.rsplit("_", 1)[-1])
        frames.append(frame)
    return pd.concat(frames, ignore_index=True)


def read_observed_rmse(run_dir: Path) -> float:
    """Read aggregate test RMSE from cv_metrics.json (NaN when absent)."""
    metrics_path = run_dir / "cv_metrics.json"
    if not metrics_path.exists():
        return float("nan")
    with metrics_path.open("r", encoding="utf-8") as handle:
        data = json.load(handle)
    return float(data.get("aggregate_metrics", {}).get("test_rmse", {}).get("mean", float("nan")))


def probe_collapse(
    sweep_dir: str | Path,
    reps: list[str],
    splits: list[str],
) -> pd.DataFrame:
    """Per (representation, split): variance decomposition of predictions, mean±std over folds."""
    sweep_dir = Path(sweep_dir)
    rows: list[dict[str, object]] = []

    for split in splits:
        for rep in reps:
            run_dir = sweep_dir / f"{rep}_{split}"
            predictions = load_run_predictions(run_dir)
            if predictions is None:
                continue

            # A cold split gives each fold different test cell lines, so compute
            # per fold and aggregate rather than pooling incomparable folds.
            fold_col = "cv_fold" if "cv_fold" in predictions.columns else None
            if fold_col is None:
                per_fold = [_fold_metrics(predictions)]
            else:
                per_fold = [
                    _fold_metrics(fold) for _, fold in predictions.groupby(fold_col, sort=True)
                ]

            fold_frame = pd.DataFrame(per_fold)
            row: dict[str, object] = {
                "representation": rep,
                "split": split,
                "n_folds": len(fold_frame),
                "observed_rmse": read_observed_rmse(run_dir),
            }
            for column in fold_frame.columns:
                row[column] = float(fold_frame[column].mean())
                row[f"{column}_std"] = float(fold_frame[column].std(ddof=0))
            rows.append(row)

    return pd.DataFrame(rows)
