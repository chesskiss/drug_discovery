from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def find_outputs_root(output_dir: str | Path) -> Path:
    output_path = Path(output_dir)
    if output_path.name == "outputs":
        return output_path

    for parent in [output_path, *output_path.parents]:
        if parent.name == "outputs":
            return parent

    if output_path.is_absolute():
        return output_path.parent

    if output_path.parts:
        return Path(output_path.parts[0])

    return output_path


def build_curve_figure_path(output_dir: str | Path, run_label: str) -> Path:
    output_path = Path(output_dir)
    outputs_root = find_outputs_root(output_path)
    relative_run_dir = output_path.relative_to(outputs_root)

    figure_dir = outputs_root / "training_curves"
    if str(relative_run_dir) == ".":
        relative_curve_dir = Path(run_label)
    else:
        relative_curve_dir = relative_run_dir

    return figure_dir / relative_curve_dir / "loss_curve.png"


def build_regression_figure_path(output_dir: str | Path, run_label: str, *, aggregate: bool = False) -> Path:
    output_path = Path(output_dir)
    outputs_root = find_outputs_root(output_path)
    relative_run_dir = output_path.relative_to(outputs_root)

    figure_dir = outputs_root / "training_curves"
    if str(relative_run_dir) == ".":
        relative_curve_dir = Path(run_label)
    else:
        relative_curve_dir = relative_run_dir

    filename = "cv_test_regression.png" if aggregate else "test_regression.png"
    return figure_dir / relative_curve_dir / filename


def save_json(path: str | Path, payload: dict[str, Any]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2)


def save_history_csv(path: str | Path, history: list[dict[str, Any]]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(history).to_csv(path, index=False)


def save_loss_curve(path: str | Path, history: list[dict[str, Any]], run_label: str) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    history_frame = pd.DataFrame(history)
    if history_frame.empty:
        return

    plt.figure(figsize=(8, 5))
    plt.plot(history_frame["epoch"], history_frame["train_mse"], label="train_mse", linewidth=2)
    plt.plot(history_frame["epoch"], history_frame["val_mse"], label="val_mse", linewidth=2)
    plt.xlabel("Epoch")
    plt.ylabel("MSE")
    plt.title(f"Train vs Val Loss: {run_label}")
    plt.legend()
    plt.tight_layout()
    plt.savefig(path, dpi=160)
    plt.close()


# Synergy_ZIP bucket edges, anchored at the biological zero (ZIP == 0 means "no
# interaction / additive"; positive == synergy, negative == antagonism).
#
# Derivation (objective, reproducible — see notebook analysis):
#   * GMM/BIC + KDE show the label distribution is UNIMODAL with fat tails, so there is
#     no data-intrinsic valley to split on; a threshold must be defined against the
#     "no interaction" null (ZIP == 0), not discovered from modality.
#   * The additive core's dispersion is estimated robustly with the MAD (resistant to
#     the synergy/antagonism tails): sigma_robust = 1.4826 * MAD ~= 4.04.
#   * Bucket edges are significance bounds of that null: |z| * sigma_robust.
#       mild   = 1.960 * sigma_robust ~= 7.9   (95% two-sided bound)
#       strong = 2.576 * sigma_robust ~= 10.4  (99% bound; ~= SynergyFinder's ZIP>10)
# Fixed constants (not recomputed per run) so buckets stay comparable across folds,
# runs, and models.
SYNERGY_ADDITIVE_BAND = 7.9
SYNERGY_STRONG_THRESHOLD = 10.4

# Ordinal synergy classes, from strongest antagonism (-2) to strongest synergy (+2).
SYNERGY_BUCKET_LABELS = {
    -2: "strong antagonism",
    -1: "mild antagonism",
    0: "additive",
    1: "mild synergy",
    2: "strong synergy",
}


def assign_synergy_bucket(
    values: "np.ndarray | list[float]",
    *,
    additive_band: float = SYNERGY_ADDITIVE_BAND,
    strong_threshold: float = SYNERGY_STRONG_THRESHOLD,
) -> np.ndarray:
    """Map continuous synergy scores onto the 5 ordinal classes in ``SYNERGY_BUCKET_LABELS``."""
    v = np.asarray(values, dtype=float)
    bucket = np.zeros(len(v), dtype=int)
    bucket[v > additive_band] = 1
    bucket[v > strong_threshold] = 2
    bucket[v < -additive_band] = -1
    bucket[v < -strong_threshold] = -2
    return bucket


def compute_synergy_bucket_metrics(
    targets: list[float],
    predictions: list[float],
    *,
    additive_band: float = SYNERGY_ADDITIVE_BAND,
    strong_threshold: float = SYNERGY_STRONG_THRESHOLD,
) -> dict[str, Any] | None:
    """Discretise true and predicted synergy into 5 ordinal classes and score agreement.

    Reports exact-bucket accuracy (predicted class == true class) and within-1-bucket
    accuracy (off by at most one ordinal step). Per-bucket, it also records the mean
    predicted score, which reveals ordinal signal even when the regression is
    miscalibrated (a well-ordered model has monotonically increasing mean prediction
    from the antagonism buckets to the synergy buckets).
    """
    frame = pd.DataFrame({"y_true": targets, "y_pred": predictions}).dropna()
    if frame.empty:
        return None

    true_bucket = assign_synergy_bucket(
        frame["y_true"].to_numpy(), additive_band=additive_band, strong_threshold=strong_threshold
    )
    pred_bucket = assign_synergy_bucket(
        frame["y_pred"].to_numpy(), additive_band=additive_band, strong_threshold=strong_threshold
    )

    exact_accuracy = float(np.mean(true_bucket == pred_bucket))
    within_one_accuracy = float(np.mean(np.abs(true_bucket - pred_bucket) <= 1))

    per_bucket: dict[int, dict[str, float]] = {}
    for code in SYNERGY_BUCKET_LABELS:
        mask = true_bucket == code
        count = int(mask.sum())
        per_bucket[code] = {
            "label": SYNERGY_BUCKET_LABELS[code],
            "count": count,
            "fraction": float(mask.mean()),
            "recall": float(np.mean(pred_bucket[mask] == code)) if count else float("nan"),
            "mean_pred": float(frame["y_pred"].to_numpy()[mask].mean()) if count else float("nan"),
        }

    return {
        "additive_band": additive_band,
        "strong_threshold": strong_threshold,
        "n_rows": len(frame),
        "exact_accuracy": exact_accuracy,
        "within_one_accuracy": within_one_accuracy,
        "per_bucket": per_bucket,
    }


def save_regression_plot(
    path: str | Path,
    *,
    targets: list[float],
    predictions: list[float],
    run_label: str,
    mse: float | None = None,
    pearson: float | None = None,
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    if not targets or not predictions:
        return

    frame = pd.DataFrame({"y_true": targets, "y_pred": predictions}).dropna()
    if frame.empty:
        return

    x = frame["y_true"].to_numpy(dtype=float)
    y = frame["y_pred"].to_numpy(dtype=float)
    fit_annotation: list[str] = []
    if len(frame) >= 2 and np.nanstd(x) > 0:
        slope, intercept = np.polyfit(x, y, deg=1)
        fitted_y = slope * x + intercept
        residual_sum_squares = float(np.sum((y - fitted_y) ** 2))
        total_sum_squares = float(np.sum((y - np.mean(y)) ** 2))
        r_squared = 1.0 - residual_sum_squares / total_sum_squares if total_sum_squares > 0 else float("nan")
        if np.isfinite(r_squared):
            fit_annotation.append(f"$R^2$ = {r_squared:.4f}")
        if np.isfinite(slope) and np.isfinite(intercept):
            sign = "+" if intercept >= 0 else "-"
            fit_annotation.append(f"$\\hat{{y}}$ = {slope:.3f}x {sign} {abs(intercept):.3f}")

    combined_min = float(min(frame["y_true"].min(), frame["y_pred"].min()))
    combined_max = float(max(frame["y_true"].max(), frame["y_pred"].max()))
    padding = max((combined_max - combined_min) * 0.05, 1e-6)
    line_min = combined_min - padding
    line_max = combined_max + padding

    title_bits = [f"Test Regression: {run_label}"]
    stats_bits: list[str] = []
    if mse is not None:
        stats_bits.append(f"MSE={mse:.4f}")
    if pearson is not None and pd.notna(pearson):
        stats_bits.append(f"Pearson={pearson:.4f}")
    if stats_bits:
        title_bits.append(" | ".join(stats_bits))

    bucket_stats = compute_synergy_bucket_metrics(frame["y_true"].tolist(), frame["y_pred"].tolist())
    if bucket_stats is not None:
        title_bits.append(
            f"5-bucket acc = {bucket_stats['exact_accuracy']:.1%} | "
            f"within ±1 = {bucket_stats['within_one_accuracy']:.1%} "
            f"(bands ±{bucket_stats['additive_band']:g} / ±{bucket_stats['strong_threshold']:g})"
        )

    plt.figure(figsize=(6.5, 6.5))
    plt.scatter(
        frame["y_true"],
        frame["y_pred"],
        s=10,
        alpha=0.35,
        edgecolors="none",
    )
    plt.plot([line_min, line_max], [line_min, line_max], linestyle="--", linewidth=2, color="black")

    # Visualise the 5 ordinal synergy classes: the additive band (+/- additive_band)
    # and the strong tails (+/- strong_threshold), applied to BOTH axes so the
    # matching diagonal cells are the exact-bucket agreements.
    if bucket_stats is not None:
        band = bucket_stats["additive_band"]
        strong = bucket_stats["strong_threshold"]
        for idx, edge in enumerate((-strong, -band, band, strong)):
            zone_label = "synergy bucket edges" if idx == 0 else None
            plt.axvline(edge, color="#9467bd", linestyle=":", linewidth=1.1, alpha=0.7, label=zone_label)
            plt.axhline(edge, color="#9467bd", linestyle=":", linewidth=1.1, alpha=0.7)
        # Name the 5 true-synergy zones along the top of the plot. The middle bands are
        # narrow, so stagger the labels across two heights to avoid overlap.
        edges = [line_min, -strong, -band, band, strong, line_max]
        names = ["str.antag", "mild antag", "additive", "mild syn", "str.syn"]
        span = line_max - line_min
        for idx, (lo, hi, name) in enumerate(zip(edges[:-1], edges[1:], names)):
            y_offset = 0.03 if idx % 2 == 0 else 0.075
            plt.text(
                (lo + hi) / 2,
                line_max - span * y_offset,
                name,
                ha="center",
                va="top",
                fontsize=7,
                color="#6a3d9a",
                alpha=0.85,
            )
        plt.legend(loc="lower right", fontsize=8, framealpha=0.85)
    if fit_annotation:
        plt.text(
            0.97,
            0.97,
            "\n".join(fit_annotation[:2]),
            transform=plt.gca().transAxes,
            ha="right",
            va="top",
            fontsize=10,
            bbox={"boxstyle": "round,pad=0.3", "facecolor": "white", "edgecolor": "0.7", "alpha": 0.85},
        )
    plt.xlim(line_min, line_max)
    plt.ylim(line_min, line_max)
    plt.xlabel("True Synergy")
    plt.ylabel("Predicted Synergy")
    plt.title("\n".join(title_bits))
    plt.tight_layout()
    plt.savefig(path, dpi=160)
    plt.close()
