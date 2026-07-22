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


HIT_SYNERGY_THRESHOLD = 10.0
HIT_TOP_K_FRACTION = 0.10


def compute_top_k_hit_rate(
    targets: list[float],
    predictions: list[float],
    *,
    hit_threshold: float = HIT_SYNERGY_THRESHOLD,
    top_k_fraction: float = HIT_TOP_K_FRACTION,
) -> dict[str, float] | None:
    """Rank rows by predicted synergy, take the top `top_k_fraction`, and report how
    many are true hits (`y_true > hit_threshold`).

    Ranking-based on purpose: the baseline's predictions are compressed toward the
    mean, so a fixed-threshold classification would never fire. Enrichment is the
    top-K hit rate over the overall hit rate (1.0 == no better than random).
    """
    frame = pd.DataFrame({"y_true": targets, "y_pred": predictions}).dropna()
    if frame.empty:
        return None

    n_rows = len(frame)
    k = max(1, int(round(n_rows * top_k_fraction)))
    is_hit = frame["y_true"] > hit_threshold
    total_hits = int(is_hit.sum())
    base_rate = float(is_hit.mean())

    top_k = frame.nlargest(k, "y_pred")
    hits_in_top_k = int((top_k["y_true"] > hit_threshold).sum())
    hit_rate = hits_in_top_k / k

    return {
        "hit_threshold": hit_threshold,
        "top_k_fraction": top_k_fraction,
        "k": k,
        "n_rows": n_rows,
        "hits_in_top_k": hits_in_top_k,
        "total_hits": total_hits,
        "hit_rate": hit_rate,
        "base_rate": base_rate,
        "enrichment": hit_rate / base_rate if base_rate > 0 else float("nan"),
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

    hit_stats = compute_top_k_hit_rate(frame["y_true"].tolist(), frame["y_pred"].tolist())
    if hit_stats is not None:
        top_pct = int(round(hit_stats["top_k_fraction"] * 100))
        enrichment = hit_stats["enrichment"]
        enrichment_text = f"{enrichment:.2f}x random" if np.isfinite(enrichment) else "n/a"
        title_bits.append(
            f"Top-{top_pct}% hit rate = {hit_stats['hit_rate']:.1%} "
            f"({hit_stats['hits_in_top_k']}/{hit_stats['k']}) | "
            f"base = {hit_stats['base_rate']:.1%} | {enrichment_text}"
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

    # Visualise the top-K hit metric: true hits lie right of the vertical line,
    # top-K ranked predictions lie above the horizontal one.
    if hit_stats is not None:
        top_k_cutoff = float(frame["y_pred"].nlargest(hit_stats["k"]).min())
        plt.axvline(
            hit_stats["hit_threshold"],
            color="#d62728",
            linestyle=":",
            linewidth=1.5,
            label=f"hit: y_true > {hit_stats['hit_threshold']:g}",
        )
        plt.axhline(
            top_k_cutoff,
            color="#2ca02c",
            linestyle=":",
            linewidth=1.5,
            label=f"top-{int(round(hit_stats['top_k_fraction'] * 100))}% pred cutoff",
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
