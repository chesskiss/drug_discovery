from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def plot_component_importance_topk(component_importance: pd.DataFrame, output_path: Path, top_k: int) -> None:
    top_df = component_importance.nsmallest(top_k, "rank").sort_values("rank")
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.bar(
        [f"C{idx}" for idx in top_df["component_idx"]],
        top_df["mean_delta_squared_error"],
        color="#2c7fb8",
    )
    ax.set_title("OCA Global Importance: Top Components")
    ax.set_xlabel("PCA Component")
    ax.set_ylabel("Mean Delta Squared Error")
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def plot_component_importance_head_tail_summary(
    component_importance: pd.DataFrame,
    output_path: Path,
    *,
    head_k: int = 10,
    tail_k: int = 5,
    zero_tolerance: float = 1e-12,
) -> None:
    ranked = component_importance.sort_values("rank").reset_index(drop=True)
    head_df = ranked.head(head_k).copy()
    tail_df = ranked.tail(tail_k).copy()

    zero_mask = ranked["mean_abs_delta_prediction"].abs() <= zero_tolerance
    zero_mask &= ranked["mean_delta_squared_error"].abs() <= zero_tolerance
    zero_mask &= ranked["mean_delta_absolute_error"].abs() <= zero_tolerance
    zero_count = int(zero_mask.sum())

    labels: list[str] = []
    values: list[float] = []
    colors: list[str] = []

    for row in head_df.itertuples(index=False):
        labels.append(f"C{int(row.component_idx)}")
        values.append(float(row.mean_delta_squared_error))
        colors.append("#2c7fb8")

    if zero_count > 0:
        labels.append(f"Zero block\n({zero_count} comps)")
        values.append(0.0)
        colors.append("#bdbdbd")

    for row in tail_df.itertuples(index=False):
        labels.append(f"C{int(row.component_idx)}")
        values.append(float(row.mean_delta_squared_error))
        colors.append("#d95f02" if row.mean_delta_squared_error < 0 else "#7570b3")

    fig, ax = plt.subplots(figsize=(14, 6))
    bars = ax.bar(labels, values, color=colors)
    ax.axhline(0.0, color="black", linewidth=1.0, alpha=0.6)
    ax.set_title("OCA Component Importance Summary: Helpful Head, Zero Block, Harmful Tail")
    ax.set_xlabel("PCA Component / Aggregated Block")
    ax.set_ylabel("Mean Delta Squared Error")
    ax.grid(axis="y", alpha=0.25)
    ax.tick_params(axis="x", rotation=35)

    for bar, value in zip(bars, values, strict=False):
        if value == 0.0:
            continue
        va = "bottom" if value > 0 else "top"
        offset = 0.003 if value > 0 else -0.003
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            value + offset,
            f"{value:.3f}",
            ha="center",
            va=va,
            fontsize=8,
        )

    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
