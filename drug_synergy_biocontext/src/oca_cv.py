from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from .oca import build_test_bundle_from_config, run_oca_analysis
from .predict import get_device, load_config, load_model
from .training_artifacts import save_json


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run fold-wise OCA for a saved cross-validation output directory and aggregate results."
    )
    parser.add_argument("--cv-output-dir", type=str, required=True, help="Path to the CV run directory.")
    parser.add_argument("--top-k", type=int, default=10, help="Top-k components used in fold plots and aggregation.")
    parser.add_argument("--local-top-n", type=int, default=5, help="Local explanation selection size per fold.")
    parser.add_argument("--mask-value", type=float, default=0.0, help="Replacement value for masked components.")
    parser.add_argument("--batch-size", type=int, default=None, help="Optional batch size override.")
    parser.add_argument("--device", type=str, default=None, help="Optional device override.")
    return parser.parse_args()


def _build_head_zero_tail_summary(aggregate_df: pd.DataFrame, *, head_k: int, tail_k: int) -> pd.DataFrame:
    ranked = aggregate_df.sort_values(["aggregate_rank", "component_idx"]).reset_index(drop=True)
    zero_mask = (
        ranked["mean_abs_delta_prediction_mean"].abs() <= 1e-12
    ) & (
        ranked["mean_delta_squared_error_mean"].abs() <= 1e-12
    ) & (
        ranked["mean_delta_absolute_error_mean"].abs() <= 1e-12
    ) & (
        ranked["mean_delta_squared_error_min"].abs() <= 1e-12
    ) & (
        ranked["mean_delta_squared_error_max"].abs() <= 1e-12
    )

    head_df = ranked.head(head_k).copy()
    tail_df = ranked.tail(tail_k).copy()
    zero_count = int(zero_mask.sum())

    rows: list[dict[str, object]] = []
    for row in head_df.itertuples(index=False):
        rows.append(
            {
                "label": f"C{int(row.component_idx)}",
                "component_idx": int(row.component_idx),
                "mean": float(row.mean_delta_squared_error_mean),
                "min": float(row.mean_delta_squared_error_min),
                "max": float(row.mean_delta_squared_error_max),
                "kind": "head",
            }
        )

    if zero_count > 0:
        rows.append(
            {
                "label": f"Zero block\n({zero_count} comps)",
                "component_idx": None,
                "mean": 0.0,
                "min": 0.0,
                "max": 0.0,
                "kind": "zero",
            }
        )

    for row in tail_df.itertuples(index=False):
        rows.append(
            {
                "label": f"C{int(row.component_idx)}",
                "component_idx": int(row.component_idx),
                "mean": float(row.mean_delta_squared_error_mean),
                "min": float(row.mean_delta_squared_error_min),
                "max": float(row.mean_delta_squared_error_max),
                "kind": "tail",
            }
        )

    return pd.DataFrame(rows)


def _plot_cv_aggregate_head_tail(aggregate_df: pd.DataFrame, output_path: Path, top_k: int) -> None:
    summary_df = _build_head_zero_tail_summary(aggregate_df, head_k=top_k, tail_k=5)
    fig, ax = plt.subplots(figsize=(15, 6))
    x_positions = list(range(len(summary_df)))

    for idx, row in enumerate(summary_df.itertuples(index=False)):
        color = "#2c7fb8" if row.kind == "head" else "#d95f02" if row.kind == "tail" else "#9e9e9e"
        ax.vlines(idx, row.min, row.max, color="black", linewidth=2.2, alpha=0.9)
        ax.scatter(idx, row.mean, color=color, s=90, zorder=3)

    ax.axhline(0.0, color="black", linewidth=1.0, alpha=0.6)
    ax.set_xticks(x_positions)
    ax.set_xticklabels(summary_df["label"].tolist(), rotation=35, ha="right")
    ax.set_title("CV OCA Aggregate: Helpful Head, Zero Block, Harmful Tail")
    ax.set_xlabel("Pathway")
    ax.set_ylabel("Mean Delta Squared Error")
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def _plot_cv_topk_frequency(aggregate_df: pd.DataFrame, output_path: Path, top_k: int) -> None:
    top_df = aggregate_df.sort_values(["top_k_count", "mean_delta_squared_error_mean"], ascending=[False, False]).head(top_k)
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.bar([f"C{int(idx)}" for idx in top_df["component_idx"]], top_df["top_k_count"], color="#2c7fb8")
    ax.set_title("CV OCA: Top-k Frequency Across Folds")
    ax.set_xlabel("Pathway")
    ax.set_ylabel("Number of Folds in Top-k")
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def _plot_cv_fold_component_heatmap(combined_df: pd.DataFrame, aggregate_df: pd.DataFrame, output_path: Path, top_k: int) -> None:
    top_components = (
        aggregate_df.sort_values(["mean_delta_squared_error_mean", "top_k_count"], ascending=[False, False])
        .head(top_k)["component_idx"]
        .tolist()
    )
    heatmap_df = combined_df[combined_df["component_idx"].isin(top_components)].copy()
    pivot = heatmap_df.pivot_table(
        index="cv_fold",
        columns="component_idx",
        values="mean_delta_squared_error",
        aggfunc="mean",
    ).reindex(columns=top_components)

    fig, ax = plt.subplots(figsize=(max(9, len(top_components) * 0.9), max(5, len(pivot) * 0.6)))
    image = ax.imshow(pivot.to_numpy(dtype=float), aspect="auto", cmap="coolwarm")
    ax.set_title("CV OCA Fold-by-Component Heatmap")
    ax.set_xlabel("Pathway")
    ax.set_ylabel("CV Fold")
    ax.set_xticks(range(len(pivot.columns)))
    ax.set_xticklabels([f"C{int(col)}" for col in pivot.columns], rotation=45, ha="right")
    ax.set_yticks(range(len(pivot.index)))
    ax.set_yticklabels([f"Fold {int(idx)}" for idx in pivot.index])
    fig.colorbar(image, ax=ax, label="Mean Delta Squared Error")
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    cv_output_dir = Path(args.cv_output_dir)
    fold_root = cv_output_dir / "fold_runs"
    if not fold_root.exists():
        raise FileNotFoundError(f"Could not find fold_runs directory under {cv_output_dir}")

    fold_dirs = sorted(path for path in fold_root.iterdir() if path.is_dir())
    if not fold_dirs:
        raise ValueError(f"No fold run directories found under {fold_root}")

    device = get_device(args.device)
    aggregate_frames: list[pd.DataFrame] = []

    for fold_dir in fold_dirs:
        model_path = fold_dir / "baseline_mlp.pt"
        config_path = fold_dir / "config.json"
        if not model_path.exists():
            raise FileNotFoundError(
                f"Missing fold checkpoint: {model_path}. "
                "Rerun CV after the checkpoint-saving patch so fold-wise OCA has models to load."
            )
        if not config_path.exists():
            raise FileNotFoundError(f"Missing fold config: {config_path}")

        config = load_config(config_path)
        if config.get("evaluation_mode") != "cross_validation":
            raise ValueError(f"{config_path} is not a CV fold config.")

        datasets, resolved = build_test_bundle_from_config(config)
        batch_size = int(args.batch_size or config.get("batch_size", 256))
        model = load_model(config, model_path, device)
        output_dir = fold_dir / "oca"

        component_df, _ = run_oca_analysis(
            model=model,
            model_path=model_path,
            config_path=config_path,
            config=config,
            datasets=datasets,
            resolved=resolved,
            output_dir=output_dir,
            batch_size=batch_size,
            device=device,
            top_k=args.top_k,
            local_top_n=args.local_top_n,
            local_row_idx=None,
            mask_value=args.mask_value,
            progress_prefix=fold_dir.name,
        )

        component_df = component_df.copy()
        component_df["fold_run"] = fold_dir.name
        component_df["cv_seed"] = int(config["cv_seed"])
        component_df["cv_fold"] = int(config["cv_fold"])
        component_df["is_top_k"] = component_df["rank"] <= args.top_k
        component_df["is_helpful"] = component_df["mean_delta_squared_error"] > 0
        component_df["is_harmful"] = component_df["mean_delta_squared_error"] < 0
        aggregate_frames.append(component_df)

    combined = pd.concat(aggregate_frames, ignore_index=True)
    aggregate_df = (
        combined.groupby("component_idx", as_index=False)
        .agg(
            mean_abs_delta_prediction_mean=("mean_abs_delta_prediction", "mean"),
            mean_abs_delta_prediction_std=("mean_abs_delta_prediction", "std"),
            mean_abs_delta_prediction_min=("mean_abs_delta_prediction", "min"),
            mean_abs_delta_prediction_max=("mean_abs_delta_prediction", "max"),
            mean_delta_squared_error_mean=("mean_delta_squared_error", "mean"),
            mean_delta_squared_error_std=("mean_delta_squared_error", "std"),
            mean_delta_squared_error_min=("mean_delta_squared_error", "min"),
            mean_delta_squared_error_max=("mean_delta_squared_error", "max"),
            mean_delta_absolute_error_mean=("mean_delta_absolute_error", "mean"),
            mean_delta_absolute_error_std=("mean_delta_absolute_error", "std"),
            mean_delta_absolute_error_min=("mean_delta_absolute_error", "min"),
            mean_delta_absolute_error_max=("mean_delta_absolute_error", "max"),
            mean_rank=("rank", "mean"),
            rank_std=("rank", "std"),
            top_k_count=("is_top_k", "sum"),
            helpful_fold_count=("is_helpful", "sum"),
            harmful_fold_count=("is_harmful", "sum"),
            fold_count=("cv_fold", "count"),
        )
        .sort_values(["mean_delta_squared_error_mean", "top_k_count", "component_idx"], ascending=[False, False, True])
        .reset_index(drop=True)
    )

    aggregate_df["aggregate_rank"] = range(1, len(aggregate_df) + 1)

    oca_root = cv_output_dir / "oca_cv"
    oca_root.mkdir(parents=True, exist_ok=True)
    combined_path = oca_root / "oca_cv_component_importance_per_fold.csv"
    aggregate_path = oca_root / "oca_cv_component_summary.csv"
    plot_path = oca_root / "oca_cv_topk_stability.png"
    heatmap_path = oca_root / "oca_cv_fold_component_heatmap.png"
    frequency_path = oca_root / "oca_cv_topk_frequency.png"
    summary_path = oca_root / "oca_cv_summary.json"

    combined.to_csv(combined_path, index=False)
    aggregate_df.to_csv(aggregate_path, index=False)
    _plot_cv_aggregate_head_tail(aggregate_df, plot_path, args.top_k)
    _plot_cv_fold_component_heatmap(combined, aggregate_df, heatmap_path, args.top_k)
    _plot_cv_topk_frequency(aggregate_df, frequency_path, args.top_k)
    save_json(
        summary_path,
        {
            "cv_output_dir": str(cv_output_dir),
            "num_folds_processed": len(fold_dirs),
            "top_k": int(args.top_k),
            "mask_value": float(args.mask_value),
            "aggregate_plot": str(plot_path),
            "heatmap_plot": str(heatmap_path),
            "topk_frequency_plot": str(frequency_path),
            "aggregate_csv": str(aggregate_path),
            "per_fold_csv": str(combined_path),
        },
    )

    print(f"[oca_cv] Saved per-fold component scores to {combined_path}")
    print(f"[oca_cv] Saved aggregate component summary to {aggregate_path}")
    print(f"[oca_cv] Saved aggregate stability plot to {plot_path}")
    print(f"[oca_cv] Saved fold/component heatmap to {heatmap_path}")
    print(f"[oca_cv] Saved top-k frequency plot to {frequency_path}")
    print(f"[oca_cv] Saved summary to {summary_path}")


if __name__ == "__main__":
    main()
