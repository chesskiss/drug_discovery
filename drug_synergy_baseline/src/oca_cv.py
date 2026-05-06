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


def _plot_cv_aggregate(aggregate_df: pd.DataFrame, output_path: Path, top_k: int) -> None:
    top_df = aggregate_df.sort_values(["mean_delta_squared_error_mean", "top_k_count"], ascending=[False, False]).head(top_k)
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.bar(
        [f"C{int(idx)}" for idx in top_df["component_idx"]],
        top_df["mean_delta_squared_error_mean"],
        yerr=top_df["mean_delta_squared_error_std"],
        capsize=4,
        color="#2c7fb8",
    )
    ax.set_title("CV OCA Aggregate: Mean Delta Squared Error Across Folds")
    ax.set_xlabel("PCA Component")
    ax.set_ylabel("Mean Delta Squared Error")
    ax.grid(axis="y", alpha=0.25)
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
            mean_delta_squared_error_mean=("mean_delta_squared_error", "mean"),
            mean_delta_squared_error_std=("mean_delta_squared_error", "std"),
            mean_delta_absolute_error_mean=("mean_delta_absolute_error", "mean"),
            mean_delta_absolute_error_std=("mean_delta_absolute_error", "std"),
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
    summary_path = oca_root / "oca_cv_summary.json"

    combined.to_csv(combined_path, index=False)
    aggregate_df.to_csv(aggregate_path, index=False)
    _plot_cv_aggregate(aggregate_df, plot_path, args.top_k)
    save_json(
        summary_path,
        {
            "cv_output_dir": str(cv_output_dir),
            "num_folds_processed": len(fold_dirs),
            "top_k": int(args.top_k),
            "mask_value": float(args.mask_value),
            "aggregate_plot": str(plot_path),
            "aggregate_csv": str(aggregate_path),
            "per_fold_csv": str(combined_path),
        },
    )

    print(f"[oca_cv] Saved per-fold component scores to {combined_path}")
    print(f"[oca_cv] Saved aggregate component summary to {aggregate_path}")
    print(f"[oca_cv] Saved aggregate stability plot to {plot_path}")
    print(f"[oca_cv] Saved summary to {summary_path}")


if __name__ == "__main__":
    main()
