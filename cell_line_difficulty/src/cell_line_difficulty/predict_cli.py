from __future__ import annotations

import argparse
from pathlib import Path

from .predictive import (
    compare_feature_views,
    load_predictive_dataset,
    run_loocv,
    save_predictive_outputs,
    save_view_comparison_outputs,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Predict cell-line synergy ease from cell-line features")
    parser.add_argument(
        "--synergy-path",
        type=str,
        default="drug_synergy_baseline/data/drugcomb.csv",
        help="Path to DrugComb synergy CSV",
    )
    parser.add_argument(
        "--pickle-path",
        type=str,
        default="data/data_compression/source_data/drugcomb.pkl",
        help="Path to DrugComb pickle with full CellLine payloads",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="cell_line_difficulty/outputs",
        help="Directory for predictive outputs",
    )
    parser.add_argument(
        "--max-cell-lines",
        type=int,
        default=None,
        help="Optional cap for smoke runs; default uses all cell lines",
    )
    parser.add_argument(
        "--feature-view-index",
        type=int,
        default=1,
        help="Which CellLine view to use: 0=23808, 1=3171, 2=627",
    )
    parser.add_argument(
        "--models",
        type=str,
        nargs="+",
        default=["ridge", "mlp"],
        help="Models to run: ridge and/or mlp",
    )
    parser.add_argument(
        "--compare-views",
        action="store_true",
        help="Compare CellLine views 0/1/2 on the same task and save summary metrics",
    )
    parser.add_argument(
        "--progress",
        action="store_true",
        help="Print LOOCV progress",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    models = tuple(args.models)

    if args.compare_views:
        comparison = compare_feature_views(
            synergy_path=args.synergy_path,
            pickle_path=args.pickle_path,
            view_indices=(0, 1, 2),
            models=models,
            max_cell_lines=args.max_cell_lines,
            progress=args.progress,
        )
        save_view_comparison_outputs(comparison, args.output_dir)
        print("Saved feature-view comparison to", Path(args.output_dir) / "view_comparison.csv")
        print(comparison.to_string(index=False))
        return

    dataset = load_predictive_dataset(
        args.synergy_path,
        args.pickle_path,
        feature_view_index=args.feature_view_index,
        max_cell_lines=args.max_cell_lines,
    )
    predictions, metrics = run_loocv(dataset, models=models, progress=args.progress)
    save_predictive_outputs(dataset, predictions, metrics, args.output_dir)

    print(f"Built predictive dataset with {len(dataset)} cell lines and {metrics['feature_dimension']} features")
    if "ridge" in metrics["models"]:
        print(f"Ridge LOOCV RMSE: {metrics['models']['ridge']['rmse']:.4f}")
    if "mlp" in metrics["models"]:
        print(f"MLP LOOCV RMSE: {metrics['models']['mlp']['rmse']:.4f}")
    print(f"Fold-mean baseline RMSE: {metrics['models']['fold_mean_baseline']['rmse']:.4f}")
    print(f"Saved predictive outputs to {Path(args.output_dir)}")


if __name__ == "__main__":
    main()
