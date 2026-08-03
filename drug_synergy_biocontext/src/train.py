from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from torch import nn
from torch.utils.data import DataLoader

from .dataset import (
    build_dataset_bundle_from_rows,
    build_cv_dataset_bundle,
    build_datasets,
    build_row_cv_folds,
    load_aligned_synergy_and_expression,
)
from .macros import DEFAULT_MACRO_FILE, DEFAULT_MACRO_PRESET, load_macro_preset
from .model import build_biocontext_model
from .training_artifacts import (
    build_curve_figure_path,
    build_regression_figure_path,
    save_history_csv,
    save_json,
    save_loss_curve,
    save_regression_plot,
)


BIO_CONTEXT_CHOICES = ["progeny", "kegg", "progeny_kegg"]


def parse_args() -> argparse.Namespace:
    macro_parser = argparse.ArgumentParser(add_help=False)
    macro_parser.add_argument("--macro-file", type=str, default=str(DEFAULT_MACRO_FILE))
    macro_parser.add_argument("--macro-preset", type=str, default=DEFAULT_MACRO_PRESET)
    macro_args, _ = macro_parser.parse_known_args()
    macro_values = load_macro_preset(
        preset=macro_args.macro_preset,
        macro_file=macro_args.macro_file,
    )

    parser = argparse.ArgumentParser(
        description="Train and evaluate a minimal DeepSynergy-style baseline on DrugComb/TDC data",
        parents=[macro_parser],
    )
    parser.add_argument(
        "--synergy-path",
        type=str,
        default="../data/data_compression/source_data/drugcomb.pkl",
        help="Synergy table. The raw pickle carries both the synergy rows and CellLine[0].",
    )
    parser.add_argument(
        "--fallback-pickle-path",
        type=str,
        default="../data/data_compression/source_data/drugcomb.pkl",
        help="Source of the raw 23808-dim CellLine[0] vectors the bio-context matrix projects.",
    )
    parser.add_argument("--output-dir", type=str, default="outputs")
    parser.add_argument("--epochs", type=int, default=int(macro_values["epochs"]))
    parser.add_argument("--batch-size", type=int, default=int(macro_values["batch_size"]))
    parser.add_argument("--lr", type=float, default=float(macro_values["lr"]))
    parser.add_argument("--smiles-dim", type=int, default=int(macro_values["smiles_dim"]))
    parser.add_argument("--hidden-dims", type=int, nargs="+", default=list(macro_values["hidden_dims"]))
    parser.add_argument("--dropout", type=float, default=float(macro_values["dropout"]))
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--max-samples", type=int, default=None)
    parser.add_argument("--train-fraction", type=float, default=float(macro_values["train_fraction"]))
    parser.add_argument("--val-fraction", type=float, default=float(macro_values["val_fraction"]))
    parser.add_argument(
        "--split-strategy",
        choices=["random", "cell_line", "drug", "drug_pair", "drug_and_cell_line"],
        default="random",
    )
    parser.add_argument("--use-gene-expression", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument(
        "--bio-context",
        choices=BIO_CONTEXT_CHOICES,
        default="progeny",
        help="Fixed pathway matrix projecting CellLine[0]: progeny=14, kegg=336, progeny_kegg=350 dims.",
    )
    parser.add_argument(
        "--mlp0-out-dim",
        type=int,
        default=8,
        help="Output width of MLP0, the trainable pathway compression head. 0 disables the cell-line branch.",
    )
    parser.add_argument(
        "--normalize-pathways",
        choices=["zscore", "none"],
        default="zscore",
        help="Z-score pathway activities using TRAIN cell lines only. Raw activities are ~1e4x the drug features.",
    )
    parser.add_argument("--cv-folds", type=int, default=1)
    parser.add_argument("--cv-seeds", type=int, nargs="*", default=None)
    parser.add_argument("--stratified-cv", action="store_true")
    parser.add_argument(
        "--holdout-test-fraction",
        type=float,
        default=0.0,
        help="Group fraction reserved as a global test set shared by all folds. 0 disables (default).",
    )
    parser.add_argument(
        "--holdout-seed",
        type=int,
        default=42,
        help="Seed for the global holdout carve. Kept separate from --cv-seeds so the holdout is identical across seeds.",
    )
    parser.add_argument(
        "--holdout-mode",
        choices=["instead", "additional"],
        default="instead",
        help=(
            "instead: the shared holdout IS each fold's test set (fold buckets become validation only). "
            "additional: keep each fold's own test set and score the shared holdout as an extra metric."
        ),
    )
    return parser.parse_args()


def resolve_mlp0_out_dim(args: argparse.Namespace) -> int:
    """Effective MLP0 width. Disabling gene expression forces the branch off."""
    if not args.use_gene_expression:
        return 0
    return max(0, int(args.mlp0_out_dim))


def make_loader(dataset, batch_size: int, shuffle: bool) -> DataLoader:
    return DataLoader(dataset, batch_size=batch_size, shuffle=shuffle, num_workers=0)


def get_device() -> torch.device:
    if torch.backends.mps.is_available():
        return torch.device("mps")
    if torch.cuda.is_available():
        return torch.device("cuda")
    return torch.device("cpu")


def run_epoch(model, loader, optimizer, criterion, device: torch.device, train: bool) -> float:
    model.train(mode=train)
    total_loss = 0.0
    total_examples = 0

    for batch in loader:
        drug_a = batch["drug_a"].to(device)
        drug_b = batch["drug_b"].to(device)
        gene_expr = batch["gene_expr"].to(device)
        target = batch["target"].to(device)

        with torch.set_grad_enabled(train):
            prediction = model(drug_a, drug_b, gene_expr)
            loss = criterion(prediction, target)

        if train:
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()

        batch_size = target.shape[0]
        total_loss += loss.item() * batch_size
        total_examples += batch_size

    return total_loss / max(total_examples, 1)


def predict_loader(model, loader, device: torch.device) -> tuple[list[float], list[float]]:
    model.eval()
    predictions: list[float] = []
    targets: list[float] = []

    with torch.no_grad():
        for batch in loader:
            drug_a = batch["drug_a"].to(device)
            drug_b = batch["drug_b"].to(device)
            gene_expr = batch["gene_expr"].to(device)
            target = batch["target"].to(device)

            prediction = model(drug_a, drug_b, gene_expr)
            predictions.extend(prediction.detach().cpu().tolist())
            targets.extend(target.detach().cpu().tolist())

    return predictions, targets


def compute_mse(predictions: list[float], targets: list[float]) -> float:
    if not targets:
        return float("nan")
    return float(sum((pred - target) ** 2 for pred, target in zip(predictions, targets, strict=False)) / len(targets))


def compute_mean_baseline_mse(targets: list[float], mean_target: float) -> float:
    if not targets:
        return float("nan")
    return float(sum((mean_target - target) ** 2 for target in targets) / len(targets))


def compute_correlations(predictions: list[float], targets: list[float]) -> tuple[float, float]:
    frame = pd.DataFrame({"pred": predictions, "target": targets})
    pearson = float(frame["pred"].corr(frame["target"], method="pearson"))
    spearman = float(frame["pred"].corr(frame["target"], method="spearman"))
    return pearson, spearman


def train_once(
    datasets,
    args: argparse.Namespace,
    *,
    run_seed: int,
    output_dir: Path | None = None,
    save_artifacts: bool = False,
    run_label: str = "single_run",
    extra_config: dict[str, object] | None = None,
) -> tuple[dict[str, object], dict[str, object]]:
    torch.manual_seed(run_seed)

    device = get_device()
    mlp0_out_dim = resolve_mlp0_out_dim(args)
    model = build_biocontext_model(
        drug_dim=datasets.drug_dim,
        pathway_dim=datasets.gene_dim,
        mlp0_out_dim=mlp0_out_dim,
        hidden_dims=tuple(args.hidden_dims),
        dropout=args.dropout,
    ).to(device)

    optimizer = torch.optim.Adam(model.parameters(), lr=args.lr)
    criterion = nn.MSELoss()

    train_loader = make_loader(datasets.train, batch_size=args.batch_size, shuffle=True)
    val_loader = make_loader(datasets.val, batch_size=args.batch_size, shuffle=False)
    test_loader = make_loader(datasets.test, batch_size=args.batch_size, shuffle=False)

    print(f"[{run_label}] Train samples: {len(datasets.train)}")
    print(f"[{run_label}] Val samples: {len(datasets.val)}")
    print(f"[{run_label}] Test samples: {len(datasets.test)}")
    print(f"[{run_label}] Drug feature dim: {datasets.drug_dim}")
    print(f"[{run_label}] Pathway dim: {datasets.gene_dim} | MLP0 out dim: {mlp0_out_dim}")
    print(f"[{run_label}] Device: {device}")

    history: list[dict[str, float]] = []
    for epoch in range(1, args.epochs + 1):
        train_loss = run_epoch(model, train_loader, optimizer, criterion, device, train=True)
        val_loss = run_epoch(model, val_loader, optimizer, criterion, device, train=False)
        history.append({"epoch": epoch, "train_mse": train_loss, "val_mse": val_loss})
        print(f"[{run_label}] Epoch {epoch:03d} | train_mse={train_loss:.6f} | val_mse={val_loss:.6f}")

    test_loss = run_epoch(model, test_loader, optimizer, criterion, device, train=False)
    print(f"[{run_label}] Test MSE: {test_loss:.6f}")

    train_mean_target = float(datasets.train_rows["target"].mean())
    val_predictions, val_targets = predict_loader(model, val_loader, device)
    test_predictions, test_targets = predict_loader(model, test_loader, device)
    val_mse = compute_mse(val_predictions, val_targets)
    test_mse = compute_mse(test_predictions, test_targets)
    val_baseline_mse = compute_mean_baseline_mse(val_targets, train_mean_target)
    test_baseline_mse = compute_mean_baseline_mse(test_targets, train_mean_target)
    val_pearson, val_spearman = compute_correlations(val_predictions, val_targets)
    test_pearson, test_spearman = compute_correlations(test_predictions, test_targets)

    metrics: dict[str, object] = {
        "run_label": run_label,
        "seed": run_seed,
        "split_strategy": args.split_strategy,
        "use_gene_expression": args.use_gene_expression,
        "bio_context": args.bio_context if args.use_gene_expression else None,
        "pathway_dim": datasets.gene_dim,
        "mlp0_out_dim": mlp0_out_dim,
        "normalize_pathways": args.normalize_pathways,
        "train_samples": len(datasets.train),
        "val_samples": len(datasets.val),
        "test_samples": len(datasets.test),
        "drug_dim": datasets.drug_dim,
        "gene_dim": datasets.gene_dim,
        "train_target_mean": train_mean_target,
        "val_mse": val_mse,
        "test_mse": test_mse,
        "val_rmse": math.sqrt(val_mse),
        "test_rmse": math.sqrt(test_mse),
        "val_mean_baseline_mse": val_baseline_mse,
        "test_mean_baseline_mse": test_baseline_mse,
        "val_pearson": val_pearson,
        "val_spearman": val_spearman,
        "test_pearson": test_pearson,
        "test_spearman": test_spearman,
        "history": history,
    }

    eval_outputs = {
        "val_predictions": val_predictions,
        "val_targets": val_targets,
        "test_predictions": test_predictions,
        "test_targets": test_targets,
    }

    holdout_dataset = getattr(datasets, "holdout", None)
    if holdout_dataset is not None:
        holdout_loader = make_loader(holdout_dataset, batch_size=args.batch_size, shuffle=False)
        holdout_predictions, holdout_targets = predict_loader(model, holdout_loader, device)
        holdout_mse = compute_mse(holdout_predictions, holdout_targets)
        holdout_baseline_mse = compute_mean_baseline_mse(holdout_targets, train_mean_target)
        holdout_pearson, holdout_spearman = compute_correlations(holdout_predictions, holdout_targets)
        metrics.update(
            {
                "holdout_samples": len(holdout_dataset),
                "holdout_mse": holdout_mse,
                "holdout_rmse": math.sqrt(holdout_mse),
                "holdout_mean_baseline_mse": holdout_baseline_mse,
                "holdout_pearson": holdout_pearson,
                "holdout_spearman": holdout_spearman,
            }
        )
        print(
            f"[{run_label}] Holdout MSE: {holdout_mse:.6f} | "
            f"Pearson/Spearman: {holdout_pearson:.4f} / {holdout_spearman:.4f}"
        )

    if output_dir is not None:
        output_dir.mkdir(parents=True, exist_ok=True)
        metrics_path = output_dir / "metrics.json"
        config_path = output_dir / "config.json"
        history_path = output_dir / "history.csv"
        curve_path = build_curve_figure_path(output_dir, run_label)
        regression_path = build_regression_figure_path(output_dir, run_label)

        config_payload = {
            "model_type": "DeepSynergyMLP",
            "synergy_path": args.synergy_path,
            "fallback_pickle_path": args.fallback_pickle_path,
            "macro_file": args.macro_file,
            "macro_preset": args.macro_preset,
            "drug_dim": datasets.drug_dim,
            "gene_dim": datasets.gene_dim,
            "hidden_dims": args.hidden_dims,
            "dropout": args.dropout,
            "seed": run_seed,
            "lr": args.lr,
            "epochs": args.epochs,
            "batch_size": args.batch_size,
            "split_strategy": args.split_strategy,
            "train_fraction": args.train_fraction,
            "val_fraction": args.val_fraction,
            "max_samples": args.max_samples,
            "use_gene_expression": args.use_gene_expression,
            "bio_context": args.bio_context if args.use_gene_expression else None,
            "pathway_dim": datasets.gene_dim,
            "mlp0_out_dim": mlp0_out_dim,
            "normalize_pathways": args.normalize_pathways,
            "cell_encoder_type": "bio_context_mlp0",
        }
        if extra_config:
            config_payload.update(extra_config)

        save_json(config_path, config_payload)
        save_json(metrics_path, metrics)
        save_history_csv(history_path, history)
        save_loss_curve(curve_path, history, run_label)
        save_regression_plot(
            regression_path,
            targets=test_targets,
            predictions=test_predictions,
            run_label=run_label,
            mse=test_mse,
            pearson=test_pearson,
        )

        print(f"[{run_label}] Saved metrics to {metrics_path}")
        print(f"[{run_label}] Saved config to {config_path}")
        print(f"[{run_label}] Saved epoch history to {history_path}")
        print(f"[{run_label}] Saved loss curve to {curve_path}")
        print(f"[{run_label}] Saved test regression plot to {regression_path}")

    if save_artifacts and output_dir is not None:
        model_path = output_dir / "baseline_mlp.pt"
        val_predictions_path = output_dir / "val_predictions.csv"
        test_predictions_path = output_dir / "test_predictions.csv"

        torch.save(model.state_dict(), model_path)

        val_predictions_df = datasets.val_rows.copy().rename(columns={"target": "y_true"})
        val_predictions_df["y_pred"] = val_predictions
        val_predictions_df.to_csv(val_predictions_path, index=False)

        test_predictions_df = datasets.test_rows.copy().rename(columns={"target": "y_true"})
        test_predictions_df["y_pred"] = test_predictions
        test_predictions_df.to_csv(test_predictions_path, index=False)

        print(f"[{run_label}] Saved model to {model_path}")
        print(f"[{run_label}] Saved validation predictions to {val_predictions_path}")
        print(f"[{run_label}] Saved test predictions to {test_predictions_path}")

    print(f"[{run_label}] Val baseline MSE: {val_baseline_mse:.6f} | Test baseline MSE: {test_baseline_mse:.6f}")
    print(f"[{run_label}] Val Pearson/Spearman: {val_pearson:.4f} / {val_spearman:.4f}")
    print(f"[{run_label}] Test Pearson/Spearman: {test_pearson:.4f} / {test_spearman:.4f}")
    return metrics, eval_outputs


def run_cross_validation(args: argparse.Namespace) -> None:
    if args.split_strategy not in {"random", "drug_and_cell_line", "cell_line"}:
        raise ValueError(
            "Cross-validation is currently supported only with --split-strategy "
            "random, cell_line, or drug_and_cell_line."
        )
    if args.split_strategy in {"drug_and_cell_line", "cell_line"} and args.stratified_cv:
        print(f"WARNING: --stratified-cv is ignored for {args.split_strategy} CV; group-aware folds are used instead.")

    cv_seeds = args.cv_seeds or [args.seed]
    synergy_df, expression_lookup, gene_dim = load_aligned_synergy_and_expression(
        synergy_path=args.synergy_path,
        fallback_pickle_path=args.fallback_pickle_path,
        bio_context=args.bio_context,
        use_gene_expression=args.use_gene_expression,
        max_samples=args.max_samples,
    )

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    per_fold_metrics: list[dict[str, object]] = []
    all_test_predictions: list[pd.DataFrame] = []

    observed_fraction = args.train_fraction + args.val_fraction
    if observed_fraction <= 0 or observed_fraction >= 1:
        val_fraction_within_train = 0.1
    else:
        val_fraction_within_train = args.val_fraction / observed_fraction

    for seed in cv_seeds:
        if args.split_strategy == "random":
            folds = build_row_cv_folds(synergy_df, num_folds=args.cv_folds, seed=seed, stratified=args.stratified_cv)
            fold_indices = range(len(folds))
        else:
            fold_indices = range(args.cv_folds)

        for fold_idx in fold_indices:
            datasets = build_cv_dataset_bundle(
                synergy_df,
                expression_lookup,
                gene_dim=gene_dim,
                smiles_dim=args.smiles_dim,
                split_strategy=args.split_strategy,
                num_folds=args.cv_folds,
                seed=seed,
                fold_idx=fold_idx,
                train_fraction=args.train_fraction,
                val_fraction=args.val_fraction,
                stratified=args.stratified_cv,
                holdout_fraction=args.holdout_test_fraction,
                holdout_seed=args.holdout_seed,
                holdout_mode=args.holdout_mode,
                normalize_pathways=args.normalize_pathways,
            )
            run_label = f"cv_seed_{seed}_fold_{fold_idx + 1}"
            fold_output_dir = output_dir / "fold_runs" / run_label
            fold_extra_config = {
                "evaluation_mode": "cross_validation",
                "cv_seed": seed,
                "cv_fold": fold_idx + 1,
                "cv_folds": args.cv_folds,
                "cv_group_strategy": "row_folds" if args.split_strategy == "random" else f"{args.split_strategy}_priority",
                "stratified_cv": args.stratified_cv if args.split_strategy == "random" else False,
                "holdout_test_fraction": args.holdout_test_fraction,
                "holdout_seed": args.holdout_seed,
                "holdout_mode": args.holdout_mode,
            }
            metrics, eval_outputs = train_once(
                datasets,
                args,
                run_seed=seed,
                output_dir=fold_output_dir,
                save_artifacts=True,
                run_label=run_label,
                extra_config=fold_extra_config,
            )
            metrics["cv_seed"] = seed
            metrics["cv_fold"] = fold_idx + 1
            per_fold_metrics.append(metrics)

            fold_predictions = datasets.test_rows.copy().rename(columns={"target": "y_true"})
            fold_predictions["y_pred"] = eval_outputs["test_predictions"]
            fold_predictions["cv_seed"] = seed
            fold_predictions["cv_fold"] = fold_idx + 1
            all_test_predictions.append(fold_predictions)

    summary_frame = pd.DataFrame(per_fold_metrics)
    metric_columns = ["val_mse", "test_mse", "val_rmse", "test_rmse", "val_pearson", "test_pearson", "val_spearman", "test_spearman"]
    holdout_metric_columns = [
        col
        for col in ("holdout_mse", "holdout_rmse", "holdout_pearson", "holdout_spearman")
        if col in summary_frame.columns
    ]
    summary = {
        "evaluation_mode": "cross_validation",
        "cv_folds": args.cv_folds,
        "cv_seeds": cv_seeds,
        "stratified_cv": args.stratified_cv if args.split_strategy == "random" else False,
        "requested_stratified_cv": args.stratified_cv,
        "cv_group_strategy": "row_folds" if args.split_strategy == "random" else f"{args.split_strategy}_priority",
        "split_strategy": args.split_strategy,
        "macro_file": args.macro_file,
        "macro_preset": args.macro_preset,
        "use_gene_expression": args.use_gene_expression,
        "bio_context": args.bio_context if args.use_gene_expression else None,
        "pathway_dim": gene_dim,
        "mlp0_out_dim": resolve_mlp0_out_dim(args),
        "normalize_pathways": args.normalize_pathways,
        "holdout_test_fraction": args.holdout_test_fraction,
        "holdout_mode": args.holdout_mode if args.holdout_test_fraction > 0 else None,
        "aggregate_metrics": {
            metric: {
                "mean": float(summary_frame[metric].mean()),
                "std": float(summary_frame[metric].std(ddof=0)),
            }
            for metric in metric_columns + holdout_metric_columns
        },
    }

    with open(output_dir / "cv_metrics.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)
    summary_frame.to_csv(output_dir / "cv_runs.csv", index=False)
    if all_test_predictions:
        combined_predictions = pd.concat(all_test_predictions, ignore_index=True)
        combined_predictions.to_csv(output_dir / "cv_test_predictions.csv", index=False)

        regression_path = build_regression_figure_path(output_dir, "cv_aggregate", aggregate=True)
        aggregate_test_mse = float(summary_frame["test_mse"].mean())
        aggregate_test_pearson = float(summary_frame["test_pearson"].mean())
        save_regression_plot(
            regression_path,
            targets=combined_predictions["y_true"].tolist(),
            predictions=combined_predictions["y_pred"].tolist(),
            run_label="cv_aggregate",
            mse=aggregate_test_mse,
            pearson=aggregate_test_pearson,
        )

    print(f"Saved CV summary to {output_dir / 'cv_metrics.json'}")
    print(f"Saved per-run metrics to {output_dir / 'cv_runs.csv'}")
    if all_test_predictions:
        print(f"Saved CV test predictions to {output_dir / 'cv_test_predictions.csv'}")
        print(f"Saved CV test regression plot to {regression_path}")


def main() -> None:
    args = parse_args()

    if args.cv_folds > 1:
        run_cross_validation(args)
        return

    datasets = build_datasets(
        synergy_path=args.synergy_path,
        fallback_pickle_path=args.fallback_pickle_path,
        use_gene_expression=args.use_gene_expression,
        bio_context=args.bio_context,
        split_strategy=args.split_strategy,
        smiles_dim=args.smiles_dim,
        train_fraction=args.train_fraction,
        val_fraction=args.val_fraction,
        random_seed=args.seed,
        max_samples=args.max_samples,
        normalize_pathways=args.normalize_pathways,
    )

    output_dir = Path(args.output_dir)
    train_once(datasets, args, run_seed=args.seed, output_dir=output_dir, save_artifacts=True)


if __name__ == "__main__":
    main()
