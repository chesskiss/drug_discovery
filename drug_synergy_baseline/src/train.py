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
    build_datasets,
    load_aligned_synergy_and_expression,
)
from .macros import DEFAULT_MACRO_FILE, DEFAULT_MACRO_PRESET, load_macro_preset
from .model import build_baseline_model


GENE_FEATURE_SET_TO_VIEW = {
    "raw": 0,
    "filtered": 1,
    "compact": 2,
}


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
    parser.add_argument("--synergy-path", type=str, default="data/drugcomb.csv")
    parser.add_argument("--cell-expression-path", type=str, default=None)
    parser.add_argument("--fallback-pickle-path", type=str, default="data/drugcomb.pkl")
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
    parser.add_argument("--split-strategy", choices=["random", "cell_line", "drug", "drug_pair"], default="random")
    parser.add_argument("--use-gene-expression", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument(
        "--gene-feature-set",
        choices=sorted(GENE_FEATURE_SET_TO_VIEW),
        default=None,
        help="Named CellLine view: raw=23808, filtered=3171, compact=627",
    )
    parser.add_argument("--cell-feature-view", type=int, default=0, help="Fallback-pickle CellLine view: 0=23808, 1=3171, 2=627")
    parser.add_argument("--cv-folds", type=int, default=1)
    parser.add_argument("--cv-seeds", type=int, nargs="*", default=None)
    parser.add_argument("--stratified-cv", action="store_true")
    return parser.parse_args()


def resolve_cell_feature_view(args: argparse.Namespace) -> int:
    if args.gene_feature_set is not None:
        return GENE_FEATURE_SET_TO_VIEW[args.gene_feature_set]
    return args.cell_feature_view


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
) -> tuple[dict[str, object], dict[str, object]]:
    torch.manual_seed(run_seed)

    device = get_device()
    model = build_baseline_model(
        drug_dim=datasets.drug_dim,
        gene_dim=datasets.gene_dim,
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
    print(f"[{run_label}] Gene expression dim: {datasets.gene_dim}")
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
        "gene_feature_set": args.gene_feature_set if args.use_gene_expression else None,
        "cell_feature_view": resolve_cell_feature_view(args) if args.use_gene_expression else None,
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

    if save_artifacts and output_dir is not None:
        output_dir.mkdir(parents=True, exist_ok=True)
        model_path = output_dir / "baseline_mlp.pt"
        metrics_path = output_dir / "metrics.json"
        config_path = output_dir / "config.json"
        val_predictions_path = output_dir / "val_predictions.csv"
        test_predictions_path = output_dir / "test_predictions.csv"

        torch.save(model.state_dict(), model_path)
        with open(config_path, "w", encoding="utf-8") as f:
            json.dump(
                {
                    "model_type": "DeepSynergyMLP",
                    "synergy_path": args.synergy_path,
                    "cell_expression_path": args.cell_expression_path,
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
                    "use_gene_expression": args.use_gene_expression,
                    "gene_feature_set": args.gene_feature_set if args.use_gene_expression else None,
                    "cell_feature_view": resolve_cell_feature_view(args) if args.use_gene_expression else None,
                    "cell_encoder_type": "identity",
                    "cell_latent_dim": datasets.gene_dim,
                },
                f,
                indent=2,
            )
        with open(metrics_path, "w", encoding="utf-8") as f:
            json.dump(metrics, f, indent=2)

        val_predictions_df = datasets.val_rows.copy().rename(columns={"target": "y_true"})
        val_predictions_df["y_pred"] = val_predictions
        val_predictions_df.to_csv(val_predictions_path, index=False)

        test_predictions_df = datasets.test_rows.copy().rename(columns={"target": "y_true"})
        test_predictions_df["y_pred"] = test_predictions
        test_predictions_df.to_csv(test_predictions_path, index=False)

        print(f"[{run_label}] Saved model to {model_path}")
        print(f"[{run_label}] Saved metrics to {metrics_path}")
        print(f"[{run_label}] Saved config to {config_path}")
        print(f"[{run_label}] Saved validation predictions to {val_predictions_path}")
        print(f"[{run_label}] Saved test predictions to {test_predictions_path}")

    print(f"[{run_label}] Val baseline MSE: {val_baseline_mse:.6f} | Test baseline MSE: {test_baseline_mse:.6f}")
    print(f"[{run_label}] Val Pearson/Spearman: {val_pearson:.4f} / {val_spearman:.4f}")
    print(f"[{run_label}] Test Pearson/Spearman: {test_pearson:.4f} / {test_spearman:.4f}")
    return metrics, eval_outputs


def build_cv_folds(frame: pd.DataFrame, num_folds: int, seed: int, stratified: bool) -> list[np.ndarray]:
    if num_folds < 2:
        raise ValueError("Cross-validation requires at least 2 folds.")

    n_rows = len(frame)
    if n_rows < num_folds:
        raise ValueError(f"Cannot build {num_folds} folds from only {n_rows} rows.")

    rng = np.random.default_rng(seed)
    fold_buckets: list[list[int]] = [[] for _ in range(num_folds)]

    if not stratified:
        indices = np.arange(n_rows)
        rng.shuffle(indices)
        for fold_idx, idx in enumerate(indices):
            fold_buckets[fold_idx % num_folds].append(int(idx))
        return [np.asarray(bucket, dtype=int) for bucket in fold_buckets]

    target = frame["target"]
    num_bins = min(num_folds, max(2, min(10, target.nunique())))
    bins = pd.qcut(target, q=num_bins, labels=False, duplicates="drop")
    strat_df = pd.DataFrame({"row_idx": np.arange(n_rows), "bin": bins})

    for _, group in strat_df.groupby("bin", dropna=False):
        indices = group["row_idx"].to_numpy(dtype=int)
        rng.shuffle(indices)
        for fold_idx, idx in enumerate(indices):
            fold_buckets[fold_idx % num_folds].append(int(idx))

    return [np.asarray(sorted(bucket), dtype=int) for bucket in fold_buckets]


def run_cross_validation(args: argparse.Namespace) -> None:
    if args.split_strategy != "random":
        raise ValueError("Cross-validation is currently only supported with --split-strategy random.")

    cv_seeds = args.cv_seeds or [args.seed]
    synergy_df, expression_lookup, gene_dim = load_aligned_synergy_and_expression(
        synergy_path=args.synergy_path,
        cell_expression_path=args.cell_expression_path,
        fallback_pickle_path=args.fallback_pickle_path,
        cell_feature_view=resolve_cell_feature_view(args),
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
        folds = build_cv_folds(synergy_df, num_folds=args.cv_folds, seed=seed, stratified=args.stratified_cv)
        for fold_idx, test_indices in enumerate(folds):
            mask = np.ones(len(synergy_df), dtype=bool)
            mask[test_indices] = False
            remaining = synergy_df.loc[mask].reset_index(drop=True)
            test_rows = synergy_df.iloc[test_indices].reset_index(drop=True)

            shuffled_indices = np.arange(len(remaining))
            fold_rng = np.random.default_rng(seed * 1000 + fold_idx)
            fold_rng.shuffle(shuffled_indices)

            val_count = max(1, int(len(shuffled_indices) * val_fraction_within_train))
            if val_count >= len(shuffled_indices):
                val_count = max(1, len(shuffled_indices) - 1)

            val_rows = remaining.iloc[shuffled_indices[:val_count]].reset_index(drop=True)
            train_rows = remaining.iloc[shuffled_indices[val_count:]].reset_index(drop=True)

            datasets = build_dataset_bundle_from_rows(
                train_rows=train_rows,
                val_rows=val_rows,
                test_rows=test_rows,
                expression_lookup=expression_lookup,
                smiles_dim=args.smiles_dim,
                gene_dim=gene_dim,
            )
            run_label = f"cv_seed_{seed}_fold_{fold_idx + 1}"
            metrics, eval_outputs = train_once(datasets, args, run_seed=seed, run_label=run_label)
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
    summary = {
        "evaluation_mode": "cross_validation",
        "cv_folds": args.cv_folds,
        "cv_seeds": cv_seeds,
        "stratified_cv": args.stratified_cv,
        "split_strategy": args.split_strategy,
        "macro_file": args.macro_file,
        "macro_preset": args.macro_preset,
        "use_gene_expression": args.use_gene_expression,
        "gene_feature_set": args.gene_feature_set if args.use_gene_expression else None,
        "cell_feature_view": resolve_cell_feature_view(args) if args.use_gene_expression else None,
        "aggregate_metrics": {
            metric: {
                "mean": float(summary_frame[metric].mean()),
                "std": float(summary_frame[metric].std(ddof=0)),
            }
            for metric in metric_columns
        },
    }

    with open(output_dir / "cv_metrics.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)
    summary_frame.to_csv(output_dir / "cv_runs.csv", index=False)
    if all_test_predictions:
        pd.concat(all_test_predictions, ignore_index=True).to_csv(output_dir / "cv_test_predictions.csv", index=False)

    print(f"Saved CV summary to {output_dir / 'cv_metrics.json'}")
    print(f"Saved per-run metrics to {output_dir / 'cv_runs.csv'}")
    if all_test_predictions:
        print(f"Saved CV test predictions to {output_dir / 'cv_test_predictions.csv'}")


def main() -> None:
    args = parse_args()

    if args.cv_folds > 1:
        run_cross_validation(args)
        return

    datasets = build_datasets(
        synergy_path=args.synergy_path,
        cell_expression_path=args.cell_expression_path,
        fallback_pickle_path=args.fallback_pickle_path,
        use_gene_expression=args.use_gene_expression,
        cell_feature_view=resolve_cell_feature_view(args),
        split_strategy=args.split_strategy,
        smiles_dim=args.smiles_dim,
        train_fraction=args.train_fraction,
        val_fraction=args.val_fraction,
        random_seed=args.seed,
        max_samples=args.max_samples,
    )

    output_dir = Path(args.output_dir)
    train_once(datasets, args, run_seed=args.seed, output_dir=output_dir, save_artifacts=True)


if __name__ == "__main__":
    main()
