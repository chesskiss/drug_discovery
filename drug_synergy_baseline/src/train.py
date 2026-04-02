from __future__ import annotations

import argparse
import json
from pathlib import Path

import torch
from torch import nn
from torch.utils.data import DataLoader

from .dataset import build_datasets
from .model import DeepSynergyMLP


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Train a minimal DeepSynergy-style baseline on DrugComb/TDC data")
    parser.add_argument("--synergy-path", type=str, default="data/drugcomb_synergy.csv")
    parser.add_argument("--cell-expression-path", type=str, default="data/cell_line_gene_expression.csv")
    parser.add_argument("--fallback-pickle-path", type=str, default="data/drugcomb.pkl")
    parser.add_argument("--output-dir", type=str, default="outputs")
    parser.add_argument("--epochs", type=int, default=10)
    parser.add_argument("--batch-size", type=int, default=64)
    parser.add_argument("--lr", type=float, default=1e-3)
    parser.add_argument("--smiles-dim", type=int, default=256)
    parser.add_argument("--hidden-dims", type=int, nargs="+", default=[512, 256])
    parser.add_argument("--dropout", type=float, default=0.2)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--max-samples", type=int, default=None)
    return parser.parse_args()


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


def main() -> None:
    args = parse_args()
    torch.manual_seed(args.seed)

    datasets = build_datasets(
        synergy_path=args.synergy_path,
        cell_expression_path=args.cell_expression_path,
        fallback_pickle_path=args.fallback_pickle_path,
        smiles_dim=args.smiles_dim,
        random_seed=args.seed,
        max_samples=args.max_samples,
    )

    device = get_device()
    model = DeepSynergyMLP(
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

    print(f"Train samples: {len(datasets.train)}")
    print(f"Val samples: {len(datasets.val)}")
    print(f"Test samples: {len(datasets.test)}")
    print(f"Drug feature dim: {datasets.drug_dim}")
    print(f"Gene expression dim: {datasets.gene_dim}")
    print(f"Device: {device}")

    for epoch in range(1, args.epochs + 1):
        train_loss = run_epoch(model, train_loader, optimizer, criterion, device, train=True)
        val_loss = run_epoch(model, val_loader, optimizer, criterion, device, train=False)
        print(f"Epoch {epoch:03d} | train_mse={train_loss:.6f} | val_mse={val_loss:.6f}")

    test_loss = run_epoch(model, test_loader, optimizer, criterion, device, train=False)
    print(f"Test MSE: {test_loss:.6f}")

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    model_path = output_dir / "baseline_mlp.pt"
    metrics_path = output_dir / "metrics.json"

    torch.save(model.state_dict(), model_path)
    with open(metrics_path, "w", encoding="utf-8") as f:
        json.dump(
            {
                "train_samples": len(datasets.train),
                "val_samples": len(datasets.val),
                "test_samples": len(datasets.test),
                "drug_dim": datasets.drug_dim,
                "gene_dim": datasets.gene_dim,
                "test_mse": test_loss,
            },
            f,
            indent=2,
        )

    print(f"Saved model to {model_path}")
    print(f"Saved metrics to {metrics_path}")


if __name__ == "__main__":
    main()
