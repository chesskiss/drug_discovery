from __future__ import annotations

import argparse
import json
from pathlib import Path

import torch

from .data_loading import load_expression_lookup, load_synergy_table
from .dataset import smiles_to_vector
from .bio_context import load_bio_context_matrix, project_expression
from .model import build_biocontext_model


BIO_CONTEXT_CHOICES = ["progeny", "kegg", "progeny_kegg"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run single-sample synergy prediction with a trained baseline model")
    parser.add_argument("--model-path", type=str, default="outputs/baseline_mlp.pt")
    parser.add_argument("--config-path", type=str, default="outputs/config.json")
    parser.add_argument("--synergy-path", type=str, default="data/drugcomb.csv")
    parser.add_argument(
        "--fallback-pickle-path",
        type=str,
        default="../data/data_compression/source_data/drugcomb.pkl",
    )
    parser.add_argument(
        "--bio-context",
        choices=BIO_CONTEXT_CHOICES,
        default=None,
        help="Optional override for the bio-context matrix (defaults to the trained config).",
    )
    parser.add_argument("--smiles-a", type=str, default=None)
    parser.add_argument("--smiles-b", type=str, default=None)
    parser.add_argument("--cell-line", type=str, default=None)
    parser.add_argument("--row-idx", type=int, default=None)
    parser.add_argument("--device", type=str, default=None)
    return parser.parse_args()


def get_device(preferred: str | None) -> torch.device:
    if preferred is not None:
        return torch.device(preferred)
    if torch.backends.mps.is_available():
        return torch.device("mps")
    if torch.cuda.is_available():
        return torch.device("cuda")
    return torch.device("cpu")


def load_config(path: str | Path) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def resolve_bio_context(args: argparse.Namespace, config: dict) -> str:
    if getattr(args, "bio_context", None):
        return str(args.bio_context)
    return str(config.get("bio_context") or "progeny")


def load_model(config: dict, model_path: str | Path, device: torch.device):
    model = build_biocontext_model(
        drug_dim=int(config["drug_dim"]),
        pathway_dim=int(config.get("pathway_dim", config.get("gene_dim", 0))),
        mlp0_out_dim=int(config.get("mlp0_out_dim", 0)),
        hidden_dims=config["hidden_dims"],
        dropout=float(config["dropout"]),
    )
    state_dict = torch.load(model_path, map_location=device)
    model.load_state_dict(state_dict)
    model.to(device)
    model.eval()
    return model


def resolve_sample(args: argparse.Namespace) -> dict[str, str | float | None]:
    if args.row_idx is not None:
        frame = load_synergy_table(args.synergy_path)
        if args.row_idx < 0 or args.row_idx >= len(frame):
            raise IndexError(f"row-idx {args.row_idx} is outside the dataset (0..{len(frame) - 1})")
        row = frame.iloc[args.row_idx]
        return {
            "smiles_a": row["smiles_a"],
            "smiles_b": row["smiles_b"],
            "cell_line": row["cell_line"],
            "target": float(row["target"]),
        }

    if args.smiles_a and args.smiles_b and args.cell_line:
        return {
            "smiles_a": args.smiles_a,
            "smiles_b": args.smiles_b,
            "cell_line": args.cell_line,
            "target": None,
        }

    raise ValueError("Provide either `--row-idx` or all of `--smiles-a`, `--smiles-b`, and `--cell-line`.")


def main() -> None:
    args = parse_args()
    config = load_config(args.config_path)
    device = get_device(args.device)
    model = load_model(config, args.model_path, device)
    sample = resolve_sample(args)

    raw_lookup = load_expression_lookup(
        cell_expression_path=None,
        fallback_pickle_path=args.fallback_pickle_path,
        feature_view_index=0,
    )
    weights, _names = load_bio_context_matrix(resolve_bio_context(args, config))
    expression_lookup = project_expression(raw_lookup, weights)
    cell_line = str(sample["cell_line"])
    if cell_line not in expression_lookup:
        raise KeyError(f"Cell line `{cell_line}` not found in expression lookup.")

    drug_a = torch.tensor(smiles_to_vector(str(sample["smiles_a"]), dim=int(config["drug_dim"])), dtype=torch.float32)
    drug_b = torch.tensor(smiles_to_vector(str(sample["smiles_b"]), dim=int(config["drug_dim"])), dtype=torch.float32)
    gene_expr = torch.tensor(expression_lookup[cell_line], dtype=torch.float32)

    with torch.no_grad():
        prediction = model(
            drug_a.unsqueeze(0).to(device),
            drug_b.unsqueeze(0).to(device),
            gene_expr.unsqueeze(0).to(device),
        ).item()

    print(f"device={device}")
    print(f"cell_line={cell_line}")
    print(f"predicted_synergy_zip={prediction:.6f}")
    if sample["target"] is not None:
        print(f"true_synergy_zip={float(sample['target']):.6f}")


if __name__ == "__main__":
    main()
