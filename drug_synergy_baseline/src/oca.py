from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch

from .dataset import build_datasets, smiles_to_vector
from .oca_plots import plot_component_importance_head_tail_summary, plot_component_importance_topk
from .predict import GENE_FEATURE_SET_TO_VIEW, get_device, load_config, load_model
from .training_artifacts import save_json


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run component-level occlusion attribution on a saved baseline checkpoint."
    )
    parser.add_argument("--model-path", type=str, required=True, help="Path to baseline_mlp.pt.")
    parser.add_argument("--config-path", type=str, required=True, help="Path to the matching config.json.")
    parser.add_argument("--synergy-path", type=str, default=None, help="Optional override for the synergy CSV.")
    parser.add_argument(
        "--cell-expression-path",
        type=str,
        default=None,
        help="Optional override for the cell-expression source.",
    )
    parser.add_argument(
        "--fallback-pickle-path",
        type=str,
        default=None,
        help="Optional override for the fallback pickle path.",
    )
    parser.add_argument(
        "--split-strategy",
        type=str,
        default=None,
        help="Optional override for the split strategy saved in config.",
    )
    parser.add_argument(
        "--gene-feature-set",
        choices=("raw", "filtered", "compact"),
        default=None,
        help="Optional override if you want to map to a specific built-in gene view.",
    )
    parser.add_argument(
        "--cell-feature-view",
        type=int,
        default=None,
        help="Optional direct override for the cell feature view index.",
    )
    parser.add_argument("--max-samples", type=int, default=None, help="Optional sample cap override.")
    parser.add_argument("--top-k", type=int, default=10, help="How many top components to plot.")
    parser.add_argument(
        "--local-row-idx",
        type=int,
        action="append",
        default=None,
        help="Explicit test-row indices for local explanations. Repeatable.",
    )
    parser.add_argument(
        "--local-top-n",
        type=int,
        default=5,
        help="How many rows to pick from each default local slice when explicit indices are not provided.",
    )
    parser.add_argument("--mask-value", type=float, default=0.0, help="Replacement value for masked components.")
    parser.add_argument(
        "--batch-size",
        type=int,
        default=None,
        help="Batch size used during attribution forward passes. Defaults to config batch size.",
    )
    parser.add_argument("--device", type=str, default=None, help="Optional device override.")
    parser.add_argument(
        "--output-dir",
        type=str,
        default=None,
        help="Where to write OCA artifacts. Defaults to <run>/oca.",
    )
    return parser.parse_args()


def _resolve_path_value(override: str | None, config: dict[str, object], key: str) -> str | None:
    if override is not None:
        return override
    value = config.get(key)
    if value in (None, ""):
        return None
    return str(value)


def _build_test_bundle(args: argparse.Namespace, config: dict[str, object]):
    synergy_path = _resolve_path_value(args.synergy_path, config, "synergy_path")
    if synergy_path is None:
        raise ValueError("No synergy path available. Pass --synergy-path or use a config with synergy_path.")

    use_gene_expression = bool(config.get("use_gene_expression", True))
    cell_expression_path = _resolve_path_value(args.cell_expression_path, config, "cell_expression_path")
    fallback_pickle_path = _resolve_path_value(args.fallback_pickle_path, config, "fallback_pickle_path")
    split_strategy = args.split_strategy or str(config.get("split_strategy", "random"))
    train_fraction = float(config.get("train_fraction", 0.8))
    val_fraction = float(config.get("val_fraction", 0.1))
    random_seed = int(config.get("seed", 42))
    smiles_dim = int(config.get("drug_dim", 256))
    max_samples = args.max_samples if args.max_samples is not None else config.get("max_samples")
    hidden_dims = config.get("hidden_dims")
    dropout = float(config.get("dropout", 0.2))

    if args.gene_feature_set is not None:
        cell_feature_view = GENE_FEATURE_SET_TO_VIEW[args.gene_feature_set]
    elif args.cell_feature_view is not None:
        cell_feature_view = int(args.cell_feature_view)
    elif config.get("cell_feature_view") is not None:
        cell_feature_view = int(config["cell_feature_view"])
    elif config.get("gene_feature_set") in GENE_FEATURE_SET_TO_VIEW:
        cell_feature_view = GENE_FEATURE_SET_TO_VIEW[str(config["gene_feature_set"])]
    else:
        cell_feature_view = 0

    datasets = build_datasets(
        synergy_path=synergy_path,
        cell_expression_path=cell_expression_path,
        fallback_pickle_path=fallback_pickle_path,
        use_gene_expression=use_gene_expression,
        cell_feature_view=cell_feature_view,
        split_strategy=split_strategy,
        smiles_dim=smiles_dim,
        train_fraction=train_fraction,
        val_fraction=val_fraction,
        random_seed=random_seed,
        max_samples=int(max_samples) if max_samples is not None else None,
    )
    return datasets, {
        "synergy_path": synergy_path,
        "cell_expression_path": cell_expression_path,
        "fallback_pickle_path": fallback_pickle_path,
        "split_strategy": split_strategy,
        "train_fraction": train_fraction,
        "val_fraction": val_fraction,
        "random_seed": random_seed,
        "smiles_dim": smiles_dim,
        "max_samples": max_samples,
        "hidden_dims": hidden_dims,
        "dropout": dropout,
        "cell_feature_view": cell_feature_view,
        "use_gene_expression": use_gene_expression,
    }


def _build_test_arrays(test_rows: pd.DataFrame, expression_lookup: dict[str, np.ndarray], smiles_dim: int):
    drug_a = np.stack([smiles_to_vector(smiles, dim=smiles_dim) for smiles in test_rows["smiles_a"]]).astype(np.float32)
    drug_b = np.stack([smiles_to_vector(smiles, dim=smiles_dim) for smiles in test_rows["smiles_b"]]).astype(np.float32)
    gene_expr = np.stack([expression_lookup[str(cell_line)] for cell_line in test_rows["cell_line"]]).astype(np.float32)
    targets = test_rows["target"].to_numpy(dtype=np.float32)
    return drug_a, drug_b, gene_expr, targets


def _predict_arrays(
    model: torch.nn.Module,
    *,
    drug_a: np.ndarray,
    drug_b: np.ndarray,
    gene_expr: np.ndarray,
    batch_size: int,
    device: torch.device,
) -> np.ndarray:
    outputs: list[np.ndarray] = []
    model.eval()
    with torch.no_grad():
        for start in range(0, len(drug_a), batch_size):
            end = start + batch_size
            batch_drug_a = torch.from_numpy(drug_a[start:end]).to(device=device, dtype=torch.float32)
            batch_drug_b = torch.from_numpy(drug_b[start:end]).to(device=device, dtype=torch.float32)
            batch_gene_expr = torch.from_numpy(gene_expr[start:end]).to(device=device, dtype=torch.float32)
            predictions = model(batch_drug_a, batch_drug_b, batch_gene_expr).detach().cpu().numpy()
            outputs.append(predictions.astype(np.float32))
    return np.concatenate(outputs, axis=0)


def _select_local_row_indices(
    test_rows: pd.DataFrame,
    base_predictions: np.ndarray,
    local_top_n: int,
    explicit_indices: list[int] | None,
) -> list[int]:
    if explicit_indices:
        return sorted({idx for idx in explicit_indices if 0 <= idx < len(test_rows)})

    targets = test_rows["target"].to_numpy(dtype=np.float32)
    residuals = np.abs(base_predictions - targets)

    high_target = np.argsort(targets)[-local_top_n:]
    low_target = np.argsort(targets)[:local_top_n]
    high_residual = np.argsort(residuals)[-local_top_n:]

    selected = sorted(set(high_target.tolist()) | set(low_target.tolist()) | set(high_residual.tolist()))
    return selected


def _plot_local_heatmap(
    local_df: pd.DataFrame,
    component_importance: pd.DataFrame,
    output_path: Path,
    top_k: int,
) -> None:
    if local_df.empty:
        return

    top_components = component_importance.nsmallest(top_k, "rank")["component_idx"].tolist()
    heatmap_df = local_df[local_df["component_idx"].isin(top_components)].copy()
    if heatmap_df.empty:
        return

    heatmap_df["sample_label"] = heatmap_df.apply(
        lambda row: f"{int(row['test_row_idx'])}:{row['cell_line']}",
        axis=1,
    )
    pivot = heatmap_df.pivot_table(
        index="sample_label",
        columns="component_idx",
        values="delta_squared_error",
        aggfunc="mean",
    ).reindex(columns=top_components, fill_value=0.0)

    fig, ax = plt.subplots(figsize=(max(8, len(top_components) * 0.8), max(4, len(pivot) * 0.45)))
    image = ax.imshow(pivot.to_numpy(dtype=np.float32), aspect="auto", cmap="viridis")
    ax.set_title("OCA Local Delta Squared Error Heatmap")
    ax.set_xlabel("PCA Component")
    ax.set_ylabel("Selected Test Rows")
    ax.set_xticks(range(len(pivot.columns)))
    ax.set_xticklabels([f"C{int(col)}" for col in pivot.columns], rotation=45, ha="right")
    ax.set_yticks(range(len(pivot.index)))
    ax.set_yticklabels(pivot.index.tolist())
    fig.colorbar(image, ax=ax, label="Delta Squared Error")
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    config = load_config(Path(args.config_path))
    datasets, resolved = _build_test_bundle(args, config)

    if not resolved["use_gene_expression"] or datasets.gene_dim <= 0:
        raise ValueError("OCA v1 requires a gene-enabled checkpoint with a positive gene dimension.")

    model_path = Path(args.model_path)
    output_dir = Path(args.output_dir) if args.output_dir else model_path.parent / "oca"
    output_dir.mkdir(parents=True, exist_ok=True)

    device = get_device(args.device)
    batch_size = int(args.batch_size or config.get("batch_size", 256))
    model = load_model(config, model_path, device)

    test_rows = datasets.test_rows.copy().reset_index(drop=True)
    expression_lookup = datasets.test.expression_lookup
    drug_a, drug_b, gene_expr, targets = _build_test_arrays(test_rows, expression_lookup, datasets.drug_dim)
    base_predictions = _predict_arrays(
        model,
        drug_a=drug_a,
        drug_b=drug_b,
        gene_expr=gene_expr,
        batch_size=batch_size,
        device=device,
    )
    base_squared_error = np.square(base_predictions - targets)
    base_absolute_error = np.abs(base_predictions - targets)
    local_indices = _select_local_row_indices(test_rows, base_predictions, args.local_top_n, args.local_row_idx)

    component_records: list[dict[str, float | int]] = []
    local_records: list[dict[str, object]] = []

    for component_idx in range(datasets.gene_dim):
        masked_gene_expr = gene_expr.copy()
        masked_gene_expr[:, component_idx] = np.float32(args.mask_value)
        masked_predictions = _predict_arrays(
            model,
            drug_a=drug_a,
            drug_b=drug_b,
            gene_expr=masked_gene_expr,
            batch_size=batch_size,
            device=device,
        )
        masked_squared_error = np.square(masked_predictions - targets)
        masked_absolute_error = np.abs(masked_predictions - targets)
        delta_prediction = masked_predictions - base_predictions
        delta_squared_error = masked_squared_error - base_squared_error
        delta_absolute_error = masked_absolute_error - base_absolute_error

        component_records.append(
            {
                "component_idx": component_idx,
                "mean_abs_delta_prediction": float(np.mean(np.abs(delta_prediction))),
                "mean_delta_squared_error": float(np.mean(delta_squared_error)),
                "mean_delta_absolute_error": float(np.mean(delta_absolute_error)),
            }
        )

        for row_idx in local_indices:
            row = test_rows.iloc[row_idx]
            local_records.append(
                {
                    "test_row_idx": row_idx,
                    "smiles_a": row["smiles_a"],
                    "smiles_b": row["smiles_b"],
                    "cell_line": row["cell_line"],
                    "component_idx": component_idx,
                    "y_true": float(targets[row_idx]),
                    "y_pred_base": float(base_predictions[row_idx]),
                    "y_pred_masked": float(masked_predictions[row_idx]),
                    "delta_prediction": float(delta_prediction[row_idx]),
                    "delta_squared_error": float(delta_squared_error[row_idx]),
                    "delta_absolute_error": float(delta_absolute_error[row_idx]),
                }
            )

    component_importance_df = pd.DataFrame(component_records).sort_values(
        by=["mean_delta_squared_error", "mean_abs_delta_prediction", "component_idx"],
        ascending=[False, False, True],
    )
    component_importance_df["rank"] = np.arange(1, len(component_importance_df) + 1)
    component_importance_df = component_importance_df[
        [
            "component_idx",
            "mean_abs_delta_prediction",
            "mean_delta_squared_error",
            "mean_delta_absolute_error",
            "rank",
        ]
    ]
    local_explanations_df = pd.DataFrame(local_records)

    component_path = output_dir / "component_importance.csv"
    local_path = output_dir / "local_explanations.csv"
    global_plot_path = output_dir / "component_importance_topk.png"
    summary_plot_path = output_dir / "component_importance_head_tail_summary.png"
    heatmap_path = output_dir / "local_explanations_heatmap.png"
    summary_path = output_dir / "oca_summary.json"

    component_importance_df.to_csv(component_path, index=False)
    local_explanations_df.to_csv(local_path, index=False)
    plot_component_importance_topk(component_importance_df, global_plot_path, args.top_k)
    plot_component_importance_head_tail_summary(component_importance_df, summary_plot_path)
    _plot_local_heatmap(local_explanations_df, component_importance_df, heatmap_path, args.top_k)
    save_json(
        summary_path,
        {
            "model_path": str(model_path),
            "config_path": str(args.config_path),
            "output_dir": str(output_dir),
            "gene_dim": int(datasets.gene_dim),
            "drug_dim": int(datasets.drug_dim),
            "test_rows": int(len(test_rows)),
            "mask_value": float(args.mask_value),
            "top_k": int(args.top_k),
            "local_row_indices": local_indices,
            "resolved": resolved,
        },
    )

    print(f"[oca] Saved component importance to {component_path}")
    print(f"[oca] Saved local explanations to {local_path}")
    print(f"[oca] Saved global importance plot to {global_plot_path}")
    print(f"[oca] Saved head/tail summary plot to {summary_plot_path}")
    print(f"[oca] Saved local heatmap to {heatmap_path}")
    print(f"[oca] Saved summary to {summary_path}")


if __name__ == "__main__":
    main()
