from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import mean_absolute_error, mean_squared_error
from sklearn.preprocessing import StandardScaler

from .analysis import analyze_cell_line_difficulty


PREDICTIVE_TARGET = "ease_score"
DEFAULT_FEATURE_VIEW_INDEX = 1
VIEW_LENGTHS = {0: 23808, 1: 3171, 2: 627}


def load_predictive_dataset(
    synergy_path: str | Path,
    pickle_path: str | Path,
    feature_view_index: int = DEFAULT_FEATURE_VIEW_INDEX,
    max_cell_lines: int | None = None,
) -> pd.DataFrame:
    targets = analyze_cell_line_difficulty(synergy_path).copy()
    targets = targets.loc[
        :,
        [
            "cell_line",
            "n_rows",
            "mean_zip",
            "median_zip",
            "std_zip",
            "max_zip",
            "q75_zip",
            "q90_zip",
            "positive_rate_zip_gt_0",
            "high_synergy_rate_zip_gt_10",
            "reliability_weight",
            "ease_score",
            "difficulty_score",
            "ease_rank",
            "difficulty_rank",
        ],
    ]

    raw = pd.read_pickle(pickle_path)
    feature_rows: dict[str, np.ndarray] = {}
    for _, row in raw[["Cell_Line_ID", "CellLine"]].dropna().iterrows():
        cell_line = str(row["Cell_Line_ID"])
        cell_features = row["CellLine"]
        if not isinstance(cell_features, (list, tuple)) or len(cell_features) <= feature_view_index:
            continue
        feature_rows.setdefault(cell_line, np.asarray(cell_features[feature_view_index], dtype=np.float32))

    if not feature_rows:
        raise ValueError("No cell-line feature vectors were recovered from the pickle file.")

    first_vector = next(iter(feature_rows.values()))
    feature_columns = [f"feature_{idx:04d}" for idx in range(len(first_vector))]
    feature_df = pd.DataFrame.from_dict(feature_rows, orient="index", columns=feature_columns)
    feature_df.index.name = "cell_line"
    feature_df = feature_df.reset_index()

    dataset = targets.merge(feature_df, on="cell_line", how="inner")
    dataset = dataset.sort_values("cell_line").reset_index(drop=True)
    dataset["feature_view_index"] = feature_view_index
    dataset["feature_view_length"] = len(first_vector)

    if max_cell_lines is not None:
        dataset = dataset.head(max_cell_lines).copy()

    if dataset["cell_line"].duplicated().any():
        raise ValueError("Predictive dataset contains duplicated cell lines after merging targets and features.")

    return dataset


def _pearson_corr(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    frame = pd.DataFrame({"y_true": y_true, "y_pred": y_pred})
    value = frame["y_true"].corr(frame["y_pred"], method="pearson")
    return float(value) if pd.notna(value) else float("nan")


def _spearman_corr(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    frame = pd.DataFrame({"y_true": y_true, "y_pred": y_pred})
    value = frame["y_true"].corr(frame["y_pred"], method="spearman")
    return float(value) if pd.notna(value) else float("nan")


def _standardize_train_test(x_train: np.ndarray, x_test: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    scaler = StandardScaler()
    x_train_scaled = scaler.fit_transform(x_train).astype(np.float64)
    x_test_scaled = scaler.transform(x_test).astype(np.float64)
    return x_train_scaled, x_test_scaled


def _fit_predict_ridge_dual(
    x_train: np.ndarray,
    y_train: np.ndarray,
    x_test: np.ndarray,
    alpha: float = 1.0,
) -> float:
    x_train_scaled = x_train.astype(np.float64)
    x_test_scaled = x_test.astype(np.float64)
    y_mean = float(np.mean(y_train))
    y_centered = y_train - y_mean
    kernel = x_train_scaled @ x_train_scaled.T
    dual_coef = np.linalg.solve(kernel + alpha * np.eye(kernel.shape[0], dtype=np.float64), y_centered)
    prediction = y_mean + float((x_test_scaled @ x_train_scaled.T) @ dual_coef)
    return prediction


def _ridge_loocv_predictions_fast(x: np.ndarray, y: np.ndarray, alpha: float = 1.0) -> np.ndarray:
    x64 = x.astype(np.float64)
    y_mean = float(np.mean(y))
    y_centered = y.astype(np.float64) - y_mean
    kernel = x64 @ x64.T
    system_inv = np.linalg.inv(kernel + alpha * np.eye(kernel.shape[0], dtype=np.float64))
    dual_coef = system_inv @ y_centered
    loo_centered = y_centered - dual_coef / np.diag(system_inv)
    return loo_centered + y_mean


def _fit_predict_tiny_mlp(x_train: np.ndarray, y_train: np.ndarray, x_test: np.ndarray) -> float:
    x_train_scaled, x_test_scaled = _standardize_train_test(x_train, x_test)

    rng = np.random.default_rng(42)
    projection = rng.normal(0.0, 1.0 / np.sqrt(x_train_scaled.shape[1]), size=(x_train_scaled.shape[1], 16))
    x_train_small = x_train_scaled @ projection
    x_test_small = x_test_scaled @ projection

    y_mean = float(np.mean(y_train))
    y_std = float(np.std(y_train))
    if y_std == 0:
        y_std = 1.0
    y_train_norm = ((y_train - y_mean) / y_std).reshape(-1, 1)

    rng = np.random.default_rng(42)
    hidden_dim = 4
    input_dim = x_train_small.shape[1]
    w1 = rng.normal(0.0, 0.1, size=(input_dim, hidden_dim))
    b1 = np.zeros((1, hidden_dim), dtype=np.float64)
    w2 = rng.normal(0.0, 0.1, size=(hidden_dim, 1))
    b2 = np.zeros((1, 1), dtype=np.float64)

    learning_rate = 0.05
    l2 = 1e-2
    sample_count = x_train_small.shape[0]

    for _ in range(20):
        hidden = np.tanh(x_train_small @ w1 + b1)
        pred = hidden @ w2 + b2
        error = pred - y_train_norm

        grad_w2 = (hidden.T @ error) / sample_count + l2 * w2
        grad_b2 = np.mean(error, axis=0, keepdims=True)

        hidden_grad = (error @ w2.T) * (1.0 - hidden**2)
        grad_w1 = (x_train_small.T @ hidden_grad) / sample_count + l2 * w1
        grad_b1 = np.mean(hidden_grad, axis=0, keepdims=True)

        w2 -= learning_rate * grad_w2
        b2 -= learning_rate * grad_b2
        w1 -= learning_rate * grad_w1
        b1 -= learning_rate * grad_b1

    hidden_test = np.tanh(x_test_small @ w1 + b1)
    pred_norm = hidden_test @ w2 + b2
    return float(pred_norm[0, 0] * y_std + y_mean)


def run_loocv(
    dataset: pd.DataFrame,
    models: tuple[str, ...] = ("ridge", "mlp"),
    progress: bool = False,
) -> tuple[pd.DataFrame, dict]:
    feature_columns = [column for column in dataset.columns if column.startswith("feature_")]
    if not feature_columns:
        raise ValueError("Predictive dataset has no feature columns.")

    x = dataset[feature_columns].to_numpy(dtype=np.float32)
    y = dataset[PREDICTIVE_TARGET].to_numpy(dtype=np.float64)
    cell_lines = dataset["cell_line"].tolist()

    run_ridge = "ridge" in models
    run_mlp = "mlp" in models
    if not run_ridge and not run_mlp:
        raise ValueError("At least one predictive model must be selected.")

    ridge_predictions: list[float] = []
    mlp_predictions: list[float] = []
    baseline_predictions: list[float] = []

    for held_out_idx in range(len(dataset)):
        if progress:
            print(f"LOOCV fold {held_out_idx + 1}/{len(dataset)}")
        train_mask = np.ones(len(dataset), dtype=bool)
        train_mask[held_out_idx] = False

        x_train = x[train_mask]
        y_train = y[train_mask]
        x_test = x[[held_out_idx]]

        if run_ridge:
            ridge_predictions.append(_fit_predict_ridge_dual(x_train, y_train, x_test))

        if run_mlp:
            mlp_predictions.append(_fit_predict_tiny_mlp(x_train, y_train, x_test))

        baseline_predictions.append(float(np.mean(y_train)))

    payload: dict[str, object] = {
        "cell_line": cell_lines,
        "y_true": y,
        "baseline_pred": baseline_predictions,
    }
    if run_ridge:
        payload["ridge_pred"] = ridge_predictions
    if run_mlp:
        payload["mlp_pred"] = mlp_predictions

    predictions = pd.DataFrame(payload)
    if run_ridge:
        predictions["ridge_abs_error"] = np.abs(predictions["ridge_pred"] - predictions["y_true"])
    if run_mlp:
        predictions["mlp_abs_error"] = np.abs(predictions["mlp_pred"] - predictions["y_true"])
    predictions["baseline_abs_error"] = np.abs(predictions["baseline_pred"] - predictions["y_true"])

    metrics = {
        "sample_count": int(len(dataset)),
        "feature_dimension": int(len(feature_columns)),
        "feature_view_index": int(dataset["feature_view_index"].iloc[0]),
        "feature_view_length": int(dataset["feature_view_length"].iloc[0]),
        "target_name": PREDICTIVE_TARGET,
        "models_requested": list(models),
        "score_definition": {
            "ease_score": "reliability_weight * mean(z_mean_zip, z_high_synergy_rate_zip_gt_10, z_positive_rate_zip_gt_0)",
            "difficulty_score": "-ease_score",
        },
        "models": {},
    }
    if run_ridge:
        metrics["models"]["ridge"] = {
            "rmse": float(np.sqrt(mean_squared_error(predictions["y_true"], predictions["ridge_pred"]))),
            "mae": float(mean_absolute_error(predictions["y_true"], predictions["ridge_pred"])),
            "pearson": _pearson_corr(predictions["y_true"].to_numpy(), predictions["ridge_pred"].to_numpy()),
            "spearman": _spearman_corr(predictions["y_true"].to_numpy(), predictions["ridge_pred"].to_numpy()),
        }
    if run_mlp:
        metrics["models"]["mlp"] = {
            "rmse": float(np.sqrt(mean_squared_error(predictions["y_true"], predictions["mlp_pred"]))),
            "mae": float(mean_absolute_error(predictions["y_true"], predictions["mlp_pred"])),
            "pearson": _pearson_corr(predictions["y_true"].to_numpy(), predictions["mlp_pred"].to_numpy()),
            "spearman": _spearman_corr(predictions["y_true"].to_numpy(), predictions["mlp_pred"].to_numpy()),
        }
    metrics["models"]["fold_mean_baseline"] = {
        "rmse": float(np.sqrt(mean_squared_error(predictions["y_true"], predictions["baseline_pred"]))),
        "mae": float(mean_absolute_error(predictions["y_true"], predictions["baseline_pred"])),
        "pearson": _pearson_corr(predictions["y_true"].to_numpy(), predictions["baseline_pred"].to_numpy()),
        "spearman": _spearman_corr(predictions["y_true"].to_numpy(), predictions["baseline_pred"].to_numpy()),
    }

    return predictions.sort_values("cell_line").reset_index(drop=True), metrics


def compare_feature_views(
    synergy_path: str | Path,
    pickle_path: str | Path,
    view_indices: tuple[int, ...] = (0, 1, 2),
    models: tuple[str, ...] = ("ridge",),
    max_cell_lines: int | None = None,
    progress: bool = False,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for view_index in view_indices:
        dataset = load_predictive_dataset(
            synergy_path=synergy_path,
            pickle_path=pickle_path,
            feature_view_index=view_index,
            max_cell_lines=max_cell_lines,
        )
        if models == ("ridge",):
            feature_columns = [column for column in dataset.columns if column.startswith("feature_")]
            x = dataset[feature_columns].to_numpy(dtype=np.float32)
            y = dataset[PREDICTIVE_TARGET].to_numpy(dtype=np.float64)
            ridge_pred = _ridge_loocv_predictions_fast(x, y)
            baseline_pred = np.array([np.mean(np.delete(y, idx)) for idx in range(len(y))], dtype=np.float64)
            metrics = {
                "sample_count": int(len(dataset)),
                "feature_view_length": int(dataset["feature_view_length"].iloc[0]),
                "target_name": PREDICTIVE_TARGET,
                "models": {
                    "ridge": {
                        "rmse": float(np.sqrt(mean_squared_error(y, ridge_pred))),
                        "mae": float(mean_absolute_error(y, ridge_pred)),
                        "pearson": _pearson_corr(y, ridge_pred),
                        "spearman": _spearman_corr(y, ridge_pred),
                    },
                    "fold_mean_baseline": {
                        "rmse": float(np.sqrt(mean_squared_error(y, baseline_pred))),
                        "mae": float(mean_absolute_error(y, baseline_pred)),
                        "pearson": _pearson_corr(y, baseline_pred),
                        "spearman": _spearman_corr(y, baseline_pred),
                    },
                },
            }
        else:
            _, metrics = run_loocv(dataset, models=models, progress=progress)
        row: dict[str, object] = {
            "feature_view_index": view_index,
            "feature_view_length": int(metrics["feature_view_length"]),
            "sample_count": int(metrics["sample_count"]),
            "target_name": metrics["target_name"],
        }
        for model_name, model_metrics in metrics["models"].items():
            for metric_name, metric_value in model_metrics.items():
                row[f"{model_name}_{metric_name}"] = metric_value
        rows.append(row)
    return pd.DataFrame(rows).sort_values("feature_view_index").reset_index(drop=True)


def save_predictive_outputs(
    dataset: pd.DataFrame,
    predictions: pd.DataFrame,
    metrics: dict,
    output_dir: str | Path,
) -> None:
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    dataset.to_csv(output_path / "predictive_dataset.csv", index=False)
    predictions.to_csv(output_path / "loocv_predictions.csv", index=False)
    with open(output_path / "predictive_metrics.json", "w", encoding="utf-8") as file:
        json.dump(metrics, file, indent=2)


def save_view_comparison_outputs(comparison: pd.DataFrame, output_dir: str | Path) -> None:
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    comparison.to_csv(output_path / "view_comparison.csv", index=False)
