from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import pandas as pd


def find_outputs_root(output_dir: str | Path) -> Path:
    output_path = Path(output_dir)
    if output_path.name == "outputs":
        return output_path

    for parent in [output_path, *output_path.parents]:
        if parent.name == "outputs":
            return parent

    if output_path.is_absolute():
        return output_path.parent

    if output_path.parts:
        return Path(output_path.parts[0])

    return output_path


def build_curve_figure_path(output_dir: str | Path, run_label: str) -> Path:
    output_path = Path(output_dir)
    outputs_root = find_outputs_root(output_path)
    relative_run_dir = output_path.relative_to(outputs_root)

    figure_dir = outputs_root / "training_curves"
    if str(relative_run_dir) == ".":
        relative_curve_dir = Path(run_label)
    else:
        relative_curve_dir = relative_run_dir

    return figure_dir / relative_curve_dir / "loss_curve.png"


def save_json(path: str | Path, payload: dict[str, Any]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2)


def save_history_csv(path: str | Path, history: list[dict[str, Any]]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(history).to_csv(path, index=False)


def save_loss_curve(path: str | Path, history: list[dict[str, Any]], run_label: str) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    history_frame = pd.DataFrame(history)
    if history_frame.empty:
        return

    plt.figure(figsize=(8, 5))
    plt.plot(history_frame["epoch"], history_frame["train_mse"], label="train_mse", linewidth=2)
    plt.plot(history_frame["epoch"], history_frame["val_mse"], label="val_mse", linewidth=2)
    plt.xlabel("Epoch")
    plt.ylabel("MSE")
    plt.title(f"Train vs Val Loss: {run_label}")
    plt.legend()
    plt.tight_layout()
    plt.savefig(path, dpi=160)
    plt.close()
