from __future__ import annotations

import argparse
import json
from pathlib import Path

from .training_artifacts import build_curve_figure_path, save_history_csv, save_loss_curve


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Backfill history.csv files and loss curves from existing metrics.json files.")
    parser.add_argument("--outputs-dir", type=str, default="outputs")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    outputs_dir = Path(args.outputs_dir)

    metrics_paths = sorted(outputs_dir.rglob("metrics.json"))
    processed = 0
    skipped = 0

    for metrics_path in metrics_paths:
        with metrics_path.open("r", encoding="utf-8") as handle:
            metrics = json.load(handle)

        history = metrics.get("history")
        if not history:
            skipped += 1
            continue

        run_dir = metrics_path.parent
        run_label = str(metrics.get("run_label", run_dir.name))
        history_path = run_dir / "history.csv"
        figure_path = build_curve_figure_path(run_dir, run_label)

        save_history_csv(history_path, history)
        save_loss_curve(figure_path, history, run_label)
        processed += 1

    print(f"Backfilled history and curves for {processed} runs.")
    if skipped:
        print(f"Skipped {skipped} runs with no history block.")


if __name__ == "__main__":
    main()
