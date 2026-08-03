from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import matplotlib.pyplot as plt

from .training_artifacts import find_outputs_root

# Splits shown as separate lines, in a stable order with fixed styling.
SPLIT_STYLE = {
    "random": {"label": "random (no exclusion)", "color": "#1f77b4", "marker": "o"},
    "cell_line": {"label": "cell-line exclusion", "color": "#2ca02c", "marker": "s"},
    "drug_and_cell_line": {"label": "drug + cell-line exclusion", "color": "#d62728", "marker": "^"},
}

# Where the "no genes" (dim 0) points sit on the reversed log x-axis: just right of dim 1.
NO_GENES_X = 0.6


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot 10-fold CV test RMSE vs compressed cell-line dimension (PCA sweep)."
    )
    parser.add_argument(
        "--sweep-dir",
        type=str,
        default="outputs/sweep",
        help="Directory holding the per-run output dirs (mlp0{K}_{split}; mlp00_* / nogenes_* = no cell-line branch).",
    )
    parser.add_argument(
        "--output",
        type=str,
        default=None,
        help="PNG path. Default: <outputs>/training_curves/compression_rmse.png",
    )
    parser.add_argument(
        "--eval-name",
        type=str,
        default="10-fold CV",
        help="Evaluation label used in the title/axis (e.g. '10-fold CV' or 'single split').",
    )
    parser.add_argument(
        "--metric",
        choices=["rmse", "mse"],
        default="rmse",
        help="Which error metric to plot (test_rmse or test_mse).",
    )
    return parser.parse_args()


def _read_metric(metrics_path: Path, metric: str = "rmse") -> tuple[float, float] | None:
    """Return (mean, std) for the given metric ('rmse' or 'mse'). Works for both
    10-fold cv_metrics.json (aggregate mean/std across folds) and single-split
    metrics.json (scalar, std 0)."""
    key = f"test_{metric}"
    with metrics_path.open("r", encoding="utf-8") as handle:
        data = json.load(handle)
    if "aggregate_metrics" in data:  # 10-fold CV
        payload = data["aggregate_metrics"].get(key, {})
        mean = payload.get("mean")
        std = payload.get("std")
        if mean is None:
            return None
        return float(mean), float(std if std is not None else 0.0)
    mean = data.get(key)  # single split
    if mean is None:
        return None
    return float(mean), 0.0


def _metrics_file(run_dir: Path) -> Path | None:
    for name in ("cv_metrics.json", "metrics.json"):
        candidate = run_dir / name
        if candidate.exists():
            return candidate
    return None


def collect_points(sweep_dir: Path, metric: str = "rmse") -> dict[str, dict[str, list]]:
    """Return {split: {"dims":[...], "rmse":[...], "err":[...], "nogenes":(rmse,err)|None}}."""
    points: dict[str, dict[str, list]] = {
        s: {"dims": [], "rmse": [], "err": [],
            "fs_dims": [], "fs_rmse": [], "fs_err": [], "nogenes": None}
        for s in SPLIT_STYLE
    }
    for run_dir in sorted(sweep_dir.iterdir()):
        metrics_file = _metrics_file(run_dir)
        if metrics_file is None:
            continue
        name = run_dir.name
        rmse = _read_metric(metrics_file, metric)
        if rmse is None:
            continue
        mlp0_match = re.match(r"mlp0(\d+)_(.+)", name)
        fs_match = None
        nog_match = re.match(r"(?:nogenes|mlp00)_(.+)", name)
        if mlp0_match and int(mlp0_match.group(1)) > 0:
            dim, split = int(mlp0_match.group(1)), mlp0_match.group(2)
            if split in points:
                points[split]["dims"].append(dim)
                points[split]["rmse"].append(rmse[0])
                points[split]["err"].append(rmse[1])
        elif fs_match:
            dim, split = int(fs_match.group(1)), fs_match.group(2)
            if split in points:
                points[split]["fs_dims"].append(dim)
                points[split]["fs_rmse"].append(rmse[0])
                points[split]["fs_err"].append(rmse[1])
        elif nog_match:
            split = nog_match.group(1)
            if split in points:
                points[split]["nogenes"] = rmse
    # Sort PCA and feature-selection points by dimension.
    for split in points:
        for dkey, rkey, ekey in (("dims", "rmse", "err"), ("fs_dims", "fs_rmse", "fs_err")):
            order = sorted(range(len(points[split][dkey])), key=lambda i: points[split][dkey][i])
            for key in (dkey, rkey, ekey):
                points[split][key] = [points[split][key][i] for i in order]
    return points


def make_plot(points: dict[str, dict[str, list]], output_path: Path,
              eval_name: str = "10-fold CV", metric: str = "rmse") -> None:
    fig, ax = plt.subplots(figsize=(8.5, 6))

    for split, style in SPLIT_STYLE.items():
        p = points[split]
        if p["dims"]:
            ax.errorbar(
                p["dims"], p["rmse"], yerr=p["err"],
                color=style["color"], marker=style["marker"], markersize=7,
                capsize=3, linewidth=1.8, label=style["label"],
            )
        # Feature-selection points (627/3171/23808): different compression method,
        # shown dashed with hollow markers in the same split colour, connected to the
        # PCA curve's high-dim end so each split reads as one continuous story.
        if p["fs_dims"]:
            # Continue the PCA curve outward: 58 -> 627 -> 3171 -> 23808 (prepend the
            # PCA endpoint so the dashed line is monotonic and never backtracks).
            link_x = ([p["dims"][-1]] if p["dims"] else []) + list(p["fs_dims"])
            link_y = ([p["rmse"][-1]] if p["dims"] else []) + list(p["fs_rmse"])
            ax.errorbar(
                link_x, link_y, yerr=([0.0] if p["dims"] else []) + list(p["fs_err"]),
                color=style["color"], marker=style["marker"], markersize=8,
                markerfacecolor="white", markeredgewidth=1.6,
                capsize=3, linewidth=1.4, linestyle="--",
            )
        if p["nogenes"] is not None:
            mean, err = p["nogenes"]
            ax.errorbar(
                [NO_GENES_X], [mean], yerr=[err],
                color=style["color"], marker="*", markersize=16,
                capsize=3, linestyle="none",
            )

    ax.set_xscale("log")
    ax.invert_xaxis()  # high dimension (more info) on the left, "none" on the right

    # Custom ticks: all real dims (PCA + feature-selection) + a dedicated "none" slot.
    all_dims = sorted({d for p in points.values() for d in (*p["dims"], *p["fs_dims"])})
    ticks = all_dims + [NO_GENES_X]
    labels = [str(d) for d in all_dims] + ["none"]
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels)
    ax.minorticks_off()

    err_note = "mean ± fold std" if eval_name == "10-fold CV" else "single split"
    metric_label = metric.upper()
    any_fs = any(points[s]["fs_dims"] for s in points)
    ax.set_xlabel("Compressed cell-line dimension (log, reversed) → less info")
    ax.set_ylabel(f"{eval_name} test {metric_label}  ({err_note})")
    ax.set_title(f"Bio-context model ({eval_name}): {metric_label} vs MLP0 width")
    ax.grid(True, which="major", axis="both", alpha=0.25)

    # Legend: split lines + method/marker key.
    extra = [plt.Line2D([], [], color="0.3", marker="*", markersize=13, linestyle="none",
                        label="no cell-line branch (MLP0=0)")]
    if any_fs:
        extra.append(plt.Line2D([], [], color="0.3", marker="o", markersize=7, linestyle="-",
                                label="bio-context (solid, filled)"))
        extra.append(plt.Line2D([], [], color="0.3", marker="o", markersize=7, linestyle="--",
                                markerfacecolor="white", label="feature-selection (dashed, hollow)"))
    handles, _ = ax.get_legend_handles_labels()
    ax.legend(handles=handles + extra, loc="best", fontsize=8)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(output_path, dpi=160)
    plt.close(fig)


def print_table(points: dict[str, dict[str, list]], metric: str = "rmse") -> None:
    print(f"{'split':<26}{'dim':>6}{metric + '_mean':>12}{metric + '_std':>10}")
    for split in SPLIT_STYLE:
        p = points[split]
        for dim, rmse, err in zip(p["fs_dims"], p["fs_rmse"], p["fs_err"]):
            print(f"{split:<26}{dim:>6}{rmse:>12.4f}{err:>10.4f}   [feat-select]")
        for dim, rmse, err in zip(p["dims"], p["rmse"], p["err"]):
            print(f"{split:<26}{dim:>6}{rmse:>12.4f}{err:>10.4f}   [mlp0]")
        if p["nogenes"] is not None:
            m, e = p["nogenes"]
            print(f"{split:<26}{'none':>6}{m:>12.4f}{e:>10.4f}")


def main() -> None:
    args = parse_args()
    sweep_dir = Path(args.sweep_dir)
    if not sweep_dir.exists():
        raise FileNotFoundError(f"Sweep directory not found: {sweep_dir}")

    points = collect_points(sweep_dir, metric=args.metric)
    if args.output is not None:
        output_path = Path(args.output)
    else:
        fname = f"compression_{args.metric}.png"
        output_path = find_outputs_root(sweep_dir) / "training_curves" / fname

    print_table(points, metric=args.metric)
    make_plot(points, output_path, eval_name=args.eval_name, metric=args.metric)
    print(f"\nSaved compression-vs-{args.metric.upper()} plot to {output_path}")


if __name__ == "__main__":
    main()
