from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import pandas as pd

from .macros import load_macro_preset


BASELINE_DIR = Path(__file__).resolve().parents[1]
RESULTS_DIR = BASELINE_DIR / "results"
CANONICAL_RESULTS_PATH = RESULTS_DIR / "baseline_experiments.csv"
SUMMARY_RESULTS_PATH = RESULTS_DIR / "baseline_summary.csv"

CANONICAL_COLUMNS = [
    "run_id",
    "stage_id",
    "status",
    "phase",
    "model_id",
    "gene_setting",
    "gene_dim",
    "split_strategy",
    "eval_mode",
    "sample_cap",
    "macro_preset",
    "seed_spec",
    "epochs",
    "lr",
    "batch_size",
    "test_mse",
    "test_rmse",
    "test_pearson",
    "test_spearman",
    "test_mean_baseline_mse",
    "improvement_vs_mean_baseline",
    "source_metrics_path",
    "notes",
]

COMPETITIVE_MODELS = {"no_genes", "pca128"}
TIE_BREAK_PRIORITY = {"no_genes": 0, "pca128": 1}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build the canonical baseline results tables.")
    parser.add_argument(
        "--results-dir",
        type=str,
        default=str(RESULTS_DIR),
        help="Directory where baseline_experiments.csv and baseline_summary.csv will be written.",
    )
    return parser.parse_args()


def _relative_to_baseline(path: str | Path | None) -> str:
    if path in (None, ""):
        return ""
    path_obj = Path(path)
    if path_obj.is_absolute():
        return str(path_obj.relative_to(BASELINE_DIR))
    return str(path_obj)


def _metric_paths(output_subdir: str, eval_mode: str) -> tuple[Path, Path]:
    run_dir = BASELINE_DIR / output_subdir
    metrics_name = "cv_metrics.json" if eval_mode == "cv10_stratified" else "metrics.json"
    return run_dir / metrics_name, run_dir / "config.json"


def _empty_row(**overrides: Any) -> dict[str, Any]:
    row = {column: "" for column in CANONICAL_COLUMNS}
    row.update(overrides)
    return row


def _legacy_specs() -> list[dict[str, Any]]:
    return [
        _empty_row(
            run_id="legacy_raw_genes_random",
            stage_id="1.1",
            status="planned",
            phase="short_single_split",
            model_id="raw_genes",
            gene_setting="raw",
            split_strategy="random",
            eval_mode="single_split",
            sample_cap="full",
            macro_preset="legacy_manual",
            seed_spec="42",
            source_metrics_path="outputs/genes_on_random/metrics.json",
            notes="Archived negative control from the pre-funnel baseline runs.",
        ),
        _empty_row(
            run_id="legacy_filtered_genes_random",
            stage_id="1.1",
            status="planned",
            phase="short_single_split",
            model_id="filtered_genes",
            gene_setting="filtered",
            split_strategy="random",
            eval_mode="single_split",
            sample_cap="full",
            macro_preset="legacy_manual",
            seed_spec="42",
            source_metrics_path="outputs/genes_filtered_random/metrics.json",
            notes="Archived negative control from the pre-funnel baseline runs.",
        ),
        _empty_row(
            run_id="legacy_no_genes_random",
            stage_id="1.1",
            status="planned",
            phase="short_single_split",
            model_id="no_genes",
            gene_setting="off",
            split_strategy="random",
            eval_mode="single_split",
            sample_cap="full",
            macro_preset="legacy_manual",
            seed_spec="42",
            source_metrics_path="outputs/genes_off_random/metrics.json",
            notes="Historical reference before the canonical short-phase matrix was frozen.",
        ),
        _empty_row(
            run_id="legacy_pca128_random",
            stage_id="1.3A",
            status="planned",
            phase="short_single_split",
            model_id="pca128",
            gene_setting="pca128",
            split_strategy="random",
            eval_mode="single_split",
            sample_cap="full",
            macro_preset="legacy_manual",
            seed_spec="42",
            source_metrics_path="outputs/pca128_random/metrics.json",
            notes="Historical reference before the canonical short-phase matrix was frozen.",
        ),
    ]


def _planned_short_phase_specs() -> list[dict[str, Any]]:
    practical = load_macro_preset("practical_research")
    planned_rows: list[dict[str, Any]] = []
    for model_id, stage_id, gene_setting, gene_dim, output_prefix in [
        ("no_genes", "1.1", "off", 0, "no_genes"),
        ("pca128", "1.3A", "pca128", 128, "pca128"),
    ]:
        for split_strategy in ("random", "drug", "cell_line", "drug_and_cell_line"):
            run_name = f"short_{output_prefix}_{split_strategy}_practical"
            planned_rows.append(
                _empty_row(
                    run_id=run_name,
                    stage_id=stage_id,
                    status="planned",
                    phase="short_single_split",
                    model_id=model_id,
                    gene_setting=gene_setting,
                    gene_dim=gene_dim,
                    split_strategy=split_strategy,
                    eval_mode="single_split",
                    sample_cap="full",
                    macro_preset="practical_research",
                    seed_spec="42",
                    epochs=int(practical["epochs"]),
                    lr=float(practical["lr"]),
                    batch_size=int(practical["batch_size"]),
                    source_metrics_path=f"outputs/{run_name}/metrics.json",
                    notes="Canonical short-phase comparison row.",
                )
            )
    return planned_rows


def _planned_long_phase_placeholder() -> dict[str, Any]:
    practical = load_macro_preset("practical_research")
    return _empty_row(
        run_id="long_cv10_10k_promoted_winner",
        stage_id="promotion_gate",
        status="planned",
        phase="long_cv_10k",
        model_id="",
        gene_setting="",
        gene_dim="",
        split_strategy="random",
        eval_mode="cv10_stratified",
        sample_cap="10000",
        macro_preset="practical_research",
        seed_spec="cv_seeds=42",
        epochs=int(practical["epochs"]),
        lr=float(practical["lr"]),
        batch_size=int(practical["batch_size"]),
        source_metrics_path="",
        notes="Placeholder row. Fill after the canonical random short-phase winner is chosen by lowest test_mse, then highest test_pearson, then simpler model.",
    )


def _explicit_long_phase_specs() -> list[dict[str, Any]]:
    practical = load_macro_preset("practical_research")
    return [
        _empty_row(
            run_id="long_cv10_10k_no_genes_practical",
            stage_id="1.1",
            status="planned",
            phase="long_cv_10k",
            model_id="no_genes",
            gene_setting="off",
            gene_dim=0,
            split_strategy="random",
            eval_mode="cv10_stratified",
            sample_cap="10000",
            macro_preset="practical_research",
            seed_spec="cv_seeds=42",
            epochs=int(practical["epochs"]),
            lr=float(practical["lr"]),
            batch_size=int(practical["batch_size"]),
            source_metrics_path="outputs/long_cv10_10k_no_genes_practical/cv_metrics.json",
            notes="Explicit long-phase CV row for the no-genes baseline.",
        ),
        _empty_row(
            run_id="long_cv10_10k_pca128_practical",
            stage_id="1.3A",
            status="planned",
            phase="long_cv_10k",
            model_id="pca128",
            gene_setting="pca128",
            gene_dim=128,
            split_strategy="random",
            eval_mode="cv10_stratified",
            sample_cap="10000",
            macro_preset="practical_research",
            seed_spec="cv_seeds=42",
            epochs=int(practical["epochs"]),
            lr=float(practical["lr"]),
            batch_size=int(practical["batch_size"]),
            source_metrics_path="outputs/long_cv10_10k_pca128_practical/cv_metrics.json",
            notes="Explicit long-phase CV row for the pca128 baseline.",
        ),
        _empty_row(
            run_id="long_cv10_10k_no_genes_drug_and_cell_line_practical",
            stage_id="1.1",
            status="planned",
            phase="long_cv_10k",
            model_id="no_genes",
            gene_setting="off",
            gene_dim=0,
            split_strategy="drug_and_cell_line",
            eval_mode="cv10_group_aware",
            sample_cap="10000",
            macro_preset="practical_research",
            seed_spec="cv_seeds=42",
            epochs=int(practical["epochs"]),
            lr=float(practical["lr"]),
            batch_size=int(practical["batch_size"]),
            source_metrics_path="outputs/long_cv10_10k_no_genes_drug_and_cell_line_practical/cv_metrics.json",
            notes="Explicit group-aware long-phase CV row for the no-genes baseline.",
        ),
        _empty_row(
            run_id="long_cv10_10k_pca128_drug_and_cell_line_practical",
            stage_id="1.3A",
            status="planned",
            phase="long_cv_10k",
            model_id="pca128",
            gene_setting="pca128",
            gene_dim=128,
            split_strategy="drug_and_cell_line",
            eval_mode="cv10_group_aware",
            sample_cap="10000",
            macro_preset="practical_research",
            seed_spec="cv_seeds=42",
            epochs=int(practical["epochs"]),
            lr=float(practical["lr"]),
            batch_size=int(practical["batch_size"]),
            source_metrics_path="outputs/long_cv10_10k_pca128_drug_and_cell_line_practical/cv_metrics.json",
            notes="Explicit group-aware long-phase CV row for the pca128 baseline.",
        ),
    ]


def _load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def _sync_single_split_row(row: dict[str, Any]) -> dict[str, Any]:
    if row["source_metrics_path"] == "":
        return row
    metrics_path = BASELINE_DIR / row["source_metrics_path"]
    config_path = metrics_path.with_name("config.json")
    if not metrics_path.exists():
        return row

    metrics = _load_json(metrics_path)
    config = _load_json(config_path) if config_path.exists() else {}

    row["status"] = "done"
    row["gene_dim"] = metrics.get("gene_dim", config.get("gene_dim", row["gene_dim"]))
    row["macro_preset"] = config.get("macro_preset", row["macro_preset"])
    row["seed_spec"] = str(config.get("seed", row["seed_spec"]))
    row["epochs"] = config.get("epochs", row["epochs"])
    row["lr"] = config.get("lr", row["lr"])
    row["batch_size"] = config.get("batch_size", row["batch_size"])
    row["test_mse"] = metrics.get("test_mse", "")
    row["test_rmse"] = metrics.get("test_rmse", "")
    row["test_pearson"] = metrics.get("test_pearson", "")
    row["test_spearman"] = metrics.get("test_spearman", "")
    row["test_mean_baseline_mse"] = metrics.get("test_mean_baseline_mse", "")
    if row["test_mse"] != "" and row["test_mean_baseline_mse"] != "":
        row["improvement_vs_mean_baseline"] = float(row["test_mean_baseline_mse"]) - float(row["test_mse"])
    return row


def _sync_cv_row(row: dict[str, Any]) -> dict[str, Any]:
    if row["source_metrics_path"] == "":
        return row
    metrics_path = BASELINE_DIR / row["source_metrics_path"]
    config_path = metrics_path.with_name("config.json")
    if not metrics_path.exists():
        return row

    metrics = _load_json(metrics_path)
    config = _load_json(config_path) if config_path.exists() else {}
    aggregate = metrics.get("aggregate_metrics", {})

    def mean_value(metric_name: str) -> Any:
        metric_payload = aggregate.get(metric_name, {})
        return metric_payload.get("mean", "")

    row["status"] = "done"
    row["macro_preset"] = config.get("macro_preset", row["macro_preset"])
    row["epochs"] = config.get("epochs", row["epochs"])
    row["lr"] = config.get("lr", row["lr"])
    row["batch_size"] = config.get("batch_size", row["batch_size"])
    row["gene_dim"] = config.get("gene_dim", row["gene_dim"])
    row["test_mse"] = mean_value("test_mse")
    row["test_rmse"] = mean_value("test_rmse")
    row["test_pearson"] = mean_value("test_pearson")
    row["test_spearman"] = mean_value("test_spearman")
    if row["test_mse"] != "":
        row["notes"] = f"{row['notes']} Aggregate metrics reflect the mean across folds."
    return row


def _canonical_random_winner(rows: list[dict[str, Any]]) -> dict[str, Any] | None:
    candidates = [
        row
        for row in rows
        if row["status"] == "done"
        and row["phase"] == "short_single_split"
        and row["split_strategy"] == "random"
        and row["macro_preset"] == "practical_research"
        and row["model_id"] in COMPETITIVE_MODELS
    ]
    if not candidates:
        return None

    return sorted(
        candidates,
        key=lambda row: (
            float(row["test_mse"]),
            -float(row["test_pearson"]),
            TIE_BREAK_PRIORITY[row["model_id"]],
        ),
    )[0]


def _promote_long_phase_row(rows: list[dict[str, Any]]) -> None:
    winner = _canonical_random_winner(rows)
    long_row = next(row for row in rows if row["run_id"] == "long_cv10_10k_promoted_winner")
    if winner is None:
        return

    long_row["stage_id"] = winner["stage_id"]
    long_row["model_id"] = winner["model_id"]
    long_row["gene_setting"] = winner["gene_setting"]
    long_row["gene_dim"] = winner["gene_dim"]
    long_row["source_metrics_path"] = f"outputs/long_cv10_10k_{winner['model_id']}_practical/cv_metrics.json"
    long_row["notes"] = (
        f"Promoted from {winner['run_id']} using the winner rule: lowest test_mse, then highest "
        f"test_pearson, then simpler model."
    )

    synced = _sync_cv_row(long_row)
    long_row.update(synced)


def build_registry() -> list[dict[str, Any]]:
    rows = (
        _legacy_specs()
        + _planned_short_phase_specs()
        + _explicit_long_phase_specs()
        + [_planned_long_phase_placeholder()]
    )

    synced_rows: list[dict[str, Any]] = []
    for row in rows:
        synced = dict(row)
        if row["eval_mode"] in {"cv10_stratified", "cv10_group_aware"}:
            synced = _sync_cv_row(synced)
        else:
            synced = _sync_single_split_row(synced)
        synced_rows.append(synced)

    _promote_long_phase_row(synced_rows)
    return synced_rows


def write_canonical_results(results_dir: Path, rows: list[dict[str, Any]]) -> Path:
    results_dir.mkdir(parents=True, exist_ok=True)
    frame = pd.DataFrame(rows, columns=CANONICAL_COLUMNS)
    frame.to_csv(results_dir / "baseline_experiments.csv", index=False)
    return results_dir / "baseline_experiments.csv"


def _preferred_row(rows: list[dict[str, Any]], model_id: str, split_strategy: str) -> dict[str, Any] | None:
    matching = [
        row
        for row in rows
        if row["status"] == "done"
        and row["model_id"] == model_id
        and row["split_strategy"] == split_strategy
        and row["eval_mode"] == "single_split"
    ]
    if not matching:
        return None

    def priority(row: dict[str, Any]) -> tuple[int, int]:
        canonical = int(
            row["phase"] == "short_single_split"
            and row["macro_preset"] == "practical_research"
            and row["run_id"].startswith("short_")
        )
        return (canonical, 1)

    return sorted(matching, key=priority, reverse=True)[0]


def build_summary(rows: list[dict[str, Any]]) -> pd.DataFrame:
    summary_rows: list[dict[str, Any]] = []
    long_row = next((row for row in rows if row["run_id"] == "long_cv10_10k_promoted_winner" and row["status"] == "done"), None)

    for model_id in ("no_genes", "pca128"):
        summary_row: dict[str, Any] = {"model_id": model_id}
        for split_strategy in ("random", "drug", "cell_line", "drug_and_cell_line"):
            preferred = _preferred_row(rows, model_id=model_id, split_strategy=split_strategy)
            summary_row[f"{split_strategy}_test_mse"] = "" if preferred is None else preferred["test_mse"]
            summary_row[f"{split_strategy}_test_pearson"] = "" if preferred is None else preferred["test_pearson"]

        summary_row["cv10_10k_mse"] = ""
        summary_row["drug_and_cell_line_cv10_10k_mse"] = ""
        summary_row["drug_and_cell_line_cv10_10k_pearson"] = ""
        if long_row is not None and long_row["model_id"] == model_id:
            summary_row["cv10_10k_mse"] = long_row["test_mse"]
        else:
            explicit_long = next(
                (
                    row
                    for row in rows
                    if row["status"] == "done"
                    and row["phase"] == "long_cv_10k"
                    and row["model_id"] == model_id
                    and row["split_strategy"] == "random"
                    and row["run_id"] != "long_cv10_10k_promoted_winner"
                ),
                None,
            )
            if explicit_long is not None:
                summary_row["cv10_10k_mse"] = explicit_long["test_mse"]

        explicit_combined_long = next(
            (
                row
                for row in rows
                if row["status"] == "done"
                and row["phase"] == "long_cv_10k"
                and row["model_id"] == model_id
                and row["split_strategy"] == "drug_and_cell_line"
            ),
            None,
        )
        if explicit_combined_long is not None:
            summary_row["drug_and_cell_line_cv10_10k_mse"] = explicit_combined_long["test_mse"]
            summary_row["drug_and_cell_line_cv10_10k_pearson"] = explicit_combined_long["test_pearson"]

        summary_rows.append(summary_row)

    return pd.DataFrame(summary_rows)


def write_summary(results_dir: Path, rows: list[dict[str, Any]]) -> Path:
    results_dir.mkdir(parents=True, exist_ok=True)
    summary = build_summary(rows)
    summary.to_csv(results_dir / "baseline_summary.csv", index=False)
    return results_dir / "baseline_summary.csv"


def main() -> None:
    args = parse_args()
    results_dir = Path(args.results_dir)
    rows = build_registry()
    canonical_path = write_canonical_results(results_dir, rows)
    summary_path = write_summary(results_dir, rows)
    print(f"Wrote canonical results table to {canonical_path}")
    print(f"Wrote summary table to {summary_path}")


if __name__ == "__main__":
    main()
