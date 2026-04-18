from __future__ import annotations

import argparse
from pathlib import Path

from .analysis import analyze_cell_line_difficulty, save_analysis_outputs


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Rank DrugComb cell lines by observed synergy outcomes")
    parser.add_argument(
        "--input",
        type=str,
        default="drug_synergy_baseline/data/drugcomb.csv",
        help="Path to DrugComb synergy CSV",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="cell_line_difficulty/outputs",
        help="Directory for ranked outputs",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    metrics = analyze_cell_line_difficulty(args.input)
    save_analysis_outputs(metrics, input_path=args.input, output_dir=args.output_dir)

    print(f"Analyzed {len(metrics)} cell lines")
    print(f"Top easiest: {metrics.iloc[0]['cell_line']} (ease_score={metrics.iloc[0]['ease_score']:.4f})")
    print(
        f"Top hardest: {metrics.sort_values('difficulty_rank').iloc[0]['cell_line']} "
        f"(difficulty_score={metrics.sort_values('difficulty_rank').iloc[0]['difficulty_score']:.4f})"
    )
    print(f"Saved outputs to {Path(args.output_dir)}")


if __name__ == "__main__":
    main()
