from __future__ import annotations

import argparse
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Entry point for baseline reverse-engineering and attribution workflows."
    )
    parser.add_argument(
        "--model-path",
        type=str,
        default="drug_synergy_baseline/outputs/pca128_random/baseline_mlp.pt",
        help="Path to a trained single-run checkpoint to analyze.",
    )
    parser.add_argument(
        "--config-path",
        type=str,
        default="drug_synergy_baseline/outputs/pca128_random/config.json",
        help="Path to the matching training config.",
    )
    parser.add_argument(
        "--predictions-path",
        type=str,
        default=None,
        help="Optional saved predictions file to align with attribution outputs later.",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="reverse_engineering/outputs",
        help="Directory for reverse-engineering outputs.",
    )
    parser.add_argument(
        "--method",
        choices=["placeholder", "occlusion"],
        default="placeholder",
        help="Attribution method to run. Only placeholder wiring exists currently.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    print("reverse_engineering module scaffold is ready.")
    print(f"method={args.method}")
    print(f"model_path={Path(args.model_path)}")
    print(f"config_path={Path(args.config_path)}")
    if args.predictions_path is not None:
        print(f"predictions_path={Path(args.predictions_path)}")
    print(f"output_dir={Path(args.output_dir)}")
    print("Next step: implement PCA-component occlusion and summary export in this module.")


if __name__ == "__main__":
    main()
