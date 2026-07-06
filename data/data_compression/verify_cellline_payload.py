from __future__ import annotations

import argparse
import csv
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Check whether a DrugComb export contains truncated `CellLine` payloads."
    )
    parser.add_argument(
        "--input",
        type=str,
        required=True,
        help="Path to a CSV file that may contain a `CellLine` column.",
    )
    parser.add_argument(
        "--sample-rows",
        type=int,
        default=20,
        help="How many rows to inspect before summarizing.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    input_path = Path(args.input)

    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    inspected_rows = 0
    ellipsis_rows = 0
    missing_column = False

    with input_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None or "CellLine" not in reader.fieldnames:
            missing_column = True
        else:
            for row in reader:
                payload = row.get("CellLine", "")
                inspected_rows += 1
                if "..." in payload:
                    ellipsis_rows += 1
                if inspected_rows >= args.sample_rows:
                    break

    print(f"file={input_path}")
    if missing_column:
        print("status=no_cellline_column")
        return

    print(f"sample_rows={inspected_rows}")
    print(f"rows_with_ellipsis={ellipsis_rows}")
    if ellipsis_rows > 0:
        print("status=truncated")
        print("message=CellLine payload contains ellipsis, so the full arrays are not preserved in this CSV.")
    else:
        print("status=not_truncated_in_sample")
        print("message=No ellipsis detected in the sampled rows.")


if __name__ == "__main__":
    main()
