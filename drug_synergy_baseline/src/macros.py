from __future__ import annotations

from pathlib import Path
from typing import Any

import tomllib


DEFAULT_MACRO_PRESET = "practical_research"
DEFAULT_MACRO_FILE = Path(__file__).resolve().parents[1] / "macros.toml"

ALLOWED_MACRO_KEYS = {
    "epochs",
    "batch_size",
    "lr",
    "smiles_dim",
    "hidden_dims",
    "dropout",
    "train_fraction",
    "val_fraction",
    "description",
}


def load_macro_preset(
    preset: str = DEFAULT_MACRO_PRESET,
    macro_file: str | Path = DEFAULT_MACRO_FILE,
) -> dict[str, Any]:
    macro_path = Path(macro_file)
    if not macro_path.exists():
        raise FileNotFoundError(f"Macro file not found: {macro_path}")

    with macro_path.open("rb") as handle:
        parsed = tomllib.load(handle)

    if preset not in parsed:
        available = ", ".join(sorted(parsed))
        raise KeyError(f"Macro preset `{preset}` not found in {macro_path}. Available presets: {available}")

    macro_values = dict(parsed[preset])
    unknown = sorted(set(macro_values).difference(ALLOWED_MACRO_KEYS))
    if unknown:
        raise ValueError(f"Unsupported macro keys in preset `{preset}`: {unknown}")

    return macro_values
