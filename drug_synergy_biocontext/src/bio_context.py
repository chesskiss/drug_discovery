"""Biological-context projection of cell-line gene expression.

Replaces the PCA compression used by the baseline. A cell line's raw 23,808-dim
expression vector `x` is projected onto pathway activities with a FIXED, biologically
derived weight matrix `W` (`[n_pathways x 23808]`):

    z = W @ x

The projection is deterministic and depends only on the cell line, so it is computed
once per unique cell line (there are only 59) rather than per row.

Unlike PCA, `W` is not fitted to this dataset and its rows are named, interpretable
pathways -- which is the whole point of the experiment: PCA at full rank is a lossless
fingerprint of the 59 cell lines (it can encode cell-line identity), whereas a small
fixed pathway basis cannot.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

# Repo-root-relative location of the aligned weight matrices.
BIO_CONTEXT_ROOT = Path(__file__).resolve().parents[2] / "data" / "bio_context"

BIO_CONTEXT_SOURCES = {
    "progeny": BIO_CONTEXT_ROOT / "progeny" / "data" / "progeny_tdc_weights.npz",
    "kegg": BIO_CONTEXT_ROOT / "kegg" / "data" / "kegg_tdc_weights.npz",
}

# Expected gene axis (TDC / NCI-60 RNA-seq composite); both matrices share it.
EXPECTED_GENE_AXIS = 23808


def _load_single(source: str, *, drop_dead_rows: bool = True) -> tuple[np.ndarray, list[str]]:
    """Load one weight matrix as (W, pathway_names).

    `drop_dead_rows` removes all-zero rows. KEGG ships 36 edgeless membership diagrams
    (Ribosome, ABC transporters, ...) whose weights are identically zero; keeping them
    would feed the model dead, always-0.0 inputs.
    """
    path = BIO_CONTEXT_SOURCES[source]
    if not path.exists():
        raise FileNotFoundError(f"Bio-context matrix not found for '{source}': {path}")

    payload = np.load(path, allow_pickle=True)
    weights = np.asarray(payload["W"], dtype=np.float32)
    if weights.ndim != 2 or weights.shape[1] != EXPECTED_GENE_AXIS:
        raise ValueError(
            f"'{source}' matrix has shape {weights.shape}; expected [n_pathways x {EXPECTED_GENE_AXIS}]."
        )

    # KEGG carries human-readable names alongside the ids; PROGENy names are the ids.
    if "pathway_names" in payload.files:
        names = [f"{source}:{n}" for n in np.asarray(payload["pathway_names"]).tolist()]
    else:
        names = [f"{source}:{n}" for n in np.asarray(payload["pathways"]).tolist()]

    if drop_dead_rows:
        keep = np.abs(weights).sum(axis=1) > 0
        dropped = int((~keep).sum())
        if dropped:
            print(f"[bio_context] {source}: dropped {dropped} all-zero pathway rows (dead inputs).")
            weights = weights[keep]
            names = [n for n, k in zip(names, keep, strict=False) if k]

    return weights, names


def load_bio_context_matrix(source: str, *, drop_dead_rows: bool = True) -> tuple[np.ndarray, list[str]]:
    """Load the projection matrix for `progeny`, `kegg`, or `progeny_kegg` (both stacked).

    Stacking is a plain row-concat because both matrices are aligned to the identical
    23,808-gene axis.
    """
    keys = source.split("_") if source == "progeny_kegg" else [source]
    unknown = [k for k in keys if k not in BIO_CONTEXT_SOURCES]
    if unknown:
        raise ValueError(
            f"Unknown bio-context source(s) {unknown}. Valid: progeny, kegg, progeny_kegg."
        )

    matrices, names = [], []
    for key in keys:
        weights, pathway_names = _load_single(key, drop_dead_rows=drop_dead_rows)
        matrices.append(weights)
        names.extend(pathway_names)

    stacked = np.vstack(matrices) if len(matrices) > 1 else matrices[0]
    print(f"[bio_context] source={source} -> pathway_dim={stacked.shape[0]}")
    return stacked, names


def project_expression(
    expression_lookup: dict[str, np.ndarray],
    weights: np.ndarray,
) -> dict[str, np.ndarray]:
    """Map {cell_line: raw 23808-dim vector} -> {cell_line: n_pathway activity vector}."""
    projected: dict[str, np.ndarray] = {}
    for cell_line, vector in expression_lookup.items():
        x = np.asarray(vector, dtype=np.float32)
        if x.shape[0] != weights.shape[1]:
            raise ValueError(
                f"Cell line '{cell_line}' has {x.shape[0]} genes; "
                f"the bio-context matrix expects {weights.shape[1]}. "
                "Raw view 0 (23808-dim) is required."
            )
        projected[cell_line] = (weights @ x).astype(np.float32)
    return projected


def fit_pathway_normalizer(
    expression_lookup: dict[str, np.ndarray],
    train_cell_lines: set[str],
) -> tuple[np.ndarray, np.ndarray]:
    """Fit per-pathway z-score stats using ONLY training cell lines.

    Fitting on train cell lines alone keeps the cold-cell-line splits honest: test
    cell lines must not influence the feature scaling. Statistics are unweighted over
    unique cell lines (each cell line is one observation of the pathway distribution),
    not row-weighted, so frequently-assayed cell lines do not dominate the centring.
    """
    vectors = [expression_lookup[c] for c in sorted(train_cell_lines) if c in expression_lookup]
    if not vectors:
        raise ValueError("No training cell lines available to fit the pathway normalizer.")

    matrix = np.vstack(vectors).astype(np.float32)
    mean = matrix.mean(axis=0)
    std = matrix.std(axis=0)
    std[std == 0] = 1.0  # constant pathway -> leave it centred at 0
    return mean.astype(np.float32), std.astype(np.float32)


def apply_pathway_normalizer(
    expression_lookup: dict[str, np.ndarray],
    mean: np.ndarray,
    std: np.ndarray,
) -> dict[str, np.ndarray]:
    """Apply fitted z-score stats to every cell line (train, val and test alike)."""
    return {
        cell_line: ((np.asarray(v, dtype=np.float32) - mean) / std).astype(np.float32)
        for cell_line, v in expression_lookup.items()
    }
