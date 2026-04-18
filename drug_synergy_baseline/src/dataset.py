from __future__ import annotations

import hashlib
from dataclasses import dataclass
from functools import lru_cache

import numpy as np
import pandas as pd
import torch
from torch.utils.data import Dataset

from .data_loading import load_expression_lookup, load_synergy_table


def _stable_bucket(token: str, dim: int) -> int:
    digest = hashlib.blake2b(token.encode("utf-8"), digest_size=8).digest()
    return int.from_bytes(digest, byteorder="big") % dim


def smiles_to_vector(smiles: str, dim: int = 256) -> np.ndarray:
    vector = np.zeros(dim, dtype=np.float32)
    if not smiles:
        return vector

    chars = list(smiles)
    for idx, ch in enumerate(chars):
        vector[(ord(ch) + idx) % dim] += 1.0
        if idx + 1 < len(chars):
            pair_token = f"{idx}:{ch}{chars[idx + 1]}"
            vector[_stable_bucket(pair_token, dim)] += 0.5

    norm = np.linalg.norm(vector)
    if norm > 0:
        vector /= norm
    return vector


@dataclass(frozen=True)
class DatasetBundle:
    train: Dataset
    val: Dataset
    test: Dataset
    gene_dim: int
    drug_dim: int
    train_rows: pd.DataFrame
    val_rows: pd.DataFrame
    test_rows: pd.DataFrame


def load_aligned_synergy_and_expression(
    synergy_path: str,
    cell_expression_path: str | None,
    fallback_pickle_path: str | None,
    *,
    cell_feature_view: int = 0,
    use_gene_expression: bool = True,
    max_samples: int | None = None,
) -> tuple[pd.DataFrame, dict[str, np.ndarray], int]:
    synergy_df = load_synergy_table(synergy_path)

    if use_gene_expression:
        expression_lookup = load_expression_lookup(
            cell_expression_path=cell_expression_path,
            fallback_pickle_path=fallback_pickle_path,
            feature_view_index=cell_feature_view,
        )
        synergy_df = synergy_df[synergy_df["cell_line"].isin(expression_lookup)].copy()
        if not expression_lookup:
            raise ValueError("Expression lookup is empty after loading gene-expression features.")
        gene_dim = len(next(iter(expression_lookup.values())))
    else:
        expression_lookup = {}
        gene_dim = 0

    if max_samples is not None:
        synergy_df = synergy_df.head(max_samples).copy()

    if synergy_df.empty:
        raise ValueError("No training rows remain after aligning cell lines with expression features.")

    return synergy_df.reset_index(drop=True), expression_lookup, gene_dim


class DrugSynergyDataset(Dataset):
    def __init__(
        self,
        rows,
        expression_lookup: dict[str, np.ndarray],
        smiles_dim: int = 256,
        gene_dim: int | None = None,
    ) -> None:
        self.rows = rows.reset_index(drop=True)
        self.expression_lookup = expression_lookup
        self.smiles_dim = smiles_dim
        if gene_dim is not None:
            self.gene_dim = gene_dim
        elif expression_lookup:
            self.gene_dim = len(next(iter(expression_lookup.values())))
        else:
            self.gene_dim = 0

    def __len__(self) -> int:
        return len(self.rows)

    @lru_cache(maxsize=16384)
    def _cached_smiles_features(self, smiles: str) -> np.ndarray:
        return smiles_to_vector(smiles, dim=self.smiles_dim)

    def __getitem__(self, index: int) -> dict[str, torch.Tensor]:
        row = self.rows.iloc[index]
        cell_line = row["cell_line"]
        if self.gene_dim == 0:
            gene_expr = np.zeros(0, dtype=np.float32)
        else:
            gene_expr = self.expression_lookup[cell_line]

        return {
            "drug_a": torch.from_numpy(self._cached_smiles_features(row["smiles_a"])).float(),
            "drug_b": torch.from_numpy(self._cached_smiles_features(row["smiles_b"])).float(),
            "gene_expr": torch.from_numpy(gene_expr).float(),
            "target": torch.tensor(row["target"], dtype=torch.float32),
        }


def _split_unique_values(
    values: list[str],
    *,
    train_fraction: float,
    val_fraction: float,
    random_seed: int,
) -> tuple[set[str], set[str], set[str]]:
    if not values:
        raise ValueError("Cannot split an empty group set.")

    unique_values = np.asarray(sorted(set(values)))
    rng = np.random.default_rng(random_seed)
    rng.shuffle(unique_values)

    n_groups = len(unique_values)
    train_count = max(1, int(n_groups * train_fraction))
    val_count = max(1, int(n_groups * val_fraction)) if n_groups >= 3 else 0

    if train_count + val_count >= n_groups:
        val_count = max(0, min(val_count, n_groups - train_count - 1))

    test_count = n_groups - train_count - val_count
    if test_count <= 0:
        if val_count > 0:
            val_count -= 1
        elif train_count > 1:
            train_count -= 1
        test_count = n_groups - train_count - val_count

    train_values = set(unique_values[:train_count].tolist())
    val_values = set(unique_values[train_count : train_count + val_count].tolist())
    test_values = set(unique_values[train_count + val_count :].tolist())
    return train_values, val_values, test_values


def split_unique_values_into_folds(
    values: list[str],
    *,
    num_folds: int,
    random_seed: int,
) -> list[set[str]]:
    if num_folds < 2:
        raise ValueError("Group-aware cross-validation requires at least 2 folds.")
    if not values:
        raise ValueError("Cannot split an empty group set.")

    unique_values = np.asarray(sorted(set(values)))
    if len(unique_values) < num_folds:
        raise ValueError(f"Cannot build {num_folds} folds from only {len(unique_values)} unique groups.")

    rng = np.random.default_rng(random_seed)
    rng.shuffle(unique_values)

    fold_buckets: list[list[str]] = [[] for _ in range(num_folds)]
    for index, value in enumerate(unique_values):
        fold_buckets[index % num_folds].append(str(value))

    return [set(bucket) for bucket in fold_buckets]


def split_drug_and_cell_line_priority_rows(
    synergy_df: pd.DataFrame,
    *,
    val_drugs: set[str],
    test_drugs: set[str],
    val_cells: set[str],
    test_cells: set[str],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    working_df = synergy_df.copy()

    def classify_joint_row(row: pd.Series) -> str:
        drug_a = str(row["smiles_a"])
        drug_b = str(row["smiles_b"])
        cell_line = str(row["cell_line"])
        if drug_a in test_drugs or drug_b in test_drugs or cell_line in test_cells:
            return "test"
        if drug_a in val_drugs or drug_b in val_drugs or cell_line in val_cells:
            return "val"
        return "train"

    working_df["split"] = working_df.apply(classify_joint_row, axis=1)
    train_rows = working_df[working_df["split"] == "train"].drop(columns=["split"]).copy()
    val_rows = working_df[working_df["split"] == "val"].drop(columns=["split"]).copy()
    test_rows = working_df[working_df["split"] == "test"].drop(columns=["split"]).copy()
    return train_rows, val_rows, test_rows


def split_synergy_rows(
    synergy_df: pd.DataFrame,
    *,
    split_strategy: str = "random",
    train_fraction: float = 0.8,
    val_fraction: float = 0.1,
    random_seed: int = 42,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    if synergy_df.empty:
        raise ValueError("Cannot split an empty synergy table.")

    if split_strategy == "random":
        rng = np.random.default_rng(random_seed)
        indices = np.arange(len(synergy_df))
        rng.shuffle(indices)

        train_end = int(len(indices) * train_fraction)
        val_end = train_end + int(len(indices) * val_fraction)

        train_rows = synergy_df.iloc[indices[:train_end]].copy()
        val_rows = synergy_df.iloc[indices[train_end:val_end]].copy()
        test_rows = synergy_df.iloc[indices[val_end:]].copy()
    elif split_strategy == "cell_line":
        train_cells, val_cells, test_cells = _split_unique_values(
            synergy_df["cell_line"].astype(str).tolist(),
            train_fraction=train_fraction,
            val_fraction=val_fraction,
            random_seed=random_seed,
        )
        train_rows = synergy_df[synergy_df["cell_line"].isin(train_cells)].copy()
        val_rows = synergy_df[synergy_df["cell_line"].isin(val_cells)].copy()
        test_rows = synergy_df[synergy_df["cell_line"].isin(test_cells)].copy()
    elif split_strategy == "drug_and_cell_line":
        all_drugs = sorted(set(synergy_df["smiles_a"].astype(str)) | set(synergy_df["smiles_b"].astype(str)))
        _, val_drugs, test_drugs = _split_unique_values(
            all_drugs,
            train_fraction=train_fraction,
            val_fraction=val_fraction,
            random_seed=random_seed,
        )
        _, val_cells, test_cells = _split_unique_values(
            synergy_df["cell_line"].astype(str).tolist(),
            train_fraction=train_fraction,
            val_fraction=val_fraction,
            random_seed=random_seed,
        )
        train_rows, val_rows, test_rows = split_drug_and_cell_line_priority_rows(
            synergy_df,
            val_drugs=val_drugs,
            test_drugs=test_drugs,
            val_cells=val_cells,
            test_cells=test_cells,
        )
    elif split_strategy in {"drug", "drug_pair"}:
        working_df = synergy_df.copy()
        if split_strategy == "drug":
            all_drugs = sorted(set(working_df["smiles_a"].astype(str)) | set(working_df["smiles_b"].astype(str)))
            train_drugs, val_drugs, test_drugs = _split_unique_values(
                all_drugs,
                train_fraction=train_fraction,
                val_fraction=val_fraction,
                random_seed=random_seed,
            )

            def classify_row(row: pd.Series) -> str | None:
                drug_a = str(row["smiles_a"])
                drug_b = str(row["smiles_b"])
                if drug_a in train_drugs and drug_b in train_drugs:
                    return "train"
                if drug_a in val_drugs and drug_b in val_drugs:
                    return "val"
                if drug_a in test_drugs and drug_b in test_drugs:
                    return "test"
                return None

            working_df["split"] = working_df.apply(classify_row, axis=1)
        else:
            pair_tokens = working_df.apply(
                lambda row: "||".join(sorted((str(row["smiles_a"]), str(row["smiles_b"])))),
                axis=1,
            )
            train_pairs, val_pairs, test_pairs = _split_unique_values(
                pair_tokens.tolist(),
                train_fraction=train_fraction,
                val_fraction=val_fraction,
                random_seed=random_seed,
            )
            working_df["pair_token"] = pair_tokens
            working_df["split"] = working_df["pair_token"].map(
                lambda token: "train"
                if token in train_pairs
                else "val"
                if token in val_pairs
                else "test"
                if token in test_pairs
                else None
            )

        train_rows = working_df[working_df["split"] == "train"].drop(
            columns=[col for col in ("split", "pair_token") if col in working_df.columns]
        ).copy()
        val_rows = working_df[working_df["split"] == "val"].drop(
            columns=[col for col in ("split", "pair_token") if col in working_df.columns]
        ).copy()
        test_rows = working_df[working_df["split"] == "test"].drop(
            columns=[col for col in ("split", "pair_token") if col in working_df.columns]
        ).copy()
    else:
        raise ValueError(f"Unsupported split strategy: {split_strategy}")

    if train_rows.empty or val_rows.empty or test_rows.empty:
        raise ValueError(
            f"Split strategy `{split_strategy}` produced an empty split "
            f"(train={len(train_rows)}, val={len(val_rows)}, test={len(test_rows)})."
        )

    return train_rows.reset_index(drop=True), val_rows.reset_index(drop=True), test_rows.reset_index(drop=True)


def build_dataset_bundle_from_rows(
    train_rows: pd.DataFrame,
    val_rows: pd.DataFrame,
    test_rows: pd.DataFrame,
    expression_lookup: dict[str, np.ndarray],
    *,
    smiles_dim: int = 256,
    gene_dim: int = 0,
) -> DatasetBundle:
    train_dataset = DrugSynergyDataset(train_rows, expression_lookup, smiles_dim=smiles_dim, gene_dim=gene_dim)
    val_dataset = DrugSynergyDataset(val_rows, expression_lookup, smiles_dim=smiles_dim, gene_dim=gene_dim)
    test_dataset = DrugSynergyDataset(test_rows, expression_lookup, smiles_dim=smiles_dim, gene_dim=gene_dim)

    return DatasetBundle(
        train=train_dataset,
        val=val_dataset,
        test=test_dataset,
        gene_dim=gene_dim,
        drug_dim=smiles_dim,
        train_rows=train_rows.reset_index(drop=True),
        val_rows=val_rows.reset_index(drop=True),
        test_rows=test_rows.reset_index(drop=True),
    )


def build_datasets(
    synergy_path: str,
    cell_expression_path: str | None,
    fallback_pickle_path: str | None,
    *,
    use_gene_expression: bool = True,
    cell_feature_view: int = 0,
    split_strategy: str = "random",
    smiles_dim: int = 256,
    train_fraction: float = 0.8,
    val_fraction: float = 0.1,
    random_seed: int = 42,
    max_samples: int | None = None,
) -> DatasetBundle:
    synergy_df, expression_lookup, gene_dim = load_aligned_synergy_and_expression(
        synergy_path=synergy_path,
        cell_expression_path=cell_expression_path,
        fallback_pickle_path=fallback_pickle_path,
        cell_feature_view=cell_feature_view,
        use_gene_expression=use_gene_expression,
        max_samples=max_samples,
    )
    train_rows, val_rows, test_rows = split_synergy_rows(
        synergy_df,
        split_strategy=split_strategy,
        train_fraction=train_fraction,
        val_fraction=val_fraction,
        random_seed=random_seed,
    )
    return build_dataset_bundle_from_rows(
        train_rows=train_rows,
        val_rows=val_rows,
        test_rows=test_rows,
        expression_lookup=expression_lookup,
        smiles_dim=smiles_dim,
        gene_dim=gene_dim,
    )
