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


class DrugSynergyDataset(Dataset):
    def __init__(
        self,
        rows,
        expression_lookup: dict[str, np.ndarray],
        smiles_dim: int = 256,
    ) -> None:
        self.rows = rows.reset_index(drop=True)
        self.expression_lookup = expression_lookup
        self.smiles_dim = smiles_dim
        self.gene_dim = len(next(iter(expression_lookup.values())))

    def __len__(self) -> int:
        return len(self.rows)

    @lru_cache(maxsize=16384)
    def _cached_smiles_features(self, smiles: str) -> np.ndarray:
        return smiles_to_vector(smiles, dim=self.smiles_dim)

    def __getitem__(self, index: int) -> dict[str, torch.Tensor]:
        row = self.rows.iloc[index]
        cell_line = row["cell_line"]
        gene_expr = self.expression_lookup[cell_line]

        return {
            "drug_a": torch.from_numpy(self._cached_smiles_features(row["smiles_a"])).float(),
            "drug_b": torch.from_numpy(self._cached_smiles_features(row["smiles_b"])).float(),
            "gene_expr": torch.from_numpy(gene_expr).float(),
            "target": torch.tensor(row["target"], dtype=torch.float32),
        }


def build_datasets(
    synergy_path: str,
    cell_expression_path: str | None,
    fallback_pickle_path: str | None,
    smiles_dim: int = 256,
    train_fraction: float = 0.8,
    val_fraction: float = 0.1,
    random_seed: int = 42,
    max_samples: int | None = None,
) -> DatasetBundle:
    synergy_df = load_synergy_table(synergy_path)
    expression_lookup = load_expression_lookup(
        cell_expression_path=cell_expression_path,
        fallback_pickle_path=fallback_pickle_path,
    )

    synergy_df = synergy_df[synergy_df["cell_line"].isin(expression_lookup)].copy()
    if max_samples is not None:
        synergy_df = synergy_df.head(max_samples).copy()

    if synergy_df.empty:
        raise ValueError("No training rows remain after aligning cell lines with expression features.")

    rng = np.random.default_rng(random_seed)
    indices = np.arange(len(synergy_df))
    rng.shuffle(indices)

    train_end = int(len(indices) * train_fraction)
    val_end = train_end + int(len(indices) * val_fraction)

    train_rows = synergy_df.iloc[indices[:train_end]].copy()
    val_rows = synergy_df.iloc[indices[train_end:val_end]].copy()
    test_rows = synergy_df.iloc[indices[val_end:]].copy()

    train_dataset = DrugSynergyDataset(train_rows, expression_lookup, smiles_dim=smiles_dim)
    val_dataset = DrugSynergyDataset(val_rows, expression_lookup, smiles_dim=smiles_dim)
    test_dataset = DrugSynergyDataset(test_rows, expression_lookup, smiles_dim=smiles_dim)

    return DatasetBundle(
        train=train_dataset,
        val=val_dataset,
        test=test_dataset,
        gene_dim=train_dataset.gene_dim,
        drug_dim=smiles_dim,
        train_rows=train_rows.reset_index(drop=True),
        val_rows=val_rows.reset_index(drop=True),
        test_rows=test_rows.reset_index(drop=True),
    )
