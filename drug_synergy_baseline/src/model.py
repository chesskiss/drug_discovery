from __future__ import annotations

import torch
from torch import nn


class DeepSynergyMLP(nn.Module):
    def __init__(
        self,
        drug_dim: int,
        gene_dim: int,
        hidden_dims: tuple[int, ...] = (512, 256),
        dropout: float = 0.2,
    ) -> None:
        super().__init__()

        input_dim = (2 * drug_dim) + gene_dim
        layers: list[nn.Module] = []
        current_dim = input_dim

        for hidden_dim in hidden_dims:
            layers.append(nn.Linear(current_dim, hidden_dim))
            layers.append(nn.ReLU())
            if dropout > 0:
                layers.append(nn.Dropout(dropout))
            current_dim = hidden_dim

        layers.append(nn.Linear(current_dim, 1))
        self.network = nn.Sequential(*layers)

    def forward(
        self,
        drug_a: torch.Tensor,
        drug_b: torch.Tensor,
        gene_expr: torch.Tensor,
    ) -> torch.Tensor:
        features = torch.cat([drug_a, drug_b, gene_expr], dim=-1)
        return self.network(features).squeeze(-1)
