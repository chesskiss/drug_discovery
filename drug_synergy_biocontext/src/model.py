from __future__ import annotations

import torch
from torch import nn


class BioContextSynergyMLP(nn.Module):
    """DeepSynergy-style MLP whose cell-line branch is a pathway-activity encoder.

    Pathway activities (`z = W @ x`, already projected and z-scored in the data layer)
    pass through a trainable compression head, MLP0, before joining the drug features:

        concat[ drug_a, drug_b, MLP0(pathways) ] -> hidden layers -> synergy

    `mlp0_out_dim == 0` disables the cell-line branch entirely (the no-genes control),
    leaving the model with only the two drug fingerprints.
    """

    def __init__(
        self,
        drug_dim: int,
        pathway_dim: int,
        mlp0_out_dim: int,
        hidden_dims: tuple[int, ...] = (512, 256),
        dropout: float = 0.2,
    ) -> None:
        super().__init__()

        self.pathway_dim = pathway_dim
        self.mlp0_out_dim = mlp0_out_dim
        self.use_pathways = mlp0_out_dim > 0 and pathway_dim > 0

        if self.use_pathways:
            mlp0_layers: list[nn.Module] = [nn.Linear(pathway_dim, mlp0_out_dim), nn.ReLU()]
            if dropout > 0:
                mlp0_layers.append(nn.Dropout(dropout))
            self.mlp0: nn.Module | None = nn.Sequential(*mlp0_layers)
        else:
            self.mlp0 = None

        input_dim = (2 * drug_dim) + (mlp0_out_dim if self.use_pathways else 0)
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
        if self.use_pathways and self.mlp0 is not None:
            features = torch.cat([drug_a, drug_b, self.mlp0(gene_expr)], dim=-1)
        else:
            features = torch.cat([drug_a, drug_b], dim=-1)
        return self.network(features).squeeze(-1)


def build_biocontext_model(
    drug_dim: int,
    pathway_dim: int,
    mlp0_out_dim: int,
    hidden_dims: list[int] | tuple[int, ...],
    dropout: float,
) -> BioContextSynergyMLP:
    return BioContextSynergyMLP(
        drug_dim=drug_dim,
        pathway_dim=pathway_dim,
        mlp0_out_dim=mlp0_out_dim,
        hidden_dims=tuple(hidden_dims),
        dropout=dropout,
    )
