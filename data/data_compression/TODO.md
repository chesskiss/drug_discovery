# Data Compression TODO

## Planned Methods

- `z-score -> Var3k -> PCA512 -> AE128`: best practical upgrade.
- `z-score -> BioFilter1k-3k -> AE128`: interpretable, more fragile.
- `raw -> AE128`: hardest to stabilize; avoid first.
