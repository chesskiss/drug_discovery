# Research Plan And Milestones

## Goal

Improve the gene expression encoder for drug synergy prediction on TDC DrugComb.

## Stage Map

- `Stage 1.0`: data grounding + initial baseline pipeline
- `Stage 1.1`: baseline variant diagnostics
- `Stage 1.2`: modular gene encoder replacement
- `Stage 1.3`: compression baseline comparisons
- `Stage 1.3A`: first explicit compression export with `z-score -> variance top-k -> PCA128`

## Strategies
- Use weak baseline and improve/replace compression methods of celline genes. Should see better improvement in compariosn to SOTA models. 
- Naive compression (627) vs my compression plan - then add to weak baseline

## Done

### Stage 1.0A. Data grounding + visualization

- Visualized the dataset 
- Confirmed dataset structure:
  - drug pair inputs from SMILES
  - cell line input from gene expression
  - target = `Synergy_ZIP`
- Confirmed `CellLine` carries 3 views:
  - `23808` raw
  - `3171` filtered
  - `627` pathway-level

### Stage 1.0B. Baseline pipeline

- Implemented a DeepSynergy-style baseline:
  - input = `[drug_a, drug_b, gene_expr]`
  - model = MLP regressor
- Added training + prediction pipeline
- Added deterministic featurization
- Added single-sample prediction CLI
- Current baseline result:
  - `test_mse = 29.0515`
  - `test_rmse = 5.3899`

### Stage 1.0C. Cell-line-only analysis

- Implemented separate `cell_line_difficulty/` module
- Ignored molecule identity and scored cell lines only from observed ZIP values
- Produced ranked outputs for “easier” vs “harder” cell lines
- Current top-level result:
  - easiest examples: `HL-60(TB)`, `MOLT-4`
  - hardest examples: `UO-31`, `HOP-92`

### Stage 1.0D. Cell-line-only predictive model

- Built a predictive dataset with:
  - `59` cell lines
  - `3171` filtered cell-line features
  - target = composite `ease_score`
- Ran LOOCV:
  - ridge `rmse = 0.7631`
  - tiny MLP `rmse = 0.8370`
  - fold-mean baseline `rmse = 0.9006`
- Brief take:
  - ridge is best so far
  - there is predictive signal in the cell-line features
  - performance is still moderate, not strong

## Next

### Stage 1.1. Baseline variant diagnostics

- Re-check baseline quality with explicit baseline variants:
  - compare genes on vs genes off
  - compare `CellLine[0]`, `CellLine[1]`, `CellLine[2]`
  - compare split strategies:
    - `random`
    - `cell_line`
    - `drug`
    - `drug_pair`
  - run `10`-fold CV, ideally stratified, on `10K` first and then full data
  - verify training stability under smaller LR / modest epoch count
- Results:
  - `________________`

### Stage 1.2. Modular gene encoder

- Replace raw `gene_expr` with:
  - `z_cell = encoder(gene_expr)`
- Keep the rest of the model fixed
- First encoder to test:
  - MLP encoder
- Results:
  - `________________`

### Stage 1.2A. Gene filtering idea

- Test filtering genes with value `> 2.3`
- Keep both:
  - expressed values
  - their original indices
- Candidate representations:
  - sparse index-value table
  - hash-based structure
  - graph-ready structure
- Goal:
  - reduce noise without losing gene identity
- Results:
  - `________________`

### Stage 1.3. Compression baselines

- Compare:
  - simple filtering + MLP
  - autoencoder bottleneck
  - graph-informed encoder
- Keep drug branch and fusion fixed
- Results:
  - `________________`

### Stage 1.3A. z-score -> variance top-k -> PCA128

- New root-level module: `data_compression/`
- Current pipeline dir:
  - `data_compression/zscore_var_pca_128/`
- Data layout:
  - raw CSV moved to `data_compression/zscore_var_pca_128/data/drugcomb_raw.csv`
  - compressed export written to `data_compression/zscore_var_pca_128/data/drugcomb.csv`
- Technical note:
  - the raw CSV is kept for provenance
  - actual raw `CellLine[0]` vectors are read from `drug_synergy_baseline/data/drugcomb.pkl`
  - because DrugComb has only `59` unique cell lines here, PCA rank is limited and the 128-d output is padded beyond the effective rank
- Result:
  - `________________`

### Stage 2. Longer-term

- Map gene indices to gene names
- Add biological priors:
  - PPI / transcription networks
- Explore self-supervised pretraining
- Aim for:
  - better performance
  - more biologically meaningful embeddings
