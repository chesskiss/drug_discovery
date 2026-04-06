# Research Plan And Milestones

## Goal

Improve the gene expression encoder for drug synergy prediction on TDC DrugComb.

## Strategies
- Use weak baseline and improve/replace compression methods of celline genes. Should see better improvement in compariosn to SOTA models.
- 

## Done

### 1. Data grounding + visualization

- Visualized the dataset 
- Confirmed dataset structure:
  - drug pair inputs from SMILES
  - cell line input from gene expression
  - target = `Synergy_ZIP`
- Confirmed `CellLine` carries 3 views:
  - `23808` raw
  - `3171` filtered
  - `627` pathway-level

### 2. Baseline pipeline

- Implemented a DeepSynergy-style baseline:
  - input = `[drug_a, drug_b, gene_expr]`
  - model = MLP regressor
- Added training + prediction pipeline
- Added deterministic featurization
- Added single-sample prediction CLI
- Current baseline result:
  - `test_mse = 29.0515`
  - `test_rmse = 5.3899`

### 3. Cell-line-only analysis

- Implemented separate `cell_line_difficulty/` module
- Ignored molecule identity and scored cell lines only from observed ZIP values
- Produced ranked outputs for “easier” vs “harder” cell lines
- Current top-level result:
  - easiest examples: `HL-60(TB)`, `MOLT-4`
  - hardest examples: `UO-31`, `HOP-92`

## Next

### 4. Baseline diagnosis

- Re-check baseline quality:
  - compare against naive mean baseline
  - inspect prediction distribution
  - verify normalization choices
- Results:
  - `________________`

### 5. Step 1B: modular gene encoder

- Replace raw `gene_expr` with:
  - `z_cell = encoder(gene_expr)`
- Keep the rest of the model fixed
- First encoder to test:
  - MLP encoder
- Results:
  - `________________`

### 6. Gene filtering idea

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

### 7. Compression baselines

- Compare:
  - simple filtering + MLP
  - autoencoder bottleneck
  - graph-informed encoder
- Keep drug branch and fusion fixed
- Results:
  - `________________`

### 8. Longer-term

- Map gene indices to gene names
- Add biological priors:
  - PPI / transcription networks
- Explore self-supervised pretraining
- Aim for:
  - better performance
  - more biologically meaningful embeddings

