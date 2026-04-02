We are working on a drug discovery research project focused on drug synergy prediction using GNNs and gene expression data (TDC DrugComb).

The long-term goal is to innovate the gene expression encoder (cell line representation), specifically improving how high-dimensional gene expression (~20k genes) is compressed into a useful latent representation for downstream synergy prediction.

We are NOT building everything from scratch. We follow a structured approach:

--------------------------------
STEP 1A — Baseline validation
--------------------------------
We use an existing implementation (MatchMaker: https://github.com/tastanlab/matchmaker) to:
- verify the full drug synergy pipeline works
- understand data flow and model structure
- identify where gene expression enters the model

--------------------------------
STEP 1B — Minimal modification
--------------------------------
We introduce a modular encoder:
    z_cell = encoder(gene_expression)

and replace raw gene input with z_cell, while keeping the rest of the pipeline unchanged.

This ensures:
- compatibility with existing models
- controlled experimentation

--------------------------------
LONG-TERM STRATEGY
--------------------------------
The research focus is improving “compression” of gene expression.

In literature, compression typically falls into 3 families:

1. Gene filtering + MLP (DeepSynergy-style)
   - Select subset of genes (variance / prior knowledge)
   - Feed into MLP → latent vector
   - Simple, reproducible, but loses structure

2. Autoencoder-style bottlenecks (AuDNNsynergy)
   - Encoder → bottleneck → decoder
   - Trained to reconstruct input
   - Explicit dimensionality reduction
   - Strong baseline for representation learning

3. Network-informed encoders (PRODeepSyn, TranSynergy)
   - Use biological graphs (PPI, gene-gene interactions)
   - Apply GNN / graph propagation
   - Inject prior biological structure
   - More expressive, harder to implement

--------------------------------
HOW COMPRESSION WORKS IN PRACTICE
--------------------------------
Typical pipeline:

Raw gene expression (~20k dims)
    ↓
Preprocessing:
    - variance filtering (remove low-signal genes)
    - z-score normalization
    ↓
Encoder:
    - MLP OR Autoencoder OR GNN-based
    ↓
Latent vector (z_cell, ~128–256 dims)

This vector is used downstream.

--------------------------------
BEST BASELINE TO START WITH
--------------------------------
Start with:
→ DeepSynergy-style MLP encoder

Why:
- simplest
- fastest to validate
- maximally reproducible on DrugComb

Then compare against:
→ Autoencoder baseline (AuDNNsynergy-style)

--------------------------------
PIPELINE OVERVIEW
--------------------------------

Drug A (SMILES) → encoder → z_A
Drug B (SMILES) → encoder → z_B
Cell line (gene expression) → encoder → z_cell

[z_A, z_B, z_cell]
        ↓
     Fusion (MLP / attention)
        ↓
   Synergy prediction (ZIP)

--------------------------------
KEY RESEARCH GOAL
--------------------------------
We are NOT optimizing the whole pipeline.

We are:
→ isolating and improving the gene expression encoder

Requirements:
- modular (plug into any model)
- improves performance on SOTA pipelines
- preserves biological structure
- avoids excessive information loss

--------------------------------
WHAT TO FOCUS ON
--------------------------------
- representation quality (z_cell)
- compression vs information retention
- leveraging biological priors (networks)
- generalization across models

Avoid:
- over-engineering the pipeline
- modifying drug encoders early
- mixing multiple innovations at once

--------------------------------
MENTAL MODEL
--------------------------------
Baseline = “does it work”
Encoder = “what we improve”
SOTA = “where we prove impact”