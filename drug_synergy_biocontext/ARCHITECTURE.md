# Baseline Architecture

## Current Model

The current baseline is a plain DeepSynergy-style MLP.

It is **not**:

- an attention model
- a graph model
- an autoencoder model
- a model with a learned cell/gene encoder

The current cell branch is effectively an **identity pass-through**:

- the cell-line vector is fed directly into the MLP
- `cell_encoder_type = "identity"` in saved `config.json` files is metadata only
- it means no learned encoder exists yet

Relevant implementation:

- [src/model.py](/Users/arnoldcheskis/Documents/Projects/drug_discovery/drug_synergy_baseline/src/model.py)
- [src/train.py](/Users/arnoldcheskis/Documents/Projects/drug_discovery/drug_synergy_baseline/src/train.py)

## Input Blocks

The model concatenates 3 inputs:

1. drug A vector
2. drug B vector
3. cell-line / gene-expression vector

Drug vectors:

- built from SMILES
- hashed/fixed-length representation
- default dim: `256`

Cell-line vector:

- depends on the run
- examples:
  - raw genes: `23808`
  - filtered genes: `3171`
  - compact view: `627`
  - PCA-compressed genes: `128`
  - no-genes model: `0`

## Block View

General form:

```text
Drug A SMILES
   -> 256-d vector

Drug B SMILES
   -> 256-d vector

Cell-line input
   -> gene vector / PCA vector / empty

[drug_a | drug_b | cell_expr]
          |
          v
   Linear(input_dim -> hidden_1)
          |
         ReLU
          |
       Dropout
          |
          v
   Linear(hidden_1 -> hidden_2)
          |
         ReLU
          |
       Dropout
          |
          v
   Linear(hidden_2 -> 1)
          |
          v
   Predicted Synergy_ZIP
```

## Practical Research Preset

For the current `practical_research` preset:

- hidden dims: `[1024, 512]`
- dropout: `0.3`
- output: scalar regression value

So the current baseline is:

```text
concatenate inputs
-> Linear(input_dim -> 1024)
-> ReLU
-> Dropout(0.3)
-> Linear(1024 -> 512)
-> ReLU
-> Dropout(0.3)
-> Linear(512 -> 1)
```

## Example Input Sizes

### No Genes

```text
Drug A: 256
Drug B: 256
Cell:     0
--------------
Total:  512
```

### PCA128

```text
Drug A: 256
Drug B: 256
Cell:   128
--------------
Total:  640
```

### Raw Genes

```text
Drug A:   256
Drug B:   256
Cell:   23808
----------------
Total: 24320
```

## Interpretation Implication

Because the model has no explicit learned cell encoder:

- attribution is currently on the **input dimensions**
- for PCA runs, attribution is first at the **component level**
- for raw/filtered runs, attribution is directly at the **feature dimension level**

For the current best gene model (`pca128`), this means:

- the model is not choosing raw genes directly
- it is using PCA components
- mapping important components back to genes requires saved PCA loading information
