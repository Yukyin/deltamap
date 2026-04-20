# EvidenceFM: Design Concept

This document describes the broader architectural vision motivating the analysis pipeline in this repository. The current codebase implements the front-end representation layer of a larger system; this document outlines the full intended design.

## Motivation

Foundation models in genomics (Geneformer, scGPT, etc.) typically operate at the cell level and learn rich representations of individual cells. But for many clinical research questions, the scientifically meaningful unit is the subject, and the meaningful signal is not a static expression snapshot but a change over time in response to disease state, treatment, or environmental perturbation.

The central question EvidenceFM is designed to answer:

> Can we learn a universal, interpretable, cross-study representation of subject-level immune evidence, built from temporal delta signals, that is aggregatable, transferable across cohorts, and useful for downstream clinical reasoning?



## Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│                        EvidenceFM System                        │
│                                                                 │
│  ┌──────────────┐    ┌──────────────┐    ┌───────────────────┐  │
│  │  Front-end   │    │  Model layer │    │  Downstream tasks │  │
│  │  (this repo) │───▶│  (planned)   │───▶│  (planned)        │  │
│  └──────────────┘    └──────────────┘    └───────────────────┘  │
│                                                                 │
│  Paired scRNA-seq    Transformer or       Disease state         │
│  → pseudobulk delta  graph model          prediction,           │
│  → evidence modules  pre-trained on       patient stratification│
│  → subject vector    multi-cohort         treatment response    │
│                      evidence vectors                           │
└─────────────────────────────────────────────────────────────────┘
```

### Layer 1: Evidence Representation (implemented here)

Input: raw paired scRNA-seq (two time points per subject)

Process:
1. Pseudobulk aggregation per sample (cells to gene-sum vector)
2. Within-subject temporal delta (D2 minus D1)
3. Evidence module scoring: project delta vector onto a small set of biologically interpretable gene programs
4. Low-dimensional subject-level evidence matrix (subjects x modules)

Output: a structured, interpretable, low-dimensional vector per subject that captures the dominant axes of immune state change.

Design principles:
- Interpretability first: module scores map directly to named biological programs (humoral, myeloid, APC/MHC-II, stress response)
- Within-subject normalization: using paired delta removes inter-subject baseline variation, isolating response signal
- Robustness checks: bootstrap stability, LOSO validation, outlier audit



### Layer 2: Foundation Model (planned)

Input: subject evidence vectors from Layer 1, across multiple cohorts and diseases

Architecture candidates:
- Transformer encoder over module-score sequences (subjects as tokens)
- Graph neural network over subject-subject similarity graphs
- Contrastive pre-training across cohort-matched pairs

Pre-training objective: learn cohort-invariant subject embeddings such that subjects with similar immune states are close in embedding space, regardless of cohort or platform.

Transfer: fine-tune on disease-specific downstream tasks with small labeled datasets.



### Layer 3: Downstream Tasks (planned)

- Binary case/control classification (requires sufficient n and cross-cohort validation)
- Disease subtype discovery (continuous immune state)
- Treatment response prediction (paired pre/post intervention)
- Cross-disease transfer (shared immune programs across ME/CFS, Long COVID, autoimmune)



## Why This Approach

| Dimension | Standard scRNA-seq analysis | EvidenceFM approach |
|-----------|---------------------------|---------------------|
| Analysis unit | Cell or sample | Subject (temporal change) |
| Representation | High-dimensional, sparse | Low-dimensional, interpretable modules |
| Cross-study use | Requires re-analysis | Standardized evidence vector |
| Clinical interpretability | Low | High (named programs) |
| Foundation model input | Cell embeddings | Subject evidence vectors |



## Current Status

| Component | Status |
|-----------|--------|
| Front-end pipeline (GSE214284) | Complete |
| Module stability validation | Complete (bootstrap, LOSO) |
| Multi-cohort generalization | Not yet implemented |
| Foundation model training | Not yet implemented |
| Cross-disease transfer | Not yet implemented |

The analysis in this repository demonstrates that a stable, interpretable 4-module evidence representation can be extracted from paired scRNA-seq in ME/CFS. Whether these modules generalize across cohorts and diseases is the key open question for future work.



## Relevant Prior Work

- Pseudobulk methods: Squair et al. (2021) Nature Communications
- scRNA-seq foundation models: Theodoris et al. Geneformer (2023), Chen et al. scGPT (2024)
- ME/CFS immune profiling: Ahmed et al. Cell Reports Medicine (2024)
- Evidence-based clinical AI: general direction, no direct prior art for this specific framing
