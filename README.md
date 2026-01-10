# Single-Cell Cancer Metabolic Modeling Using Flux Balance Cellular Automata

---

## 🧬 Project Overview

This project implements a **Flux Balance Cellular Automata (FBCA)** framework integrated with **single-cell RNA sequencing (scRNA-seq)** data to study **metabolic heterogeneity and phenotypic diversity in cancer**.

Cancer tissues exhibit strong intra-tumoral heterogeneity, where individual cells differ in:
- Metabolic programs  
- Proliferation capacity  
- Survival and dormancy states  
- Drug response and resistance  

FBCA provides a **systems biology approach** that combines:
- **Single-cell transcriptomics** (to infer metabolic activity)
- **Flux-based reasoning** (to approximate cellular energetic states)
- **Cellular automata** (to simulate spatial and temporal tumor evolution)

---

## 🎯 Objectives

- Quantify **metabolic heterogeneity** at single-cell resolution  
- Map scRNA-seq data to **metabolic pathway activities**
- Simulate **cell state transitions** (proliferation, survival, dormancy)
- Model **tumor evolution over time** using FBCA rules
- Provide a **reproducible and extensible framework** for cancer systems biology

---

## 🧪 Dataset

- **Type:** Single-cell RNA sequencing  
- **Cancer:** Breast cancer  
- **Source:** GSE75688
- **Technology:** Smart-seq / full-length scRNA-seq  


---

## 🔬 Methodology

### 1. scRNA-seq Analysis
- Quality control and filtering
- Normalization and scaling
- Dimensionality reduction (PCA, UMAP)
- Cell clustering using Seurat

### 2. Single-Cell Metabolic Scoring
Metabolic pathway activities are inferred using **GSVA (ssGSEA)** for:
- Glycolysis
- TCA cycle
- Oxidative phosphorylation (OxPhos)
- Glutamine metabolism
- Fatty acid synthesis

Each cell is assigned a **metabolic profile** based on gene expression.

### 3. Flux-Informed Cellular Automata (FBCA)
- Each cell is mapped onto a 2D spatial grid
- Cellular ATP proxy is estimated from metabolic scores
- Cells transition between states using biologically inspired rules:

| Metabolic Signal | Cell State |
|------------------|------------|
| High ATP         | Proliferative |
| Moderate ATP     | Survival |
| Low ATP          | Dormant |

- The system evolves over discrete time steps, simulating tumor dynamics.

---

## 🧠 Key Concepts

- **Single-cell systems biology**
- **Metabolic heterogeneity**
- **Tumor evolution modeling**
- **Flux balance reasoning (conceptual)**
- **Cellular automata**

---

## 📊 Outputs

- UMAP plots of metabolic pathway activity
- Spatial tumor maps across time steps
- Classification of metabolic cell states
- CSV file containing FBCA simulation results
- Publication-quality figures

---

## 📂 Repository Structure

```text
FBCA-scRNAseq/
│
├── data/
│   └── counts_matrix.csv        # scRNA-seq input data
│
├── scripts/
│   └── fbca_scRNAseq.R          # Main analysis & simulation script
│
├── results/
│   ├── figures/                 # Generated plots
│   └── fbca_simulation.csv      # FBCA time-series results
│
├── README.md
└── LICENSE
