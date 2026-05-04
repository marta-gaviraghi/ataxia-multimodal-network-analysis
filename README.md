# Ataxia Multimodal Network Analysis

## Overview

This repository contains the **statistical and machine learning analysis pipeline** used for the study:

**Multimodal imaging and virtual brain modelling reveal the diversity of network changes in cerebellar ataxias**  
submitted to *Brain Communications*.

The repository provides a **fully reproducible workflow** integrating:

- statistical analysis of brain volumes across functional networks
- network-based regressions
- dimensionality reduction using Principal Component Analysis (PCA)
- unsupervised clustering analyses

The analyses aim to characterise **network-level alterations** and identify data-driven patterns of heterogeneity across cerebellar ataxia patients.

---

## Repository Structure
├── volume_analysis.R
├── regression_analysis.R
├── pca_clustering.R
└── README.md


### Scripts

#### volume_analysis.R
Example implementation of the statistical workflow applied to brain volume measures.

This script serves as a **template** for analysing additional feature sets (eg. graph metrics on SC and FC, TVB parameters).

---

#### regression_analysis.R
Performs **network-specific stepwise linear regression** analyses to investigate associations between imaging features and clinical variables.

Each functional network is analysed independently.

---

#### pca_clustering.R
Implements a data-driven workflow including:

- feature normalisation
- Principal Component Analysis (PCA)
- dimensionality reduction
- unsupervised clustering for patient stratification

---

## Data Requirements

The repository does **not** include raw data.

Input datasets should contain:

- subject identifiers
- group labels
- regional or network-derived metrics
- clinical variables (when required)

Data paths must be adapted locally within each script.

---

## Dependencies

Analyses were developed in **R**.

Required packages typically include:

tidyverse
dplyr
ggplot2
stats
FactoMineR
factoextra
cluster

## Notes

- The provided scripts illustrate the **analysis strategy**, not a fixed pipeline.
- The same analytical framework can be reused across **multiple datasets and imaging modalities**.
- The volume analysis is included as an **example implementation** of the general workflow.

---

## Citation

If you use this code, please cite:

**Gaviraghi et al.**  
*Multimodal imaging and virtual brain modelling reveal the diversity of network changes in cerebellar ataxias*  
*Brain Communications* (under review)
