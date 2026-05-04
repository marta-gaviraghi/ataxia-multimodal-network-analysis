# Ataxia Multimodal Network Analysis – Volumetric Pipeline

## Overview

This repository contains the **statistical analysis pipeline** used for the volumetric network analyses reported in the manuscript:

**Multimodal imaging and virtual brain modelling reveal the diversity of network changes in cerebellar ataxias**  
submitted to *Brain Communications*.

The code implements a **fully reproducible workflow** for group-level comparisons of brain volumes across functional networks.

---

## Study Groups

- **HC** – Healthy Controls  
- **JS** – Joubert Syndrome  
- **SP** – Slowly progressive ataxia 

---

## Data Types

- Regional brain volumes  
- Network-level volumetric measures  

---

## Networks Analyzed

- Somatomotor Network  
- Ventral Attention Network  

The pipeline is **generalizable to other networks**.

---

## Statistical Pipeline

For each region / feature:

- **Shapiro–Wilk test** for normality  
- **One-way ANOVA** (normal data) or **Kruskal–Wallis** (non-normal data)  
- **Post-hoc tests**:
  - Tukey HSD (ANOVA)  
  - Dunn test with Holm correction (Kruskal–Wallis)  
- Significance threshold: **p < 0.05**  

Only **significant regions** are visualized.

---

## Reproducibility DA SISTEMARE!!!

Raw clinical data cannot be shared due to ethical restrictions.

An **example dataset** is provided in:

`example_data/example_volumes.xlsx`

---

## Relationship to Full Study

The same statistical pipeline was applied across all modalities in the study:

- structural connectivity (SC) graph metrics  
- functional connectivity (FC) graph metrics  
- The Virtual Brain (TVB) parameters  
- volumetric analyses  

This repository provides a **representative implementation** of the full analytical framework.

---

## Requirements

R (≥ 4.2)

Packages:
- ggplot2  
- dplyr  
- dunn.test  
- readxl  
- tidyr  

---

## Run Analysis

```r
source("volume_analysis.R")
