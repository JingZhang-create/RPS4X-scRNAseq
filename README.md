# 🧬 RPS4X Single-Cell Transcriptomics

**Single-cell RNA-seq analysis of control and *Rps4x*-deficient MB49 tumors**

---

## 📖 Overview

This repository contains the R scripts used for single-cell RNA sequencing (scRNA-seq) analysis of subcutaneous tumors from C57BL/6J mice implanted with control or *Rps4x*-deficient MB49 cells. We provide the computational workflow for major cell-type annotation, CD45-positive immune-cell subclustering, and downstream gene-expression analyses.

---

## ⚙️ Computational Workflow

Raw scRNA-seq data were processed using the publicly accessible **OmicStudio cloud platform**. Downstream single-cell analyses were performed in **R** using **Seurat (v5.2.1)**. The R scripts provided in this repository document the analytical procedures and key parameters used for the reported scRNA-seq analyses.

---

## 🔬 Analysis Scripts

| Script | Description |
|---|---|
| `01_major_celltype_annotation.R` | Major cell-type annotation based on gene markers |
| `02_CD45pos_immune_subclustering.R` | Reclustering of CD45-positive immune cells |
| `03_downstream_analysis.R` | Visualization of selected gene expression in epithelial cells and average-expression profiling of immune-cell clusters |

---

## 📂 Input Data

Raw sequencing data are available through the data repository specified in the associated manuscript.

---

## 📌 Repository Structure

```text
RPS4X-single-cell-transcriptomics/
│
├── README.md
├── 01_major_celltype_annotation.R
├── 02_CD45pos_immune_subclustering.R
└── 03_downstream_analysis.R
```

