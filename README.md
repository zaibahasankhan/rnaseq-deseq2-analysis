# RNA-seq Differential Expression Analysis using DESeq2

## Overview
This repository contains a complete and reproducible RNA-seq differential expression analysis workflow implemented in R using the DESeq2 package. 
The pipeline performs normalization, statistical testing, and downstream visualization to identify genes differentially expressed between experimental conditions (cancer vs. normal).

## Input Data
- Raw gene count matrix (genes × samples)
- Sample metadata describing experimental conditions
Example input files are provided for demonstration purposes

## Tools & Packages
- R (>= 4.0)
- DESeq2
- EnhancedVolcano
- ggplot2
- pheatmap

## Workflow
1. Import count data and metadata
2. Pre-processing and filtering of low-count genes
3. Differential expression analysis using DESeq2
5. Visualization (volcano plot, heatmap)
6. Visualization and quality control:
   Volcano plot
   MA plot
   Principal Component Analysis (PCA)
   Heatmap of top differentially expressed genes

## Output
- Differentially expressed genes (DEGs):
  upregulated and downregulated 
- Volcano plot
- MA plot
- Heatmap of significant genes (top 100)

## How to Run
```r
source("scripts/deseq2_analysis.R")
