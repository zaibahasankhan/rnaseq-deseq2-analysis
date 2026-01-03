# RNA-seq Differential Expression Analysis using DESeq2

## Overview
This repository provides a complete RNA-seq differential expression
analysis workflow using the DESeq2 package in R.

## Input Data
- Raw gene count matrix (genes × samples)
- Sample metadata describing experimental conditions

## Tools & Packages
- R (>= 4.0)
- DESeq2
- ggplot2
- pheatmap
- dplyr

## Workflow
1. Import count data and metadata
2. Quality control and normalization
3. Differential expression analysis
4. Visualization (volcano plot, heatmap)

## Output
- Differentially expressed genes (DEGs)
- Volcano plot
- Heatmap of significant genes

## How to Run
```r
source("scripts/deseq2_analysis.R")
```r
source("scripts/deseq2_analysis.R")
