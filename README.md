# COVID-19 vs Influenza Differential Expression Analysis

Group final project for BINF 6310 (Introduction to Bioinformatics) at Northeastern 
University. This project is a reproducibility study of transcriptomic signatures in 
COVID-19 using publicly available RNA-seq data.

**Dataset:** GSE161731 — RNA-seq data from COVID-19 and influenza patients

**My contribution:** Differential expression analysis using the limma/voom pipeline in R

## Key Findings

194 significantly differentially expressed genes identified between COVID-19 and 
influenza patients (adj. p < 0.05) out of 15,556 genes tested.

## Methods

- TMM normalization using edgeR
- Voom transformation to model mean-variance relationship
- Linear modeling with limma and empirical Bayes moderation
- Benjamini-Hochberg correction for multiple testing
- Visualizations: volcano plot, MA plot, heatmap of top 50 DEGs

## Files

- `differential_expression.R` — full analysis script
- `DE_results_COVID_vs_Influenza_full.tsv` — complete results table
- `volcano_plot.pdf` — volcano plot of all tested genes
- `heatmap_top50_DEGs.pdf` — heatmap of top 50 DEGs
- `MA_plot.pdf` — MA plot

## Tools

R, limma, edgeR, ggplot2, pheatmap

## Note
This was a 4-person group project. This repository contains my individual 
contribution — the differential expression analysis component. The full group 
project repository including pathway enrichment analysis and final report can 
be found here: https://github.com/NakamuraNamie/Group1-FinalProject