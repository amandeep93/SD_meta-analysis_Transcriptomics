Sleep Deprivation Meta-analysis

The study integrates independent transcriptomic datasets to identify conserved molecular signatures associated with acute sleep deprivation across multiple rodent brain regions.

## Overview

This repository contains the complete computational workflow used for:

- Microarray preprocessing and differential expression analysis
- RNA-seq preprocessing and DESeq2 analysis
- Random-effects meta-analysis using metafor
- Cross-platform integration
- Forest plot generation
- Volcano plot generation
- Heatmap visualization
- Brain-region specific analysis
- Functional enrichment analysis
- PPI network preparation
- Reproducibility analysis

## Software Requirements

- R >= 4.4.1
- Bioconductor
- DESeq2
- limma
- metafor
- tidyverse
- pheatmap

## External software
  FastQC v0.12.1
  Trimmomatic v0.39
  STAR v2.7.10a
  featureCounts v2.0.6
  Cytoscape v3.10
  STRING database v12.0

## Repository Structure

- data/: Input datasets
- script/: All analysis scripts
- additional information: software versions used in this project


Workflow Summary
1. Microarray Processing

Raw CEL files were normalized using Robust Multi-array Average (RMA) normalization through the affy package. Differential expression analysis was performed using limma with empirical Bayes moderation.

2. RNA-seq Processing

RNA-seq datasets were processed using FastQC, Trimmomatic, STAR alignment, and featureCounts quantification. Differential expression analysis was performed using DESeq2.

3. Meta-analysis

Random-effects meta-analysis with REML estimation and Knapp–Hartung adjustment was performed using the metafor package.

4. Cross-platform Integration

Microarray and RNA-seq meta-analysis results were integrated using effect-size based approaches and Fisher’s combined probability method.

5. Brain Region-specific Analysis

Transcriptomic signatures were mapped to brain regions using metadata-guided tissue annotation and cross-platform concordance analysis. 

