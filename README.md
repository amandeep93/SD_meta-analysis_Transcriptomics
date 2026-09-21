# Sleep Deprivation Meta-analysis Framework

An integrated, cross-platform computational workflow designed to synthesize independent multi-cohort transcriptomic datasets (microarray and RNA-seq). The framework identifies conserved and directionally concordant molecular signatures associated with acute sleep deprivation across rodent brain regions.

---

## Pipeline Architecture Overview

This repository hosts the complete end-to-end computational workflow for:

- Automated upstream RNA-seq quality control, adapter trimming, splice-aware genomic alignment, and gene-level quantification.
- Standalone microarray preprocessing using Robust Multi-array Average (RMA) normalization and `limma` differential-expression modeling.
- Individual RNA-seq differential-expression analysis using the negative-binomial generalized linear model implemented in `DESeq2`.
- Platform-specific random-effects meta-analysis using Restricted Maximum Likelihood (REML) with Knapp–Hartung inference via `metafor`.
- Treatment-independent microarray probe selection based on mean expression rather than treatment-associated effect size.
- Cross-platform integration using precision-weighted synthesis of platform-level effect estimates and their standard errors.
- Directional concordance assessment between microarray and RNA-seq estimates.
- Brain-region-stratified random-effects meta-analysis and cross-platform synthesis where region-level metadata permit independent regional assignment.
- Automated publication-ready visualization, including forest plots, transcriptomic heatmaps, regional heatmaps, and regional gene-count plots.

---

## Technical Specifications & Requirements

### Computational Environment

- **Operating System:** Linux / macOS terminal environment for upstream RNA-seq processing.
- **R Environment:** `R >= 4.4.1`

### Required R / Bioconductor Packages

| Package | Source | Core Pipeline Utility |
|---|---|---|
| `affy` | Bioconductor | Microarray preprocessing and RMA normalization |
| `limma` | Bioconductor | Linear modeling and empirical Bayes moderation |
| `DESeq2` | Bioconductor | Negative-binomial RNA-seq differential-expression modeling |
| `metafor` | CRAN | REML random-effects meta-analysis and Knapp–Hartung inference |
| `dplyr` | CRAN | Data manipulation |
| `tidyr` | CRAN | Data restructuring |
| `purrr` | CRAN | Iterative analysis across datasets |
| `pheatmap` | CRAN | Hierarchical clustering heatmaps |
| `ggplot2` | CRAN | Programmatic statistical visualization |

### Upstream Bash Binary Dependencies

- **FastQC** `v0.12.1` — raw sequencing-read quality assessment.
- **Trimmomatic** `v0.39` — adapter removal and quality trimming.
- **STAR** `v2.7.10a` — splice-aware genomic alignment.
- **featureCounts** (`subread`) `v2.0.6` — gene-level read quantification.

---

## Reference Genome and Annotation

RNA-seq preprocessing uses a mouse reference framework consisting of:

- **Genome:** GRCm39
- **Gene annotation:** GENCODE vM32
- **STAR index:** GRCm39-GencodeM32

The same reference framework is maintained throughout STAR alignment and gene-level quantification.

---

## Repository Directory Layout

```text
├── data/
│   ├── CEL_files/
│   │   └── <microarray dataset folders>
│   ├── RNAseq_raw/
│   │   └── Raw paired-end FASTQ files
│   ├── RNAseq_counts/
│   │   └── featureCounts gene-level count matrices
│   └── metadata/
│       └── Cohort-specific phenotype and experimental annotations
│
├── script/
│   ├── 01_rnaseq_preprocessing.sh
│   │   └── FastQC, Trimmomatic, STAR, featureCounts
│   │
│   ├── 02_microarray_preprocessing.R
│   │   └── RMA normalization and limma differential expression
│   │
│   ├── 03_microarray_meta_analysis.R
│   │   └── Microarray REML random-effects meta-analysis
│   │
│   ├── 04_microarray_visualization.R
│   │   └── Forest plots and microarray heatmap
│   │
│   ├── 05_rnaseq_differential_expression.R
│   │   └── Individual-cohort DESeq2 analysis
│   │
│   ├── 06_rnaseq_meta_analysis.R
│   │   └── RNA-seq REML random-effects meta-analysis
│   │
│   ├── 07_cross_platform_meta_analysis.R
│   │   └── Precision-weighted cross-platform synthesis
│   │
│   └── 08_brain_region_analysis.R
│       └── Brain-region-stratified meta-analysis and
│           cross-platform regional synthesis
│
├── results/
│   ├── rna_individual/
│   └── brain_region/
│
└── figures/
    ├── microarray_forests/
    └── brain_region/
