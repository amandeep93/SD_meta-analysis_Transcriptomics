# Sleep Deprivation Meta-analysis Framework

An integrated, cross-platform computational workflow designed to synthesize independent multi-cohort transcriptomic datasets (Microarray and RNA-seq). This framework identifies conserved, directional molecular signatures associated with acute sleep deprivation across distinct rodent brain regions.

## Pipeline Architecture Overview
This repository hosts the complete end-to-end reproducible pipeline for:
- Automated upstream raw sequencing quality control, adapter clipping, and genomic alignment.
- Standalone Microarray preprocessing, Robust Multi-array Average (RMA) normalization, and `limma` modeling.
- Individual cohort negative-binomial generalized linear modeling via `DESeq2`.
- Rigorous Random-Effects Meta-Analysis utilizing Restricted Maximum Likelihood (REML) estimations with Knapp–Hartung adjustments via `metafor`.
- Cross-platform directional harmony validation and joint Fisher-to-FDR probability integration.
- Automated publication-ready figure generation (Forest plots, spatial clustering heatmaps, and tissue-intersect network graphs).

---

## Technical Specifications & Requirements

### Computational Environment
- **Operating System:** Linux / macOS terminal environment (for upstream bash sequencing alignment scripts)
- **R Environment Version:** `R >= 4.4.1`

### Required R / Bioconductor Packages
| Package Name | Minimum Version | Core Pipeline Utility |
| :--- | :--- | :--- |
| `affy` | Bioconductor | Binary CEL array background correction & normalization |
| `limma` | Bioconductor | General Linear Modeling & empirical Bayes moderation |
| `DESeq2` | Bioconductor | Negative Binomial GLM raw read counts moderation |
| `metafor` | CRAN | REML Random-effects modeling & Knapp-Hartung calculation |
| `dplyr` / `tidyr` / `purrr` | CRAN | Tidyverse data frame manipulations & parallel iterations |
| `pheatmap` | CRAN | High-resolution hierarchical clustering heatmaps |
| `ggplot2` | CRAN | Base multi-layered programmatic charting engine |
| `igraph` / `ggraph` | CRAN | Relational network matrix node & edge structural plotting |

### Upstream Bash Binary Dependencies
- **FastQC** `v0.12.1` (Baseline raw read quality analytics tracking)
- **Trimmomatic** `v0.39` (Phred-score guided adapter and low-quality base clipping)
- **STAR Alignment Engine** `v2.7.10a` (Splice-aware genomic index mapping)
- **featureCounts (subread)** `v2.0.6` (Digital transcript feature counts quantification)

---

## Repository Directory Layout

```text
├── data/
│   ├── CEL_files/             <- Raw microarray binary file subfolders (GSE6514, GSE33302)
│   ├── RNAseq_raw/            <- Raw paired-end sequencing inputs (*_R1.fastq.gz, *_R2.fastq.gz)
│   ├── RNAseq_counts/         <- Output destination for featureCounts digital matrices
│   └── metadata/              <- Cohort-specific phenotype annotations (*_metadata.csv)
├── script/
│   ├── 01_rnaseq_preprocessing.sh        <- Quality check, adapter trimming, & STAR alignments
│   ├── 02_microarray_preprocessing.R     <- RMA processing & exact standard error extraction
│   ├── 03_microarray_meta_analysis.R     <- Microarray REML meta-analysis & probe collapse
│   ├── 04_microarray_visualization.R      <- Microarray forest profiling & clustering heatmaps
│   ├── 05_rnaseq_differential_expression.R <- DESeq2 individual study reference group loops
│   ├── 06_rnaseq_meta_analysis.R         <- RNA-seq REML profiling utilizing native lfcSE weights
│   ├── 07_cross_platform_meta_analysis.R  <- Directional harmony validation & Fisher integration
│   └── 08_brain_region_analysis.R        <- Tissue-specific Fisher loops, networks, & heatmaps
├── results/                   <- Auto-generated analytical tables (*_DEG.csv, meta assets)
└── figures/                   <- High-resolution publication-ready vector/raster charts (.png)
```
Sequential Execution Roadmap
To completely reproduce the findings detailed in the manuscript, execute the analysis scripts sequentially from your terminal:
Phase 1: Upstream Processing & Individual Study Modeling
# 1. Execute raw RNA-seq quality trimming, alignment, and gene counting loops
bash script/01_rnaseq_preprocessing.sh

# 2. Process raw Microarray arrays, generate limma models, and extract true SE
Rscript script/02_microarray_preprocessing.R

# 3. Model individual raw RNA-seq counts through DESeq2 with strict control contrasts
Rscript script/05_rnaseq_differential_expression.R

Phase 2: Standalone Platform Meta-Analyses
# 4. Run REML meta-analysis across consolidated microarray cohorts
Rscript script/03_microarray_meta_analysis.R

# 5. Generate high-resolution forest charts and microarray clustering arrays
Rscript script/04_microarray_visualization.R

# 6. Run REML meta-analysis across DESeq2 pipelines using native lfcSE variance weights
Rscript script/06_rnaseq_meta_analysis.R

Phase 3: Cross-Platform Global & Regional Integration
# 7. Execute Fisher's integration combined with strict directional harmony constraints
Rscript script/07_cross_platform_meta_analysis.R

# 8. Conduct localized tissue Fisher integration, spatial heatmaps, and intersection networks
Rscript script/08_brain_region_analysis.R


Core Analytical Workflow Details:
Microarray Normalization & Moderation: Raw arrays are background corrected, quantile normalized, and log2 transformed using Robust Multi-array Average (RMA). True moderated standard errors are captured straight from the unscaled standard error fields of limma models (SE = abs(logFC) / abs(t_stat)).
RNA-seq Dispersion Scaling: Read allocations are modeled using a negative binomial generalized linear model inside DESeq2. Native log fold change standard errors (lfcSE) are preserved to prevent sample weight distortions during downstream pooling steps.
Robust Synthesis Estimation: Studies are synthesized via a Random-Effects Meta-Analysis with Restricted Maximum Likelihood (REML) and Knapp–Hartung adjustments via the metafor package. Multi-probe/transcript overlapping instances are filtered by max-responsiveness to strictly preserve statistical independence parameters.
Fisher Combined Probability Transition: Merging cross-platform data relies on joint Chi-squared calculations (df=4). It implements a strict directional enforcement script (sign(meta_logFC.x) == sign(meta_logFC.y)) to drop anti-correlated regulatory variations, running Benjamini-Hochberg False Discovery Rate (FDR) corrections globally.
Spatial Structural Intersections: Target markers are segmented based on their associated localized brain regions. Concordance scatters, horizontal distribution bars, and stress-layout interaction networks are programmatically compiled to illustrate overlapping transcriptomic networks across distinct anatomical boundaries.

