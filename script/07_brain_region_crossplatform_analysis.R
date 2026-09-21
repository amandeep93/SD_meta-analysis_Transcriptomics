############################################################
# PRECISION-WEIGHTED CROSS-PLATFORM META-ANALYSIS
# Script: script/07_cross_platform_meta_analysis.R
#
# Microarray + RNA-seq integration using inverse-variance
# weighting of platform-level random-effects estimates.
############################################################

############################################################
# SETTINGS & CONFIGURATIONS
############################################################

options(stringsAsFactors = FALSE)
set.seed(123)

############################################################
# LIBRARIES
############################################################

library(dplyr)

############################################################
# INPUT FILES
############################################################

micro_path <- "results/Microarray_Meta_AllGenes.csv"
rna_path   <- "results/RNAseq_Meta_AllGenes.csv"

if(!file.exists(micro_path)) {
  stop(
    "Critical Error: Microarray meta-analysis file not found: ",
    micro_path,
    "\nRun Script 03 first."
  )
}

if(!file.exists(rna_path)) {
  stop(
    "Critical Error: RNA-seq meta-analysis file not found: ",
    rna_path,
    "\nRun Script 06 first."
  )
}

############################################################
# LOAD PLATFORM-SPECIFIC META-ANALYSIS RESULTS
############################################################

cat("Loading platform-specific meta-analysis results...\n")

micro <- read.csv(
  micro_path,
  stringsAsFactors = FALSE
)

rna <- read.csv(
  rna_path,
  stringsAsFactors = FALSE
)

############################################################
# VERIFY REQUIRED COLUMNS
############################################################

required_micro <- c(
  "Gene",
  "meta_logFC",
  "meta_SE",
  "meta_pval"
)

required_rna <- c(
  "Gene",
  "meta_logFC",
  "meta_SE",
  "meta_pval"
)

missing_micro <- setdiff(required_micro, colnames(micro))
missing_rna   <- setdiff(required_rna, colnames(rna))

if(length(missing_micro) > 0) {
  stop(
    "Microarray meta-analysis file is missing required columns: ",
    paste(missing_micro, collapse = ", "),
    "\nModify Script 03 to export meta_SE."
  )
}

if(length(missing_rna) > 0) {
  stop(
    "RNA-seq meta-analysis file is missing required columns: ",
    paste(missing_rna, collapse = ", "),
    "\nModify Script 06 to export meta_SE."
  )
}

############################################################
# CLEAN PLATFORM DATA
############################################################

micro_clean <- micro %>%
  mutate(
    Gene = toupper(trimws(Gene))
  ) %>%
  filter(
    !is.na(Gene),
    Gene != "",
    is.finite(meta_logFC),
    is.finite(meta_SE),
    meta_SE > 0,
    is.finite(meta_pval)
  ) %>%
  select(
    Gene,
    meta_logFC,
    meta_SE,
    meta_pval,
    CI_lb,
    CI_ub,
    I2,
    tau2,
    n_studies
  ) %>%
  rename(
    micro_logFC = meta_logFC,
    micro_SE    = meta_SE,
    micro_pval  = meta_pval,
    micro_CI_lb = CI_lb,
    micro_CI_ub = CI_ub,
    micro_I2    = I2,
    micro_tau2  = tau2,
    micro_n     = n_studies
  )

rna_clean <- rna %>%
  mutate(
    Gene = toupper(trimws(Gene))
  ) %>%
  filter(
    !is.na(Gene),
    Gene != "",
    is.finite(meta_logFC),
    is.finite(meta_SE),
    meta_SE > 0,
    is.finite(meta_pval)
  ) %>%
  select(
    Gene,
    meta_logFC,
    meta_SE,
    meta_pval,
    CI_lb,
    CI_ub,
    I2,
    tau2,
    n_studies
  ) %>%
  rename(
    rna_logFC = meta_logFC,
    rna_SE    = meta_SE,
    rna_pval  = meta_pval,
    rna_CI_lb = CI_lb,
    rna_CI_ub = CI_ub,
    rna_I2    = I2,
    rna_tau2  = tau2,
    rna_n     = n_studies
  )

############################################################
# MERGE MICROARRAY AND RNA-SEQ RESULTS
############################################################

cat("Identifying genes represented by both platforms...\n")

common <- inner_join(
  micro_clean,
  rna_clean,
  by = "Gene"
)

cat(
  "Genes shared between microarray and RNA-seq: ",
  nrow(common),
  "\n",
  sep = ""
)

if(nrow(common) == 0) {
  stop(
    "Critical Error: No genes were shared between microarray and RNA-seq results."
  )
}

############################################################
# CHECK DIRECTIONAL CONSISTENCY
############################################################
#
# We retain genes for which the platform-level effects have
# the same direction.
#
# This preserves the original manuscript's requirement for
# cross-platform directional concordance.
############################################################

common_harmonized <- common %>%
  filter(
    sign(micro_logFC) == sign(rna_logFC)
  )

cat(
  "Genes with concordant direction across platforms: ",
  nrow(common_harmonized),
  "\n",
  sep = ""
)

if(nrow(common_harmonized) == 0) {
  stop(
    "Critical Error: No genes showed concordant direction across platforms."
  )
}

############################################################
# PRECISION-WEIGHTED CROSS-PLATFORM INTEGRATION
############################################################
#
# Inverse-variance weighting:
#
# weight = 1 / SE^2
#
# combined effect =
#   sum(weight_i * effect_i) / sum(weight_i)
#
# combined SE =
#   sqrt(1 / sum(weight_i))
#
############################################################

cat("Performing inverse-variance weighted cross-platform integration...\n")

common_meta <- common_harmonized %>%
  mutate(

    # ------------------------------------------------------
    # Platform-specific inverse-variance weights
    # ------------------------------------------------------

    weight_microarray = 1 / (micro_SE^2),

    weight_rnaseq = 1 / (rna_SE^2),

    # ------------------------------------------------------
    # Precision-weighted combined log2 fold change
    # ------------------------------------------------------

    combined_logFC =
      (
        weight_microarray * micro_logFC +
        weight_rnaseq * rna_logFC
      ) /
      (
        weight_microarray +
        weight_rnaseq
      ),

    # ------------------------------------------------------
    # Standard error of combined effect
    # ------------------------------------------------------

    combined_SE =
      sqrt(
        1 /
        (
          weight_microarray +
          weight_rnaseq
        )
      ),

    # ------------------------------------------------------
    # 95% confidence interval
    # ------------------------------------------------------

    combined_CI_lb =
      combined_logFC -
      1.96 * combined_SE,

    combined_CI_ub =
      combined_logFC +
      1.96 * combined_SE,

    # ------------------------------------------------------
    # Z statistic
    # ------------------------------------------------------

    combined_z =
      combined_logFC / combined_SE,

    # ------------------------------------------------------
    # Two-sided p-value
    # ------------------------------------------------------

    combined_p =
      2 * pnorm(
        -abs(combined_z)
      )
  )

############################################################
# BENJAMINI-HOCHBERG MULTIPLE TESTING CORRECTION
############################################################

cat("Applying Benjamini-Hochberg FDR correction...\n")

common_meta$FDR <- p.adjust(
  common_meta$combined_p,
  method = "BH"
)

############################################################
# CALCULATE RELATIVE PLATFORM CONTRIBUTIONS
############################################################

common_meta <- common_meta %>%
  mutate(

    total_weight =
      weight_microarray +
      weight_rnaseq,

    microarray_weight_fraction =
      weight_microarray / total_weight,

    rnaseq_weight_fraction =
      weight_rnaseq / total_weight

  )

############################################################
# SAVE COMPLETE CROSS-PLATFORM RESULTS
############################################################

write.csv(
  common_meta,
  "results/CrossPlatform_AllGenes_Weighted.csv",
  row.names = FALSE
)

cat(
  "Complete precision-weighted cross-platform results saved to:\n",
  "results/CrossPlatform_AllGenes_Weighted.csv\n"
)

############################################################
# SELECT SIGNIFICANT CROSS-PLATFORM GENES
############################################################
#
# Same threshold used previously:
# FDR < 0.05
# absolute combined log2FC >= 1
#
############################################################

final_genes <- common_meta %>%
  filter(
    FDR < 0.05,
    abs(combined_logFC) >= 1
  ) %>%
  arrange(FDR)

############################################################
# SAVE FINAL CROSS-PLATFORM GENE LIST
############################################################

write.csv(
  final_genes,
  "results/CrossPlatform_Weighted_genes.csv",
  row.names = FALSE
)

############################################################
# OPTIONAL: SAVE UP/DOWN REGULATED LISTS
############################################################

up_genes <- final_genes %>%
  filter(combined_logFC > 0)

down_genes <- final_genes %>%
  filter(combined_logFC < 0)

write.csv(
  up_genes,
  "results/CrossPlatform_Weighted_UP.csv",
  row.names = FALSE
)

write.csv(
  down_genes,
  "results/CrossPlatform_Weighted_DOWN.csv",
  row.names = FALSE
)

############################################################
# SUMMARY
############################################################

cat("============================================================\n")
cat("PRECISION-WEIGHTED CROSS-PLATFORM META-ANALYSIS COMPLETE\n")
cat("============================================================\n")

cat(
  "Genes shared between platforms: ",
  nrow(common),
  "\n",
  sep = ""
)

cat(
  "Genes with concordant direction: ",
  nrow(common_harmonized),
  "\n",
  sep = ""
)

cat(
  "Significant cross-platform genes (FDR < 0.05 and |log2FC| >= 1): ",
  nrow(final_genes),
  "\n",
  sep = ""
)

cat(
  "Upregulated genes: ",
  nrow(up_genes),
  "\n",
  sep = ""
)

cat(
  "Downregulated genes: ",
  nrow(down_genes),
  "\n",
  sep = ""
)

cat("============================================================\n")
cat("Output files:\n")
cat("1. results/CrossPlatform_AllGenes_Weighted.csv\n")
cat("2. results/CrossPlatform_Weighted_genes.csv\n")
cat("3. results/CrossPlatform_Weighted_UP.csv\n")
cat("4. results/CrossPlatform_Weighted_DOWN.csv\n")
cat("============================================================\n")
