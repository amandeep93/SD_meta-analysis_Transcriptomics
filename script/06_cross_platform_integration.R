############################################################
# PRECISION-WEIGHTED CROSS-PLATFORM META-ANALYSIS
#
# Script: script/07_cross_platform_meta_analysis.R
#
# Microarray + RNA-seq platform-level synthesis
#
# Method:
#   - Match genes present in both platforms
#   - Require concordant direction
#   - Precision-weight platform-level estimates using 1/SE^2
#   - Calculate combined log2FC and SE
#   - Derive 95% CI and two-sided p-value
#   - Apply Benjamini-Hochberg FDR
############################################################

############################################################
# SETTINGS
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

if (!file.exists(micro_path)) {
  stop(
    "Critical Error: Microarray meta-analysis file missing. ",
    "Run Script 03 first."
  )
}

if (!file.exists(rna_path)) {
  stop(
    "Critical Error: RNA-seq meta-analysis file missing. ",
    "Run Script 06 first."
  )
}

############################################################
# LOAD PLATFORM-LEVEL META-ANALYSIS RESULTS
############################################################

micro <- read.csv(
  micro_path,
  stringsAsFactors = FALSE
)

rna <- read.csv(
  rna_path,
  stringsAsFactors = FALSE
)

############################################################
# CHECK REQUIRED COLUMNS
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

missing_micro <- setdiff(
  required_micro,
  colnames(micro)
)

missing_rna <- setdiff(
  required_rna,
  colnames(rna)
)

if (length(missing_micro) > 0) {
  stop(
    "Missing columns in microarray meta-analysis results: ",
    paste(
      missing_micro,
      collapse = ", "
    )
  )
}

if (length(missing_rna) > 0) {
  stop(
    "Missing columns in RNA-seq meta-analysis results: ",
    paste(
      missing_rna,
      collapse = ", "
    )
  )
}

############################################################
# CLEAN PLATFORM RESULTS
############################################################

micro_clean <- micro %>%
  
  mutate(
    Gene = toupper(
      trimws(Gene)
    )
  ) %>%
  
  filter(
    !is.na(Gene),
    Gene != "",
    !is.na(meta_logFC),
    !is.na(meta_SE),
    is.finite(meta_logFC),
    is.finite(meta_SE),
    meta_SE > 0
  ) %>%
  
  select(
    Gene,
    meta_logFC,
    meta_SE,
    meta_pval
  )

rna_clean <- rna %>%
  
  mutate(
    Gene = toupper(
      trimws(Gene)
    )
  ) %>%
  
  filter(
    !is.na(Gene),
    Gene != "",
    !is.na(meta_logFC),
    !is.na(meta_SE),
    is.finite(meta_logFC),
    is.finite(meta_SE),
    meta_SE > 0
  ) %>%
  
  select(
    Gene,
    meta_logFC,
    meta_SE,
    meta_pval
  )

############################################################
# CHECK FOR DUPLICATE GENES
############################################################

if (anyDuplicated(micro_clean$Gene) > 0) {
  
  stop(
    "Duplicate Gene entries detected in microarray ",
    "meta-analysis results."
  )
}

if (anyDuplicated(rna_clean$Gene) > 0) {
  
  stop(
    "Duplicate Gene entries detected in RNA-seq ",
    "meta-analysis results."
  )
}

############################################################
# MATCH GENES ACROSS PLATFORMS
############################################################

cat(
  "Synchronizing genes across microarray and RNA-seq ",
  "platforms...\n"
)

common <- inner_join(
  micro_clean,
  rna_clean,
  by = "Gene",
  suffix = c(
    "_micro",
    "_rna"
  )
)

cat(
  "Genes shared between platforms: ",
  nrow(common),
  "\n",
  sep = ""
)

if (nrow(common) == 0) {
  
  stop(
    "No overlapping genes were found between the two ",
    "platform-level meta-analysis results."
  )
}

############################################################
# REQUIRE CONCORDANT DIRECTION
#
# Positive = upregulated in SD
# Negative = downregulated in SD
############################################################

cat(
  "Applying cross-platform directional concordance filter...\n"
)

common_harmonized <- common %>%
  
  filter(
    sign(meta_logFC_micro) ==
      sign(meta_logFC_rna)
  )

cat(
  "Concordant genes retained: ",
  nrow(common_harmonized),
  "\n",
  sep = ""
)

############################################################
# PRECISION-WEIGHTED CROSS-PLATFORM SYNTHESIS
#
# Weight = 1 / SE^2
#
# Combined effect:
#   sum(weight_i * effect_i) / sum(weight_i)
#
# Combined SE:
#   sqrt(1 / sum(weight_i))
############################################################

cat(
  "Performing precision-weighted cross-platform synthesis...\n"
)

common_meta <- common_harmonized %>%
  
  mutate(
    
    ########################################################
    # Inverse-variance weights
    ########################################################
    
    weight_micro =
      1 / (meta_SE_micro^2),
    
    weight_rna =
      1 / (meta_SE_rna^2),
    
    ########################################################
    # Total precision
    ########################################################
    
    total_weight =
      weight_micro +
      weight_rna,
    
    ########################################################
    # Precision-weighted pooled effect
    ########################################################
    
    combined_logFC =
      (
        weight_micro * meta_logFC_micro +
        weight_rna * meta_logFC_rna
      ) /
      total_weight,
    
    ########################################################
    # Standard error of pooled estimate
    ########################################################
    
    combined_SE =
      sqrt(
        1 / total_weight
      ),
    
    ########################################################
    # 95% confidence interval
    ########################################################
    
    CI_lb =
      combined_logFC -
      1.96 * combined_SE,
    
    CI_ub =
      combined_logFC +
      1.96 * combined_SE,
    
    ########################################################
    # Wald z statistic
    ########################################################
    
    z =
      combined_logFC /
      combined_SE,
    
    ########################################################
    # Two-sided p-value
    ########################################################
    
    combined_pval =
      2 * pnorm(
        -abs(z)
      )
  )

############################################################
# BENJAMINI-HOCHBERG FDR
############################################################

cat(
  "Applying Benjamini-Hochberg FDR correction...\n"
)

common_meta <- common_meta %>%
  
  mutate(
    FDR = p.adjust(
      combined_pval,
      method = "BH"
    )
  )

############################################################
# EXPORT ALL CROSS-PLATFORM RESULTS
############################################################

cross_platform_all <- common_meta %>%
  
  select(
    Gene,
    
    meta_logFC_micro,
    meta_SE_micro,
    meta_pval_micro,
    
    meta_logFC_rna,
    meta_SE_rna,
    meta_pval_rna,
    
    weight_micro,
    weight_rna,
    
    combined_logFC,
    combined_SE,
    
    CI_lb,
    CI_ub,
    
    z,
    combined_pval,
    FDR
  )

write.csv(
  cross_platform_all,
  "results/CrossPlatform_AllGenes_Weighted.csv",
  row.names = FALSE
)

############################################################
# SIGNIFICANT CROSS-PLATFORM GENES
#
# Criteria:
#   FDR < 0.05
#   |combined log2FC| >= 1
############################################################

final_genes <- cross_platform_all %>%
  
  filter(
    FDR < 0.05,
    abs(combined_logFC) >= 1
  ) %>%
  
  arrange(
    FDR
  )

############################################################
# UPREGULATED GENES
############################################################

final_up <- final_genes %>%
  
  filter(
    combined_logFC > 0
  )

############################################################
# DOWNREGULATED GENES
############################################################

final_down <- final_genes %>%
  
  filter(
    combined_logFC < 0
  )

############################################################
# EXPORT RESULTS
############################################################

write.csv(
  final_genes,
  "results/CrossPlatform_Weighted_genes.csv",
  row.names = FALSE
)

write.csv(
  final_up,
  "results/CrossPlatform_Weighted_UP.csv",
  row.names = FALSE
)

write.csv(
  final_down,
  "results/CrossPlatform_Weighted_DOWN.csv",
  row.names = FALSE
)

############################################################
# SUMMARY
############################################################

cat("\n")
cat("============================================================\n")
cat("CROSS-PLATFORM META-ANALYSIS COMPLETED\n")
cat("============================================================\n")

cat(
  "Genes shared between platforms: ",
  nrow(common),
  "\n",
  sep = ""
)

cat(
  "Concordant genes: ",
  nrow(common_harmonized),
  "\n",
  sep = ""
)

cat(
  "Significant genes (FDR < 0.05, |combined log2FC| >= 1): ",
  nrow(final_genes),
  "\n",
  sep = ""
)

cat(
  "Upregulated genes: ",
  nrow(final_up),
  "\n",
  sep = ""
)

cat(
  "Downregulated genes: ",
  nrow(final_down),
  "\n",
  sep = ""
)

cat("\nOutput files:\n")
cat("  results/CrossPlatform_AllGenes_Weighted.csv\n")
cat("  results/CrossPlatform_Weighted_genes.csv\n")
cat("  results/CrossPlatform_Weighted_UP.csv\n")
cat("  results/CrossPlatform_Weighted_DOWN.csv\n")

cat("============================================================\n")
