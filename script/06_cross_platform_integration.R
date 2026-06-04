############################################################
# CROSS-PLATFORM INTEGRATION VIA FISHER'S COMBINED TEST
# Script: script/07_cross_platform_meta_analysis.R
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
# LOAD THE INDEPENDENT PLATFORM-LEVEL REML TARGETS
############################################################
micro_path <- "results/Microarray_Meta_AllGenes.csv"
rna_path   <- "results/RNAseq_Meta_AllGenes.csv"

if(!file.exists(micro_path) | !file.exists(rna_path)) {
  stop("Critical Error: Platform-specific meta-analysis input assets missing. Execute scripts 03 and 06 first.")
}

micro <- read.csv(micro_path)
rna   <- read.csv(rna_path)

############################################################
# NOMENCLATURE CASE NORMALIZATION & INITIAL JOIN
############################################################
cat("Synchronizing nomenclature patterns across platforms...\n")

micro_clean <- micro %>% 
  mutate(Gene = toupper(trimws(Gene))) %>% 
  filter(!is.na(Gene) & Gene != "")

rna_clean <- rna %>% 
  mutate(Gene = toupper(trimws(Gene))) %>% 
  filter(!is.na(Gene) & Gene != "")

# Merge both individual platforms using a strict Inner Join
common <- inner_join(micro_clean, rna_clean, by = "Gene")

cat(paste("Intersection analysis discovered", nrow(common), "overlapping genes shared between platforms.\n"))

############################################################
# ENFORCE DIRECTIONAL HARMONY CONSTRAINT
# (CRITICAL: Prevents opposing biological regulation anomalies)
############################################################
cat("Applying strict directional harmony filters...\n")

common_harmonized <- common %>%
  # Retain only genes moving in the identical direction (both positive or both negative logFC)
  filter(sign(meta_logFC.x) == sign(meta_logFC.y))

cat(paste("Directional check complete. Retained", nrow(common_harmonized), 
          "genes showing cross-platform biological consensus.\n"))

############################################################
# COMPUTE GLOBAL FISHER STATISTICAL PORTFOLIO
############################################################
cat("Executing Fisher's Combined Probability calculations...\n")

common_meta <- common_harmonized %>%
  mutate(
    # Unweighted integration average of effect sizes across platforms
    combined_logFC = (meta_logFC.x + meta_logFC.y) / 2,
    
    # Fisher's test statistic calculation formula: -2 * sum(ln(p))
    # Employs a low-end baseline ceiling to prevent infinite value errors in p = 0 anomalies
    p_micro_safe = ifelse(meta_pval.x == 0, 1e-300, meta_pval.x),
    p_rna_safe   = ifelse(meta_pval.y == 0, 1e-300, meta_pval.y),
    
    fisher_stat = -2 * (log(p_micro_safe) + log(p_rna_safe)),
    
    # Map the statistic directly to the Chi-Square distribution 
    # Degrees of Freedom (df) = 2 * k meta-analyses (k = 2 distinct platforms, df = 4)
    fisher_p = pchisq(fisher_stat, df = 4, lower.tail = FALSE)
  )

############################################################
# INTERMEDIATE RE-APPLICATION OF GLOBAL FDR CORRECTION
# (The Subsequent Fisher-to-FDR Transition step)
############################################################
cat("Applying Benjamini-Hochberg multi-testing correction on combined probabilities...\n")
common_meta$FDR <- p.adjust(common_meta$fisher_p, method = "BH")

# Save the un-filtered combined platform asset for global repository record tracking
write.csv(
  common_meta %>% select(Gene, meta_logFC.x, meta_pval.x, meta_logFC.y, meta_pval.y, combined_logFC, fisher_stat, fisher_p, FDR),
  "results/CrossPlatform_AllGenes_Combined.csv",
  row.names = FALSE
)

############################################################
# ISOLATE STRICT CROSS-PLATFORM SIGNIFICANT DISCOVERY LIST
############################################################
final_genes <- common_meta %>%
  filter(FDR < 0.05 & abs(combined_logFC) >= 1)

# Export candidate cohort to matching file name referenced in your manuscript text
write.csv(
  final_genes,
  "results/CrossPlatform_901_genes.csv",
  row.names = FALSE
)

cat("========================================================================\n")
cat(paste("Success! Final dataset generated: results/CrossPlatform_901_genes.csv\n"))
cat(paste("Identified exactly", nrow(final_genes), "robust cross-platform candidates.\n"))
cat("========================================================================\n")
