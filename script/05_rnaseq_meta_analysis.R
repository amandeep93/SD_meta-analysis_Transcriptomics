############################################################
# ROBUST RANDOM-EFFECTS RNA-SEQ META-ANALYSIS (REML)
# Script: script/06_rnaseq_meta_analysis.R
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
library(purrr)
library(metafor)

############################################################
# CONSOLIDATE DYNAMIC DESEQ2 EXTRACTS
############################################################
cat("Gathers standalone DESeq2 analytical files...\n")
rna_files <- list.files(
  "results/rna_individual", 
  pattern = "_DESeq2_DEG.csv", 
  full.names = TRUE
)

if(length(rna_files) == 0) {
  stop("Critical Error: No standalone individual DESeq2 files discovered in results/rna_individual/. Execute script 05 first.")
}

# Combine files across cohorts into a single operational dataset dataframe
rna_data <- map_df(rna_files, read.csv)

############################################################
# SANITIZATION AND STRUCTURAL INDEPENDENCE ENFORCEMENT
############################################################
cat("Applying strict multi-transcript/isoform independence controls...\n")
rna_data_cleaned <- rna_data %>%
  mutate(Gene = toupper(trimws(Gene))) %>%
  filter(!is.na(Gene) & Gene != "" & !is.na(log2FoldChange) & !is.na(lfcSE)) %>%
  # If any overlapping transcript profiles remain for a gene within a single cohort,
  # subset to retain the single variant with maximum absolute change to safeguard independence.
  group_by(dataset, Gene) %>%
  filter(abs(log2FoldChange) == max(abs(log2FoldChange))) %>%
  slice(1) %>%
  ungroup()

############################################################
# EXECUTE PARALLELIZED REML META-ANALYSIS LOOP
############################################################
cat("Executing random-effects meta-analysis via native DESeq2 standard errors...\n")

rna_meta_res <- rna_data_cleaned %>%
  group_by(Gene) %>%
  group_split() %>%
  map_df(function(df){
    
    current_gene <- df$Gene[1]
    study_count  <- nrow(df)
    
    # Must be present across at least 2 separate cohorts to run cross-study pooling
    if(study_count < 2) return(NULL)
    
    # Execute random-effects engine with Knapp-Hartung adjustments using native lfcSE parameters
    model <- tryCatch(
      rma(
        yi = df$log2FoldChange,
        sei = df$lfcSE, # CRITICAL: Passes true native standard errors
        method = "REML",
        test = "knha"
      ),
      error = function(e) NULL
    )
    
    if(is.null(model)) return(NULL)
    
    data.frame(
      Gene        = current_gene,
      meta_logFC  = as.numeric(model$b),
      CI_lb       = model$ci.lb,
      CI_ub       = model$ci.ub,
      meta_pval   = model$pval,
      I2          = model$I2,
      tau2        = model$tau2,
      n_studies   = study_count
    )
  })

if(nrow(rna_meta_res) == 0) {
  stop("Critical Error: RNA-seq meta-analysis returned empty matrix. Check individual file values.")
}

############################################################
# MULTIPLE TESTING BENJAMINI-HOCHBERG ADJUSTMENTS
############################################################
cat("Applying Benjamini-Hochberg FDR correction parameters...\n")
rna_meta_res$FDR <- p.adjust(
  rna_meta_res$meta_pval,
  method = "BH"
)

############################################################
# EXPORT POOLED TRANSCRIPTOMIC PORTFOLIOS
############################################################
# This file serves as the core entry required by your upcoming cross-platform scripts!
write.csv(
  rna_meta_res,
  "results/RNAseq_Meta_AllGenes.csv",
  row.names = FALSE
)

# Extract significant candidate subset to separate matrix for independent inspection
sig_rna <- rna_meta_res %>%
  filter(FDR < 0.05 & abs(meta_logFC) >= 1)

write.csv(
  sig_rna,
  "results/RNAseq_FDR_FC1.csv",
  row.names = FALSE
)

cat(paste("Success! Completed RNA-seq meta-analysis. Identified", nrow(sig_rna), "highly significant candidates.\n"))
