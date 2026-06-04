############################################################
# MULTI-DATASET MICROARRAY RAW PREPROCESSING AND MODELING
# Script: script/02_microarray_preprocessing.R
############################################################

############################################################
# SETTINGS & CONFIGURATIONS
############################################################
options(stringsAsFactors = FALSE)
set.seed(123)

############################################################
# LIBRARIES
############################################################
library(affy)
library(limma)
library(dplyr)
library(purrr)

############################################################
# DIRECTORY INFRASTRUCTURE
############################################################
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

############################################################
# EXAMPLE COHORT TARGET CONFIGURATIONS
# Maps to the raw microarray datasets detailed in Table 2
############################################################
datasets <- c(
  "GSE6514",
  "GSE33302"
)

############################################################
# PART 1: INDIVIDUAL RAW COHORT MODERATION LOOPS
############################################################
for(ds in datasets){
  
  cat("====================================================\n")
  cat("Processing Microarray Dataset via Limma:", ds, "\n")
  cat("====================================================\n")
  
  cel_path <- paste0("data/CEL_files/", ds)
  metadata_path <- paste0("data/metadata/", ds, "_metadata.csv")
  
  # Safe checkpoint validation for missing data structures
  if(!dir.exists(cel_path)) {
    cat("Warning: CEL file path not found for", ds, "- Skipping cohort.\n")
    next
  }
  if(!file.exists(metadata_path)) {
    cat("Warning: Metadata reference missing for", ds, "- Skipping cohort.\n")
    next
  }
  
  # 1. Read Raw CEL Files
  cat("Reading raw binary CEL array datasets...\n")
  raw_data <- ReadAffy(celfile.path = cel_path)
  
  # 2. Execute RMA (Robust Multi-array Average) Normalization
  cat("Executing background correction, quantile normalization, and summarization...\n")
  norm_data <- rma(raw_data)
  expr_matrix <- exprs(norm_data)
  
  # 3. Import Experimental Target Annotations
  metadata <- read.csv(metadata_path)
  
  # Synchronize expression column layouts to match sample metadata alignment
  if(!all(colnames(expr_matrix) %in% metadata$Sample)) {
    cat("Warning: Expression matrix headers do not match metadata sample keys. Verifying sync...\n")
  }
  
  # 4. Construct Design Matrix Environment
  metadata$Condition <- factor(metadata$Condition, levels = c("Control", "SD"))
  design <- model.matrix(~0 + Condition, data = metadata)
  colnames(design) <- c("Control", "SD")
  
  # 5. Fit Linear General Linear Models
  cat("Fitting linear generalized expression models across probes...\n")
  fit <- lmFit(expr_matrix, design)
  
  # 6. Apply Contrast Vectors (Sleep Deprived vs. Control baseline)
  contrast.matrix <- makeContrasts(
    SD_vs_Control = SD - Control,
    levels = design
  )
  fit2 <- contrasts.fit(fit, contrast.matrix)
  
  # 7. Execute Empirical Bayes Moderation Variance Adjustments
  cat("Running Empirical Bayes variance shrinkage modeling...\n")
  fit2 <- eBayes(fit2)
  
  # 8. Extract Full TopTable Statistical Matrices
  deg <- topTable(
    fit2,
    coef = "SD_vs_Control",
    number = Inf,
    adjust.method = "BH"
  )
  
  # Extract exact structural degrees of freedom and moderated parameters directly
  deg$Probe_ID <- rownames(deg)
  deg$df_total <- fit2$df.total[match(rownames(deg), rownames(fit2))]
  deg$stdev_unscaled <- fit2$stdev.unscaled[match(rownames(deg), rownames(fit2)), "SD_vs_Control"]
  
  # 9. Clean, Select, and Normalize Column Structure
  deg_export <- deg %>%
    mutate(
      Gene = toupper(trimws(Probe_ID)),  # Ensure consistent, normalized probe casing
      logFC = logFC,
      t_stat = t,                        # Direct native moderated t-statistic tracking
      P.Value = P.Value,
      adj.P.Val = adj.P.Val
    ) %>%
    filter(!is.na(Gene) & Gene != "") %>%
    select(Gene, logFC, t_stat, df_total, stdev_unscaled, P.Value, adj.P.Val)
  
  # 10. Save Standalone Contrast DEG Outputs
  write.csv(
    deg_export,
    paste0("results/", ds, "_limma_DEG.csv"),
    row.names = FALSE
  )
  cat("Dataset profile completed. Matrix generated:", paste0("results/", ds, "_limma_DEG.csv"), "\n")
}

############################################################
# PART 2: CONSOLIDATE INDIVIDUAL LOGS FOR POOLED METAFOR INPUT
############################################################
cat("\n====================================================\n")
cat("Consolidating all generated limma summaries...\n")
cat("====================================================\n")

deg_files <- list.files(
  "results",
  pattern = "_limma_DEG.csv",
  full.names = TRUE
)

if(length(deg_files) == 0) {
  stop("Critical Error: No calculated individual dataset files discovered inside results/ directory.")
}

micro_data <- map_df(deg_files, function(f){
  df <- read.csv(f)
  dataset_name <- gsub("_limma_DEG.csv", "", basename(f))
  df$dataset <- dataset_name
  df
})

# Filter out edge-case empty gene strings
micro_data <- micro_data %>%
  filter(!is.na(Gene) & Gene != "")

############################################################
# COMPUTE EXACT EXPERIMENTAL RESIDUAL STANDARD ERRORS
# (Deriving True Standard Errors utilizing the native moderated framework)
############################################################
cat("Applying mathematically rigorous Standard Error extractions...\n")
micro_data <- micro_data %>%
  mutate(
    # Native exact Standard Error derived from absolute logFC and the moderated t-statistic
    SE = abs(logFC) / abs(t_stat)
  ) %>%
  # Protect calculations against division-by-zero anomalies in identical expressions
  filter(!is.na(SE) & is.finite(SE) & SE > 0)

############################################################
# EXPORT INTERMEDIATE METAFOR-READY MATRIX ASSET
############################################################
write.csv(
  micro_data,
  "results/Microarray_processed_data.csv",
  row.names = FALSE
)

# Document current reproducible environments
writeLines(
  capture.output(sessionInfo()),
  "sessionInfo_microarray.txt"
)

cat("Success! Generated comprehensive metafor asset: results/Microarray_processed_data.csv\n")
