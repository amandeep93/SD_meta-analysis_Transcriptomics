############################################################
# ROBUST RANDOM-EFFECTS RNA-SEQ META-ANALYSIS
# REML + Knapp-Hartung
#
# Script: script/06_rnaseq_meta_analysis.R
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
library(purrr)
library(metafor)

############################################################
# DIRECTORY
############################################################

dir.create(
  "results",
  showWarnings = FALSE,
  recursive = TRUE
)

############################################################
# 1. COLLECT INDIVIDUAL DESEQ2 RESULTS
############################################################

cat(
  "Collecting standalone DESeq2 analytical files...\n"
)

rna_files <- list.files(
  "results/rna_individual",
  pattern = "_DESeq2_DEG\\.csv$",
  full.names = TRUE
)

if (length(rna_files) == 0) {
  
  stop(
    "Critical Error: No individual DESeq2 files found in ",
    "results/rna_individual/. Run Script 05 first."
  )
}

cat(
  "Number of individual RNA-seq result files found: ",
  length(rna_files),
  "\n",
  sep = ""
)

############################################################
# Read and combine all datasets
############################################################

rna_data <- map_df(
  rna_files,
  ~ read.csv(
    .x,
    stringsAsFactors = FALSE
  )
)

############################################################
# 2. CHECK REQUIRED COLUMNS
############################################################

required_columns <- c(
  "Gene",
  "log2FoldChange",
  "lfcSE",
  "dataset"
)

missing_columns <- setdiff(
  required_columns,
  colnames(rna_data)
)

if (length(missing_columns) > 0) {
  
  stop(
    "Required columns missing from RNA-seq DESeq2 results: ",
    paste(
      missing_columns,
      collapse = ", "
    )
  )
}

############################################################
# 3. CLEAN DATA
############################################################

cat(
  "Cleaning and validating study-level RNA-seq estimates...\n"
)

rna_data_cleaned <- rna_data %>%
  
  mutate(
    Gene = toupper(
      trimws(Gene)
    ),
    
    dataset = trimws(
      dataset
    )
  ) %>%
  
  filter(
    !is.na(Gene),
    Gene != "",
    !is.na(log2FoldChange),
    !is.na(lfcSE),
    is.finite(log2FoldChange),
    is.finite(lfcSE),
    lfcSE > 0,
    !is.na(dataset),
    dataset != ""
  )

############################################################
# 4. ENFORCE ONE GENE ESTIMATE PER DATASET
#
# DESeq2/featureCounts should already provide gene-level
# estimates. Therefore, duplicate dataset-Gene records are
# not resolved by selecting the largest effect.
#
# If duplicates occur, retain only a unique record when
# the duplicated entries are identical.
# Otherwise stop and investigate.
############################################################

duplicate_check <- rna_data_cleaned %>%
  
  count(
    dataset,
    Gene,
    name = "n"
  ) %>%
  
  filter(
    n > 1
  )

if (nrow(duplicate_check) > 0) {
  
  cat(
    "Duplicate dataset-Gene combinations detected.\n"
  )
  
  print(
    duplicate_check
  )
  
  stop(
    "Duplicate dataset-Gene estimates detected. ",
    "Do not select the largest effect automatically. ",
    "Investigate the underlying DESeq2 files."
  )
}

############################################################
# 5. RANDOM-EFFECTS META-ANALYSIS
#
# One estimate per independent RNA-seq dataset.
#
# Effect size:
#   DESeq2 log2FoldChange
#
# Standard error:
#   Native DESeq2 lfcSE
#
# Model:
#   REML random effects
#
# Inference:
#   Knapp-Hartung adjustment
############################################################

cat(
  "Executing REML random-effects meta-analysis...\n"
)

rna_meta_res <- rna_data_cleaned %>%
  
  group_by(Gene) %>%
  
  group_split() %>%
  
  map_df(
    
    function(df) {
      
      current_gene <- df$Gene[1]
      
      ######################################################
      # Number of independent datasets
      ######################################################
      
      study_count <- n_distinct(
        df$dataset
      )
      
      ######################################################
      # Require at least two independent datasets
      ######################################################
      
      if (study_count < 2) {
        return(NULL)
      }
      
      ######################################################
      # Random-effects model
      ######################################################
      
      model <- tryCatch(
        
        rma(
          yi = df$log2FoldChange,
          sei = df$lfcSE,
          method = "REML",
          test = "knha"
        ),
        
        error = function(e) {
          
          cat(
            "Meta-analysis failed for ",
            current_gene,
            ": ",
            e$message,
            "\n",
            sep = ""
          )
          
          return(NULL)
        }
      )
      
      if (is.null(model)) {
        return(NULL)
      }
      
      ######################################################
      # Prediction interval
      ######################################################
      
      pred <- tryCatch(
        
        predict(model),
        
        error = function(e) {
          return(NULL)
        }
      )
      
      ######################################################
      # Extract prediction interval
      ######################################################
      
      if (!is.null(pred)) {
        
        PI_lb <- as.numeric(
          pred$pi.lb
        )
        
        PI_ub <- as.numeric(
          pred$pi.ub
        )
        
      } else {
        
        PI_lb <- NA_real_
        PI_ub <- NA_real_
      }
      
      ######################################################
      # Return meta-analysis statistics
      ######################################################
      
      data.frame(
        
        Gene = current_gene,
        
        meta_logFC = as.numeric(
          model$b
        ),
        
        meta_SE = as.numeric(
          model$se
        ),
        
        CI_lb = as.numeric(
          model$ci.lb
        ),
        
        CI_ub = as.numeric(
          model$ci.ub
        ),
        
        PI_lb = PI_lb,
        
        PI_ub = PI_ub,
        
        meta_pval = as.numeric(
          model$pval
        ),
        
        I2 = as.numeric(
          model$I2
        ),
        
        tau2 = as.numeric(
          model$tau2
        ),
        
        n_studies = study_count,
        
        stringsAsFactors = FALSE
      )
    }
  )

############################################################
# 6. CHECK META-ANALYSIS OUTPUT
############################################################

if (
  nrow(rna_meta_res) == 0
) {
  
  stop(
    "Critical Error: RNA-seq meta-analysis returned ",
    "no results. Check individual DESeq2 files."
  )
}

cat(
  "RNA-seq meta-analysis generated estimates for ",
  nrow(rna_meta_res),
  " genes.\n",
  sep = ""
)

############################################################
# 7. BENJAMINI-HOCHBERG FDR
############################################################

cat(
  "Applying Benjamini-Hochberg FDR correction...\n"
)

rna_meta_res <- rna_meta_res %>%
  
  mutate(
    FDR = p.adjust(
      meta_pval,
      method = "BH"
    )
  )

############################################################
# 8. EXPORT ALL META-ANALYSIS RESULTS
############################################################

write.csv(
  rna_meta_res,
  "results/RNAseq_Meta_AllGenes.csv",
  row.names = FALSE
)

############################################################
# 9. IDENTIFY SIGNIFICANT GENES
#
# Criteria:
#   FDR < 0.05
#   absolute pooled log2FC >= 1
############################################################

sig_rna <- rna_meta_res %>%
  
  filter(
    FDR < 0.05,
    abs(meta_logFC) >= 1
  ) %>%
  
  arrange(
    FDR
  )

############################################################
# 10. EXPORT SIGNIFICANT RESULTS
############################################################

write.csv(
  sig_rna,
  "results/RNAseq_FDR_FC1.csv",
  row.names = FALSE
)

############################################################
# 11. EXPORT UP- AND DOWN-REGULATED GENES
############################################################

sig_rna_up <- sig_rna %>%
  
  filter(
    meta_logFC > 0
  )

sig_rna_down <- sig_rna %>%
  
  filter(
    meta_logFC < 0
  )

write.csv(
  sig_rna_up,
  "results/RNAseq_FDR_FC1_UP.csv",
  row.names = FALSE
)

write.csv(
  sig_rna_down,
  "results/RNAseq_FDR_FC1_DOWN.csv",
  row.names = FALSE
)

############################################################
# 12. SUMMARY
############################################################

cat("\n")
cat("====================================================\n")
cat("RNA-SEQ META-ANALYSIS COMPLETED\n")
cat("====================================================\n")

cat(
  "Individual DESeq2 files: ",
  length(rna_files),
  "\n",
  sep = ""
)

cat(
  "Genes with >=2 independent datasets: ",
  nrow(rna_meta_res),
  "\n",
  sep = ""
)

cat(
  "Significant genes (FDR < 0.05, |log2FC| >= 1): ",
  nrow(sig_rna),
  "\n",
  sep = ""
)

cat(
  "Upregulated genes: ",
  nrow(sig_rna_up),
  "\n",
  sep = ""
)

cat(
  "Downregulated genes: ",
  nrow(sig_rna_down),
  "\n",
  sep = ""
)

cat("\nOutput files:\n")
cat("  results/RNAseq_Meta_AllGenes.csv\n")
cat("  results/RNAseq_FDR_FC1.csv\n")
cat("  results/RNAseq_FDR_FC1_UP.csv\n")
cat("  results/RNAseq_FDR_FC1_DOWN.csv\n")

cat("====================================================\n")
