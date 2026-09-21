############################################################
# ROBUST RANDOM-EFFECTS MICROARRAY META-ANALYSIS (REML)
# Script: script/03_microarray_meta_analysis.R
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
# LOAD PROCESSED MICROARRAY DATA
############################################################

input_path <- "results/Microarray_processed_data.csv"

if(!file.exists(input_path)) {
  stop(
    "Critical Error: Intermediate file 'results/Microarray_processed_data.csv' ",
    "missing. Please execute Script 02 first."
  )
}

micro_data <- read.csv(
  input_path,
  stringsAsFactors = FALSE
)

############################################################
# CHECK REQUIRED COLUMNS
############################################################

required_columns <- c(
  "Gene",
  "logFC",
  "SE",
  "mean_expression",
  "dataset"
)

missing_columns <- setdiff(
  required_columns,
  colnames(micro_data)
)

if(length(missing_columns) > 0) {
  stop(
    "Critical Error: The following required columns are missing from ",
    "Microarray_processed_data.csv: ",
    paste(missing_columns, collapse = ", "),
    "\nPlease rerun Script 02 after adding mean_expression."
  )
}

############################################################
# CLEAN DATA AND REMOVE DUPLICATE PROBES
############################################################
#
# IMPORTANT:
# Probe selection is now independent of the observed
# sleep-deprivation effect.
#
# For genes represented by multiple probes within a dataset,
# the probe with the highest mean expression across samples
# is retained.
#
############################################################

cat(
  "Applying treatment-independent probe selection...\n"
)

micro_data_cleaned <- micro_data %>%
  
  mutate(
    Gene = toupper(trimws(Gene))
  ) %>%
  
  filter(
    !is.na(Gene),
    Gene != "",
    !is.na(logFC),
    !is.na(SE),
    is.finite(SE),
    SE > 0,
    !is.na(mean_expression),
    is.finite(mean_expression)
  ) %>%
  
  group_by(dataset, Gene) %>%
  
  slice_max(
    order_by = mean_expression,
    n = 1,
    with_ties = FALSE
  ) %>%
  
  ungroup()

cat(
  "Probe-independent dataset-level records retained: ",
  nrow(micro_data_cleaned),
  "\n",
  sep = ""
)

############################################################
# RANDOM-EFFECTS META-ANALYSIS
############################################################

cat(
  "Executing random-effects meta-analysis using REML ",
  "with Knapp-Hartung adjustment...\n"
)

micro_meta_res <- micro_data_cleaned %>%
  
  group_by(Gene) %>%
  
  group_split() %>%
  
  map_df(function(df) {
    
    current_gene <- df$Gene[1]
    
    # Number of independent datasets contributing
    study_count <- nrow(df)
    
    ########################################################
    # REQUIRE AT LEAST TWO INDEPENDENT DATASETS
    ########################################################
    
    if(study_count < 2) {
      return(NULL)
    }
    
    ########################################################
    # RANDOM-EFFECTS MODEL
    ########################################################
    
    model <- tryCatch(
      
      rma(
        yi   = df$logFC,
        sei  = df$SE,
        method = "REML",
        test = "knha"
      ),
      
      error = function(e) {
        return(NULL)
      }
    )
    
    ########################################################
    # SKIP FAILED MODELS
    ########################################################
    
    if(is.null(model)) {
      return(NULL)
    }
    
    ########################################################
    # PREDICTION INTERVAL
    ########################################################
    
    pred <- tryCatch(
      
      predict(model),
      
      error = function(e) {
        return(NULL)
      }
    )
    
    ########################################################
    # IF PREDICTION INTERVAL CANNOT BE CALCULATED
    ########################################################
    
    if(is.null(pred)) {
      
      PI_lb <- NA_real_
      PI_ub <- NA_real_
      
    } else {
      
      PI_lb <- as.numeric(pred$pi.lb)
      PI_ub <- as.numeric(pred$pi.ub)
      
    }
    
    ########################################################
    # RETURN META-ANALYSIS RESULTS
    ########################################################
    
    data.frame(
      
      Gene = current_gene,
      
      # Pooled effect
      meta_logFC = as.numeric(model$b),
      
      # Standard error of pooled effect
      meta_SE = as.numeric(model$se),
      
      # 95% confidence interval
      CI_lb = as.numeric(model$ci.lb),
      CI_ub = as.numeric(model$ci.ub),
      
      # 95% prediction interval
      PI_lb = PI_lb,
      PI_ub = PI_ub,
      
      # Statistical significance
      meta_pval = as.numeric(model$pval),
      
      # Heterogeneity
      I2 = as.numeric(model$I2),
      tau2 = as.numeric(model$tau2),
      
      # Number of contributing datasets
      n_studies = study_count
      
    )
    
  })

############################################################
# VALIDATE META-ANALYSIS OUTPUT
############################################################

if(nrow(micro_meta_res) == 0) {
  
  stop(
    "Critical Error: Meta-analysis yielded zero output rows. ",
    "Check source data and SE values."
  )
}

cat(
  "Meta-analysis completed for ",
  nrow(micro_meta_res),
  " genes.\n",
  sep = ""
)

############################################################
# BENJAMINI-HOCHBERG FDR CORRECTION
############################################################

cat(
  "Applying Benjamini-Hochberg FDR correction...\n"
)

micro_meta_res$FDR <- p.adjust(
  micro_meta_res$meta_pval,
  method = "BH"
)

############################################################
# SAVE COMPLETE META-ANALYSIS RESULTS
############################################################

write.csv(
  micro_meta_res,
  "results/Microarray_Meta_AllGenes.csv",
  row.names = FALSE
)

cat(
  "Complete microarray meta-analysis saved to:\n",
  "results/Microarray_Meta_AllGenes.csv\n"
)

############################################################
# SIGNIFICANT MICROARRAY GENES
############################################################
#
# Manuscript threshold:
# FDR < 0.05
# absolute meta-analysis log2FC >= 1
#
############################################################

sig_micro <- micro_meta_res %>%
  
  filter(
    FDR < 0.05,
    abs(meta_logFC) >= 1
  ) %>%
  
  arrange(FDR)

############################################################
# SAVE SIGNIFICANT RESULTS
############################################################

write.csv(
  sig_micro,
  "results/Microarray_FDR_FC1.csv",
  row.names = FALSE
)

############################################################
# SUMMARY
############################################################

cat(
  "============================================================\n"
)

cat(
  "Microarray random-effects meta-analysis completed.\n"
)

cat(
  "Total genes meta-analyzed: ",
  nrow(micro_meta_res),
  "\n",
  sep = ""
)

cat(
  "Significant genes (FDR < 0.05 and |log2FC| >= 1): ",
  nrow(sig_micro),
  "\n",
  sep = ""
)

cat(
  "============================================================\n"
)
