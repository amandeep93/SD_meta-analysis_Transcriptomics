############################################################
# BRAIN-REGION STRATIFIED CROSS-PLATFORM META-ANALYSIS
#
# Script:
# script/07_brain_region_crossplatform_analysis.R
#
# Workflow:
#
# Individual study estimates
#          ↓
# Region-specific microarray REML meta-analysis
#          ↓
# Region-specific RNA-seq REML meta-analysis
#          ↓
# Cross-platform regional REML synthesis
#          ↓
# BH FDR correction
#
# Statistical framework:
#   - REML random-effects meta-analysis
#   - Knapp-Hartung inference
#   - Native study-level SE
#   - Region-specific analysis
#   - Concordant direction across platforms
#
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
library(tidyr)
library(purrr)
library(metafor)

############################################################
# DIRECTORIES
############################################################

dir.create(
  "results/brain_region",
  showWarnings = FALSE,
  recursive = TRUE
)

dir.create(
  "figures/brain_region",
  showWarnings = FALSE,
  recursive = TRUE
)

############################################################
# INPUT FILES
############################################################

micro_path <-
  "results/Microarray_processed_data.csv"

############################################################
# RNA-SEQ INDIVIDUAL RESULTS
############################################################

rna_dir <-
  "results/rna_individual"

############################################################
# METADATA DIRECTORY
############################################################

metadata_dir <-
  "data/metadata"

############################################################
# CHECK INPUTS
############################################################

if (!file.exists(micro_path)) {
  
  stop(
    "Microarray processed data not found: ",
    micro_path
  )
}

if (!dir.exists(rna_dir)) {
  
  stop(
    "RNA-seq results directory not found: ",
    rna_dir
  )
}

if (!dir.exists(metadata_dir)) {
  
  stop(
    "Metadata directory not found: ",
    metadata_dir
  )
}

############################################################
# 1. LOAD MICROARRAY STUDY-LEVEL DATA
############################################################

cat(
  "Loading microarray study-level differential-expression data...\n"
)

micro <- read.csv(
  micro_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

############################################################
# REQUIRED MICROARRAY COLUMNS
############################################################

required_micro <- c(
  "Gene",
  "dataset",
  "logFC",
  "SE"
)

missing_micro <-
  setdiff(
    required_micro,
    colnames(micro)
  )

if (length(missing_micro) > 0) {
  
  stop(
    "Microarray file is missing required columns: ",
    paste(
      missing_micro,
      collapse = ", "
    )
  )
}

############################################################
# 2. LOAD MICROARRAY METADATA
############################################################

cat(
  "Loading brain-region metadata...\n"
)

metadata_files <- list.files(
  metadata_dir,
  pattern = "_metadata\\.csv$",
  full.names = TRUE
)

if (length(metadata_files) == 0) {
  
  stop(
    "No dataset-specific metadata files found in ",
    metadata_dir
  )
}

metadata_list <- lapply(
  metadata_files,
  function(x) {
    
    df <- read.csv(
      x,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    
    df$metadata_file <-
      basename(x)
    
    df
  }
)

metadata <- bind_rows(
  metadata_list
)

############################################################
# CHECK BRAIN REGION FIELD
############################################################

if (!"brain_region" %in% colnames(metadata)) {
  
  stop(
    "The metadata files do not contain a 'brain_region' column.\n",
    "Regional analysis cannot be performed reliably without ",
    "region-level metadata."
  )
}

############################################################
# CHECK DATASET FIELD
############################################################

if (!"dataset" %in% colnames(metadata)) {
  
  stop(
    "The metadata files do not contain a 'dataset' column."
  )
}

############################################################
# CLEAN METADATA
############################################################

metadata <- metadata %>%
  
  mutate(
    dataset = trimws(dataset),
    brain_region = trimws(brain_region)
  ) %>%
  
  filter(
    !is.na(dataset),
    dataset != "",
    !is.na(brain_region),
    brain_region != ""
  ) %>%
  
  distinct(
    dataset,
    brain_region,
    .keep_all = TRUE
  )

############################################################
# 3. STANDARDIZE MICROARRAY DATA
############################################################

micro <- micro %>%
  
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
    !is.na(dataset),
    dataset != "",
    is.finite(logFC),
    is.finite(SE),
    SE > 0
  )

############################################################
# 4. SELECT ONE MICROARRAY PROBE PER
#    DATASET-GENE
#
# IMPORTANT:
# Use mean expression rather than treatment-associated
# effect size for probe selection.
############################################################

if (!"mean_expression" %in% colnames(micro)) {
  
  stop(
    "Microarray data do not contain 'mean_expression'. ",
    "Use the revised Script 02 before running this analysis."
  )
}

micro <- micro %>%
  
  group_by(
    dataset,
    Gene
  ) %>%
  
  arrange(
    desc(mean_expression),
    .by_group = TRUE
  ) %>%
  
  slice(1) %>%
  
  ungroup()

############################################################
# 5. ADD BRAIN REGION TO MICROARRAY DATA
############################################################

micro_region <- micro %>%
  
  left_join(
    metadata %>%
      select(
        dataset,
        brain_region
      ) %>%
      distinct(),
    
    by = "dataset"
  )

############################################################
# CHECK UNMATCHED MICROARRAY DATASETS
############################################################

unmatched_micro <- micro_region %>%
  
  filter(
    is.na(brain_region)
  ) %>%
  
  distinct(
    dataset
  )

if (nrow(unmatched_micro) > 0) {
  
  cat(
    "\nWARNING: The following microarray datasets have no ",
    "brain-region annotation:\n"
  )
  
  print(
    unmatched_micro
  )
}

micro_region <- micro_region %>%
  
  filter(
    !is.na(brain_region)
  )

############################################################
# 6. LOAD RNA-SEQ INDIVIDUAL DESEQ2 RESULTS
############################################################

cat(
  "\nLoading individual RNA-seq DESeq2 results...\n"
)

rna_files <- list.files(
  rna_dir,
  pattern = "_DESeq2_DEG\\.csv$",
  full.names = TRUE
)

if (length(rna_files) == 0) {
  
  stop(
    "No RNA-seq DESeq2 result files found."
  )
}

rna <- map_df(
  rna_files,
  ~ read.csv(
    .x,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
)

############################################################
# CHECK RNA-SEQ COLUMNS
############################################################

required_rna <- c(
  "Gene",
  "dataset",
  "log2FoldChange",
  "lfcSE"
)

missing_rna <-
  setdiff(
    required_rna,
    colnames(rna)
  )

if (length(missing_rna) > 0) {
  
  stop(
    "RNA-seq results are missing required columns: ",
    paste(
      missing_rna,
      collapse = ", "
    )
  )
}

############################################################
# STANDARDIZE RNA-SEQ DATA
############################################################

rna <- rna %>%
  
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
    !is.na(dataset),
    dataset != "",
    is.finite(log2FoldChange),
    is.finite(lfcSE),
    lfcSE > 0
  )

############################################################
# CHECK DUPLICATE DATASET-GENE COMBINATIONS
############################################################

rna_duplicates <- rna %>%
  
  count(
    dataset,
    Gene
  ) %>%
  
  filter(
    n > 1
  )

if (nrow(rna_duplicates) > 0) {
  
  stop(
    "Duplicate dataset-Gene combinations detected in ",
    "RNA-seq results. Do not select the largest effect ",
    "automatically; investigate the underlying results."
  )
}

############################################################
# 7. ADD BRAIN REGION TO RNA-SEQ DATA
############################################################

rna_region <- rna %>%
  
  left_join(
    metadata %>%
      select(
        dataset,
        brain_region
      ) %>%
      distinct(),
    
    by = "dataset"
  )

############################################################
# CHECK UNMATCHED RNA-SEQ DATASETS
############################################################

unmatched_rna <- rna_region %>%
  
  filter(
    is.na(brain_region)
  ) %>%
  
  distinct(
    dataset
  )

if (nrow(unmatched_rna) > 0) {
  
  cat(
    "\nWARNING: The following RNA-seq datasets have no ",
    "brain-region annotation:\n"
  )
  
  print(
    unmatched_rna
  )
}

rna_region <- rna_region %>%
  
  filter(
    !is.na(brain_region)
  )

############################################################
# 8. REGION-SPECIFIC MICROARRAY META-ANALYSIS
############################################################

cat(
  "\nRunning region-specific microarray random-effects ",
  "meta-analysis...\n"
)

micro_region_meta <- micro_region %>%
  
  group_by(
    brain_region,
    Gene
  ) %>%
  
  group_split() %>%
  
  map_df(
    
    function(df) {
      
      region <- df$brain_region[1]
      gene <- df$Gene[1]
      
      n_studies <- n_distinct(
        df$dataset
      )
      
      ######################################################
      # Require >=2 independent datasets
      ######################################################
      
      if (n_studies < 2) {
        return(NULL)
      }
      
      ######################################################
      # REML RANDOM-EFFECTS MODEL
      ######################################################
      
      model <- tryCatch(
        
        rma(
          yi = df$logFC,
          sei = df$SE,
          method = "REML",
          test = "knha"
        ),
        
        error = function(e) {
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
        error = function(e) NULL
      )
      
      if (!is.null(pred)) {
        PI_lb <- as.numeric(pred$pi.lb)
        PI_ub <- as.numeric(pred$pi.ub)
      } else {
        PI_lb <- NA_real_
        PI_ub <- NA_real_
      }
      
      data.frame(
        
        brain_region = region,
        Gene = gene,
        
        micro_logFC = as.numeric(model$b),
        micro_SE = as.numeric(model$se),
        
        micro_CI_lb = as.numeric(model$ci.lb),
        micro_CI_ub = as.numeric(model$ci.ub),
        
        micro_PI_lb = PI_lb,
        micro_PI_ub = PI_ub,
        
        micro_p = as.numeric(model$pval),
        
        micro_I2 = as.numeric(model$I2),
        micro_tau2 = as.numeric(model$tau2),
        
        micro_n_studies = n_studies,
        
        stringsAsFactors = FALSE
      )
    }
  )

############################################################
# 9. REGION-SPECIFIC RNA-SEQ META-ANALYSIS
############################################################

cat(
  "Running region-specific RNA-seq random-effects ",
  "meta-analysis...\n"
)

rna_region_meta <- rna_region %>%
  
  group_by(
    brain_region,
    Gene
  ) %>%
  
  group_split() %>%
  
  map_df(
    
    function(df) {
      
      region <- df$brain_region[1]
      gene <- df$Gene[1]
      
      n_studies <- n_distinct(
        df$dataset
      )
      
      ######################################################
      # Require >=2 independent datasets
      ######################################################
      
      if (n_studies < 2) {
        return(NULL)
      }
      
      ######################################################
      # REML RANDOM-EFFECTS MODEL
      ######################################################
      
      model <- tryCatch(
        
        rma(
          yi = df$log2FoldChange,
          sei = df$lfcSE,
          method = "REML",
          test = "knha"
        ),
        
        error = function(e) {
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
        error = function(e) NULL
      )
      
      if (!is.null(pred)) {
        PI_lb <- as.numeric(pred$pi.lb)
        PI_ub <- as.numeric(pred$pi.ub)
      } else {
        PI_lb <- NA_real_
        PI_ub <- NA_real_
      }
      
      data.frame(
        
        brain_region = region,
        Gene = gene,
        
        rna_logFC = as.numeric(model$b),
        rna_SE = as.numeric(model$se),
        
        rna_CI_lb = as.numeric(model$ci.lb),
        rna_CI_ub = as.numeric(model$ci.ub),
        
        rna_PI_lb = PI_lb,
        rna_PI_ub = PI_ub,
        
        rna_p = as.numeric(model$pval),
        
        rna_I2 = as.numeric(model$I2),
        rna_tau2 = as.numeric(model$tau2),
        
        rna_n_studies = n_studies,
        
        stringsAsFactors = FALSE
      )
    }
  )

############################################################
# CHECK REGION-SPECIFIC RESULTS
############################################################

if (nrow(micro_region_meta) == 0) {
  
  warning(
    "No region-specific microarray meta-analysis results."
  )
}

if (nrow(rna_region_meta) == 0) {
  
  warning(
    "No region-specific RNA-seq meta-analysis results."
  )
}

############################################################
# 10. CROSS-PLATFORM REGIONAL MERGE
############################################################

cat(
  "\nMerging regional microarray and RNA-seq estimates...\n"
)

regional_common <- inner_join(
  
  micro_region_meta,
  
  rna_region_meta,
  
  by = c(
    "brain_region",
    "Gene"
  )
)

cat(
  "Region-gene combinations represented by both platforms: ",
  nrow(regional_common),
  "\n",
  sep = ""
)

if (nrow(regional_common) == 0) {
  
  stop(
    "No region-gene combinations were shared between ",
    "microarray and RNA-seq."
  )
}

############################################################
# 11. CROSS-PLATFORM DIRECTIONAL CONCORDANCE
############################################################

regional_concordant <- regional_common %>%
  
  filter(
    sign(micro_logFC) ==
      sign(rna_logFC)
  )

cat(
  "Concordant region-gene combinations: ",
  nrow(regional_concordant),
  "\n",
  sep = ""
)

############################################################
# 12. REGION-SPECIFIC CROSS-PLATFORM REML SYNTHESIS
############################################################
#
# Each region-gene combination contains:
#
#   microarray pooled effect + SE
#   RNA-seq pooled effect + SE
#
# These two platform estimates are then synthesized using
# a REML random-effects model.
#
############################################################

cat(
  "\nPerforming regional cross-platform REML synthesis...\n"
)

regional_crossplatform <- regional_concordant %>%
  
  group_by(
    brain_region,
    Gene
  ) %>%
  
  group_split() %>%
  
  map_df(
    
    function(df) {
      
      region <- df$brain_region[1]
      gene <- df$Gene[1]
      
      effects <- c(
        df$micro_logFC[1],
        df$rna_logFC[1]
      )
      
      ses <- c(
        df$micro_SE[1],
        df$rna_SE[1]
      )
      
      ######################################################
      # REML model across the two platforms
      ######################################################
      
      model <- tryCatch(
        
        rma(
          yi = effects,
          sei = ses,
          method = "REML",
          test = "knha"
        ),
        
        error = function(e) NULL
      )
      
      if (is.null(model)) {
        return(NULL)
      }
      
      ######################################################
      # Prediction interval
      ######################################################
      
      pred <- tryCatch(
        predict(model),
        error = function(e) NULL
      )
      
      if (!is.null(pred)) {
        PI_lb <- as.numeric(pred$pi.lb)
        PI_ub <- as.numeric(pred$pi.ub)
      } else {
        PI_lb <- NA_real_
        PI_ub <- NA_real_
      }
      
      data.frame(
        
        brain_region = region,
        Gene = gene,
        
        combined_logFC =
          as.numeric(model$b),
        
        combined_SE =
          as.numeric(model$se),
        
        combined_CI_lb =
          as.numeric(model$ci.lb),
        
        combined_CI_ub =
          as.numeric(model$ci.ub),
        
        combined_PI_lb =
          PI_lb,
        
        combined_PI_ub =
          PI_ub,
        
        combined_p =
          as.numeric(model$pval),
        
        combined_I2 =
          as.numeric(model$I2),
        
        combined_tau2 =
          as.numeric(model$tau2),
        
        platform_n = 2,
        
        micro_n_studies =
          df$micro_n_studies[1],
        
        rna_n_studies =
          df$rna_n_studies[1],
        
        stringsAsFactors = FALSE
      )
    }
  )

############################################################
# 13. BENJAMINI-HOCHBERG CORRECTION
#
# Correction is performed across all tested
# region-gene combinations.
############################################################

regional_crossplatform <- regional_crossplatform %>%
  
  mutate(
    FDR = p.adjust(
      combined_p,
      method = "BH"
    )
  )

############################################################
# 14. EXPORT COMPLETE REGIONAL RESULTS
############################################################

write.csv(
  regional_crossplatform,
  "results/brain_region/BrainRegion_CrossPlatform_AllGenes.csv",
  row.names = FALSE
)

############################################################
# 15. SIGNIFICANT REGIONAL GENES
#
# Same threshold as the main analysis:
#
# FDR < 0.05
# |combined log2FC| >= 1
############################################################

regional_sig <- regional_crossplatform %>%
  
  filter(
    FDR < 0.05,
    abs(combined_logFC) >= 1
  ) %>%
  
  arrange(
    brain_region,
    FDR
  )

############################################################
# 16. UP/DOWN REGULATED REGIONAL GENES
############################################################

regional_up <- regional_sig %>%
  
  filter(
    combined_logFC > 0
  )

regional_down <- regional_sig %>%
  
  filter(
    combined_logFC < 0
  )

############################################################
# 17. EXPORT REGIONAL SIGNIFICANT RESULTS
############################################################

write.csv(
  regional_sig,
  "results/brain_region/BrainRegion_CrossPlatform_Significant.csv",
  row.names = FALSE
)

write.csv(
  regional_up,
  "results/brain_region/BrainRegion_CrossPlatform_UP.csv",
  row.names = FALSE
)

write.csv(
  regional_down,
  "results/brain_region/BrainRegion_CrossPlatform_DOWN.csv",
  row.names = FALSE
)

############################################################
# 18. REGIONAL GENE COUNTS
############################################################

regional_counts <- regional_sig %>%
  
  group_by(
    brain_region
  ) %>%
  
  summarise(
    significant_genes = n(),
    .groups = "drop"
  ) %>%
  
  arrange(
    desc(significant_genes)
  )

write.csv(
  regional_counts,
  "results/brain_region/BrainRegion_SignificantGeneCounts.csv",
  row.names = FALSE
)

############################################################
# 19. SUMMARY
############################################################

cat("\n")
cat("============================================================\n")
cat("BRAIN-REGION CROSS-PLATFORM META-ANALYSIS COMPLETED\n")
cat("============================================================\n")

cat(
  "Microarray region-gene meta-analysis results: ",
  nrow(micro_region_meta),
  "\n",
  sep = ""
)

cat(
  "RNA-seq region-gene meta-analysis results: ",
  nrow(rna_region_meta),
  "\n",
  sep = ""
)

cat(
  "Shared region-gene combinations: ",
  nrow(regional_common),
  "\n",
  sep = ""
)

cat(
  "Concordant region-gene combinations: ",
  nrow(regional_concordant),
  "\n",
  sep = ""
)

cat(
  "Significant regional genes (FDR < 0.05, |log2FC| >= 1): ",
  nrow(regional_sig),
  "\n",
  sep = ""
)

cat(
  "Upregulated regional genes: ",
  nrow(regional_up),
  "\n",
  sep = ""
)

cat(
  "Downregulated regional genes: ",
  nrow(regional_down),
  "\n",
  sep = ""
)

cat("\nRegional gene counts:\n")
print(regional_counts)

cat("\nOutput files:\n")

cat(
  "  results/brain_region/BrainRegion_CrossPlatform_AllGenes.csv\n"
)

cat(
  "  results/brain_region/BrainRegion_CrossPlatform_Significant.csv\n"
)

cat(
  "  results/brain_region/BrainRegion_CrossPlatform_UP.csv\n"
)

cat(
  "  results/brain_region/BrainRegion_CrossPlatform_DOWN.csv\n"
)

cat(
  "  results/brain_region/BrainRegion_SignificantGeneCounts.csv\n"
)

cat("============================================================\n")
