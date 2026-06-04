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
# LOAD PROCESSED HARD VARIANCE INTERMEDIATE DATA
############################################################
input_path <- "results/Microarray_processed_data.csv"

if(!file.exists(input_path)) {
  stop("Critical Error: Intermediate file 'results/Microarray_processed_data.csv' missing. Please execute script 02 first.")
}

micro_data <- read.csv(input_path)

############################################################
# CLEAN UP & ENFORCE MULTI-PROBE INDEPENDENCE SELECTION
# (Maintains statistical independence as required by peer review)
############################################################
cat("Applying probe independence controls and structural cleaning...\n")

micro_data_cleaned <- micro_data %>%
  mutate(Gene = toupper(trimws(Gene))) %>%
  filter(!is.na(Gene) & Gene != "" & !is.na(logFC) & !is.na(SE)) %>%
  # If a dataset features multiple probe variations for a single gene symbol,
  # retain only the single most responsive probe per dataset to protect independence assumptions
  group_by(dataset, Gene) %>%
  filter(abs(logFC) == max(abs(logFC))) %>%
  slice(1) %>% # Eliminate absolute value matching ties safely
  ungroup()

############################################################
# EXECUTE PARALLELIZED REML META-ANALYSIS LOOP
############################################################
cat("Executing random-effects meta-analysis via Restricted Maximum Likelihood (REML)...\n")

micro_meta_res <- micro_data_cleaned %>%
  group_by(Gene) %>%
  group_split() %>%
  map_df(function(df){
    
    current_gene <- df$Gene[1]
    study_count  <- nrow(df)
    
    # Validation constraint: Must be present in at least 2 independent cohorts
    # to run cross-study synthesis modeling
    if(study_count < 2) return(NULL)
    
    # Run random-effects model with Knapp-Hartung adjustments for robust Type I error control
    model <- tryCatch(
      rma(
        yi = df$logFC,
        sei = df$SE,
        method = "REML",
        test = "knha"
      ),
      error = function(e) NULL
    )
    
    if(is.null(model)) return(NULL)
    
    # Construct unified statistical matrix row safely
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

if(nrow(micro_meta_res) == 0) {
  stop("Critical Error: Meta-analysis yielded zero output rows. Verify source variant structures.")
}

############################################################
# APPLY MULTIPLE TESTING CORRECTIONS (BENJAMINI-HOCHBERG)
############################################################
cat("Applying False Discovery Rate (FDR) control parameters...\n")
micro_meta_res$FDR <- p.adjust(
  micro_meta_res$meta_pval,
  method = "BH"
)

############################################################
# SAVE COMPLETE SYNTHESIZED MATRICES
############################################################
write.csv(
  micro_meta_res,
  "results/Microarray_Meta_AllGenes.csv",
  row.names = FALSE
)
cat("Full portfolio matrix successfully saved to: results/Microarray_Meta_AllGenes.csv\n")

############################################################
# EXTRACT CANDIDATE PORTFOLIOS MATCHING MANUSCRIPT CRITERIA
############################################################
sig_micro <- micro_meta_res %>%
  filter(FDR < 0.05 & abs(meta_logFC) >= 1)

write.csv(
  sig_micro,
  "results/Microarray_FDR_FC1.csv",
  row.names = FALSE
)

cat(paste("Success! Identified", nrow(sig_micro), "statistically significant cross-study microarray candidates.\n"))
