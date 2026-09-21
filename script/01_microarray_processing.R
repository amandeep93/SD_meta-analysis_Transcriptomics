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
# MICROARRAY DATASETS
############################################################
#
# IMPORTANT:
# This script currently contains only the two datasets below.
# If the manuscript states that 11 microarray datasets were
# analyzed, the complete dataset list must be restored here
# before the final analysis.
#
############################################################

datasets <- c(
  "GSE6514",
  "GSE33302"
)

############################################################
# PART 1: INDIVIDUAL MICROARRAY DATASET PROCESSING
############################################################

for(ds in datasets){

  cat("====================================================\n")
  cat("Processing Microarray Dataset via Limma:", ds, "\n")
  cat("====================================================\n")

  ##########################################################
  # FILE PATHS
  ##########################################################

  cel_path <- paste0(
    "data/CEL_files/",
    ds
  )

  metadata_path <- paste0(
    "data/metadata/",
    ds,
    "_metadata.csv"
  )

  ##########################################################
  # CHECK REQUIRED FILES
  ##########################################################

  if(!dir.exists(cel_path)) {

    cat(
      "Warning: CEL file path not found for ",
      ds,
      ". Skipping cohort.\n",
      sep = ""
    )

    next
  }

  if(!file.exists(metadata_path)) {

    cat(
      "Warning: Metadata file not found for ",
      ds,
      ". Skipping cohort.\n",
      sep = ""
    )

    next
  }

  ##########################################################
  # 1. READ RAW CEL FILES
  ##########################################################

  cat(
    "Reading raw CEL files...\n"
  )

  raw_data <- ReadAffy(
    celfile.path = cel_path
  )

  ##########################################################
  # 2. RMA NORMALIZATION
  ##########################################################

  cat(
    "Performing RMA normalization...\n"
  )

  norm_data <- rma(
    raw_data
  )

  expr_matrix <- exprs(
    norm_data
  )

  ##########################################################
  # 3. READ SAMPLE METADATA
  ##########################################################

  metadata <- read.csv(
    metadata_path,
    stringsAsFactors = FALSE
  )

  ##########################################################
  # 4. VERIFY SAMPLE IDENTIFIERS
  ##########################################################

  if(!"Sample" %in% colnames(metadata)) {

    stop(
      "Critical Error: Metadata for ",
      ds,
      " does not contain a 'Sample' column."
    )

  }

  if(!"Condition" %in% colnames(metadata)) {

    stop(
      "Critical Error: Metadata for ",
      ds,
      " does not contain a 'Condition' column."
    )

  }

  ##########################################################
  # CHECK THAT ALL EXPRESSION SAMPLES HAVE METADATA
  ##########################################################

  missing_metadata <- setdiff(
    colnames(expr_matrix),
    metadata$Sample
  )

  if(length(missing_metadata) > 0) {

    stop(
      "Critical Error: The following expression samples ",
      "are missing from metadata for ",
      ds,
      ": ",
      paste(missing_metadata, collapse = ", ")
    )

  }

  ##########################################################
  # REORDER METADATA TO EXACTLY MATCH EXPRESSION MATRIX
  ##########################################################

  metadata <- metadata[
    match(
      colnames(expr_matrix),
      metadata$Sample
    ),
    ,
    drop = FALSE
  ]

  ##########################################################
  # FINAL SAMPLE ORDER CHECK
  ##########################################################

  if(!all(
    colnames(expr_matrix) == metadata$Sample
  )) {

    stop(
      "Critical Error: Expression matrix and metadata ",
      "sample order could not be synchronized for ",
      ds
    )

  }

  ##########################################################
  # 5. DEFINE EXPERIMENTAL CONDITION
  ##########################################################

  metadata$Condition <- factor(
    metadata$Condition,
    levels = c(
      "Control",
      "SD"
    )
  )

  ##########################################################
  # CHECK THAT BOTH CONDITIONS ARE PRESENT
  ##########################################################

  if(any(
    is.na(metadata$Condition)
  )) {

    stop(
      "Critical Error: Metadata contains conditions other ",
      "than 'Control' or 'SD' for ",
      ds
    )

  }

  if(length(
    unique(metadata$Condition)
  ) < 2) {

    stop(
      "Critical Error: Both Control and SD groups are not ",
      "present in dataset ",
      ds
    )

  }

  ##########################################################
  # 6. DESIGN MATRIX
  ##########################################################

  design <- model.matrix(
    ~0 + Condition,
    data = metadata
  )

  colnames(design) <- c(
    "Control",
    "SD"
  )

  ##########################################################
  # 7. LIMMA MODEL FITTING
  ##########################################################

  cat(
    "Fitting limma models...\n"
  )

  fit <- lmFit(
    expr_matrix,
    design
  )

  ##########################################################
  # 8. SD VS CONTROL CONTRAST
  ##########################################################

  contrast.matrix <- makeContrasts(
    SD_vs_Control = SD - Control,
    levels = design
  )

  fit2 <- contrasts.fit(
    fit,
    contrast.matrix
  )

  ##########################################################
  # 9. EMPIRICAL BAYES MODERATION
  ##########################################################

  fit2 <- eBayes(
    fit2
  )

  ##########################################################
  # 10. EXTRACT DIFFERENTIAL EXPRESSION RESULTS
  ##########################################################

  deg <- topTable(
    fit2,
    coef = "SD_vs_Control",
    number = Inf,
    adjust.method = "BH"
  )

  ##########################################################
  # 11. ADD PROBE IDENTIFIERS
  ##########################################################

  deg$Probe_ID <- rownames(
    deg
  )

  ##########################################################
  # 12. CALCULATE TREATMENT-INDEPENDENT MEAN EXPRESSION
  ##########################################################
  #
  # This is used later for probe selection.
  #
  # IMPORTANT:
  # Probe selection is therefore based on overall expression,
  # NOT on the magnitude of the SD effect.
  #
  ##########################################################

  deg$mean_expression <- rowMeans(
    expr_matrix[
      deg$Probe_ID,
      ,
      drop = FALSE
    ],
    na.rm = TRUE
  )

  ##########################################################
  # 13. EXTRACT LIMMA PARAMETERS
  ##########################################################

  deg$df_total <- fit2$df.total[
    match(
      rownames(deg),
      rownames(fit2)
    )
  ]

  deg$stdev_unscaled <- fit2$stdev.unscaled[
    match(
      rownames(deg),
      rownames(fit2)
    ),
    "SD_vs_Control"
  ]

  ##########################################################
  # 14. CREATE STANDARD ERROR
  ##########################################################
  #
  # The moderated t-statistic is:
  #
  # t = logFC / SE
  #
  # Therefore:
  #
  # SE = |logFC| / |t|
  #
  ##########################################################

  deg$SE <- abs(
    deg$logFC
  ) / abs(
    deg$t
  )

  ##########################################################
  # 15. PREPARE EXPORT TABLE
  ##########################################################

  deg_export <- deg %>%
    
    mutate(
      Gene = toupper(
        trimws(
          Probe_ID
        )
      ),
      
      logFC = as.numeric(
        logFC
      ),
      
      t_stat = as.numeric(
        t
      ),
      
      P.Value = as.numeric(
        P.Value
      ),
      
      adj.P.Val = as.numeric(
        adj.P.Val
      ),
      
      mean_expression = as.numeric(
        mean_expression
      ),
      
      SE = as.numeric(
        SE
      )
    ) %>%
    
    filter(
      !is.na(Gene),
      Gene != "",
      !is.na(logFC),
      is.finite(logFC),
      !is.na(SE),
      is.finite(SE),
      SE > 0,
      !is.na(mean_expression),
      is.finite(mean_expression)
    ) %>%
    
    select(
      Gene,
      Probe_ID,
      logFC,
      t_stat,
      SE,
      df_total,
      stdev_unscaled,
      mean_expression,
      P.Value,
      adj.P.Val
    )

  ##########################################################
  # 16. SAVE INDIVIDUAL DATASET RESULTS
  ##########################################################

  output_file <- paste0(
    "results/",
    ds,
    "_limma_DEG.csv"
  )

  write.csv(
    deg_export,
    output_file,
    row.names = FALSE
  )

  cat(
    "Dataset completed: ",
    output_file,
    "\n",
    sep = ""
  )

}

############################################################
# PART 2: CONSOLIDATE DATASETS FOR META-ANALYSIS
############################################################

cat(
  "\n====================================================\n"
)

cat(
  "Consolidating microarray results...\n"
)

cat(
  "====================================================\n"
)

############################################################
# FIND GENERATED DATASET FILES
############################################################

deg_files <- list.files(
  "results",
  pattern = "_limma_DEG.csv$",
  full.names = TRUE
)

if(length(deg_files) == 0) {

  stop(
    "Critical Error: No microarray DEG files found in results/."
  )

}

############################################################
# COMBINE DATASETS
############################################################

micro_data <- map_df(
  deg_files,
  function(f) {

    df <- read.csv(
      f,
      stringsAsFactors = FALSE
    )

    dataset_name <- gsub(
      "_limma_DEG.csv",
      "",
      basename(f)
    )

    df$dataset <- dataset_name

    df

  }
)

############################################################
# FINAL DATA CLEANING
############################################################

micro_data <- micro_data %>%
  
  filter(
    !is.na(Gene),
    Gene != "",
    !is.na(logFC),
    !is.na(SE),
    is.finite(SE),
    SE > 0,
    !is.na(mean_expression),
    is.finite(mean_expression)
  )

############################################################
# SAVE META-ANALYSIS INPUT
############################################################

write.csv(
  micro_data,
  "results/Microarray_processed_data.csv",
  row.names = FALSE
)

############################################################
# SAVE SESSION INFORMATION
############################################################

writeLines(
  capture.output(
    sessionInfo()
  ),
  "sessionInfo_microarray.txt"
)

############################################################
# FINAL SUMMARY
############################################################

cat(
  "====================================================\n"
)

cat(
  "Microarray preprocessing completed.\n"
)

cat(
  "Datasets processed: ",
  length(unique(micro_data$dataset)),
  "\n",
  sep = ""
)

cat(
  "Total dataset-gene records: ",
  nrow(micro_data),
  "\n",
  sep = ""
)

cat(
  "Output: results/Microarray_processed_data.csv\n"
)

cat(
  "====================================================\n"
)
