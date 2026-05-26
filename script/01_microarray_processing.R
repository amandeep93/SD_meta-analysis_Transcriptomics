# Description:
# This script performs:
# 1. Reading microarray CEL files
# 2. RMA normalization
# 3. Probe summarization
# 4. Differential expression analysis using limma
# 5. Export of DEG matrices
# 6. Preparation for meta-analysis
############################################################

############################################################
# SETTINGS
############################################################

options(stringsAsFactors = FALSE)
set.seed(123)

############################################################
# LIBRARIES
############################################################

library(affy)
library(limma)
library(readxl)
library(dplyr)
library(purrr)

############################################################
# OUTPUT DIRECTORIES
############################################################

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

############################################################
# PART 1:
# RAW MICROARRAY PREPROCESSING
############################################################

# NOTE:
# This section assumes raw CEL files are available
# inside:
#
# data/CEL_files/
#
# with separate folders for each dataset

############################################################
# EXAMPLE DATASET LIST
############################################################

datasets <- c(
  "GSE6514",
  "GSE33302"
)

############################################################
# LOOP THROUGH DATASETS
############################################################

for(ds in datasets){

  cat("Processing dataset:", ds, "\n")

  ##########################################################
  # READ CEL FILES
  ##########################################################

  cel_path <- paste0(
    "data/CEL_files/",
    ds
  )

  raw_data <- ReadAffy(
    celfile.path = cel_path
  )

  ##########################################################
  # RMA NORMALIZATION
  ##########################################################

  norm_data <- rma(raw_data)

  ##########################################################
  # EXTRACT EXPRESSION MATRIX
  ##########################################################

  expr_matrix <- exprs(norm_data)

  ##########################################################
  # SAMPLE METADATA
  ##########################################################

  # Example:
  # metadata should contain:
  #
  # Sample | Condition
  #
  # SD1    | SD
  # SD2    | SD
  # C1     | Control

  metadata <- read.csv(
    paste0(
      "data/metadata/",
      ds,
      "_metadata.csv"
    )
  )

  ##########################################################
  # LIMMA DESIGN MATRIX
  ##########################################################

  metadata$Condition <- factor(
    metadata$Condition
  )

  design <- model.matrix(
    ~0 + Condition,
    data = metadata
  )

  colnames(design) <- levels(
    metadata$Condition
  )

  ##########################################################
  # LINEAR MODEL FITTING
  ##########################################################

  fit <- lmFit(
    expr_matrix,
    design
  )

  ##########################################################
  # CONTRAST MATRIX
  ##########################################################

  contrast.matrix <- makeContrasts(
    SD_vs_Control = SD - Control,
    levels = design
  )

  ##########################################################
  # FIT CONTRASTS
  ##########################################################

  fit2 <- contrasts.fit(
    fit,
    contrast.matrix
  )

  ##########################################################
  # EMPIRICAL BAYES MODERATION
  ##########################################################

  fit2 <- eBayes(fit2)

  ##########################################################
  # DIFFERENTIAL EXPRESSION RESULTS
  ##########################################################

  deg <- topTable(
    fit2,
    coef = "SD_vs_Control",
    number = Inf,
    adjust.method = "BH"
  )

  ##########################################################
  # ADD GENE COLUMN
  ##########################################################

  deg$Gene <- rownames(deg)

  ##########################################################
  # KEEP REQUIRED COLUMNS
  ##########################################################

  deg_export <- deg %>%
    select(
      Gene,
      logFC,
      P.Value,
      adj.P.Val
    )

  ##########################################################
  # SAVE DEG RESULTS
  ##########################################################

  write.csv(
    deg_export,
    paste0(
      "results/",
      ds,
      "_limma_DEG.csv"
    ),
    row.names = FALSE
  )
}

############################################################
# PART 2:
# COMBINE ALL LIMMA OUTPUTS
############################################################

deg_files <- list.files(
  "results",
  pattern = "_limma_DEG.csv",
  full.names = TRUE
)

micro_data <- map_df(deg_files, function(f){

  df <- read.csv(f)

  dataset_name <- gsub(
    "_limma_DEG.csv",
    "",
    basename(f)
  )

  df$dataset <- dataset_name

  df
})

############################################################
# REMOVE MISSING GENES
############################################################

micro_data <- micro_data %>%
  filter(!is.na(Gene))

############################################################
# CALCULATE STANDARD ERRORS
############################################################

micro_data <- micro_data %>%
  mutate(

    pval = ifelse(
      P.Value < 1e-300,
      1e-300,
      P.Value
    ),

    Z = qnorm(1 - pval/2),

    SE = abs(logFC)/Z
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
# SESSION INFO
############################################################

writeLines(
  capture.output(sessionInfo()),
  "sessionInfo.txt"
)

############################################################
# MICROARRAY DATA PROCESSING
############################################################
############################################################
# NOTE:
# This script assumes differential expression analysis
# using limma has already been performed separately
# for each microarray dataset and that processed
# DEG matrices containing Gene, logFC, and P.Value
# columns are available.
############################################################

options(stringsAsFactors = FALSE)
set.seed(123)

library(readxl)
library(dplyr)
library(purrr)

file <- "data/microarray_only_matrix.xlsx"
sheets <- excel_sheets(file)

micro_data <- map_df(sheets, function(s){
  
  df <- read_excel(file, sheet = s)
  
  df <- df[,1:3]
  
  colnames(df) <- c("Gene","logFC","P.Value")
  
  df$dataset <- s
  
  df
})

micro_data <- micro_data %>%
  filter(!is.na(Gene))

micro_data <- micro_data %>%
  mutate(
    pval = ifelse(P.Value < 1e-300, 1e-300, P.Value),
    Z = qnorm(1 - pval/2),
    SE = abs(logFC)/Z
  )

write.csv(micro_data,
          "results/Microarray_processed_data.csv",
          row.names = FALSE)
