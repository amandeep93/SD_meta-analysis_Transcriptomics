############################################################
# MULTI-DATASET RNA-SEQ DESEQ2 DIFFERENTIAL EXPRESSION
#
# Script: script/05_rnaseq_differential_expression.R
#
# Workflow:
# featureCounts → QC/filtering → DESeq2 → SD vs Control
############################################################

############################################################
# SETTINGS
############################################################

options(stringsAsFactors = FALSE)
set.seed(123)

############################################################
# LIBRARIES
############################################################

library(DESeq2)
library(dplyr)
library(tibble)

############################################################
# DIRECTORY INFRASTRUCTURE
############################################################

dir.create(
  "results/rna_individual",
  showWarnings = FALSE,
  recursive = TRUE
)

############################################################
# RNA-SEQ DATASETS
############################################################

rna_datasets <- c(
  "GSE132074",
  "GSE144957",
  "GSE222410",
  "GSE113754",
  "GSE211088",
  "GSE211301",
  "GSE156925",
  "GSE166831"
)

############################################################
# PROCESS EACH DATASET
############################################################

for (ds in rna_datasets) {
  
  cat("\n")
  cat("====================================================\n")
  cat("Processing RNA-seq dataset:", ds, "\n")
  cat("====================================================\n")
  
  ##########################################################
  # FILE PATHS
  ##########################################################
  
  counts_path <- paste0(
    "data/RNAseq_counts/",
    ds,
    "_raw_counts.txt"
  )
  
  metadata_path <- paste0(
    "data/metadata/",
    ds,
    "_metadata.csv"
  )
  
  ##########################################################
  # CHECK INPUT FILES
  ##########################################################
  
  if (!file.exists(counts_path)) {
    
    warning(
      "Count matrix missing for ",
      ds,
      ". Skipping dataset."
    )
    
    next
  }
  
  if (!file.exists(metadata_path)) {
    
    warning(
      "Metadata missing for ",
      ds,
      ". Skipping dataset."
    )
    
    next
  }
  
  ##########################################################
  # 1. IMPORT FEATURECOUNTS MATRIX
  ##########################################################
  
  cat("Reading featureCounts matrix...\n")
  
  counts_raw <- read.table(
    counts_path,
    header = TRUE,
    row.names = 1,
    sep = "\t",
    skip = 1,
    check.names = FALSE
  )
  
  ##########################################################
  # featureCounts structure:
  #
  # Geneid
  # Chr
  # Start
  # End
  # Strand
  # Length
  # sample1
  # sample2
  # ...
  #
  # Therefore columns 6 onward are sample counts.
  ##########################################################
  
  if (ncol(counts_raw) < 7) {
    
    warning(
      "Unexpected featureCounts structure for ",
      ds,
      ". Skipping dataset."
    )
    
    next
  }
  
  counts <- counts_raw[
    ,
    6:ncol(counts_raw),
    drop = FALSE
  ]
  
  ##########################################################
  # Convert count matrix to numeric matrix
  ##########################################################
  
  counts <- as.matrix(counts)
  mode(counts) <- "numeric"
  
  ##########################################################
  # Clean sample names
  ##########################################################
  
  colnames(counts) <- gsub(
    "\\.",
    "-",
    colnames(counts)
  )
  
  ##########################################################
  # 2. IMPORT SAMPLE METADATA
  ##########################################################
  
  cat("Reading sample metadata...\n")
  
  coldata <- read.csv(
    metadata_path,
    row.names = 1,
    check.names = FALSE
  )
  
  ##########################################################
  # Clean metadata sample names
  ##########################################################
  
  rownames(coldata) <- gsub(
    "\\.",
    "-",
    rownames(coldata)
  )
  
  ##########################################################
  # CHECK CONDITION VARIABLE
  ##########################################################
  
  if (!"condition" %in% colnames(coldata)) {
    
    warning(
      "No 'condition' column found in metadata for ",
      ds,
      ". Skipping dataset."
    )
    
    next
  }
  
  ##########################################################
  # Clean condition labels
  ##########################################################
  
  coldata$condition <- trimws(
    as.character(coldata$condition)
  )
  
  coldata$condition <- factor(
    coldata$condition,
    levels = c(
      "Control",
      "SD"
    )
  )
  
  ##########################################################
  # Check for invalid condition values
  ##########################################################
  
  if (any(is.na(coldata$condition))) {
    
    warning(
      "Metadata for ",
      ds,
      " contains condition values other than ",
      "'Control' or 'SD'.\n",
      "Check the metadata file before analysis."
    )
    
    print(
      unique(
        read.csv(
          metadata_path,
          stringsAsFactors = FALSE
        )$condition
      )
    )
    
    next
  }
  
  ##########################################################
  # 3. CRITICAL SAMPLE-METADATA MATCHING
  ##########################################################
  
  count_samples <- colnames(counts)
  metadata_samples <- rownames(coldata)
  
  ##########################################################
  # Identify samples missing from either file
  ##########################################################
  
  missing_metadata <- setdiff(
    count_samples,
    metadata_samples
  )
  
  missing_counts <- setdiff(
    metadata_samples,
    count_samples
  )
  
  if (length(missing_metadata) > 0) {
    
    stop(
      "Samples present in count matrix but absent from metadata for ",
      ds,
      ": ",
      paste(missing_metadata, collapse = ", ")
    )
  }
  
  if (length(missing_counts) > 0) {
    
    stop(
      "Samples present in metadata but absent from count matrix for ",
      ds,
      ": ",
      paste(missing_counts, collapse = ", ")
    )
  }
  
  ##########################################################
  # REORDER METADATA TO EXACTLY MATCH COUNT MATRIX
  ##########################################################
  
  coldata <- coldata[
    count_samples,
    ,
    drop = FALSE
  ]
  
  ##########################################################
  # Final identity check
  ##########################################################
  
  if (!identical(
    colnames(counts),
    rownames(coldata)
  )) {
    
    stop(
      "Count matrix and metadata are not correctly aligned for ",
      ds
    )
  }
  
  cat(
    "Sample-metadata matching confirmed for ",
    ds,
    ".\n",
    sep = ""
  )
  
  ##########################################################
  # 4. CHECK EXPERIMENTAL GROUPS
  ##########################################################
  
  condition_table <- table(
    coldata$condition
  )
  
  print(condition_table)
  
  if (!all(
    c("Control", "SD") %in%
      names(condition_table)
  )) {
    
    warning(
      ds,
      " does not contain both Control and SD samples. ",
      "Skipping dataset."
    )
    
    next
  }
  
  ##########################################################
  # 5. CREATE DESEQDATASET
  ##########################################################
  
  dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = coldata,
    design = ~ condition
  )
  
  ##########################################################
  # 6. LOW-COUNT FILTERING
  #
  # Retain genes with at least 10 total reads across
  # all samples in the dataset.
  ##########################################################
  
  keep <- rowSums(
    counts(dds)
  ) >= 10
  
  dds <- dds[
    keep,
    ,
    drop = FALSE
  ]
  
  cat(
    "Genes retained after low-count filtering: ",
    nrow(dds),
    "\n",
    sep = ""
  )
  
  ##########################################################
  # 7. DESEQ2 MODEL FITTING
  ##########################################################
  
  cat(
    "Running DESeq2 size-factor, dispersion and GLM ",
    "estimation...\n"
  )
  
  dds <- DESeq(
    dds
  )
  
  ##########################################################
  # 8. SD VS CONTROL CONTRAST
  ##########################################################
  
  res <- results(
    dds,
    contrast = c(
      "condition",
      "SD",
      "Control"
    )
  )
  
  ##########################################################
  # 9. EXTRACT RESULTS
  ##########################################################
  
  res_df <- as.data.frame(
    res
  ) %>%
    
    rownames_to_column(
      var = "Gene"
    ) %>%
    
    mutate(
      Gene = toupper(
        trimws(Gene)
      ),
      dataset = ds
    ) %>%
    
    filter(
      !is.na(Gene),
      Gene != "",
      !is.na(log2FoldChange),
      !is.na(lfcSE)
    ) %>%
    
    select(
      Gene,
      log2FoldChange,
      lfcSE,
      pvalue,
      padj,
      dataset
    )
  
  ##########################################################
  # 10. SAVE INDIVIDUAL DATASET RESULTS
  ##########################################################
  
  output_file <- paste0(
    "results/rna_individual/",
    ds,
    "_DESeq2_DEG.csv"
  )
  
  write.csv(
    res_df,
    output_file,
    row.names = FALSE
  )
  
  cat(
    "Dataset completed successfully.\n",
    "Results saved to: ",
    output_file,
    "\n",
    sep = ""
  )
}

############################################################
# COMPLETION
############################################################

cat("\n")
cat("====================================================\n")
cat(
  "All RNA-seq DESeq2 differential-expression models ",
  "have been processed.\n"
)
cat("====================================================\n")
