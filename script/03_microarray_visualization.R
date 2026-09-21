############################################################
# MICROARRAY META-ANALYSIS PLOTS
# FOREST PLOTS & HEATMAP
#
# Script: script/04_microarray_visualization.R
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
library(pheatmap)

############################################################
# LOAD DATA
############################################################

meta_all_path <- "results/Microarray_Meta_AllGenes.csv"
sig_genes_path <- "results/Microarray_FDR_FC1.csv"
processed_data_path <- "results/Microarray_processed_data.csv"

if (!file.exists(meta_all_path) ||
    !file.exists(sig_genes_path) ||
    !file.exists(processed_data_path)) {
  
  stop(
    "Critical Error: Required microarray results are missing. ",
    "Run Scripts 02 and 03 first."
  )
}

micro_meta_res <- read.csv(
  meta_all_path,
  stringsAsFactors = FALSE
)

sig_micro <- read.csv(
  sig_genes_path,
  stringsAsFactors = FALSE
)

micro_data <- read.csv(
  processed_data_path,
  stringsAsFactors = FALSE
)

############################################################
# OUTPUT DIRECTORIES
############################################################

dir.create(
  "figures/microarray_forests",
  showWarnings = FALSE,
  recursive = TRUE
)

############################################################
# CHECK REQUIRED COLUMNS
############################################################

required_meta <- c(
  "Gene",
  "FDR"
)

required_data <- c(
  "Gene",
  "dataset",
  "logFC",
  "SE",
  "mean_expression"
)

missing_meta <- setdiff(
  required_meta,
  colnames(micro_meta_res)
)

missing_data <- setdiff(
  required_data,
  colnames(micro_data)
)

if (length(missing_meta) > 0) {
  stop(
    "Missing columns in Microarray_Meta_AllGenes.csv: ",
    paste(missing_meta, collapse = ", ")
  )
}

if (length(missing_data) > 0) {
  stop(
    "Missing columns in Microarray_processed_data.csv: ",
    paste(missing_data, collapse = ", ")
  )
}

############################################################
# CLEAN AND SELECT ONE PROBE PER GENE PER DATASET
#
# IMPORTANT:
# Probe selection is based on mean expression rather than
# treatment-associated effect size.
#
# This avoids outcome-dependent probe selection.
############################################################

micro_data_cleaned <- micro_data %>%
  
  mutate(
    Gene = toupper(trimws(Gene)),
    dataset = trimws(dataset)
  ) %>%
  
  filter(
    !is.na(Gene),
    Gene != "",
    !is.na(logFC),
    !is.na(SE),
    !is.na(mean_expression),
    SE > 0
  ) %>%
  
  group_by(dataset, Gene) %>%
  
  arrange(
    desc(mean_expression),
    Probe_ID,
    .by_group = TRUE
  ) %>%
  
  slice(1) %>%
  
  ungroup()

cat(
  "Unique dataset-gene observations after expression-based ",
  "probe selection: ",
  nrow(micro_data_cleaned),
  "\n"
)

############################################################
# PART 1
# FOREST PLOTS FOR TOP 10 SIGNIFICANT GENES
############################################################

cat(
  "Generating forest plots for the top 10 significant ",
  "microarray genes...\n"
)

############################################################
# Select top genes by FDR
############################################################

top_forest_genes <- sig_micro %>%
  
  mutate(
    Gene = toupper(trimws(Gene))
  ) %>%
  
  arrange(FDR) %>%
  
  slice_head(n = 10) %>%
  
  pull(Gene)

############################################################
# Generate forest plot for each gene
############################################################

for (g in top_forest_genes) {
  
  df <- micro_data_cleaned %>%
    filter(Gene == g) %>%
    arrange(dataset)
  
  ##########################################################
  # Require at least two independent datasets
  ##########################################################
  
  if (nrow(df) < 2) {
    cat(
      "Skipping ",
      g,
      ": fewer than two contributing datasets.\n",
      sep = ""
    )
    next
  }
  
  ##########################################################
  # Random-effects meta-analysis
  ##########################################################
  
  model <- tryCatch(
    
    rma(
      yi = logFC,
      sei = SE,
      method = "REML",
      test = "knha",
      data = df
    ),
    
    error = function(e) {
      cat(
        "Meta-analysis failed for ",
        g,
        ": ",
        e$message,
        "\n",
        sep = ""
      )
      return(NULL)
    }
  )
  
  if (is.null(model)) {
    next
  }
  
  ##########################################################
  # Prediction interval
  ##########################################################
  
  pred <- tryCatch(
    predict(model),
    error = function(e) NULL
  )
  
  ##########################################################
  # Output file
  ##########################################################
  
  output_file <- paste0(
    "figures/microarray_forests/Forest_",
    g,
    ".png"
  )
  
  png(
    filename = output_file,
    width = 2000,
    height = 1500,
    res = 300
  )
  
  ##########################################################
  # Forest plot
  ##########################################################
  
  forest(
    model,
    slab = df$dataset,
    xlab = "Log2 Fold Change",
    mlab = "Random-effects model",
    cex = 0.8,
    font = 1
  )
  
  ##########################################################
  # Title
  ##########################################################
  
  pooled <- round(
    as.numeric(model$b),
    2
  )
  
  ci.lb <- round(
    model$ci.lb,
    2
  )
  
  ci.ub <- round(
    model$ci.ub,
    2
  )
  
  if (!is.null(pred)) {
    
    pi.lb <- round(
      pred$pi.lb,
      2
    )
    
    pi.ub <- round(
      pred$pi.ub,
      2
    )
    
    title(
      main = paste0(
        g,
        " | pooled log2FC = ",
        pooled,
        " | 95% CI [",
        ci.lb,
        ", ",
        ci.ub,
        "] | PI [",
        pi.lb,
        ", ",
        pi.ub,
        "]"
      ),
      cex.main = 0.85,
      line = 1
    )
    
  } else {
    
    title(
      main = paste0(
        g,
        " | pooled log2FC = ",
        pooled,
        " | 95% CI [",
        ci.lb,
        ", ",
        ci.ub,
        "]"
      ),
      cex.main = 0.9,
      line = 1
    )
  }
  
  dev.off()
  
  cat(
    "Generated forest plot: ",
    output_file,
    "\n",
    sep = ""
  )
}

cat(
  "Forest plots completed.\n"
)

############################################################
# PART 2
# MICROARRAY HEATMAP
############################################################

cat(
  "Generating microarray meta-analysis heatmap...\n"
)

############################################################
# Select top 100 genes according to FDR
############################################################

top_heatmap_genes <- micro_meta_res %>%
  
  mutate(
    Gene = toupper(trimws(Gene))
  ) %>%
  
  filter(
    !is.na(FDR)
  ) %>%
  
  arrange(FDR) %>%
  
  slice_head(n = 100) %>%
  
  pull(Gene)

############################################################
# Prepare heatmap data
#
# IMPORTANT:
# Missing observations remain NA rather than being converted
# to log2FC = 0.
############################################################

heatmap_data <- micro_data_cleaned %>%
  
  filter(
    Gene %in% top_heatmap_genes
  ) %>%
  
  select(
    Gene,
    dataset,
    logFC
  ) %>%
  
  distinct(
    Gene,
    dataset,
    .keep_all = TRUE
  ) %>%
  
  pivot_wider(
    names_from = dataset,
    values_from = logFC
  )

############################################################
# Convert to matrix
############################################################

mat <- as.matrix(
  heatmap_data[, -1, drop = FALSE]
)

rownames(mat) <- heatmap_data$Gene

############################################################
# Ensure numeric matrix
############################################################

mode(mat) <- "numeric"

############################################################
# Remove genes with insufficient information
#
# At least two non-missing study estimates are required.
############################################################

n_observed <- rowSums(
  !is.na(mat)
)

mat <- mat[
  n_observed >= 2,
  ,
  drop = FALSE
]

############################################################
# Remove genes with zero variance among observed values
############################################################

row_variances <- apply(
  mat,
  1,
  function(x) {
    
    x <- x[!is.na(x)]
    
    if (length(x) < 2) {
      return(NA_real_)
    }
    
    var(x)
  }
)

mat <- mat[
  !is.na(row_variances) &
    row_variances > 0,
  ,
  drop = FALSE
]

############################################################
# Row-wise Z-score normalization
############################################################

mat_scaled <- t(
  apply(
    mat,
    1,
    function(x) {
      
      m <- mean(
        x,
        na.rm = TRUE
      )
      
      s <- sd(
        x,
        na.rm = TRUE
      )
      
      if (is.na(s) || s == 0) {
        return(rep(0, length(x)))
      }
      
      (x - m) / s
    }
  )
)

colnames(mat_scaled) <- colnames(mat)
rownames(mat_scaled) <- rownames(mat)

############################################################
# Heatmap
############################################################

pheatmap(
  
  mat_scaled,
  
  color = colorRampPalette(
    c(
      "blue",
      "white",
      "red"
    )
  )(100),
  
  na_col = "grey90",
  
  clustering_distance_rows = "euclidean",
  
  clustering_distance_cols = "euclidean",
  
  clustering_method = "ward.D2",
  
  fontsize_row = 5,
  
  fontsize_col = 10,
  
  border_color = NA,
  
  filename = "figures/Microarray_Heatmap.png",
  
  width = 6,
  
  height = 10
)

cat(
  "Heatmap successfully generated:\n",
  "figures/Microarray_Heatmap.png\n"
)

############################################################
# COMPLETION MESSAGE
############################################################

cat(
  "\n============================================================\n",
  "SCRIPT 04 COMPLETED SUCCESSFULLY\n",
  "============================================================\n"
)
