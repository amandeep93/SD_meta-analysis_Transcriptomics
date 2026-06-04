############################################################
# MICROARRAY META-ANALYSIS PLOTS (FOREST & HEATMAPS)
# Script: script/04_microarray_visualization.R
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
library(tidyr)
library(purrr)
library(metafor)
library(pheatmap)

############################################################
# LOAD DATA ASSETS
############################################################
meta_all_path <- "results/Microarray_Meta_AllGenes.csv"
sig_genes_path <- "results/Microarray_FDR_FC1.csv"
processed_data_path <- "results/Microarray_processed_data.csv"

if(!file.exists(meta_all_path) | !file.exists(sig_genes_path)) {
  stop("Critical Error: Standalone meta-analysis files missing. Run script 03 first.")
}

micro_meta_res <- read.csv(meta_all_path)
sig_micro      <- read.csv(sig_genes_path)
micro_data     <- read.csv(processed_data_path)

# Ensure output figures directory structure is available
dir.create("figures/microarray_forests", showWarnings = FALSE, recursive = TRUE)

############################################################
# PART 1: HIGHEST SIGNIFICANT GENE FOREST PLOT AUTOMATION
############################################################
cat("Generating individual study forest plots for top significant genes...\n")

# Target top 10 most highly altered significant discovery genes for explicit forest plotting
top_forest_genes <- sig_micro %>%
  arrange(FDR) %>%
  slice(1:10) %>%
  pull(Gene)

# Enforce identical probe independence constraints as the meta-analysis step
micro_data_cleaned <- micro_data %>%
  mutate(Gene = toupper(trimws(Gene))) %>%
  group_by(dataset, Gene) %>%
  filter(abs(logFC) == max(abs(logFC))) %>%
  slice(1) %>%
  ungroup()

for (g in top_forest_genes) {
  df <- micro_data_cleaned %>% filter(Gene == g)
  
  if (nrow(df) < 2) next
  
  # Recalculate model parameters to pass to the forest engine
  model <- tryCatch(
    rma(yi = df$logFC, sei = df$SE, method = "REML", test = "knha"),
    error = function(e) NULL
  )
  
  if (is.null(model)) next
  
  png(
    filename = paste0("figures/microarray_forests/Forest_", g, ".png"),
    width = 2000,
    height = 1500,
    res = 300
  )
  
  # Execute base forest engine mapping configurations
  forest(
    model,
    slab = df$dataset,
    xlab = "Log2 Fold Change",
    mlab = "Random-effects model",
    cex = 0.8,
    font = 1
  )
  
  pooled <- round(as.numeric(model$b), 2)
  ci.lb  <- round(model$ci.lb, 2)
  ci.ub  <- round(model$ci.ub, 2)
  
  title(
    main = paste0(g, " (meta log2FC = ", pooled, ", 95% CI [", ci.lb, "; ", ci.ub, "])"),
    cex.main = 0.9,
    line = 1
  )
  
  dev.off()
}
cat("Forest plots successfully compiled in figures/microarray_forests/\n")

############################################################
# PART 2: DUAL-PLATFORM COMPREHENSIVE CLUSTERING HEATMAP
############################################################
cat("Compiling global microarray transcriptomic heatmap...\n")

# Extract top 100 discovery features ranked by False Discovery Rate control
top_heatmap_genes <- micro_meta_res %>%
  arrange(FDR) %>%
  slice(1:100) %>%
  pull(Gene)

# Pivot data into a structured layout matching Gene vs. Dataset
heatmap_data <- micro_data_cleaned %>%
  filter(Gene %in% top_heatmap_genes) %>%
  select(Gene, dataset, logFC) %>%
  pivot_wider(names_from = dataset, values_from = logFC, values_fill = 0)

mat <- as.matrix(heatmap_data[, -1])
rownames(mat) <- heatmap_data$Gene

# Safeguard feature validation: Remove features with zero variation across studies to prevent Z-score scaling crashes
row_variances <- apply(mat, 1, var)
mat <- mat[row_variances > 0, , drop = FALSE]

# Perform row-wise normalization scaling (Z-scores across cohorts)
mat_scaled <- t(scale(t(mat)))

# Export high-resolution transcriptomic plot matrix
pheatmap(
  mat_scaled,
  color = colorRampPalette(c("blue", "white", "red"))(100), # Standard peer-review color scheme
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

cat("Success! Generated comprehensive heatmap plot asset: figures/Microarray_Heatmap.png\n")
