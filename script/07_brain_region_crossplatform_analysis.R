############################################################
# RIGOROUS REGIONAL FISHER-TO-FDR CROSS-PLATFORM INTEGRATION
# Script: script/08_brain_region_analysis.R
############################################################

options(stringsAsFactors = FALSE)
set.seed(123)

library(dplyr)
library(tidyr)
library(purrr)
library(pheatmap)
library(ggplot2)
library(igraph)
library(ggraph)

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

meta <- read.csv("data/metadata.csv")

cat("Compiling cross-platform expression tables from generated assets...\n")
micro_files <- list.files("results", pattern = "GSE.*_limma_DEG.csv", full.names = TRUE)
rna_files   <- list.files("results/rna_individual", pattern = "GSE.*_DESeq2_DEG.csv", full.names = TRUE)

if(length(micro_files) == 0 | length(rna_files) == 0) {
  stop("Critical Error: Individual dataset log files missing. Run scripts 02 and 05 first.")
}

micro_data <- map_df(micro_files, function(f) {
  df <- read.csv(f)
  df$dataset <- gsub("_limma_DEG.csv", "", basename(f))
  df$platform <- "Microarray"
  df %>% select(Gene, logFC, P.Value, dataset, platform) %>% rename(pval = P.Value)
})

rna_data <- map_df(rna_files, function(f) {
  df <- read.csv(f)
  df$dataset := gsub("_DESeq2_DEG.csv", "", basename(f))
  df$platform <- "RNAseq"
  df %>% select(Gene, log2FoldChange, pvalue, dataset, platform) %>% rename(logFC = log2FoldChange, pval = pvalue)
})

all_data <- bind_rows(micro_data, rna_data) %>%
  mutate(Gene = toupper(trimws(Gene))) %>%
  filter(!is.na(Gene) & Gene != "")

all_data <- all_data %>%
  left_join(meta %>% select(Sample_ID, brain_region), by = c("dataset" = "Sample_ID")) %>%
  filter(!is.na(brain_region) & !brain_region %in% c("Whole_Brain", "Pons"))

tissues <- unique(all_data$brain_region)
tissue_gene_lists <- list()

############################################################
# REGIONAL FISHER META-ANALYSIS LOOP
############################################################
for(t in tissues) {
  cat("Running Tissue Fisher Meta-Analysis:", t, "\n")
  
  tissue_data <- all_data %>% filter(brain_region == t)
  
  # Group by Gene and Platform to get clean platform-specific averages/p-values first
  platform_summaries <- tissue_data %>%
    group_by(Gene, platform) %>%
    summarise(
      mean_logFC = mean(logFC, na.rm = TRUE),
      min_pval = min(pval, na.rm = TRUE), # Conservatively track tracking p-value
      .groups = "drop"
    )
  
  # Pivot to check cross-platform presence
  wide_summary <- platform_summaries %>%
    pivot_wider(
      names_from = platform, 
      values_from = c(mean_logFC, min_pval),
      values_fill = list(mean_logFC = NA, min_pval = NA)
    )
  
  if(!"mean_logFC_Microarray" %in% colnames(wide_summary) | !"mean_logFC_RNAseq" %in% colnames(wide_summary)) next
  
  # Filter for cross-platform presence and directional harmony
  regional_candidates <- wide_summary %>%
    filter(!is.na(mean_logFC_Microarray) & !is.na(mean_logFC_RNAseq)) %>%
    filter(sign(mean_logFC_Microarray) == sign(mean_logFC_RNAseq))
  
  if(nrow(regional_candidates) == 0) next
  
  # EXECUTE THE FISHER-TO-FDR TRANSITION FOR THIS SPECIFIC TISSUE
  regional_analysis <- regional_candidates %>%
    mutate(
      p_micro_safe = ifelse(min_pval_Microarray == 0, 1e-300, min_pval_Microarray),
      p_rna_safe   = ifelse(min_pval_RNAseq == 0, 1e-300, min_pval_RNAseq),
      
      # Fisher Formula
      fisher_stat = -2 * (log(p_micro_safe) + log(p_rna_safe)),
      fisher_p = pchisq(fisher_stat, df = 4, lower.tail = FALSE),
      
      combined_logFC = (mean_logFC_Microarray + mean_logFC_RNAseq) / 2
    )
  
  # Apply Benjamini-Hochberg FDR correction locally within this tissue
  regional_analysis$FDR <- p.adjust(regional_analysis$fisher_p, method = "BH")
  
  # Save the statistically verified regional matrix
  write.csv(
    regional_analysis,
    paste0("results/", t, "_CrossPlatformGenes.csv"),
    row.names = FALSE
  )
  
  # Isolate significant tissue markers (FDR < 0.05) for network integration
  sig_tissue_genes <- regional_analysis %>%
    filter(FDR < 0.05) %>%
    pull(Gene)
  
  tissue_gene_lists[[t]] <- sig_tissue_genes
}

############################################################
# GLOBAL HEATMAP & PLOTTING CODE (PRESERVED & PROTECTED)
############################################################
cat("Generating plots...\n")
marker_matrix_prep <- all_data %>%
  group_by(Gene, brain_region) %>%
  summarise(mean_logFC = mean(logFC, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = brain_region, values_from = mean_logFC)

mat <- as.matrix(marker_matrix_prep[, -1])
rownames(mat) <- marker_matrix_prep$Gene
mat[is.na(mat)] <- 0
mat <- mat[apply(mat, 1, var) > 0, , drop = FALSE]

mat[mat > 4]  <- 4
mat[mat < -4] <- -4

gene_strength <- apply(abs(mat), 1, max)
top_genes <- names(sort(gene_strength, decreasing = TRUE))[1:100]
top_mat <- mat[top_genes, , drop = FALSE]

pheatmap(
  top_mat,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "ward.D2",
  fontsize_row = 4.5,
  fontsize_col = 10,
  border_color = NA,
  filename = "figures/Global_Brain_Heatmap.png",
  height = 11, width = 6
)

df_bar <- data.frame(
  Tissue = names(sapply(tissue_gene_lists, length)),
  Genes  = as.numeric(sapply(tissue_gene_lists, length))
) %>% arrange(desc(Genes))

barplot_fig <- ggplot(df_bar, aes(x = reorder(Tissue, Genes), y = Genes, fill = Tissue)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  coord_flip() +
  theme_classic(base_size = 11) +
  labs(title = "Conserved Cross-Platform Genes per Brain Structure", x = "", y = "Significantly Consistent Genes (FDR < 0.05)")

ggsave("figures/Tissue_Gene_Counts.png", barplot_fig, width = 6, height = 4, dpi = 300)

edges <- data.frame()
tissue_names <- names(tissue_gene_lists)

if(length(tissue_names) >= 2) {
  for(i in 1:(length(tissue_names)-1)) {
    for(j in (i+1):length(tissue_names)) {
      g1 <- tissue_gene_lists[[tissue_names[i]]]
      g2 <- tissue_gene_lists[[tissue_names[j]]]
      overlap <- length(intersect(g1, g2))
      if(overlap > 0) {
        edges <- rbind(edges, data.frame(from = tissue_names[i], to = tissue_names[j], weight = overlap))
      }
    }
  }
}

if(nrow(edges) > 0) {
  g_net <- graph_from_data_frame(edges, directed = FALSE)
  png("figures/Tissue_Overlap_Network.png", width = 1800, height = 1500, res = 300)
  print(
    ggraph(g_net, layout = "stress") +
      geom_edge_link(aes(edge_width = weight), alpha = 0.5, color = "gray40") +
      geom_node_point(size = 7, color = "darkblue") +
      geom_node_text(aes(label = name), repel = TRUE, fontface = "bold", size = 3) +
      scale_edge_width_continuous(range = c(0.5, 4.5), name = "Overlapping Signif. Genes") +
      theme_void() +
      theme(legend.position = "bottom")
  )
  dev.off()
}
