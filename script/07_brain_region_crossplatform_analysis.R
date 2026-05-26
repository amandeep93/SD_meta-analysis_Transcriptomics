############################################################
# BRAIN REGION-SPECIFIC CROSS-PLATFORM ANALYSIS
#
# Description:
# Tissue-level integration of microarray and
# RNA-seq datasets using metadata-guided
# regional mapping and cross-platform filtering.
############################################################

options(stringsAsFactors = FALSE)
set.seed(123)

############################################################
# LIBRARIES
############################################################

library(readxl)
library(dplyr)
library(tidyr)
library(pheatmap)
library(ggplot2)
library(ggrepel)
library(igraph)
library(ggraph)

############################################################
# OUTPUT DIRECTORIES
############################################################

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

############################################################
# READ METADATA
############################################################

meta <- read.csv("data/metadata.csv")

############################################################
# READ MICROARRAY DATA
############################################################

micro_file <- "data/microarray_only_matrix.xlsx"

micro_sheets <- excel_sheets(micro_file)

micro_list <- lapply(micro_sheets, function(s){
  
  df <- read_excel(micro_file, sheet = s)
  
  colnames(df)[1:3] <- c(
    "Gene",
    "logFC",
    "pval"
  )
  
  df$dataset <- s
  
  df$platform <- "Microarray"
  
  df
})

micro_data <- bind_rows(micro_list)

############################################################
# READ RNA-SEQ DATA
############################################################

rna_file <- "data/RNA_only_matrix.xlsx"

rna_sheets <- excel_sheets(rna_file)

rna_list <- lapply(rna_sheets, function(s){
  
  df <- read_excel(rna_file, sheet = s)
  
  colnames(df)[1:3] <- c(
    "Gene",
    "logFC",
    "pval"
  )
  
  df$dataset <- s
  
  df$platform <- "RNAseq"
  
  df
})

rna_data <- bind_rows(rna_list)

############################################################
# COMBINE DATASETS
############################################################

all_data <- bind_rows(
  micro_data,
  rna_data
)

############################################################
# MAP BRAIN REGIONS USING METADATA
############################################################

all_data <- all_data %>%
  left_join(
    meta,
    by = c("dataset" = "Sample_ID")
  )

############################################################
# FIX DUPLICATED PLATFORM COLUMNS
############################################################

if("platform.x" %in% colnames(all_data)){
  
  all_data <- all_data %>%
    mutate(platform = platform.x) %>%
    select(-platform.x, -platform.y)
}

############################################################
# REMOVE WHOLE BRAIN / PONS
############################################################

all_data <- all_data %>%
  filter(
    !brain_region %in% c(
      "Whole_Brain",
      "Pons"
    )
  )

############################################################
# TISSUE LIST
############################################################

tissues <- unique(all_data$brain_region)

tissues <- tissues[!is.na(tissues)]

############################################################
# STORAGE
############################################################

tissue_gene_lists <- list()

############################################################
# LOOP THROUGH TISSUES
############################################################

for(t in tissues){
  
  cat("Processing:", t, "\n")
  
  tissue_data <- all_data %>%
    filter(brain_region == t)
  
  ##########################################################
  # GENE SUPPORT ACROSS PLATFORMS
  ##########################################################
  
  support <- tissue_data %>%
    
    distinct(
      Gene,
      dataset,
      platform
    ) %>%
    
    group_by(
      Gene,
      platform
    ) %>%
    
    summarise(
      n = n(),
      .groups = "drop"
    )
  
  support_wide <- support %>%
    
    pivot_wider(
      names_from = platform,
      values_from = n,
      values_fill = 0
    )
  
  if(!"Microarray" %in% colnames(support_wide)){
    support_wide$Microarray <- 0
  }
  
  if(!"RNAseq" %in% colnames(support_wide)){
    support_wide$RNAseq <- 0
  }
  
  ##########################################################
  # KEEP GENES PRESENT IN BOTH PLATFORMS
  ##########################################################
  
  genes_keep <- support_wide %>%
    
    filter(
      Microarray >= 1 &
        RNAseq >= 1
    ) %>%
    
    pull(Gene)
  
  tissue_data <- tissue_data %>%
    filter(Gene %in% genes_keep)
  
  ##########################################################
  # PLATFORM EFFECTS
  ##########################################################
  
  platform_effect <- tissue_data %>%
    
    group_by(
      Gene,
      platform
    ) %>%
    
    summarise(
      mean_logFC = mean(logFC),
      .groups = "drop"
    )
  
  platform_effect <- platform_effect %>%
    
    pivot_wider(
      names_from = platform,
      values_from = mean_logFC
    )
  
  ##########################################################
  # SAVE GENE LISTS
  ##########################################################
  
  write.csv(
    platform_effect,
    paste0(
      "results/",
      t,
      "_CrossPlatformGenes.csv"
    ),
    row.names = FALSE
  )
  
  tissue_gene_lists[[t]] <- platform_effect$Gene
  
  ##########################################################
  # CROSS-PLATFORM CONCORDANCE
  ##########################################################
  
  scatter <- ggplot(
    platform_effect,
    aes(Microarray, RNAseq)
  ) +
    
    geom_point(
      color = "black",
      size = 2
    ) +
    
    geom_smooth(
      method = "lm",
      color = "red"
    ) +
    
    theme_classic(base_size = 14) +
    
    ggtitle(
      paste(
        "Cross-platform concordance:",
        t
      )
    )
  
  ggsave(
    paste0(
      "figures/",
      t,
      "_ConcordanceScatter.png"
    ),
    scatter,
    width = 6,
    height = 6,
    dpi = 300
  )
}

############################################################
# GLOBAL BRAIN REGION MATRIX
############################################################

marker_data <- all_data %>%
  
  group_by(
    Gene,
    brain_region
  ) %>%
  
  summarise(
    mean_logFC = mean(logFC),
    .groups = "drop"
  )

marker_matrix <- marker_data %>%
  
  pivot_wider(
    names_from = brain_region,
    values_from = mean_logFC
  )

############################################################
# MATRIX CONVERSION
############################################################

mat <- as.matrix(marker_matrix[,-1])

rownames(mat) <- marker_matrix$Gene

############################################################
# HANDLE MISSING VALUES
############################################################

mat <- mat[
  rowSums(is.na(mat)) <= 2,
]

############################################################
# REMOVE EXTREME VALUES
############################################################

mat[mat > 4]  <- 4
mat[mat < -4] <- -4

############################################################
# SELECT TOP GENES
############################################################

gene_strength <- apply(
  abs(mat),
  1,
  max
)

top_genes <- names(
  sort(
    gene_strength,
    decreasing = TRUE
  )
)[1:120]

top_genes <- intersect(
  top_genes,
  rownames(mat)
)

top_mat <- mat[
  top_genes,
  ,
  drop = FALSE
]

############################################################
# GLOBAL HEATMAP
############################################################

pheatmap(
  
  top_mat,
  
  color = colorRampPalette(
    c("green","black","red")
  )(100),
  
  clustering_distance_rows = "euclidean",
  
  clustering_distance_cols = "euclidean",
  
  clustering_method = "ward.D2",
  
  fontsize_row = 4,
  
  fontsize_col = 12,
  
  border_color = NA,
  
  filename = "figures/Global_Brain_Heatmap.png",
  
  height = 14,
  
  width = 7
)

############################################################
# TISSUE GENE COUNTS
############################################################

tissue_counts <- sapply(
  tissue_gene_lists,
  length
)

df_bar <- data.frame(
  Tissue = names(tissue_counts),
  Genes = tissue_counts
)

barplot_fig <- ggplot(
  df_bar,
  aes(Tissue, Genes, fill = Tissue)
) +
  
  geom_bar(stat = "identity") +
  
  theme_classic(base_size = 14)

ggsave(
  "figures/Tissue_Gene_Counts.png",
  barplot_fig,
  width = 7,
  height = 5,
  dpi = 300
)

############################################################
# GENE OVERLAP NETWORK
############################################################

edges <- data.frame()

for(i in 1:(length(tissues)-1)){
  
  for(j in (i+1):length(tissues)){
    
    g1 <- tissue_gene_lists[[tissues[i]]]
    
    g2 <- tissue_gene_lists[[tissues[j]]]
    
    overlap <- length(
      intersect(g1, g2)
    )
    
    edges <- rbind(
      edges,
      data.frame(
        from = tissues[i],
        to = tissues[j],
        weight = overlap
      )
    )
  }
}

edges <- edges %>%
  filter(weight > 0)

g <- graph_from_data_frame(edges)

png(
  "figures/Tissue_Overlap_Network.png",
  width = 1800,
  height = 1400,
  res = 300
)

ggraph(g, layout = "fr") +
  
  geom_edge_link(
    aes(width = weight),
    alpha = 0.8
  ) +
  
  geom_node_point(
    size = 8,
    color = "red"
  ) +
  
  geom_node_text(
    aes(label = name),
    repel = TRUE
  ) +
  
  theme_void()

dev.off()

############################################################
# SAVE SESSION INFO
############################################################

writeLines(
  capture.output(sessionInfo()),
  "sessionInfo.txt"
)

############################################################
# FINAL MESSAGE
############################################################

cat(
  "Brain region-specific cross-platform analysis completed successfully\n"
)
