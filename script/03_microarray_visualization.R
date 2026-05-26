############################################################
)

png(
  paste0("figures/Forest_", g, ".png"),
  width = 1800,
  height = 1400,
  res = 300
)

forest(
  model,
  slab = df$dataset,
  xlab = "Log2 Fold Change",
  mlab = "Random-effects model",
  cex = 0.8
)

pooled <- round(model$b,2)
ci.lb <- round(model$ci.lb,2)
ci.ub <- round(model$ci.ub,2)

title(
  paste0(
    g,
    " (meta log2FC = ",
    pooled,
    ", 95% CI [",
    ci.lb,
    "; ",
    ci.ub,
    "])"
  )
)

dev.off()
}

############################################################
# Heatmap
############################################################

top_heatmap_genes <- micro_meta_res %>%
  arrange(FDR) %>%
  slice(1:250) %>%
  pull(Gene)

micro_data_unique <- micro_data %>%
  group_by(Gene, dataset) %>%
  summarise(logFC = mean(logFC), .groups = "drop")

heatmap_data <- micro_data_unique %>%
  filter(Gene %in% top_heatmap_genes) %>%
  pivot_wider(names_from = dataset,
              values_from = logFC)

mat <- as.matrix(heatmap_data[,-1])
rownames(mat) <- heatmap_data$Gene

mat[is.na(mat)] <- 0

mat_scaled <- t(scale(t(mat)))

pheatmap(
  mat_scaled,
  fontsize_row = 6,
  fontsize_col = 8,
  filename = "figures/Microarray_Heatmap.png"
)
