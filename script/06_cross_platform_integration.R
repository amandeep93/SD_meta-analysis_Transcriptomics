############################################################
# CROSS-PLATFORM META-ANALYSIS
############################################################

options(stringsAsFactors = FALSE)
set.seed(123)

library(dplyr)

micro <- read.csv("results/Microarray_Meta_AllGenes.csv")
rna <- read.csv("results/RNAseq_Meta_AllGenes.csv")

common <- inner_join(micro, rna, by = "Gene")

common <- common %>%
  mutate(
    combined_logFC = (meta_logFC.x + meta_logFC.y)/2
  )

common <- common %>%
  mutate(
    fisher_stat = -2*(
      log(meta_pval.x) + log(meta_pval.y)
    ),
    fisher_p = pchisq(
      fisher_stat,
      df = 4,
      lower.tail = FALSE
    )
  )

common$FDR <- p.adjust(
  common$fisher_p,
  method = "BH"
)

final_genes <- common %>%
  filter(FDR < 0.05 & abs(combined_logFC) >= 1)

write.csv(
  final_genes,
  "results/CrossPlatform_901_genes.csv",
  row.names = FALSE
)
