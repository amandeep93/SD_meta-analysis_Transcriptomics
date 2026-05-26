############################################################
# RNA-SEQ META-ANALYSIS
############################################################

options(stringsAsFactors = FALSE)
set.seed(123)

library(dplyr)
library(purrr)
library(metafor)

rna_data <- read.csv("data/RNAseq_all_datasets.csv")

rna_data <- rna_data %>%
  mutate(
    pval = ifelse(P.Value < 1e-300, 1e-300, P.Value),
    Z = qnorm(1 - pval/2),
    SE = abs(logFC)/Z
  )

rna_meta_res <- rna_data %>%
  group_by(Gene) %>%
  group_split() %>%
  map_df(function(df){
    
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
    
    data.frame(
      Gene = unique(df$Gene),
      meta_logFC = as.numeric(model$b),
      CI_lb = model$ci.lb,
      CI_ub = model$ci.ub,
      meta_pval = model$pval,
      I2 = model$I2,
      tau2 = model$tau2,
      n_studies = nrow(df)
    )
  })

rna_meta_res$FDR <- p.adjust(
  rna_meta_res$meta_pval,
  method = "BH"
)

write.csv(
  rna_meta_res,
  "results/RNAseq_Meta_AllGenes.csv",
  row.names = FALSE
)