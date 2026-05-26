############################################################

library(dplyr)
library(purrr)
library(metafor)

micro_data <- read.csv("results/Microarray_processed_data.csv")

micro_meta_res <- micro_data %>%
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

micro_meta_res$FDR <- p.adjust(
  micro_meta_res$meta_pval,
  method = "BH"
)

write.csv(
  micro_meta_res,
  "results/Microarray_Meta_AllGenes.csv",
  row.names = FALSE
)

sig_micro <- micro_meta_res %>%
  filter(FDR < 0.05 & abs(meta_logFC) >= 1)

write.csv(
  sig_micro,
  "results/Microarray_FDR_FC1.csv",
  row.names = FALSE
)