############################################################
# MICROARRAY DATA PROCESSING
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