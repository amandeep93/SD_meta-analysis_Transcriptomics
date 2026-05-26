############################################################
# RNA-SEQ DIFFERENTIAL EXPRESSION ANALYSIS
############################################################

options(stringsAsFactors = FALSE)
set.seed(123)

library(DESeq2)

counts <- read.csv("data/counts.csv", row.names = 1)
coldata <- read.csv("data/metadata.csv")

coldata$condition <- factor(coldata$condition)

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = coldata,
  design = ~ condition
)

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

dds <- DESeq(dds)

res <- results(dds)

res_df <- as.data.frame(res)
res_df$Gene <- rownames(res_df)

write.csv(
  res_df,
  "results/RNAseq_DEG.csv",
  row.names = FALSE
)
