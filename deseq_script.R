library(data.table)
library(tidyverse)
library(DESeq2)

args = commandArgs(trailingOnly=TRUE)
setwd(args[1])
counts <- fread("raw_counts.csv")
counts <- counts %>% remove_rownames() %>% column_to_rownames(var = "V1")
coldata <- fread("sample_table.csv",header = T)
coldata <- coldata %>% remove_rownames() %>% column_to_rownames(var = "V1")

dds <- DESeqDataSetFromMatrix(countData = counts,
                             colData = coldata,
                             design = ~ Type)

#Preserve all rows that have more than 10% expression in the cohort
keep <- rowSums(counts(dds)) >= length(rownames(coldata))%/%10
dds <- dds[keep,]
dds <- DESeq(dds); res <- results(dds)
write.csv(as.data.frame(res),file = "../deseq_output/diff_exp.csv")
rld <- vst(dds)
transformed_counts <- assay(rld)
write.csv(transformed_counts,file = "../deseq_output/vst_counts.csv")
write.csv(counts(dds,normalized=TRUE),file="../deseq_output/normalized_counts.csv")