rm(list = ls())
library("affyPLM")
library("affy")

memory.limit(16000)

setwd("/home/minzhang/workspace/meta-Analysis/BRCA_metastasis/cel/GSE102484_GPL570")
raw.set <- ReadAffy()
rma.data <- rma(raw.set)
GSE102484 <- exprs(rma.data)
write.table(GSE102484, "/home/minzhang/workspace/meta-Analysis/BRCA_metastasis/expr/GSE102484_GPL570_expr.tsv", sep="\t")
rm(raw.set, rma.data)
gc()


