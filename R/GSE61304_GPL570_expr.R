rm(list = ls())
library("affyPLM")
library("affy")

memory.limit(16000)

setwd("/home/minzhang/workspace/meta-Analysis/BRCA_metastasis/cel/GSE61304_GPL570")
raw.set <- ReadAffy()
rma.data <- rma(raw.set)
GSE61304 <- exprs(rma.data)
write.table(GSE61304, "/home/minzhang/workspace/meta-Analysis/BRCA_metastasis/expr/GSE61304_GPL570_expr.tsv", sep="\t")
rm(raw.set, rma.data)
gc()


