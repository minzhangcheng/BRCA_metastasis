rm(list = ls())
library("affyPLM")
library("affy")

memory.limit(16000)

setwd("/home/minzhang/workspace/meta-Analysis/BRCA_metastasis/cel/GSE20685_GPL570")
raw.set <- ReadAffy()
rma.data <- rma(raw.set)
GSE20685 <- exprs(rma.data)
write.table(GSE20685, "/home/minzhang/workspace/meta-Analysis/BRCA_metastasis/expr/GSE20685_GPL570_expr.tsv", sep="\t")
rm(raw.set, rma.data)
gc()


