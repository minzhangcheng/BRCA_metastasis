rm(list = ls())
library("affyPLM")
library("affy")

memory.limit(16000)

setwd("/home/minzhang/workspace/meta-Analysis/BRCA_metastasis/cel/GSE5327_GPL96")
raw.set <- ReadAffy()
rma.data <- rma(raw.set)
GSE5327 <- exprs(rma.data)
write.table(GSE5327, "/home/minzhang/workspace/meta-Analysis/BRCA_metastasis/expr/GSE5327_GPL96_expr.tsv", sep="\t")
rm(raw.set, rma.data)
gc()


