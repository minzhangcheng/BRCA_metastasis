rm(list = ls())
library("MAMA")
library("metaMA")
library("affyPLM")
library("affy")
library("CONOR")
library("DEDS")

memory.limit(16000)

a <- read.table("/home/minzhang/workspace/meta-Analysis/BRCA_metastasis/expr/GSE102484_GPL570_expr.tsv")
b <- read.table("/home/minzhang/workspace/meta-Analysis/BRCA_metastasis/expr/GSE20685_GPL570_expr.tsv")
m <- xpn(a, b, iterations=10)
a <- m$x
b <- m$y
m <- merge(a, b, by="row.names")
row.names(m) <- m$Row.names
a <- m[, -1]
a <- a[, order(colnames(a))]
a <- as.matrix(a)
rm(b, m)
gc()

b <- read.table("/home/minzhang/workspace/meta-Analysis/BRCA_metastasis/expr/GSE58984_GPL570_expr.tsv")
m <- xpn(a, b, iterations=10)
a <- m$x
b <- m$y
m <- merge(a, b, by="row.names")
row.names(m) <- m$Row.names
a <- m[, -1]
a <- a[, order(colnames(a))]
a <- as.matrix(a)
rm(b, m)
gc()

b <- read.table("/home/minzhang/workspace/meta-Analysis/BRCA_metastasis/expr/GSE29431_GPL570_expr.tsv")
m <- xpn(a, b, iterations=10)
a <- m$x
b <- m$y
m <- merge(a, b, by="row.names")
row.names(m) <- m$Row.names
a <- m[, -1]
a <- a[, order(colnames(a))]
a <- as.matrix(a)
rm(b, m)
gc()

b <- read.table("/home/minzhang/workspace/meta-Analysis/BRCA_metastasis/expr/GSE61304_GPL570_expr.tsv")
m <- xpn(a, b, iterations=10)
a <- m$x
b <- m$y
m <- merge(a, b, by="row.names")
row.names(m) <- m$Row.names
a <- m[, -1]
a <- a[, order(colnames(a))]
a <- as.matrix(a)
rm(b, m)
gc()

merged.expr <- a
save(merged.expr, file="/home/minzhang/workspace/meta-Analysis/BRCA_metastasis/meta/merged.expr.RData")
rm(a)
gc()

annot <- read.csv("/home/minzhang/workspace/meta-Analysis/BRCA_metastasis/meta/clinic.csv", row.names=1, header=T)
merged <- as.data.frame(merged.expr)
annot <- annot[order(row.names(annot)), ]
merged <- merged[, order(colnames(merged))]
merged <- as.matrix(merged)
rm(merged.expr)
gc()
deds <- deds.stat.linkC(merged, annot, B=1300)
save(deds, file="/home/minzhang/workspace/meta-Analysis/BRCA_metastasis/meta/deds.RData")
gc()

r <- topgenes(deds, number=25, sort.by="fc")
write.csv(r, file="/home/minzhang/workspace/meta-Analysis/BRCA_metastasis/meta/mDEDS_25.csv")
r <- topgenes(deds, number=50, sort.by="fc")
write.csv(r, file="/home/minzhang/workspace/meta-Analysis/BRCA_metastasis/meta/mDEDS_50.csv")
r <- topgenes(deds, number=100, sort.by="fc")
write.csv(r, file="/home/minzhang/workspace/meta-Analysis/BRCA_metastasis/meta/mDEDS_100.csv")
r <- topgenes(deds, number=250, sort.by="fc")
write.csv(r, file="/home/minzhang/workspace/meta-Analysis/BRCA_metastasis/meta/mDEDS_250.csv")
r <- topgenes(deds, number=500, sort.by="fc")
write.csv(r, file="/home/minzhang/workspace/meta-Analysis/BRCA_metastasis/meta/mDEDS_500.csv")
r <- topgenes(deds, number=1000, sort.by="fc")
write.csv(r, file="/home/minzhang/workspace/meta-Analysis/BRCA_metastasis/meta/mDEDS_1000.csv")
r <- topgenes(deds, number=10000, sort.by="fc")
write.csv(r, file="/home/minzhang/workspace/meta-Analysis/BRCA_metastasis/meta/mDEDS_10000.csv")
write.table(rownames(merged), file="/home/minzhang/workspace/meta-Analysis/BRCA_metastasis/meta/geneOrder.tsv", sep="	")
