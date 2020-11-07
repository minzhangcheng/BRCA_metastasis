rm(list = ls())
library("MetaDE")
library("MetaQC")

setwd("/home/minzhang/workspace/meta-Analysis/BRCA_metastasis/metaqc")
memory.limit(16000)

study.names <- c("GSE102484_GPL570", "GSE20685_GPL570", "GSE58984_GPL570", "GSE29431_GPL570", "GSE61304_GPL570")
raw <- MetaDE.Read(study.names, skip=rep(1, 5), via="txt", matched=T, log=T)
gc()
merged <- MetaDE.merge(raw)
gc()
Data.QC<-list()
for(i in 1:5){
  colnames(merged[[i]][[1]])<-merged[[i]][[2]]
  Data.QC[[i]]<-impute.knn(merged[[i]][[1]])$data
}
names(Data.QC)<-names(merged)
QC <- MetaQC(Data.QC, "/home/minzhang/workspace/meta-Analysis/c2.all.v6.2.symbols.gmt", filterGenes=F,verbose=TRUE, isParallel=T, resp.type="Twoclass")
gc()
runQC(QC, B=1e4, fileForCQCp="/home/minzhang/workspace/meta-Analysis/c2.all.v6.2.symbols.gmt")
png(file="/home/minzhang/workspace/meta-Analysis/BRCA_metastasis/metaqc/metaQC.png", width=2048, height=2048, bg="transparent")
plot(QC)
dev.off()
sink("/home/minzhang/workspace/meta-Analysis/BRCA_metastasis/metaqc/metaQC.txt")
print(QC)
sink()
save.image("/home/minzhang/workspace/meta-Analysis/BRCA_metastasis/metaqc/metaQC.RData")