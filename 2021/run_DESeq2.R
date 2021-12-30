library(DESeq2)
library(ggplot2)

ctx <- as.matrix(read.csv("RNAseq_counts_matrix.csv",row.names=1))
#ctx <- ctx[,which(!colnames(ctx) %in% c("GSM4403289","GSM4636681"))]
#l <- unlist(lapply(colnames(ctx),function(x){x %in% c("GSM4403300","GSM4403307","GSM4403291","GSM4403296")}))
#names(l) <- colnames(ctx)
#coldata <- data.frame(condition=factor(l))

coldata <- read.csv("RNAseq_tpm_matrix.patientCorrCluster.csv",row.names=1,header=T)
coldata <- coldata[coldata$glia14 %in% c(1,3),]
#coldata <- coldata[coldata$tissue=="LPS",]
coldata$condition <- factor(coldata$glia14)
print(coldata)
ctx <- ctx[, rownames(coldata)]
print(head(ctx))

dds <- DESeqDataSetFromMatrix(countData = ctx,
                              colData = coldata,
                              design = ~ condition)

dds <- dds[ rowSums(counts(dds)) > 100, ]
dds$sample <- rownames(coldata)
dds$donor <- coldata[,"donorID"]
dds <- collapseReplicates(dds,dds$donor,dds$sample)

dds <- DESeq(dds)
res <- results(dds)
res <- res[!is.na(res$padj),]
write.table(res[order(res$padj,decreasing=FALSE),], file="DESeq2.cluster.tsv", quote=FALSE, sep="\t")
