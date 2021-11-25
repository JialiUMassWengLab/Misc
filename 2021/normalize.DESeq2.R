library(DESeq2)
library(limma)
library(pracma)

counts <- read.csv("../Output/ROSMAP_geneLevel_counts.csv",row.names=1)
#sample 791_130530 doesn't have RIN score
counts <- counts[,-which(colnames(counts)=="X791_130530")]
counts2 <- apply(counts,2,as.integer)
rownames(counts2) <- rownames(counts)

meta <- read.csv("../metadata/ROSMAP_assay_rnaSeq_metadata.csv",row.names=1)
rownames(meta) <- unlist(lapply(rownames(meta),function(x){if (grepl("^[0-9]",x)) {paste0("X",x)} else {x}}))
meta <- meta[which(row.names(meta) %in% colnames(counts2)),]
levels(meta$libraryBatch)[levels(meta$libraryBatch)=="0, 6, 7"] <- "0_6_7"

demo <- read.csv("../Output/demo_clinical.csv",row.names=1)
rownames(demo) <- unlist(lapply(rownames(demo),function(x){if (grepl("^[0-9]",x)) {paste0("X",x)} else {x}}))
meta <- merge(meta,demo[,-3],by=0)
rownames(meta) <- meta$Row.names
print(dim(meta))

coldata <- meta[,c("libraryBatch","RIN","diagnosis")]
counts2 <- counts2[,rownames(coldata)]

dds <- DESeqDataSetFromMatrix(countData = counts2,colData = coldata,design = ~ diagnosis)
vsd <- vst(dds,blind=FALSE)
mat <- limma::removeBatchEffect(assay(vsd),batch=vsd$libraryBatch,covariates=vsd$RIN)

write.csv(2^mat,"ROSMAP_geneLevel_vstBlindFalse_normed_counts.csv",quote=FALSE)
#write.csv(2^assay(vsd),"uncorrected_vst_counts.csv",quote=FALSE)
