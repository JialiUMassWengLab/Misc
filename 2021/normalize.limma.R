library(edgeR)
library(pracma)

counts <- read.csv("../Output/ROSMAP_geneLevel_counts.csv",row.names=1)
#sample 791_130530 doesn't have RIN score
counts <- counts[,-which(colnames(counts)=="X791_130530")]
d0 <- DGEList(counts)
d0 <- calcNormFactors(d0)

meta <- read.csv("../metadata/ROSMAP_assay_rnaSeq_metadata.csv",row.names=1)
rownames(meta) <- unlist(lapply(rownames(meta),function(x){if (grepl("^[0-9]",x)) {paste0("X",x)} else {x}}))
meta <- meta[which(row.names(meta) %in% colnames(counts)),]
levels(meta$libraryBatch)[levels(meta$libraryBatch)=="0, 6, 7"] <- "0_6_7"

demo <- read.csv("../Output/demo_clinical.csv",row.names=1)
rownames(demo) <- unlist(lapply(rownames(demo),function(x){if (grepl("^[0-9]",x)) {paste0("X",x)} else {x}}))
meta <- merge(meta,demo[,-3],by=0)
rownames(meta) <- meta$Row.names
print(dim(meta))

mm <- model.matrix(~1 + meta$libraryBatch + meta$RIN + meta$sex + meta$race + meta$apoeGenotype + meta$diagnosis)
counts <- counts[,rownames(meta)]
norm <- voom(counts,mm)
#fit <- lmFit(norm$E,design=norm$design,weights=norm$weights)
#residualExpr <- residuals.MArrayLM(fit,norm$E)
#coef <- fit$coefficients[,-1]
#d <- norm$design[,-1]
#noise <- apply(coef, 1, function(x) dot(x, t(d)))
#residual <- norm$E - t(noise)
residual <- limma::removeBatchEffect(norm$E,batch=meta$libraryBatch,covariates=meta$RIN)
#write.csv(2^residual,"ROSMAP_geneLevel_limma_normed_counts.csv",quote=FALSE)
write.csv(2^norm$E,"uncorrected_limma_counts.csv",quote=FALSE)

