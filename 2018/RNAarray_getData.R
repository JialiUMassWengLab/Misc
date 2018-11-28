# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
# R scripts generated  Fri Jul 21 17:17:44 EDT 2017

################################################################
#   Differential expression analysis with limma
library(Biobase)
library(GEOquery)
library(limma)

# load series and platform data from GEO

gset <- getGEO("GSE110226")
#gset <- getGEO("GSE110226", GSEMatrix=TRUE)
#if (length(gset) > 1) idx <- grep("GPL10379", attr(gset, "names")) else idx <- 1
#gset <- gset[[idx]]

ex <- exprs(gset[[1]])
#feature <- pData(featureData(gset[[1]]))[,c("RefSeq_ID","Symbol")]
#feature <- pData(featureData(gset[[1]]))[,c("ENTREZ_GENE_ID","Gene Symbol")]
feature <- pData(featureData(gset[[1]]))[,c("EntrezGeneID","GeneSymbol")]
pheno <- pData(phenoData(gset[[1]]))[,c(1,6,8)]
rownames(ex) <- feature[rownames(ex),"GeneSymbol"]
colnames(ex) <- pheno[colnames(ex),"title"]
#write.table(ex, file=stdout(), row.names=T, sep="\t", quote=F)
write.table(ex, "ChoroidPlexus.expr.tsv", row.names=T, sep="\t", quote=F)
