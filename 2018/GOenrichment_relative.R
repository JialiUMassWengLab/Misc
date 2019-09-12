library(topGO)
library(limma)
library(org.Hs.eg.db)

args <- commandArgs(trailingOnly=TRUE)
d0 = read.table("~/NASH/2018AprReval2/1804Revalidation2.congregated.tsv",sep="\t",header=T,row.names=1)
d0 <- d0[!grepl("^ERCC-",row.names(d0)),]
row.names(d0) <- gsub(".[0-9]*$","",row.names(d0))
data = read.table("H1_TPM_p11_a0.mat.tsv",sep="\t",header=T,row.names=1)
#data = read.table("H_rlog_p6_a0.mat.tsv",sep="\t",header=T,row.names=1)
#data = read.table("H_DESeq_p7_a0.mat.tsv",sep="\t",header=T,row.names=1)
data <- t(apply(data,2,function(x){x/max(x)}))
data <- data[apply(data,1,function(x){length(x[x>0.66])<=1}),]

for (i in (1:dim(data)[2]))
#for (i in (7:7))
{
    thred <- 0.3
    #geneList <- factor(as.integer(row.names(data) %in% row.names(d)))
    #names(geneList) <- row.names(data)
    #GOdata <- new("topGOdata",description="list",ontology="BP",allGenes=geneList,annot=annFUN.org,nodeSize=5,ID="ensembl",mapping = "org.Hs.eg.db")
    #resultFisher <- runTest(GOdata,algorithm="classic",statistic="fisher")
    #table <- GenTable(GOdata,classicFisher=resultFisher,orderBy="classicFisher",topNodes=length(score(resultFisher)))
    
    #geneList2 <- as.integer(row.names(d0) %in% row.names(data))
    print(dim(data[data[,i]==1,]))
    geneList2 <- as.numeric()
    for (gene in row.names(d0)) {if (gene %in% row.names(data)) {geneList2<-append(geneList2,data[gene,i])} else {geneList2<-append(geneList2,0.0)}}
    names(geneList2) <- row.names(d0)
    topDiffGenes <- function(allScores) {return(allScores==1)}
    GOdata2 <- new("topGOdata",description="p-value",ontology="BP",allGenes=geneList2,geneSel=topDiffGenes,annot=annFUN.org,nodeSize=10,ID="ensembl",mapping = "org.Hs.eg.db")
    resultFisher <- runTest(GOdata2,algorithm="classic",statistic="fisher")
    #resultKS <- runTest(GOdata2,algorithm="classic",statistic="ks")
    #table <- GenTable(GOdata2,classicFisher=resultFisher,classicKS=resultKS,orderBy="classicFisher",topNodes=length(score(resultFisher)))
    table <- GenTable(GOdata2,classicFisher=resultFisher,orderBy="classicFisher",topNodes=length(score(resultFisher)))
    table$adjustP1 = p.adjust(table$classicFisher,method='fdr')
    #table$adjustP2 = p.adjust(table$classicKS,method='fdr')
    write.table(table,paste0('H1_TPM_p11_a0.comp',i-1,'.GOenrich.txt'),sep=",",quote=FALSE)
    #write.table(table,paste0('H_rlog_p6_a0.comp',i-1,'.GOenrich.txt'),sep=",",quote=FALSE)
    #write.table(table,paste0('H_DESeq_p7_a0.comp',i-1,'.GOenrich.txt'),sep=",",quote=FALSE)
}
