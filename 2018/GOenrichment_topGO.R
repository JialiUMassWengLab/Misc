library(topGO)
library(limma)
library(org.Hs.eg.db)

args <- commandArgs(trailingOnly=TRUE)
data = read.table("~/NASH/2018AprReval2/1804Revalidation2.congregated.tsv",sep="\t",header=T,row.names=1)
data <- data[!grepl("^ERCC-",row.names(data)),]
row.names(data) <- gsub(".[0-9]*$","",row.names(data))

#files <- Sys.glob('center.*.txt')
#files <- Sys.glob('Kentucky.*.txt')
#files <- Sys.glob('split321.training.*.txt')
#files <- Sys.glob('nonKen.*.txt')
files <- Sys.glob('nonKenxUCSD20.*.txt')
for (fn in files)
{
    print(fn)
    thred <- 0.05
    d <- read.table(fn,sep=",",header=T,row.names=1)
    geneList <- factor(as.integer(row.names(data) %in% row.names(d)))
    names(geneList) <- row.names(data)
    GOdata <- new("topGOdata",description="list",ontology="BP",allGenes=geneList,annot=annFUN.org,nodeSize=5,ID="ensembl",mapping = "org.Hs.eg.db")
    resultFisher <- runTest(GOdata,algorithm="classic",statistic="fisher")
    table <- GenTable(GOdata,classicFisher=resultFisher,orderBy="classicFisher",topNodes=length(score(resultFisher)))
    
    #d1 <- d[d$cat<0 | d$cat_fdr > thred,]
    #geneList2 <- as.vector(d1$cat_fdr)
    #names(geneList2) <- row.names(d1)
    #topDiffGenes <- function(allScores) {return(allScores<thred)}
    #GOdata2 <- new("topGOdata",description="p-value",ontology="BP",allGenes=geneList2,geneSel=topDiffGenes,annot=annFUN.org,nodeSize=10,ID="ensembl",mapping = "org.Hs.eg.db")
    #resultFisher <- runTest(GOdata2,algorithm="classic",statistic="fisher")
    #resultKS <- runTest(GOdata2,algorithm="classic",statistic="ks")
    #table <- GenTable(GOdata2,classicFisher=resultFisher,classicKS=resultKS,orderBy="classicFisher",topNodes=length(score(resultFisher)))
    table$adjustP = p.adjust(table$classicFisher,method='fdr')
    write.table(table,gsub('.txt','.GOenrich.txt',fn),sep=",",quote=FALSE)

}
