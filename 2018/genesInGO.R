library(topGO)
library(limma)
library(org.Hs.eg.db)

args <- commandArgs(trailingOnly=TRUE)
data = read.table("~/NASH/2018AprReval2/1804Revalidation2.congregated.tsv",sep="\t",header=T,row.names=1)
data <- data[!grepl("^ERCC-",row.names(data)),]
row.names(data) <- gsub(".[0-9]*$","",row.names(data))

#fn <- 'center.downGenes.txt'
fn <- 'center.upGenes.txt'
print(fn)
thred <- 0.05
d <- read.table(fn,sep=",",header=T,row.names=1)
geneList <- factor(as.integer(row.names(data) %in% row.names(d)))
names(geneList) <- row.names(data)
GOdata <- new("topGOdata",description="list",ontology="BP",allGenes=geneList,annot=annFUN.org,nodeSize=5,ID="ensembl",mapping = "org.Hs.eg.db")
genes <- genesInTerm(GOdata,args[1])[[1]]
genes <- genes[genes %in% row.names(d)]
print(genes)
write(genes,paste0(args[1],".genes.txt"),ncolumns=1)
