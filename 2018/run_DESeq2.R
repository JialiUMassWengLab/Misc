library(tximport)
library(DESeq2)
args <- commandArgs(trailingOnly=TRUE)
nonblood <- read.table("~/TCGA/NonBlood.txt",header=FALSE,sep="\t")
nonblood$V1 <- gsub("\\.[0-9]+$","",nonblood$V1)
tx2gene <- read.table("/mnt/shares2/annotations/hg38/gencode_hg38.transcript2gene.mapping.csv",header=FALSE,sep=",")[,c(2,1)]
removeList <- c('11218','11238','11226','11222','11202','11206','11182','11230','11210','11194','17906','11310','11338','11382','11298','11326','11318','11330','11342','11378','11402','11290','11270','11334','11302','11358','11294','11350','11322','11394','11374')
d <- read.table("../sampleInfo.csv",sep=",",header=TRUE)
d <- d[d$Center!="University of Kentucky",]
d <- d[!(d$X %in% removeList),]

f <- as.character()
c <- as.character()
for (i in d[d$Disease=="AD","X"]) {
    a <- Sys.glob(paste0("../RSEM_output/",i,"*.dedup.genes.results"))
    f <- append(a,f)
    c <- append(rep(d[d$X==i,"Center"],length(a)),c)
}
ADf <- f[which(file.exists(f))]
ADc <- c

f <- as.character()
c <- as.character()
for (i in d[d$Disease=="NCI","X"]) {
    a <- Sys.glob(paste0("../RSEM_output/",i,"*.dedup.genes.results"))
    f <- append(a,f)
    c <- append(rep(d[d$X==i,"Center"],length(a)),c)
}
CTRLf <- f[which(file.exists(f))]
CTRLc <- c

#Run DESeq2 on Salmon genes output
fd <- ADf
fn <- CTRLf
print(length(fn))
print(length(fd))
files <- append(fn,fd)
centers <- append(CTRLc,ADc)
print(length(centers))
names(files) <- paste0("sample",1:length(files))
txi0 <- tximport(files,type="rsem")
coldata <- data.frame(condition=factor( c(rep("Normal",length(fn)),rep("DISEASE",length(fd))) ),center=factor(centers))
rownames(coldata) <- colnames(txi0$counts)

txi <- list()
txi$length <- txi0$length[apply(txi0$length,1,function(x){all(x>0)}),]
txi$counts <- txi0$counts[apply(txi0$length,1,function(x){all(x>0)}),]
txi$abundance <- txi0$abundance[apply(txi0$length,1,function(x){all(x>0)}),]
txi$countsFromAbundance <- txi0$countsFromAbundance
ddsTxi <- DESeqDataSetFromTximport(txi, colData=coldata, design=~ center + condition)
#ddsTxi <- DESeqDataSetFromTximport(txi, colData=coldata, design=~ condition)
dim(ddsTxi)
#ddsTxi <- ddsTxi[ rowSums(counts(ddsTxi)) > 25, ]
#Kentucky only
ddsTxi <- ddsTxi[ rowSums(counts(ddsTxi)) > 5, ]
dim(ddsTxi)

ddsTxi$run <- unlist(lapply(files,function(x){vec=strsplit(basename(x),"-")[[1]];paste(vec[1],vec[length(vec)-1],sep="-")}))
ddsTxi$sample <- factor(unlist(lapply(files,function(x){vec=strsplit(basename(x),"-")[[1]];vec[1]})))
ddsColl <- collapseReplicates(ddsTxi,ddsTxi$sample,ddsTxi$run)
#ddsColl <- ddsColl[ ,unlist(lapply(colData(ddsColl)$runsCollapsed,function(x){length(strsplit(x,",")[[1]])>1})) ]

dds <- DESeq(ddsColl)
res <- results(dds)
res <- res[!is.na(res$padj),]
#write.table(res[order(res$padj,decreasing=FALSE),], file="DESeq2.center.tsv", quote=FALSE, sep="\t")
#write.table(res[order(res$padj,decreasing=FALSE),], file="DESeq2.Kentucky.tsv", quote=FALSE, sep="\t")
#write.table(res[order(res$padj,decreasing=FALSE),], file="DESeq2.ADvsNCI.tsv", quote=FALSE, sep="\t")
write.table(res[order(res$padj,decreasing=FALSE),], file="DESeq2.nonKenxUCSD20.tsv", quote=FALSE, sep="\t")
