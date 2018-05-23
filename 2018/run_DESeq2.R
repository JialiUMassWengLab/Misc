library(tximport)
library(DESeq2)
args <- commandArgs(trailingOnly=TRUE)
nonblood <- read.table("~/TCGA/NonBlood.txt",header=FALSE,sep="\t")
nonblood$V1 <- gsub("\\.[0-9]+$","",nonblood$V1)
tx2gene <- read.table("/mnt/shares2/annotations/hg38/gencode_hg38.transcript2gene.mapping.csv",header=FALSE,sep=",")[,c(2,1)]
d <- read.table("../sampleInfo.csv",sep=",",header=TRUE)

f <- as.character()
for (i in d[((d$Disease=="NASH" | d$Disease=="NAFLD") & d$PathFibrosis=="F0"),"X"]) {f <- append(Sys.glob(paste0("../RSEM_output/",i,"*.dedup.genes.results")),f)}
F0f <- f[which(file.exists(f))]

f <- as.character()
for (i in d[d$Disease=="Normal Control","X"]) {f <- append(Sys.glob(paste0("../RSEM_output/",i,"*.dedup.genes.results")),f)}
CTRLf <- f[which(file.exists(f))]

f <- as.character()
for (i in d[d$Disease=="NAFLD" & (d$PathFibrosis=="F0" | d$PathFibrosis=="F1"),"X"]) {f <- append(Sys.glob(paste0("../RSEM_output/",i,"*.dedup.genes.results")),f)}
NAFLD01f <- f[which(file.exists(f))]

f <- as.character()
for (i in d[d$Disease=="NASH" & (d$PathFibrosis=="F0" | d$PathFibrosis=="F1"),"X"]) {f <- append(Sys.glob(paste0("../RSEM_output/",i,"*.dedup.genes.results")),f)}
NASH01f <- f[which(file.exists(f))]

f <- as.character()
for (i in d[d$PathFibrosis=="F3" | d$PathFibrosis=="F4","X"]) {f <- append(Sys.glob(paste0("../RSEM_output/",i,"*.dedup.genes.results")),f)}
F34f <- f[which(file.exists(f))]

f <- as.character()
for (i in d[d$PathFibrosis=="F0" | d$PathFibrosis=="F1","X"]) {f <- append(Sys.glob(paste0("../RSEM_output/",i,"*.dedup.genes.results")),f)}
F01f <- f[which(file.exists(f))]

f <- as.character()
for (i in d[d$PathFibrosis=="F3" | d$PathFibrosis=="F2","X"]) {f <- append(Sys.glob(paste0("../RSEM_output/",i,"*.dedup.genes.results")),f)}
F23f <- f[which(file.exists(f))]


#Run DESeq2 on Salmon genes output
fd <- F34f
fn <- F01f
print(length(fn))
print(length(fd))
files <- append(fn,fd)
names(files) <- paste0("sample",1:length(files))
txi0 <- tximport(files,type="rsem")
coldata <- data.frame(condition=factor( c(rep("Normal",length(fn)),rep("DISEASE",length(fd))) ))
rownames(coldata) <- colnames(txi0$counts)

txi <- list()
txi$length <- txi0$length[apply(txi0$length,1,function(x){all(x>0)}),]
txi$counts <- txi0$counts[apply(txi0$length,1,function(x){all(x>0)}),]
txi$abundance <- txi0$abundance[apply(txi0$length,1,function(x){all(x>0)}),]
txi$countsFromAbundance <- txi0$countsFromAbundance
ddsTxi <- DESeqDataSetFromTximport(txi, colData=coldata, design=~ condition)
dim(ddsTxi)
ddsTxi <- ddsTxi[ rowSums(counts(ddsTxi)) > 50, ]
dim(ddsTxi)

ddsTxi$run <- unlist(lapply(files,function(x){vec=strsplit(basename(x),"-")[[1]];paste(vec[1],vec[2],sep="-")}))
ddsTxi$sample <- factor(unlist(lapply(files,function(x){vec=strsplit(basename(x),"-")[[1]];vec[1]})))
ddsColl <- collapseReplicates(ddsTxi,ddsTxi$sample,ddsTxi$run)
ddsColl <- ddsColl[ ,unlist(lapply(colData(ddsColl)$runsCollapsed,function(x){length(strsplit(x,",")[[1]])>1})) ]

dds <- DESeq(ddsColl)
res <- results(dds)
res <- res[!is.na(res$padj),]
write.table(res[order(res$padj,decreasing=FALSE),], file="DESeq2.F01vsF34.tsv", quote=FALSE, sep="\t")
