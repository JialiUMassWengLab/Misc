library(tximport)
library(limma)
library(DESeq2)
args <- commandArgs(trailingOnly=TRUE)
nonblood <- read.table("~/TCGA/NonBlood.txt",header=FALSE,sep="\t")
nonblood$V1 <- gsub("\\.[0-9]+$","",nonblood$V1)
tx2gene <- read.table("/mnt/shares2/annotations/hg38/gencode_hg38.transcript2gene.mapping.csv",header=FALSE,sep=",")[,c(2,1)]
removeList <- c('11218','11238','11226','11222','11202','11206','11182','11230','11210','11194','17906','11310','11338','11382','11298','11326','11318','11330','11342','11378','11402','11290','11270','11334','11302','11358','11294','11350','11322','11394','11374')
d <- read.table("sampleInfo.csv",sep=",",header=TRUE)
d <- d[d$X != '11374',]
#d <- d[!(d$X %in% removeList),]
#d <- d[d$Center=="University of Kentucky",]
#d <- d[d$CDR=="0.50",]

f <- as.character()
c <- as.character()
for (i in d[d$Disease=="AD","X"]) {
    a <- Sys.glob(paste0("RSEM_output/",i,"*.dedup.genes.results"))
    f <- append(a,f)
    c <- append(rep(d[d$X==i,"Center"],length(a)),c)
}
ADf <- f[which(file.exists(f))]
ADc <- c

f <- as.character()
c <- as.character()
for (i in d[d$Disease=="NCI","X"]) {
    a <- Sys.glob(paste0("RSEM_output/",i,"*.dedup.genes.results"))
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
dim(ddsTxi)
ddsTxi <- ddsTxi[ rowSums(counts(ddsTxi)) > 5, ]
#Kentucky only
#ddsTxi <- ddsTxi[ rowSums(counts(ddsTxi)) > 100, ]
dim(ddsTxi)

ddsTxi$run <- unlist(lapply(files,function(x){vec=strsplit(basename(x),"-")[[1]];paste(vec[1],vec[length(vec)],sep="-")}))
ddsTxi$sample <- factor(unlist(lapply(files,function(x){vec=strsplit(basename(x),"-")[[1]];vec[1]})))
ddsColl <- collapseReplicates(ddsTxi,ddsTxi$sample,ddsTxi$run)
#ddsColl <- ddsColl[ ,unlist(lapply(colData(ddsColl)$runsCollapsed,function(x){length(strsplit(x,",")[[1]])>1})) ]

dds <- DESeq(ddsColl)
vsd <- vst(dds,blind=FALSE)
assay(vsd) <- limma::removeBatchEffect(assay(vsd),vsd$center)
head(assay(vsd),20)
write.table(assay(vsd),file="1819GatesADcombined.batchCorrected.DESeq2VST.tsv",sep="\t",quote=FALSE)
