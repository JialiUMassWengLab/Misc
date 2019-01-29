library(tximport)
library(DESeq2)
args <- commandArgs(trailingOnly=TRUE)
tx2gene <- read.table("/mnt/shares2/annotations/hg38/gencode_hg38.transcript2gene.mapping.csv",header=FALSE,sep=",")[,c(2,1)]

f <- as.character()
f <- Sys.glob(paste0("/home/jzhuang@ms.local/Neuro/GatesADcombined/RSEM_output/*.dedup.genes.results"))
fn <- f[which(file.exists(f))]


#Run DESeq2 on Salmon genes output
files <- fn
names(files) <- gsub(".dedup.genes.results","",basename(files))
#names(files) <- paste0("sample",1:length(files))
txi0 <- tximport(files,type="rsem")
coldata <- data.frame(condition=factor( rep("Unknown",length(fn)) ))
rownames(coldata) <- colnames(txi0$counts)

txi <- list()
txi$length <- txi0$length[apply(txi0$length,1,function(x){all(x>0)}),]
txi$counts <- txi0$counts[apply(txi0$length,1,function(x){all(x>0)}),]
txi$abundance <- txi0$abundance[apply(txi0$length,1,function(x){all(x>0)}),]
txi$countsFromAbundance <- txi0$countsFromAbundance
ddsTxi <- DESeqDataSetFromTximport(txi, colData=coldata, design=~ 1)
dim(ddsTxi)
ddsTxi <- ddsTxi[ rowSums(counts(ddsTxi)) > 10, ]
dim(ddsTxi)

dds <- DESeq(ddsTxi)
rld <- rlog(dds,blind=FALSE)
#vst <- vst(dds,blind=FALSE)
head(assay(rld),20)
write.table(assay(rld),file="18GatesADcombined.DESeq2rlog.tsv",sep="\t",quote=FALSE)

