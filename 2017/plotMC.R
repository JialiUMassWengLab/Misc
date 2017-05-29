require(rtracklayer)
library(hexbin)
args <- commandArgs(trailingOnly=TRUE)
gtf <- readGFF(args[2],version=2L,tags=c("gene_id","transcript_id"))
ids <- unique(gtf[,c("gene_id","transcript_id")])

counts1 <- read.table(paste0(args[1],".htseq.out"),sep="\t",header=FALSE)
counts2 <- read.table(paste0(args[1],".dedup.htseq.out"),sep="\t",header=FALSE)
covered <- read.table("tmp.bed",sep="\t",header=FALSE)
covered <- covered[covered$V5>0,]
ids <- ids[which(ids$transcript_id %in% covered$V1),]
covered["gene_id"] <- apply(covered,1,function(x){return(ids[which(ids$transcript_id==x["V1"]),"gene_id"])})
covFracs <- aggregate(covered$V8,by=list(covered$gene_id),FUN=max)
mat <- merge(counts1,counts2,by="V1")
mat <- mat[mat$V2.x>0 & mat$V2.y>0,]
mat2 <- merge(mat,covFracs,by.x="V1",by.y="Group.1")
names(mat2) <- c("gene_id","read_counts","uniq_counts","covered_fraction")
mat2["multiplicity"] <- mat2["read_counts"]/mat2["uniq_counts"]
write.table(mat2,file=paste0(args[1],".MC.tsv"),sep="\t",row.names=FALSE,quote=FALSE)

pdf(paste0(args[1],".MCplot.pdf"))
zones <- matrix(c(2,0,1,3), ncol=2, byrow=TRUE)
layout(zones, widths=c(4/5,1/5), heights=c(1/5,4/5))
xhist <- hist(mat2$covered_fraction, plot=FALSE, breaks=50)
yhist <- hist(mat2[mat2$multiplicity<100,]$multiplicity, plot=FALSE, breaks=50)
par(mar=c(5,5,1,1))
plot(mat2$covered_fraction,mat2$multiplicity,ylab="Multiplicity",xlab="Fraction of genes covered",pch=19,ylim=c(0,100),xlim=c(0,1))
par(mar=c(0,5,1,1))
barplot(xhist$counts,axes=FALSE,space=0)
par(mar=c(5,0,1,1))
barplot(yhist$counts,axes=FALSE,space=0,horiz=TRUE)
dev.off()

library(RColorBrewer)
mat2 <- mat2[mat2$multiplicity < 100,]
rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
pdf(paste0(args[1],".MC.hexbin.pdf"))
bin <- hexbin(mat2$covered_fraction,mat2$multiplicity,xbins=100,ylab="Multiplicity",xlab="Fraction of genes covered")
#plot(bin, colramp= function(n){LinOCS(n,beg=15,end=225)})
plot(bin, colramp=rf)
dev.off()
