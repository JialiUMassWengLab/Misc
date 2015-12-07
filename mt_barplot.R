args<-commandArgs(trailingOnly=TRUE)
name<-paste(args[1],".vects",sep="")
data<-read.table(name, header=FALSE, sep="\t")
dm<-data[order(data$V3,-data$V2,decreasing=FALSE),]
ds<-data[order(data$V2,-data$V3,decreasing=FALSE),]

d1<-mat.or.vec(16,2)
d1[,2]<-c(rev(tail(ds$V3,n=8)),tail(dm$V3,n=8))
d1[,1]<-c(rev(tail(ds$V2,n=8)),tail(dm$V2,n=8))
print(d1)
d2<-c(as.character(rev(tail(ds$V1,n=8))),as.character(tail(dm$V1,n=8)))
print(d2)

name1<-paste(args[1],".barplot.pdf",sep="")
pdf(name1,9,8)
par(mar=c(5,11,3,2)+0.1, las=1)
#a<-barplot(t(d1), beside=TRUE, width=0.5, horiz=TRUE, col=c("grey","black"), main=args[1], names.arg=d2, ylim=c(0,20), xlab="Percentage of samples", xlim=c(0,0.6), font=3)
a<-barplot(t(d1), beside=TRUE, width=0.4, horiz=TRUE, col=c("orange","purple"), main=args[1], yaxt="n", ylim=c(0,20), xlab="Percentage of samples", xlim=c(0,0.6))
text(cex=1.3, font=4, x=0, y=seq(1,length.out=16,by=1.2)-0.2, d2, xpd=TRUE, pos=2, col=c(rep("orange",8),rep("purple",8)))
dev.off()
