data6<-readLines("UCEC.homology.len")
d6<-strsplit(data6,"\t")
for (x in c(1:length(d6))) d6[[x]]<-type.convert(d6[[x]])
med6<-mat.or.vec(1,length(d6))
for (x in c(1:length(d6))) med6[x]<-median(d6[[x]])
d6<-d6[order(med6,decreasing=TRUE)]

data1<-readLines("GBM.homology.len")
d1<-strsplit(data1,"\t")
for (x in c(1:length(d1))) d1[[x]]<-type.convert(d1[[x]])
med1<-mat.or.vec(1,length(d1))
for (x in c(1:length(d1))) med1[x]<-median(d1[[x]])
d1<-d1[order(med1,decreasing=TRUE)]

data2<-readLines("CESC.homology.len")
d2<-strsplit(data2,"\t")
for (x in c(1:length(d2))) d2[[x]]<-type.convert(d2[[x]])
med2<-mat.or.vec(1,length(d2))
for (x in c(1:length(d2))) med2[x]<-median(d2[[x]])
d2<-d2[order(med2,decreasing=TRUE)]

data3<-readLines("BRCA.homology.len")
d3<-strsplit(data3,"\t")
for (x in c(1:length(d3))) d3[[x]]<-type.convert(d3[[x]])
med3<-mat.or.vec(1,length(d3))
for (x in c(1:length(d3))) med3[x]<-median(d3[[x]])
d3<-d3[order(med3,decreasing=TRUE)]

data4<-readLines("LAML.homology.len")
d4<-strsplit(data4,"\t")
for (x in c(1:length(d4))) d4[[x]]<-type.convert(d4[[x]])
med4<-mat.or.vec(1,length(d4))
for (x in c(1:length(d4))) med4[x]<-median(d4[[x]])
d4<-d4[order(med4,decreasing=TRUE)]

data5<-readLines("SARC.homology.len")
d5<-strsplit(data5,"\t")
for (x in c(1:length(d5))) d5[[x]]<-type.convert(d5[[x]])
med5<-mat.or.vec(1,length(d5))
for (x in c(1:length(d5))) med5[x]<-median(d5[[x]])
d5<-d5[order(med5,decreasing=TRUE)]

all<-append(d6,d1)
all<-append(all,d5)
all<-append(all,d2)
all<-append(all,d3)
all<-append(all,d4)
pdf("Homology_length_boxplot.pdf", 18, 7)
boxplot(all, xlab=NULL, ylab="Homology length (bp)", col=c(rep("grey",15),rep("blue",13),rep("orange",16),rep("purple",12),rep("green",20),rep("red",21)), cex.lab=1.5, cex.axis=1.5)
dev.off()
