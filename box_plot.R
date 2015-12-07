genes<-c("BACH1","BARD1","BLM","BRCA1","BRCA2","BRIP1","DNA2","EXO1","FAM175A","FAM175B","MRE11A","NBN","RAD50","RAD51","RBBP8","TOPBP1","UIMC1","WRN")
#genes<-c("TP53BP1","MDC1","RIF1","PAXIP1","H2AFX","XRCC6","XRCC5","LIG4","XRCC4","NHEJ1","PRKDC","DCLRE1C","LIG1","LIG3","APTX","APLF","PNKP","POLM","POLL","TDP1","TDP2")

pvalues<-numeric()
for (g in c(1:length(genes))) {
  name<-paste(genes[g],".tmp",sep="")
  data<-readLines(name)
  d<-strsplit(data,"\t")
  for (x in c(1:length(d))) d[[x]]<-type.convert(d[[x]])
  pvalues<-append(pvalues,wilcox.test(d[[1]],d[[5]])$p.value)
  pvalues<-append(pvalues,wilcox.test(d[[2]],d[[5]])$p.value)
  pvalues<-append(pvalues,wilcox.test(d[[3]],d[[5]])$p.value)
  pvalues<-append(pvalues,wilcox.test(d[[4]],d[[5]])$p.value)
}
qvalues<-p.adjust(pvalues, method="fdr")
print(pvalues)
print(qvalues)

pdf("HR_expr_boxplot.pdf")
#pdf("NHEJ_expr_boxplot.pdf")
par(mfrow=c(2,3),mar=c(4,4,4,3)+0.1)
i<-1
for (g in c(1:length(genes))) {
  name<-paste(genes[g],".tmp",sep="")
  data<-readLines(name)
  d<-strsplit(data,"\t")
  for (x in c(1:length(d))) d[[x]]<-type.convert(d[[x]])
  pvalue1<-format(qvalues[i], scientific=TRUE, digits=3)
  pvalue2<-format(qvalues[i+1], scientific=TRUE, digits=3)
  pvalue3<-format(qvalues[i+2], scientific=TRUE, digits=3)
  pvalue4<-format(qvalues[i+3], scientific=TRUE, digits=3)
  i<-i+4
  
  d1<-append(d[[1]],d[[2]])
  d1<-append(d1,d[[3]])
  d1<-append(d1,d[[4]])
  d1<-append(d1,d[[5]])
  d2<-c(rep("UCEC",length(d[[1]])),rep("SARC",length(d[[2]])),rep("CESC",length(d[[3]])),rep("BRCA",length(d[[4]])),rep("LAML",length(d[[5]])))
  df<-data.frame(d1,d2)
  test<-aov(d1 ~ d2, data=df)
  pvalue0<-format(summary(test)[[1]][["Pr(>F)"]][[1]], scientific=TRUE, digits=3)
#  print(pvalue0)

  mt<-paste(genes[g],pvalue1,pvalue2,pvalue3,pvalue4, sep="\n")
#  mt<-paste(genes[g],pvalue0, sep="\n")
#  boxplot(d, main=mt, names=c("BRCA","CESC","UCEC","SARC","LAML"), boxwex=0.5, cex.main=0.7, col="orange", outpch=19, outcol="grey", ylab="TPM (Transcripts per Million)")
#  boxplot(d, main=mt, names=c("BRCA","CESC","UCEC","SARC","LAML"), boxwex=0.5, col="orange", outpch=19, outcol="grey", ylab="TPM (Transcripts per Million)")
  sample<-c("UCEC","SARC","CESC","BRCA","AML")
  boxplot(d, main=mt, xaxt="n", boxwex=0.5, cex.main=0.7, col=c("grey","orange","purple","green","red"), outpch=19, outcol="grey", ylab="Expression Level (TPM)")
  axis(1, labels = FALSE)
  text(cex=0.9, x=seq_along(sample), y=par("usr")[3]-0.04*(par("usr")[4]-par("usr")[3]), adj=c(1,1), sample, xpd=TRUE, srt=45)
}
dev.off()
