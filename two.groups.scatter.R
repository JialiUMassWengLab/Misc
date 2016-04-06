data<-read.table("TARGET-RT.recurrent.expression", header=TRUE, sep="\t")
perturbation<-rnorm(31,0,0.05)
pdf("Anton_recurrent_genes.pdf")
par(mfrow=c(2,2))
for (x in 1:nrow(data)) {
    title<-paste("expr_",data[x,1],".pdf",sep="")
    plot(as.numeric(data[x,33:63])+perturbation,as.numeric(data[x,2:32]),type="n",ylab="RPKM",xaxt="n",xlab="",main=data[x,1],xlim=c(-0.5,1.5))
#    axis(1, labels = FALSE)
    text(cex=0.9, x=c(0.15,1.05), y=par("usr")[3]-0.04*(par("usr")[4]-par("usr")[3]), adj=c(1,1.2), c("No SV","SV"), xpd=TRUE)
    for (y in 2:32) lines(data[x,y+31]+perturbation[y-1],data[x,y],type="p",pch=19,col=ifelse(data[x,y+31]==1,"red","black"))
}
dev.off()
