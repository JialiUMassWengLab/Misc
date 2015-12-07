args<-commandArgs(trailingOnly=TRUE)

#Read data
data<-read.table(args[1], header=TRUE, sep="\t")
N <- nrow(data)
M <- ncol(data)-1

#Normalization
m <- as.matrix(data[,-1])
norm <- apply(m, 1, max)
m <- m/norm

#Make individual line plots
pdf("individual_line_plots.pdf")
par(mfrow=c(3,2), cex.main=1.0, mar=c(3,4,6,1.5)+0.1)

names <- c("00dpp","02dpp","04dpp","07dpp","10dpp","12dpp","14dpp","17dpp","20dpp","26dpp","42dpp")
for (x in 1:N) {
    title<-paste(data[x,1],norm[x],sep="\n");
    plot(1:M, m[x,], type="l", lwd=2.5, xlab="Time point", ylab="Expression (normalized to max)", main=title, ylim=c(0,1), xaxt="n")
    axis(1, at=1:M, labels = FALSE)
    text(x=seq_along(names), y=par("usr")[3]-0.08*(par("usr")[4]-par("usr")[3]), adj=c(1,1), names, xpd=TRUE, srt=45)
}

dev.off()

#Plot clusters
k <- as.numeric(args[2])
cl <- kmeans(m, k, iter.max=200000)
pdf("k-means_cluster.pdf")
for (y in 1:k) {
    plot(1:M, cl$centers[y,], main=paste("Cluster",y,sep=" "), type="n", ylim=c(0,1), xlab="Time point", ylab="Relative abundance", xaxt="n")
    axis(1, at=1:M, cex=0.9, label=names, las=2)
    num <- 0
    for (x in 1:N)
        if (cl$cluster[x]==y) {
            lines(1:M, m[x,], col="grey")
            num <- num+1
        }
    lines(1:M, cl$centers[y,], lwd=3, col=colors()[y*10+2], main=paste("Cluster",k," (n=",num,")",sep=""))
}
dev.off()

write(cl$cluster, file="cluster.txt", sep="\n")
