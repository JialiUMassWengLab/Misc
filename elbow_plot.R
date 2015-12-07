data <- read.table("mouse_miRNA_expr.table", header=FALSE, sep="\t")
N <- nrow(data)
matr <- as.matrix(data[,-1])
norm <- apply(matr, 1, max)
matr <- matr/norm

wss <- vector(length=15)
wss[1] <- (N-1)*sum(apply(matr,2,var))
for (k in 2:20) wss[k] <- kmeans(matr, k, iter.max=200000)$tot.withinss
pdf("elbow.pdf")
plot(1:20, wss, type="b", xlab="Number of Clusters",ylab="Within groups sum of squares")
dev.off()
