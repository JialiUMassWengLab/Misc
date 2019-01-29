library(RUVSeq)
data <- read.table("1812ADcohort.congregated.RC.tsv",sep="\t",row.names=1,header=TRUE)
data <- data[,1:(dim(data)[2]-2)]
ercc <- data[grep("^ERCC-",rownames(data)),]
ercc1 <- apply(ercc, 1, function(x) all(x>=2))
#spikes <- rownames(data)[grep("^ERCC-", rownames(data))]
set1 <- RUVg(as.matrix(data), names(ercc1[ercc1]), k=1)
write.table(set1$normalizedCounts,"1812ADcohort.RUVnorm.tsv",sep="\t",quote=FALSE)
