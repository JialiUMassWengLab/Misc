library(oligo)
library(huex10sttranscriptcluster.db)

celFiles <- list.celfiles()
affyRaw <- read.celfiles(celFiles)
eset <- rma(affyRaw)

my_frame <- data.frame(exprs(eset))
Annot <- data.frame(ACCNUM=sapply(contents(huex10sttranscriptclusterACCNUM), paste, collapse=", "), SYMBOL=sapply(contents(huex10sttranscriptclusterSYMBOL), paste, collapse=", "), DESC=sapply(contents(huex10sttranscriptclusterGENENAME), paste, collapse=", "))
all <- merge(Annot, my_frame, by.x=0, by.y=0)
write.table(all,file="data.ann.txt",sep="\t",row.names = F)
