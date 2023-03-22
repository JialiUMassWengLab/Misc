library(edgeR)

counts <- read.table("gene.results.merged.count.txt",sep="\t",header=TRUE,row.names=1)
#head(counts)
d0 <- DGEList(counts)
d0 <- calcNormFactors(d0)
cutoff <- 1
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,]

snames <- colnames(counts)
trt <- gsub("[0-9]+$","",snames)
trt[which(trt=="il")] <- "IL12"

mm <- model.matrix(~0 + trt)
png("M-A.plot.png")
y <- voom(d, mm, plot = T)
dev.off()

fit <- lmFit(y, mm)
#head(coef(fit))
contr <- makeContrasts(trtLPAM - trtIgG, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 30)
write.table(top.table,file="DEA.limma_voom.LPAM_vs_IgG.csv",sep=",",quote=FALSE)

contr <- makeContrasts(trtIL12 - trtIgG, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 30)
write.table(top.table,file="DEA.limma_voom.IL12_vs_IgG.csv",sep=",",quote=FALSE)
