library(STdeconvolve)
library(SpatialExperiment)

#For single sample
#se <- read10xVisium(samples="JY97A/outs/",type="sparse",data="filtered")
#cd <- se@assays@data@listData$counts
#pos <- spatialCoords(se)
#colnames(pos) <- c("y", "x")
#counts <- cleanCounts(cd, min.lib.size = 100, min.reads = 10)
#corpus <- restrictCorpus(counts, removeAbove=1.0, removeBelow = 0.05, nTopOD = 1000)

paths <- c("JY97A","JY98A","JY101A","JY102A")
SEs <- lapply(paths,function(x){
  se <- read10xVisium(samples=paste0(x,"/outs"),type="sparse",data="filtered")
  return(se)
  })
names(SEs) <- paths

Pos <- lapply(SEs, function(x){
  pos<-spatialCoords(x)
  colnames(pos)<-c("y","x")
  return(pos)
  })

Counts <- lapply(seq_along(SEs), function(x){
  cd <- SEs[[x]]@assays@data@listData$counts
  colnames(cd) <- gsub("-1",paste0("-",paths[x]),colnames(cd))
  return(cd)
  })
names(Counts) <- names(SEs)
comb <- do.call(cbind, Counts)

pdf("STdeconvolve.graphs.pdf")
counts <- cleanCounts(comb, min.lib.size = 2000, min.detected = 2, plot=TRUE)
corpus <- restrictCorpus(counts, removeAbove=1.0, 
                         removeBelow = 0.05, nTopOD = 2000, plot=TRUE)
ldas <- fitLDA(t(as.matrix(corpus)), Ks = seq(5,20), ncores=8, plot=TRUE)
dev.off()
saveRDS(ldas,"STdeconvolve.rds")

optLDA <- optimalModel(models = ldas, opt = 12)
results <- getBetaTheta(optLDA, perc.filt = 0.02, betaScale = 1000)
deconProp <- results$theta
deconGexp <- results$beta
write.table(deconProp,"STdeconvolve.K12.pixelProp.csv",quote=FALSE,sep=",")
write.table(deconGexp,"STdeconvolve.K12.geneExpr.csv",quote=FALSE,sep=",")
