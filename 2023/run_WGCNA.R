library(WGCNA)

plotPower <- function(d1,prefix) {
  powers = c(c(1:10), seq(from = 12, to=30, by=2))
  sft = pickSoftThreshold(d1, powerVector = powers, verbose = 5)
  cex1 = 0.9
  pdf(file = paste0(prefix,"_scaleFreePlot.pdf"))
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",main = paste("Scale independence"))
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red")
  dev.off()
}

run <- function(d1,sftPower,MEDissThres,prefix) {
  adjacency = adjacency(d1, power = sftPower)
  TOM = TOMsimilarity(adjacency)
  dissTOM = 1-TOM
  geneTree = hclust(as.dist(dissTOM), method = "average")
  minModuleSize = 10
  #dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2,
  dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 4,
                              pamRespectsDendro = FALSE,
                              minClusterSize = minModuleSize)
  dynamicColors = labels2colors(dynamicMods)
  print(table(dynamicColors))
  
  #sizeGrWindow(8,6)
  #plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05,main = "Gene dendrogram and module colors")
  #dev.off()
  
  MEList = moduleEigengenes(d1, colors = dynamicMods)
  MEs = MEList$eigengenes
  colnames(MEs) = unlist(lapply(colnames(MEs),function(x){gsub("^ME",paste0(prefix,"_module"),x)}))
  
  MEDiss = 1-cor(MEs,use="pairwise.complete.obs")
  METree = hclust(as.dist(MEDiss), method = "average")
  pdf(file = paste0(prefix,"_cut.pdf"), wi = 7, he = 6)
  plot(METree, main = "Clustering of module eigengenes",xlab = "", sub = "")
  abline(h=MEDissThres, col = "red")
  dev.off()
  
  merge = mergeCloseModules(d1, dynamicMods, cutHeight = MEDissThres, iterate = FALSE)
  mergedMEs = merge$newMEs
  mergedMods = merge$colors
  mergedColors = labels2colors(mergedMods)
  print(table(mergedMods))
  #mergedMods = dynamicMods
  
  pdf(file = paste0(prefix,"_dendrogram.pdf"), wi = 9, he = 6)
  #plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
  #			c("Dynamic Tree Cut", "Merged dynamic"),
  #    		dendroLabels = FALSE, hang = 0.03,
  #			addGuide = TRUE, guideHang = 0.05)
  plotDendroAndColors(geneTree, mergedColors,
                      "Dynamic Tree Cut",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
  dev.off()
  
  df <- data.frame(matrix(0,ncol=length(table(mergedMods)),nrow=length(mergedMods)))
  colnames(df) <- names(table(mergedMods))
  rownames(df) <- colnames(d1)
  for (module in colnames(df))
  {
    df[which(mergedMods==module),module] <- 1
  }
  #df = df[,-which(colnames(df)=="0")]
  colnames(df) <- unlist(lapply(colnames(df),function(x){paste0(prefix,"_module",x)}))
  write.table(df,paste0(prefix,".modules.genes.txt"))
  mergedMEs = mergedMEs[,-which(colnames(mergedMEs)=="ME0")]
  colnames(mergedMEs) = unlist(lapply(colnames(mergedMEs),function(x){gsub("^ME",paste0(prefix,"_module"),x)}))
  return(mergedMEs)
}

calCor <- function(d1,prefix) {
  data <- read.csv(paste0(prefix,".modules.genes.txt"),sep=" ",row.names=1)
  data <- data[data[,1]==0,]
  mat <- d1[,colnames(d1) %in% rownames(data)]
  corMat <- cor(x=mat,use="pairwise.complete.obs",method="pearson")
  write.table(corMat,file=paste0(prefix,".corr.csv"),sep=",",quote=FALSE)
}

args <- commandArgs(trailingOnly=TRUE)
mat <- read.csv("combined_pheno_forconsortium_v1_NPX.tsv",sep='\t',row.names=1)
colnames(mat) <- unlist(lapply(colnames(mat),function(x){l<-unlist(strsplit(x,"\\.")); return(paste(l[1],l[3],sep='_'))}))

#df <- run(mat,8,0.05,args[1])
#df <- run(mat,6,0.32,args[1])
#write.table(df,paste0(args[1],".moduleLevels.csv"),quote=FALSE,sep=",")
calCor(mat,args[1])
