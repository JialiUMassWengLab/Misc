library(WGCNA)

plotPower <- function(d1,prefix) {
    powers = c(c(1:10), seq(from = 12, to=50, by=2))
    sft = pickSoftThreshold(d1, powerVector = powers, verbose = 5)
    cex1 = 0.9
    pdf(file = paste0(prefix,"_scaleFreePlot.pdf"))
    plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",main = paste("Scale independence"))
    text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red")
    dev.off()
}


run <- function(geneList,datExpr,sftPower,MEDissThres,prefix) {
    d1 <- t(datExpr[which(row.names(datExpr) %in% geneList),])

    adjacency = adjacency(d1, power = sftPower)
    TOM = TOMsimilarity(adjacency)
    dissTOM = 1-TOM
    geneTree = hclust(as.dist(dissTOM), method = "average")
    minModuleSize = 30
    dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,deepSplit = 2, pamRespectsDendro = FALSE,minClusterSize = minModuleSize)
    dynamicColors = labels2colors(dynamicMods)
    print(table(dynamicColors))

    #sizeGrWindow(8,6)
    #plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05,main = "Gene dendrogram and module colors")
    #dev.off()

    MEList = moduleEigengenes(d1, colors = dynamicMods)
    MEs = MEList$eigengenes
    colnames(MEs) = unlist(lapply(colnames(MEs),function(x){gsub("^ME",prefix,x)}))

    MEDiss = 1-cor(MEs)
    METree = hclust(as.dist(MEDiss), method = "average")
    pdf(file = paste0(prefix,"_cut.pdf"), wi = 7, he = 6)
    plot(METree, main = "Clustering of module eigengenes",xlab = "", sub = "")
    abline(h=MEDissThres, col = "red")
    dev.off()
    
    merge = mergeCloseModules(d1, dynamicMods, cutHeight = MEDissThres, iterate = FALSE)
    mergedMEs = merge$newMEs
    mergedMods = merge$colors
    print(table(mergedMods))
    #mergedMods = dynamicMods

    df <- data.frame(matrix(0,ncol=length(table(mergedMods)),nrow=length(mergedMods)))
    colnames(df) <- names(table(mergedMods))
    rownames(df) <- colnames(d1)
    for (module in colnames(df))
    {
       df[which(mergedMods==module),module] <- 1
    }
    df = df[,-which(colnames(df)=="0")]
    colnames(df) <- unlist(lapply(colnames(df),function(x){paste0(prefix,x)}))
    write.table(df,paste0("modules.",prefix,".genes.txt"))
    mergedMEs = mergedMEs[,-which(colnames(mergedMEs)=="ME0")]
    colnames(mergedMEs) = unlist(lapply(colnames(mergedMEs),function(x){gsub("^ME",prefix,x)}))
    return(mergedMEs)
}

args <- commandArgs(trailingOnly=TRUE)
H <- t(read.csv("H_limma_p8_a0.mat.tsv",sep="\t",row.names=1))
H <- t(apply(H,1,function(x){x/sum(x)}))
glialGenes <- row.names(H[which(H[,1]>0.25 | H[,3]>0.25 | H[,4]>0.25 | H[,5]>0.25 | H[,7]>0.25 | H[,8]>0.25),])
neuronGenes <- row.names(H[which((H[,2] > 0.25 | H[,6] > 0.25) & !row.names(H) %in% glialGenes),])
#houseGenes <- row.names(H[which(!row.names(H) %in% glialGenes & !row.names(H) %in% neuronGenes),])
houseGenes <- row.names(H[apply(H,1,function(x){all(x<=0.25)}),])
assignedGenes <- row.names(H[which(!row.names(H) %in% houseGenes),])
print(c(length(neuronGenes),length(glialGenes),length(houseGenes)))

datExpr <- read.csv(paste0(args[1],".csv"),row.names=1)
row.names(datExpr) <- unlist(lapply(row.names(datExpr),function(x){gsub(".[0-9]*$","",x)}))

plotPower(t(datExpr[which(row.names(datExpr) %in% assignedGenes),]),"assigned")
#df1 = run(neuronGenes,datExpr,20,0.12,"neuron")
#df2 = run(glialGenes,datExpr,10,0.26,"glia")
#df3 = run(houseGenes,datExpr,8,0.25,"unassigned")
#df = run(assignedGenes,datExpr,10,"assigned")

#df = cbind(df1,df2)
#write.table(df,paste0(args[1],".moduleLevels.csv"),quote=FALSE,sep=",")
