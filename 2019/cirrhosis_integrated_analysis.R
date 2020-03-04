library(Seurat)
library(ggplot2)
library(cowplot)

dirs <- Sys.glob("GSM40411*")
liver.list <- list()
for (file in dirs) {
    name <- strsplit(file,"_")[[1]][2]
    mat <- Read10X(data.dir=file)
    liver.list[name] <- CreateSeuratObject(counts=mat,project=name)
}
liver.integrated <- AddMetaData(liver.integrated,unlist(lapply(liver.integrated@meta.data$orig.ident, function(x){substr(x,1,nchar(x)-1)})),col.name="disease")

for (i in 1:length(liver.list)) {
    liver.list[[i]] <- NormalizeData(liver.list[[i]], verbose = FALSE)
    liver.list[[i]] <- FindVariableFeatures(liver.list[[i]], selection.method = "vst",
                                            nfeatures = 2000, verbose = FALSE)
}

liver.anchors <- FindIntegrationAnchors(object.list = liver.list, dims = 1:30)
liver.integrated <- IntegrateData(anchorset = liver.anchors, dims = 1:30)
DefaultAssay(liver.integrated) <- "integrated"
liver.integrated <- ScaleData(liver.integrated)
liver.integrated <- RunPCA(liver.integrated, npcs=30)
liver.integrated <- RunUMAP(liver.integrated, reduction = "pca", dims = 1:20)

liver.integrated <- FindNeighbors(liver.integrated,dims=1:20)
liver.integrated <- FindClusters(liver.integrated,resolution=0.5)
markers <- FindAllMarkers(liver.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(markers,"cirrhosis_markers_res05.csv",sep=",",quote=FALSE,row.names=FALSE)

top20 <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
pdf("heatmap.pdf",width=10,height=45)
DoHeatmap(liver.integrated, features = as.character(top20$gene)) + NoLegend()
dev.off()

#find cell type specific DE genes and plot scatter plots
stellet <- subset(liver.integrated, idents = 7)
Idents(stellet) <- "disease"
avg.stellet <- log1p(AverageExpression(stellet)$RNA)
stellet.de.genes <- FindMarkers(stellet,ident.1="healthy",ident.2="cirrhotic",min.pct=0.2,logfc.threshold = 0.4)
pdf("stellet.cell.scatter.pdf")
p1 <- ggplot(avg.stellet,aes(healthy,cirrhotic)) + geom <- point() + ggtitle("Stellet cells")
LabelPoints(plot=p1, points = rownames(head(stellet.de.genes,n=20)),repel=TRUE)
dev.off()

