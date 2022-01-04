library(Seurat)
library(BUSpaRse)
library(tidyverse)
library(DropletUtils)
library(ggplot2)
library(schex)

#Read from BUStools output
#res_mat <- read_count_output("2019-010/",name="counts",tcc=FALSE)

#Read from 10X output
res_mat <- Read10X("/home/jzhuang/AMP/Source/brain_scRNAseq/")
tot_counts <- Matrix::colSums(res_mat)
summary(tot_counts)

bc_rank <- barcodeRanks(res_mat)
pdf("UMI_cdf.pdf")
qplot(bc_rank$total, bc_rank$rank, geom = "line") +
  geom_vline(xintercept = metadata(bc_rank)$knee, color = "blue", linetype = 2) +
  geom_vline(xintercept = metadata(bc_rank)$inflection, color = "green", linetype = 2) +
  annotate("text", y = 1000, x = 1.5 * c(metadata(bc_rank)$knee, metadata(bc_rank)$inflection),
           label = c("knee", "inflection"), color = c("blue", "green")) +
  scale_x_log10() +
  scale_y_log10() +
  labs(y = "Barcode rank", x = "Total UMI count")
dev.off()

res_mat <- res_mat[, tot_counts > metadata(bc_rank)$inflection]
obj <- CreateSeuratObject(res_mat, project="brain", min.cells = 3, min.features=200) %>% 
  NormalizeData(verbose = FALSE) %>% 
  ScaleData(verbose = FALSE) %>% 
  FindVariableFeatures(verbose = FALSE)

obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
#obj <- AddMetaData(obj, metadata = PercentageFeatureSet(obj, pattern = "^MT-"), col.name = "percent.mt")
pdf("UMI_genes_corr.pdf")
FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()

obj <- RunPCA(obj, verbose = FALSE, npcs = 30)
pdf("elbowPlot.pdf")
ElbowPlot(obj, ndims = 30)
dev.off()

pdf("pca.pdf")
DimPlot(obj,reduction="pca")
dev.off()

obj <- RunUMAP(obj, dims = 1:13)
pdf("umap.pdf")
DimPlot(obj,reduction="umap")
dev.off()

obj <- FindNeighbors(obj,dims=1:13)
obj <- FindClusters(obj,resolution=0.5)
pdf("clustering_r0.5.pdf")
DimPlot(obj,reduction="umap",label=TRUE,pt.size=0.5) + NoLegend()
dev.off()

markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
write.table(markers,file="marker.genes.csv",sep=",")
options(tibble.print_max = Inf)
markers %>% group_by(cluster) %>% slice_max(n = 8, order_by = avg_logFC)

saveRDS(obj,"seuratOjb.rds")

#cells = obj@active.ident[obj@active.ident==15 | obj@active.ident==19]
#pdf("gene_expr_plot.pdf")
#FeaturePlot(obj,features=c("C1QA","C1QB","C1QC","FCGBP","FCGR3A","CD14","C3","FYB"),cells=names(cells))
#dev.off()

obj <- make_hexbin(obj, nbins = 100, dimension_reduction = "UMAP")
genes <- read.csv("H_limma_p8_a0.comp6.genes.csv",row.names=1,header=FALSE)
pdf(file="comp6_topGenes_hexbin.pdf",onefile=TRUE)
for (g in genes[,1])
{
   if (g %in% obj@assays[[1]]@data@Dimnames[[1]])
   {
      print(plot_hexbin_gene(obj,type="scale.data",gene=g,action="mean",xlab = "UMAP1", ylab = "UMAP2"))
   }
}
dev.off()

new.cluster.ids <- c("Oligodendrocyte","Excitatory Neu","Excitatory Neu","Oligodendrocyte","Astrocyte","Astrocyte","OPC","Inhibitory Neu","Excitatory Neu","Excitatory Neu","Inhibitory Neu","Inhibitory Neu","Excitatory Neu","Excitatory Neu","Inhibitory Neu","Microglia","Astrocyte","Excitatory Neu","Other","Microglia","Excitatory Neu","Endothelial cell","Other")
names(new.cluster.ids) <- levels(obj)
obj <- RenameIdents(obj,new.cluster.ids)
avg <- AverageExpression(obj,slot="data")[[1]]
write.table(avg,file="average_celltype_expr.csv",sep=",",quote=FALSE)
