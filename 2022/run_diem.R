library(diem)
library(Seurat)

args <- commandArgs(trailingOnly=TRUE)
sample <- args[1]

counts <- Read10X_h5(paste0(sample,"_filtered_feature_bc_matrix.h5"))
data <- create_SCE(counts, name=sample)
mt_genes <- grep(pattern = "^mt-", x = rownames(data@gene_data),
                 ignore.case = TRUE, value = TRUE)
data <- get_gene_pct(x = data, genes = mt_genes, name = "pct.mt")
genes <- grep(pattern = "^malat1$", x = rownames(data@gene_data),
              ignore.case = TRUE, value = TRUE)
data <- get_gene_pct(x = data, genes = genes, name = "MALAT1")

pdf(paste0(sample,".diem.QC_plots.pdf"))
require(gridExtra)
p1 <- plot_data(mb_small, feat_x = "total_counts", feat_y = "n_genes", 
                log_x = TRUE, log_y = TRUE, ret = TRUE, data_type = "all")
p2 <- plot_data(mb_small, feat_x = "n_genes", feat_y = "pct.mt", 
                log_x = TRUE, log_y = FALSE, ret = TRUE, data_type = "all")
#p3 <- plot_data(mb_small, feat_x = "n_genes", feat_y = "MALAT1", 
#                log_x = TRUE, log_y = FALSE, ret = TRUE, data_type = "all")
#p4 <- plot_data(mb_small, feat_x = "pct.mt", feat_y = "MALAT1", 
#                log_x = FALSE, log_y = FALSE, ret = TRUE, data_type = "all")
grid.arrange(p1, p2, ncol = 2)

barcode_rank_plot(data, title = sample)
#data <- set_debris_test_set(data, min_counts = 500)
data <- set_debris_test_set(data, min_counts = 700)
length(data@test_set)
length(data@bg_set)

data <- filter_genes(data, cpm_thresh = 0)
genes <- gene_data(data)
summary(genes)

data <- get_pcs(data,n_var_genes = 1500,n_pcs = 25,min_genes=200)
data <- init(data,k_init = 15,nstart_init = 30,min_size_init = 10)
data <- run_em(data,threads=8)
data <- assign_clusters(data)

drop_data <- droplet_data(data, type="test")
table(drop_data[,"Cluster"])
tapply(drop_data[,"n_genes"],
       drop_data[,"Cluster"],
       mean)
tapply(drop_data[,"pct.mt"],
       drop_data[,"Cluster"],
       mean)
tapply(drop_data[,"MALAT1"],
       drop_data[,"Cluster"],
       mean)

data <- estimate_dbr_score(data,thresh_genes = 100,thresh_p = 0.05)
de_genes <- debris_genes(data, p_adj = 0.05)

require(gridExtra)
p1 <- plot_clust(data, feat_x = "n_genes", feat_y = "score.debris", 
                 log_x = TRUE, ret = TRUE)
p2 <- plot_clust(data, feat_x = "total_counts", feat_y = "score.debris",
                 log_x = TRUE, ret = TRUE)
p3 <- plot_clust(data, feat_x = "pct.mt", feat_y = "score.debris", 
                 log_x = FALSE, ret = TRUE)
p4 <- plot_clust(data, feat_x = "MALAT1", feat_y = "score.debris", 
                 log_x = FALSE, ret = TRUE)
grid.arrange(p1, p2, p3, p4, ncol = 2)

sm <- summarize_clusters(data, top_n = 50)
write.table(sm,paste0(sample,".diem.clusters.csv"),quote=FALSE,sep=",")
par(mfrow=c(1,2))
plot(sm[,"avg_n_genes"], sm[,"avg_dbr_score"], pch= NA,
     xlab="Average number Genes", ylab="Avergae debris score")
text(sm[,"avg_n_genes"], sm[,"avg_dbr_score"], sm[,"Cluster"])
plot(sm[,"avg_n_counts"], sm[,"avg_dbr_score"], pch= NA,
     xlab="Average total counts", ylab="Avergae debris score")
text(sm[,"avg_n_counts"], sm[,"avg_dbr_score"], sm[,"Cluster"])
dev.off()

data <- call_targets(data,thresh_score = 0.5,min_genes = 0)
drop_data <- droplet_data(data, type="test")
write.table(drop_data,paste0(sample,".diem.droplets.csv"),quote=FALSE,sep=",")

tapply(drop_data[,"n_genes"],
       drop_data[,"Call"],
       mean)
tapply(drop_data[,"score.debris"],
       drop_data[,"Call"],
       mean)
tapply(drop_data[,"pct.mt"],
       drop_data[,"Call"],
       mean)
tapply(drop_data[,"MALAT1"],
       drop_data[,"Call"],
       mean)
tapply(drop_data[,"Cluster"],
       drop_data[,"Call"],
       table)

