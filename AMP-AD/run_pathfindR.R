library(pathfindR)

sig <- read.table("/mnt/DESeq2.cluster.sig",sep=",",header=F)
sig <- sig[,c("V4","V2","V3")]

#output_df <- run_pathfindR(sig, gene_sets = "KEGG", output_dir = "/mnt/results", max_to_plot=20, iterations = 25, plot_enrichment_chart=TRUE, min_gset_size=50, max_gset_size=500)
#output_df <- run_pathfindR(sig, gene_sets = "GO-BP", output_dir = "/mnt/GO_results", max_to_plot=20, iterations = 25, plot_enrichment_chart=TRUE, min_gset_size=50, max_gset_size=500)
output_df <- run_pathfindR(sig, gene_sets = "Reactome", output_dir = "/mnt/Reactome_results", max_to_plot=20, iterations = 25, plot_enrichment_chart=TRUE, min_gset_size=50, max_gset_size=500)
