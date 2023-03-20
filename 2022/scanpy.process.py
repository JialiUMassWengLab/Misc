#! /usr/bin/python

import re
import os
import sys
import glob
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

def main():
    sc.settings.verbosity = 3
    sc.logging.print_header()

    files = glob.glob('*_filtered_feature_bc_matrix.h5')
    dataList = []
    for fn in files:
        adata = sc.read_10x_h5(fn)
        adata.var_names_make_unique()
        adata.obs['sample'] = [fn.split('_')[0]] * adata.obs.shape[0]
        if fn.split('_')[0] in ['Ov11','Ov113','Ov117']:
            sc.pp.filter_cells(adata, min_counts=1400)
        else:
            sc.pp.filter_cells(adata, min_counts=1000)
        dataList.append(adata)

    from matplotlib.backends.backend_pdf import PdfPages
    pdf = PdfPages('scanpy.plots.pdf')
    #create anndata and processing
    adata = dataList[0].concatenate(dataList[1:])
    adata.obs = adata.obs.drop(['n_counts'],axis=1)
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=10)
    #remove 2 dominant genes that are unannotated
    keep_genes = [name for name in adata.var_names if not name in ['Rn60_1_2212.2','Rn60_1_2212.4']]
    adata = adata[:, keep_genes]
    #calculate percent mitochondrial reads
    mito_gene_list = list(sc.queries.mitochondrial_genes('rnorvegicus')['external_gene_name'])
    adata.var['mt'] = np.in1d(adata.var_names, mito_gene_list)
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], groupby='sample',stripplot=False, multi_panel=True, rotation=40)
    plt.savefig(pdf,format='pdf')
    plt.close()

    fig,axes = plt.subplots(1,2,figsize=(11,6))
    sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt', ax=axes[0])
    sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', ax=axes[1])
    plt.savefig(pdf,format='pdf')
    plt.close()

    #normalization & variable gene selection
    adata = adata[adata.obs.n_genes_by_counts < 6000, :]
    adata = adata[adata.obs.pct_counts_mt < 20, :]
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, min_mean=0.025, max_mean=3, min_disp=0.5)
    sc.pl.highly_variable_genes(adata)
    plt.savefig(pdf,format='pdf')
    plt.close()

    adata.raw = adata
    adata = adata[:, adata.var.highly_variable]
    sc.pp.regress_out(adata, ['total_counts','pct_counts_mt'])

    #dimension reduction
    sc.pp.scale(adata, max_value=8)
    sc.tl.pca(adata, svd_solver='arpack')
    sc.pl.pca(adata,color='sample')
    sc.pl.pca_variance_ratio(adata, log=False)
    plt.savefig(pdf,format='pdf')
    plt.close()

    #used 25 components and Leiden
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=25)
    sc.tl.umap(adata)
    sc.tl.leiden(adata,resolution=.45)
    sc.pl.umap(adata,color='leiden', legend_loc='on data', frameon=False)
    plt.savefig(pdf,format='pdf')
    plt.close()
    pdf.close()

    #find marker genes
    sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon', pts=True)
    array = []
    for i in range(len(adata.uns['rank_genes_groups']['names'][0])):
        index = list(map(lambda x: x[i],adata.uns['rank_genes_groups']['names']))
        score = pd.Series(map(lambda x: x[i],adata.uns['rank_genes_groups']['scores']),index=index,name='scores')
        logfc = pd.Series(map(lambda x: x[i],adata.uns['rank_genes_groups']['logfoldchanges']),index=index,name='logFC')
        pv = pd.Series(map(lambda x: x[i],adata.uns['rank_genes_groups']['pvals_adj']),index=index,name='pvals_adj')
        markers1 = pd.concat([score,logfc,pv],join='inner',axis=1)
        markers1 = markers1[markers1['logFC']>0.8].head(n=200)
        markers1['cluster'] = [i] * markers1.shape[0]
        array.append(markers1)
        
    markers = pd.concat(array,join='inner',axis=0)
    markers.to_csv('scanpy.marker.genes.csv')
    adata.uns['rank_genes_groups']['pts'].to_csv('scanpy.marker.pts.csv')

    meta = pd.read_csv('metadata.csv',index_col=0)
    meta.index = list(map(lambda x: 'Ov'+str(x),meta.index))
    new_obs = adata.obs.merge(meta[['Time','Drug','Treatment']], left_on='sample', right_index=True, how='inner')
    adata.obs = new_obs
    adata.write('scanpy.h5ad')


if __name__=='__main__':
    main()
