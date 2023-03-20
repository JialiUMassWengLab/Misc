#! /usr/bin/python

import re
import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

def refMarkersPlot(h5adFile,refMarkers,prefix):
    sc.settings.verbosity = 3
    sc.logging.print_header()
    adata = sc.read_h5ad(h5adFile)
    h5cluster = 'leiden'
    optigon_clusters = {'0':'Empty_droplets','1':'Fshr+ GC','2':'Lhcgr+ GC','3':'Native GC','4':'Empty_droplets','5':'Theca','6':'Erythroid','7':'Amh+ GC','8':'Stroma','9':'Epithelial','10':'BEC','11':'Doublets','12':'Oocyte','13':'Immune','14':'LEC'}

    from pybiomart import Dataset
    dataset = Dataset(name='rnorvegicus_gene_ensembl',host='http://www.ensembl.org')
    names = dataset.query(attributes=['ensembl_gene_id','external_gene_name',\
                                      'mmusculus_homolog_ensembl_gene',\
                                      'mmusculus_homolog_associated_gene_name'])
    names.index = names['Gene name']
    #names = names[~names.index.duplicated()]
    names = names.loc[(~names.duplicated(subset='Gene name',keep=False)) & \
                       (~names.duplicated(subset='Mouse gene name',keep=False)),:]
    print(names)
        
    ref = pd.read_csv(refMarkers,index_col=0)
    sns.set(font_scale=1.5)
    from matplotlib.backends.backend_pdf import PdfPages
    pdf = PdfPages(prefix+'.celltype.anno.pdf')
    for cluster in list(ref['cluster'].unique()):
        ref1 = ref[ref['cluster']==cluster].copy()
        ref1 = pd.concat([ref1,names['Mouse gene name']],join='inner',axis=1)
        ref1 = ref1.dropna(subset=['Mouse gene name']).drop_duplicates(subset=['Mouse gene name'])
        geneList = ref1['Mouse gene name']
        geneList = geneList[map(lambda x: x in adata.raw.var_names,geneList)].head(n=25)
        print(cluster)
        print(geneList)
        '''
        fig,axes = plt.subplots(3,4,figsize=(30,15))
        for n,gene in enumerate(list(geneList)):
            i = int(n/4)
            j = n % 4
            sc.pl.violin(adata,gene,groupby='seurat_clusters',ax=axes[i][j],size=1.2)
        ''' 
        fig,ax = plt.subplots(figsize=(16,12))
        vi = sc.pl.stacked_violin(adata,geneList,groupby=h5cluster,title=optigon_clusters[str(cluster)],standard_scale='var',return_fig=True,ax=ax)
        ax1 = vi.get_axes()['mainplot_ax']
        labels = ax1.get_xticklabels()
        ax1.set_xticklabels(labels,rotation=40,ha='right')
        plt.tight_layout()
        plt.savefig(pdf,format='pdf')
        plt.close()

        fig,ax = plt.subplots(figsize=(16,12))
        dp = sc.pl.dotplot(adata,geneList,groupby=h5cluster,title=optigon_clusters[str(cluster)],standard_scale='var',return_fig=True,ax=ax)
        ax1 = dp.get_axes()['mainplot_ax']
        labels = ax1.get_xticklabels()
        ax1.set_xticklabels(labels,rotation=40,ha='right')
        ax1.set_facecolor('white')
        plt.tight_layout()
        plt.savefig(pdf,format='pdf')
        plt.close()

    pdf.close()
    

def checkClustDist(h5adFile,field):
    sc.settings.verbosity = 3
    sc.logging.print_header()
    adata = sc.read_h5ad(h5adFile)
    h5sample = field

    totalCategories = len(adata.obs[h5sample].unique())
    nrows = int(totalCategories/4) if (totalCategories % 4)==0 else int(totalCategories/4.0)+1
    sns.set(font_scale=1.5)
    fig,axes = plt.subplots(nrows,4,figsize=(20,nrows*4))
    for n,sample in enumerate(sorted(adata.obs[h5sample].unique())):
        i = int(n/4)
        j = n % 4
        print(sample,i,j)
        sc.pl.umap(adata,color=h5sample,frameon=False,groups=[sample],title=sample,legend_loc=None,na_in_legend=False,palette=['blue'],ax=axes[i][j])
    plt.savefig(h5adFile.replace('.h5ad','.%sDist.png' % field))
    plt.close()
    

def plotUMAPwithNewCluster(h5adFile):
    sc.settings.verbosity = 3
    sc.logging.print_header()
    adata = sc.read_h5ad(h5adFile)
    new_cluster_names = {'0':'Proliferating GC','1':'Antral GC','2':'Luteinizing GC','3':'Atretic GC','4':'Preantral GC','5':'Empty','6':'Theca','7':'Luteinizing GC','8':'Unannotated','9':'Stroma','10':'Luteinizing GC','11':'Empty','12':'Proliferating GC','13':'Empty','14':'Luteinizing TC?','15':'Epithelium','16':'BEC','17':'Immune','18':'Reticulocyte','19':'LEC'}
    adata.obs['ann_cluster'] = adata.obs['leiden'].map(new_cluster_names)
    shown_groups = list(set(new_cluster_names.values()))
    shown_groups.remove('Empty')
    print(shown_groups)

    sns.set(font_scale=1.5)
    fig,ax = plt.subplots(figsize=(10,10))
    sc.pl.umap(adata,color='ann_cluster',frameon=False,groups=shown_groups,legend_loc='on data',na_in_legend=False,s=10,ax=ax)
    plt.savefig(h5adFile.replace('.h5ad','.clusterDist.png'))
    plt.close()
    

def main():
    #refMarkersPlot(sys.argv[1],'../Optigon_snRNAseq_RatOvaries20_June2022/scanpy.marker.genes.csv','ratOvary')
    #checkClustDist(sys.argv[1],'Drug')
    #checkClustDist(sys.argv[1],'Treatment')
    plotUMAPwithNewCluster(sys.argv[1])
    

if __name__=='__main__':
    main()
