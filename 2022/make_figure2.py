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

def loadMetaData(fn):
    df = pd.read_csv(fn,sep='\t',header=None,index_col=0)
    df.index = list(map(lambda x: 'CB'+str(x),df.index))
    df.columns = ['FSH_dose','FSH','hCG_dose','hCG']
    df['FSH_dose'] = df['FSH_dose'].astype(int)
    df['treatment'] = df.apply(lambda x: '049:'+str(x['FSH_dose'])+'_'+'302:'+str(x['hCG_dose']),axis=1)
    mapping = {'049:0_302:0.0':'Vehicle','049:0_302:0.0375':'hCG_m','049:0_302:1.2':'hCG_hi',
               '049:1_302:0.0':'FSH_1','049:3_302:0.0':'FSH_3',
               '049:3_302:0.0375':'FSH_3 + hCG_m','049:3_302:1.2':'FSH_3 + hCG_hi',
               '049:10_302:0.0':'FSH_10','049:10_302:0.0375':'FSH_10 + hCG_m',
               '049:10_302:1.2':'FSH_10 + hCG_hi'}
    df['treatment1'] = df['treatment'].map(mapping)
    return df
    

def main():
    #meta = loadMetaData('metadata.tsv')
    h5adFile = sys.argv[1]
    sc.settings.verbosity = 3
    sc.logging.print_header()
    adata = sc.read_h5ad(h5adFile)
    shown_groups = ['GC','Theca','Stroma','Epithelial','BEC','Oocyte','Immune','LEC','Luteinizing GC']
    adata = adata[list(map(lambda x: x in shown_groups,adata.obs['ann_cluster'])),:].copy()
    
    fig,ax = plt.subplots(figsize=(10,12))
    sc.pl.umap(adata,color=['Fshr','Cyp19a1','Amh','Lhcgr'],ncols=2,s=12,frameon=False)
    plt.savefig('genes.png')
    plt.close()

    from matplotlib.backends.backend_pdf import PdfPages
    pdf = PdfPages('Figure2.pdf')
    fig,ax = plt.subplots(figsize=(10,10))
    sc.pl.umap(adata,color='new_cluster',legend_loc='on data',legend_fontsize='xx-large',s=12,frameon=False,ax=ax)
    plt.savefig(pdf,format='pdf')
    plt.close()

    #ncell = adata.obs['new_cluster'].value_counts()
    ncell = adata.obs['ann_cluster'].value_counts()
    fig,ax = plt.subplots(figsize=(12,10))
    ncell.plot(kind='bar',xlabel='Cell type',ylabel='Number of cells',edgecolor='k',color='darkgrey',linewidth=1.8,fontsize=20,ax=ax)
    ax.set_xticklabels(ax.get_xticklabels(),rotation=40,ha='right')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.tight_layout()
    plt.savefig(pdf,format='pdf')
    plt.close()

    discard_groups = ['Granulosa_4','Granulosa_5','Theca_2']
    adata = adata[list(map(lambda x: not x in discard_groups,adata.obs['new_cluster'])),:].copy()
    marker_genes_dict = {'BEC': ['Tek','Emcn','Vwf'],
                         'EpiC': ['Fras1','Ppp2r2b','Shroom3'],
                         'Granulosa': ['Fshr','Cyp19a1','Grem2','Nr5a2'],
                         'Immune': ['Ptprc','Csf1r','Fyb'],
                         'LEC': ['Ccl21','Prox1','Reln'],
                         'Luteinizing': ['Prss35','Parm1','Cemip'],
                         'Oocyte': ['Ooep','Dazl','Bmp15'],
                         'Stroma': ['Dcn','Pdgfra','Col6a3'],
                         'Theca': ['Cyp17a1','Ano1','Smoc2'],
    }
    plt.rcParams.update({'font.size': 15})
    #g = sc.pl.matrixplot(adata,marker_genes_dict,'new_cluster',standard_scale='var',var_group_rotation=40,figsize=(20,8),return_fig=True)
    fig,ax = plt.subplots(figsize=(18,7))
    g = sc.pl.matrixplot(adata,marker_genes_dict,'ann_cluster',standard_scale='var',var_group_rotation=40,ax=ax,return_fig=True)
    ax = g.get_axes()['mainplot_ax']
    ax.set_xticklabels(ax.get_xticklabels(),rotation=40,ha='right',fontsize=15,fontstyle='italic')
    plt.tight_layout()
    plt.savefig(pdf,format='pdf')
    plt.close()

    pdf.close()
    
    
if __name__=='__main__':
    main()
