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
    meta = loadMetaData('../metadata.tsv')
    h5adFile = sys.argv[1]
    sc.settings.verbosity = 3
    sc.logging.print_header()
    adata = sc.read_h5ad(h5adFile)
    shown_groups = ['GC','Theca','Stroma','Epithelial','BEC','Oocyte','Immune','LEC']
    adata = adata[list(map(lambda x: x in shown_groups,adata.obs['ann_cluster'])),:].copy()
    adatas = {}
    adatas.update({'GC': adata[adata.obs['ann_cluster']=='GC',:].copy()})
    adatas.update({'TC': adata[adata.obs['ann_cluster']=='Theca',:].copy()})
    colors = {}
    colors['GC'] = {'Proliferating':'blue','Antral':'red','Atretic':'green','Preantral':'orange'}
    colors['TC'] = {'Native':'green','Lhcgr-':'orange','2':'grey'}
    '''
    #barplot
    from matplotlib.backends.backend_pdf import PdfPages
    pdf = PdfPages('Figure3.pdf')
    sns.set(context='talk',font_scale=1.2)
    for celltype in ['GC']:
        adata1 = adatas[celltype]
        totalCells = adata1.obs['sample'].value_counts()
        array = []
        for cluster in adata1.obs[celltype+'_cluster'].unique():
            cells = adata1.obs[adata1.obs[celltype+'_cluster']==cluster]['sample'].value_counts()
            df1 = pd.concat([cells,totalCells],join='outer',axis=1).fillna(0)
            df1[cluster] = df1.apply(lambda x: x[0]/float(x[1]),axis=1)
            array.append(df1[cluster])
        df = pd.concat(array,join='outer',axis=1).fillna(0)
        df = pd.concat([df,meta[['FSH_dose','hCG_dose']]],join='inner',axis=1)
        df = df.sort_values(['FSH_dose','hCG_dose'],ascending=True)
        if celltype == 'GC':
            df = df[['Atretic','Preantral','Proliferating','Antral']]
        else:
            df = df[['Native','Lhcgr-','2']]
        print(df)

        fig,ax = plt.subplots(figsize=(16,6))
        df.plot(kind='bar',stacked=True,xlabel='Sample',ylabel='Fraction of cells',edgecolor='white',color=colors[celltype],linewidth=1.5,width=.85,ax=ax)
        ax.set_xticklabels(df.index,rotation=40,ha='right')
        ax.set_facecolor('white')
        ax.set_title(celltype)
        ax.legend(bbox_to_anchor=(1.01,1),loc='upper left')
        plt.tight_layout()
        plt.savefig(pdf,format='pdf')
        plt.close()
    pdf.close()
    '''
    #stacked violin plot or matrixplot
    geneLists = {'GC': ['Pawr','Lox','Ctgf','Mitf','Ablim1','Dock9','Amhr2','Amh',\
                        'Top2a','Mki67','Prc1','Smc4','Lhcgr','Pappa','Alas1','Cyp19a1'],
                 'TC': ['Ano1','Smoc2','Lhcgr','Mark1','Timp1','Fam126a','Mob4','Star']}
    widths = {'GC': 16, 'TC': 12}
    excluded = {'GC': ['4','5'], 'TC': ['2']}
    orders = {'GC': ['Atretic','Preantral','Proliferating','Antral'],
              'TC': ['Native','Lhcgr-']}
    marker_genes_dict = {'Antral': ['Cyp19a1','Lhcgr','Alas1','Inhba','Hsd3b2'],
                         'Atretic': ['Pawr','Mitf','Pik3ip1','Ghr','Lox'],
                         'Preantral': ['Amh','Amhr2','Pcsk6','Igfbp5','Dock9'],
                         'Proliferating': ['Top2a','Mki67','Ccnd2','Prc1','Smc4'],
    }
    
    from matplotlib.backends.backend_pdf import PdfPages
    pdf = PdfPages('Figure3d.pdf')
    #sns.set(context='talk',font_scale=1.2)
    for celltype in ['GC']:
        adata1 = adatas[celltype]
        adata1 = adata1[list(map(lambda x: not x in excluded[celltype],adata1.obs[celltype+'_cluster'])),:].copy()

        #fig,ax = plt.subplots(figsize=(widths[celltype],8))
        #vi = sc.pl.stacked_violin(adata1,geneLists[celltype],groupby=celltype+'_cluster',standard_scale='var',categories_order=orders[celltype],title=celltype,return_fig=True,ax=ax)
        plt.rcParams.update({'font.size': 15})
        fig,ax = plt.subplots(figsize=(18,6))
        vi = sc.pl.matrixplot(adata,marker_genes_dict,celltype+'_cluster',standard_scale='var',var_group_rotation=40,ax=ax,return_fig=True)
        ax1 = vi.get_axes()['mainplot_ax']
        labels = ax1.get_xticklabels()
        ax1.set_xticklabels(labels,rotation=40,ha='right',fontsize=15,fontstyle='italic')
        plt.tight_layout()
        plt.savefig(pdf,format='pdf')
        plt.close()
    pdf.close()
    '''
    #dot plot
    geneLists = ['Pawr','Lox','Amhr2','Amh','Top2a','Mki67','Lhcgr','Cyp19a1','Alas1']
    excluded = {'GC': ['4','5'], 'TC': ['2']}
    mapping = {'Vehicle':'Vehicle','hCG_m':'hCG (0.0375)','hCG_hi':'hCG (1.2)',
               'FSH_1':'FSH (1)','FSH_3':'FSH (3)',
               'FSH_3 + hCG_m':'FSH (3) + hCG (0.0375)','FSH_3 + hCG_hi':'FSH (3) + hCG (1.2)',
               'FSH_10':'FSH (10)','FSH_10 + hCG_m':'FSH (10) + hCG (0.0375)',
               'FSH_10 + hCG_hi':'FSH (10) + hCG (1.2)'}
    orders = ['Vehicle','hCG (0.0375)','hCG (1.2)','FSH (1)','FSH (3)','FSH (3) + hCG (0.0375)',
              'FSH (3) + hCG (1.2)','FSH (10)','FSH (10) + hCG (0.0375)','FSH (10) + hCG (1.2)']
    from matplotlib.backends.backend_pdf import PdfPages
    pdf = PdfPages('Figure3c.pdf')
    sns.set(context='talk',font_scale=1.2)
    for celltype in ['GC']:
        adata1 = adatas[celltype]
        adata1 = adata1[list(map(lambda x: not x in excluded[celltype],adata1.obs[celltype+'_cluster'])),:].copy()
        adata1.obs['treatment1'] = adata1.obs['treatment'].map(mapping)

        fig,ax = plt.subplots(figsize=(16,10))
        vi = sc.pl.dotplot(adata1,geneLists,groupby='treatment1',standard_scale='var',title=celltype,categories_order=orders,swap_axes=True,return_fig=True,ax=ax)
        ax1 = vi.get_axes()['mainplot_ax']
        ax1.set_facecolor('white')
        labels = ax1.get_xticklabels()
        ax1.set_xticklabels(labels,rotation=30,ha='right')
        ax1.set_yticklabels(ax1.get_yticklabels(),fontstyle='italic')
        ax2 = vi.get_axes()['size_legend_ax']
        ax2.set_facecolor('white')
        plt.tight_layout()
        plt.savefig(pdf,format='pdf')
        plt.close()
    pdf.close()
    '''
    
    
if __name__=='__main__':
    main()
