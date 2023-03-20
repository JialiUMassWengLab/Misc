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

def plotViolin(adata,prefix,geneList,cluster,meta):
    h5sample = 'sample'
    h5cluster = 'leiden'
    new_cluster_names = {'0':'GC','1':'GC','2':'GC','3':'GC','4':'GC','5':'Empty','6':'TC','7':'GC','8':'unknown','9':'Stromal','10':'GC','11':'Empty','12':'GC','13':'Empty','14':'unknown','15':'EpiC','16':'BEC','17':'Immune','18':'Erythroid','19':'LEC'}
    adata.obs['new_cluster'] = adata.obs[h5cluster].map(new_cluster_names)
    h5cluster = 'new_cluster'

    adata = adata[adata.obs[h5cluster]==cluster,:].copy()
    adata.obs[h5cluster] = adata.obs[h5cluster].astype('category')
    meta = meta.sort_values(['Time','Drug','Treatment','sample'],ascending=True)
    print(meta)

    sns.set(font_scale=1.2)
    palette = []
    colors = sns.color_palette('tab10')
    palette += [colors[0]]*2 + [colors[1]]*2 + [colors[2]] + [colors[3]]*2 + colors[4:7]
    palette2 = []
    for color in palette:
        palette2 += [color]*2
    palette += palette2

    fig,axes = plt.subplots(len(geneList),1,figsize=(12,3*len(geneList)))
    for n,gene in enumerate(list(geneList)):
        sc.pl.violin(adata,gene,groupby=h5sample,ax=axes[n],stripplot=False,rotation=40,order=list(meta.index),palette=palette)
        axes[n].set_xticklabels(meta['Treatment'],rotation=40,ha='right')
        for xtick,color in zip(axes[n].get_xticklabels(),['Blue']*10+['Red']*20):
            xtick.set_color(color)
    plt.tight_layout()
    plt.savefig(prefix+'.violin.geneBySample.%s.pdf' % cluster)
    plt.close()
    

def plotClustDist(adata,prefix,meta):
    new_cluster_names = {'0':'Proliferating GC','1':'Antral mural GC','2':'Luteinizing GC','3':'Atretic GC','4':'Preantral GC','5':'Empty','6':'Theca','7':'Luteinizing GC','8':'Unannotated','9':'Stroma','10':'Luteinizing GC','11':'Empty','12':'Proliferating GC','13':'Empty','14':'Luteinizing TC?','15':'EpiC','16':'BEC','17':'Immune','18':'RBC','19':'LEC'}
    adata.obs['ann_cluster'] = adata.obs['leiden'].map(new_cluster_names)
    h5sample = 'sample'
    h5cluster = 'ann_cluster'

    exclude_groups = ['Empty']
    adata = adata[list(map(lambda x: not x in exclude_groups,adata.obs[h5cluster])),:].copy()
    adata.obs[h5cluster] = adata.obs[h5cluster].astype('category')
    clusterList = sorted(adata.obs[h5cluster].unique())
    totalCells = adata.obs[h5sample].value_counts()
    print(sum(totalCells))

    sns.set(font_scale=1.2)
    fig,axes = plt.subplots(len(clusterList),1,figsize=(12,3*len(clusterList)))
    for n,cluster in enumerate(clusterList):
        cells = adata.obs[adata.obs[h5cluster]==cluster][h5sample].value_counts()
        df = pd.concat([cells,totalCells,meta],join='outer',axis=1)
        df['pct_cells'] = df.apply(lambda x: x[0]/float(x[1]),axis=1)
        df = df.sort_values(['Time','Drug','Treatment'],ascending=True)
        #print(cluster)
        #print(df)
        b = sns.barplot(x=df.index,y='pct_cells',hue='Drug',order=list(df.index),data=df,dodge=False,ax=axes[n])
        b.legend_.remove()
        #plt.setp(b.get_legend().get_texts(),fontsize='8')
        axes[n].set_xticklabels(df['Treatment'],rotation=40,ha='right')
        for xtick,color in zip(axes[n].get_xticklabels(),['Blue']*10+['Red']*20):
            xtick.set_color(color)
        axes[n].set_facecolor('white')
        axes[n].set_title(cluster)        
    plt.tight_layout()
    plt.savefig(prefix+'.Bar.ClusterDistBySample.pdf')
    plt.close()


def plotTreatmentDEgenes(adata,meta,celltype,control,treatment):
    new_cluster_names = {'0':'GC','1':'GC','2':'GC','3':'GC','4':'GC','5':'Empty','6':'TC','7':'GC','8':'unknown','9':'Stromal','10':'GC','11':'Empty','12':'GC','13':'Empty','14':'unknown','15':'EpiC','16':'BEC','17':'Immune','18':'Erythroid','19':'LEC'}
    adata.obs['new_cluster'] = adata.obs['leiden'].map(new_cluster_names)
    adata = adata[adata.obs['new_cluster']==celltype,:].copy()
    #remove luteinized samples for T57
    if not re.search(r'_T57$',control)==None:
        removeSampleList = ['Ov12','Ov18','Ov52','Ov76','Ov92']
        adata = adata[list(map(lambda x: not x in removeSampleList,adata.obs['sample'])),:].copy()
    if not 'base' in adata.uns['log1p']:
        adata.uns['log1p']['base'] = None

    comparison = treatment + '_vs_' + control
    sc.tl.rank_genes_groups(adata,'Treatment',groups=[treatment],reference=control,method='wilcoxon',pts=True,key_added=celltype)

    index = list(map(lambda x: x[0],adata.uns[celltype]['names']))
    logfc = pd.Series(map(lambda x: x[0],adata.uns[celltype]['logfoldchanges']),index=index,name='logFC')
    pv = pd.Series(map(lambda x: x[0],adata.uns[celltype]['pvals_adj']),index=index,name='pvals_adj')
    markers = pd.concat([logfc,pv],join='inner',axis=1)
    #pcts = pd.read_csv('scanpy.%s.treatment.pts.csv' % celltype,index_col=0)
    pcts = adata.uns[celltype]['pts']
    #only keep genes that are expressed in >10% of the cells in either condition
    geneList = list(pcts[(pcts[control]>.1) | (pcts[treatment]>.1)].index)
    markers = markers.loc[map(lambda x: x in geneList,markers.index),:]
    from pybiomart import Dataset
    dataset = Dataset(name='rnorvegicus_gene_ensembl',host='http://www.ensembl.org')
    names = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name'])
    names.index = names['Gene name']
    names = names[~names.index.duplicated()]
    markers = pd.concat([markers,names.drop('Gene name',axis=1)],join='inner',axis=1)
    markers['-logP'] = markers.apply(lambda x: -np.log10(x['pvals_adj']) if x['pvals_adj']>0 else 320,axis=1)
    markers['significant'] = markers.apply(lambda x: x['-logP']>2 and abs(x['logFC'])>1,axis=1)
    markers.drop('-logP',axis=1).to_csv('scanpy.%s.DE.%s.csv' % (celltype,comparison))
    markers = markers[(markers['logFC']>-5) & (markers['logFC']<5)]
    sig = markers[(markers['-logP']>10) & ((markers['logFC']>1) | (markers['logFC']<-1))]
    '''
    from adjustText import adjust_text
    sns.set(font_scale=1.2,context='talk')
    fig,ax = plt.subplots(figsize=(15,12))
    sns.scatterplot(x='logFC',y='-logP',hue='significant',data=markers,s=150,linewidth=1.5,palette=['darkgrey','red'],ax=ax)
    texts = []
    for label,x,y in zip(sig.index,sig['logFC'],sig['-logP']):
        texts.append(plt.text(x,y,label,ha='center',size=18))

    adjust_text(texts,autoalign='xy',arrowprops=dict(arrowstyle="simple, head_width=0.25, tail_width=0.05", color='r', lw=0.5, alpha=0.8))
    plt.savefig('scanpy.%s.volcano.%s.png' % (celltype,comparison))
    plt.close()
    '''

def calCelltypeAvgExpr(adata,meta):
    h5sample = 'sample'
    new_cluster_names = {'0':'Ambient RNA','1':'Lhcgr+ GC','2':'Lhcgr+ GC','3':'Lhcgr- GC','4':'Ambient RNA','5':'Theca','6':'Erythroid','7':'Lhcgr- GC','8':'Stroma','9':'Epithelial','10':'BEC','11':'Lhcgr+ GC','12':'Oocyte','13':'Immune','14':'LEC'}
    adata.obs['ann_cluster'] = adata.obs['leiden'].map(new_cluster_names)
    shown_groups = ['Lhcgr+ GC','Lhcgr- GC','Theca','Stroma','Epithelial','BEC','Oocyte','Immune','LEC']
    h5cluster = 'ann_cluster'

    adata = adata[list(map(lambda x: x in shown_groups,adata.obs['ann_cluster'])),:]
    adata.obs[h5cluster] = adata.obs[h5cluster].astype('category')
    df = sc.get.obs_df(adata,keys=[h5cluster,*list(adata.raw.var_names)],use_raw=True)
    df.iloc[:,1:] = df.iloc[:,1:].apply(lambda x: 2**x-1,axis=0)
    df = df.groupby(h5cluster).agg(np.mean).transpose()
    df.to_csv('scanpy.%s.avgExpr.csv' % h5cluster)
    

def main():
    h5adFile = sys.argv[1]
    sc.settings.verbosity = 3
    sc.logging.print_header()
    adata = sc.read_h5ad(h5adFile)
    meta = adata.obs[['sample','Time','Drug','Treatment']].copy()
    meta.index = meta['sample']
    meta = meta.drop_duplicates().drop('sample',axis=1)
    print(meta)

    #calCelltypeAvgExpr(adata,meta)
    plotClustDist(adata,h5adFile.replace('.h5ad',''),meta)
    #plotViolin(adata,h5adFile.replace('.h5ad',''),['Inha','Inhba','Inhbb','Lhcgr','Cyp19a1','Plxnc1','Fshr','Ccnd2','Pparg','Foxo1','Foxp2','Ctgf','Gls','Amh','Amhr2','Ablim1','Hdac4','Nav2','Ereg','Areg','Ptgs2','Pgr'],'GC',meta)
    #plotViolin(adata,h5adFile.replace('.h5ad',''),['Lhcgr','Cyp17a1','Cyp11a1','Star','Fdx1','Fdxr','Por','Ldlr','Scarb1','Ldb3'],'TC',meta)
    #plotViolin(adata,'5690StreeGenes',['Zeb2','Glg1','Dpp7','Ubc','Ddx17','Rbm6'],'GC',meta)

    #findTreatmentDEgenes(adata,meta,'GC')
    #plotTreatmentDEgenes(adata,meta,'GC','vehicle_T28','FSH_T28')
    #plotTreatmentDEgenes(adata,meta,'GC','vehicle_T28','4221:6_T28')
    #plotTreatmentDEgenes(adata,meta,'GC','vehicle_T28','2599:0.5_T28')
    #plotTreatmentDEgenes(adata,meta,'GC','vehicle_T28','2599:6_T28')
    #plotTreatmentDEgenes(adata,meta,'GC','vehicle_T28','0353:10_T28')
    #plotTreatmentDEgenes(adata,meta,'GC','vehicle_T28','0353:30_T28')
    #plotTreatmentDEgenes(adata,meta,'GC','vehicle_T28','5690:1_T28')
    #plotTreatmentDEgenes(adata,meta,'GC','vehicle_T28','5690:0.3_T28')
    
    #plotTreatmentDEgenes(adata,meta,'GC','vehicle_T57','FSH_T57')
    #plotTreatmentDEgenes(adata,meta,'GC','vehicle_T57','4221:6_T57')
    #plotTreatmentDEgenes(adata,meta,'GC','vehicle_T57','2599:0.5_T57')
    #plotTreatmentDEgenes(adata,meta,'GC','vehicle_T57','2599:6_T57')
    #plotTreatmentDEgenes(adata,meta,'GC','vehicle_T57','0353:10_T57')
    #plotTreatmentDEgenes(adata,meta,'GC','vehicle_T57','0353:30_T57')
    #plotTreatmentDEgenes(adata,meta,'GC','vehicle_T57','5690:1_T57')
    #plotTreatmentDEgenes(adata,meta,'GC','vehicle_T57','5690:0.3_T57')
    

if __name__=='__main__':
    main()
