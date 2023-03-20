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
    new_cluster_names = {'0':'GC','1':'GC','2':'GC','3':'GC','4':'GC','5':'Empty','6':'TC','7':'GC','8':'Unannotated','9':'Stromal','10':'GC','11':'Empty','12':'GC','13':'Empty','14':'TC','15':'EpiC','16':'BEC','17':'Immune','18':'RBC','19':'LEC'}
    adata.obs['new_cluster'] = adata.obs[h5cluster].map(new_cluster_names)
    h5cluster = 'new_cluster'
    treatments = ['vehicle','FSH','FSH+hCG','4221:6','0353:10','0353:30','2599:0.5','2599:6','5690:0.3','5690:1']
    treatment_order = [i+'_T28' for i in treatments] + [i+'_T57' for i in treatments]

    adata = adata[adata.obs[h5cluster]==cluster,:].copy()
    adata.obs[h5cluster] = adata.obs[h5cluster].astype('category')
    meta['rank'] = meta.apply(lambda x: treatments.index(x['Treatment'].split('_')[0]) \
                              if x['Time']=='T28' else \
                              treatments.index(x['Treatment'].split('_')[0])+10,axis=1)
    meta = meta.sort_values(['rank','sample'],ascending=True)
    print(meta)

    sns.set(font_scale=1.2)
    colors = sns.color_palette('tab10')
    palette = colors[:4] + [colors[4]]*2 + [colors[5]]*2 + [colors[6]]*2
    #palette *= 2
    palette2 = []
    for color in palette:
        palette2 += [color]*2
    palette += palette2

    fig,axes = plt.subplots(len(geneList),1,figsize=(12,3*len(geneList)))
    for n,gene in enumerate(list(geneList)):
        #sc.pl.violin(adata,gene,groupby='Treatment',ax=axes[n],stripplot=False,rotation=40,order=treatment_order,palette=palette)
        #axes[n].set_xticklabels(list(map(lambda x: x.split('_')[0],treatment_order)),rotation=40,ha='right')
        sc.pl.violin(adata,gene,groupby='sample',ax=axes[n],stripplot=False,rotation=40,order=list(meta.index),palette=palette)
        axes[n].set_xticklabels(list(map(lambda x: x.split('_')[0],meta['Treatment'])),rotation=40,ha='right')
        for xtick,color in zip(axes[n].get_xticklabels(),['Blue']*10+['Red']*20):
            xtick.set_color(color)
    plt.tight_layout()
    plt.savefig(prefix+'.violin.geneBySample.%s.pdf' % cluster)
    plt.close()
    

def plotViolin2(adata,prefix,geneList,cluster,meta):
    h5sample = 'sample'
    h5cluster = 'leiden'
    new_cluster_names = {'0':'GC','1':'GC','2':'GC','3':'GC','4':'GC','5':'Empty','6':'TC','7':'GC','8':'Unannotated','9':'Stromal','10':'GC','11':'Empty','12':'GC','13':'Empty','14':'TC','15':'EpiC','16':'BEC','17':'Immune','18':'RBC','19':'LEC'}
    adata.obs['new_cluster'] = adata.obs[h5cluster].map(new_cluster_names)
    h5cluster = 'new_cluster'
    treatments = ['vehicle','FSH','FSH+hCG','4221:6','5690:0.3','5690:1','0353:10','0353:30','2599:0.5','2599:6']
    treatment_order = [i+'_T28' for i in treatments] + [i+'_T57' for i in treatments]

    adata = adata[adata.obs[h5cluster]==cluster,:].copy()
    adata.obs[h5cluster] = adata.obs[h5cluster].astype('category')
    meta['rank'] = meta.apply(lambda x: treatments.index(x['Treatment'].split('_')[0]) \
                              if x['Time']=='T28' else \
                              treatments.index(x['Treatment'].split('_')[0])+10,axis=1)
    meta = meta.sort_values(['rank','sample'],ascending=True)
    print(meta)

    sns.set(font_scale=1.2)
    colors = sns.color_palette('tab10')
    palettes = [[],[]]
    palettes[0] = colors[:4] + [colors[4]]*2 + [colors[5]]*2 + [colors[6]]*2
    for color in palettes[0]:
        palettes[1] += [color]*2

    for i,time in enumerate(['T28','T57']):
        adata1 = adata[adata.obs['Time']==time,:].copy()
        meta1 = meta[meta['Time']==time].copy()
        fig,axes = plt.subplots(len(geneList),1,figsize=((i+1)*5,3*len(geneList)))
        for n,gene in enumerate(list(geneList)):
            sc.pl.violin(adata1,gene,groupby='sample',ax=axes[n],stripplot=False,rotation=40,order=list(meta1.index),palette=palettes[i])
            xlabels = [i.split('_')[0] for i in meta1['Treatment']]
            xlabels = [i.split(':')[0]+' (%.1f)' % float(i.split(':')[1]) if not re.search(r':',i)==None else i for i in xlabels]
            axes[n].set_xticklabels(xlabels,rotation=40,ha='right')
        plt.tight_layout()
        plt.savefig(prefix+'.violin.geneBySample2.%s.%s.pdf' % (cluster,time))
        plt.close()
    

def plotClustDist(adata,prefix,meta):
    new_cluster_names = {'0':'Proliferating GC','1':'Antral mural GC','2':'Luteinizing GC','3':'Atretic GC','4':'Preantral GC','5':'Empty','6':'Theca','7':'Luteinizing GC','8':'Unannotated','9':'Stroma','10':'Luteinizing GC','11':'Empty','12':'Proliferating GC','13':'Empty','14':'Luteinizing TC?','15':'EpiC','16':'BEC','17':'Immune','18':'RBC','19':'LEC'}
    treatments = ['vehicle','FSH','FSH+hCG','4221:6','0353:10','0353:30','2599:0.5','2599:6','5690:0.3','5690:1']
    treatment_order = [i+'_T28' for i in treatments] + [i+'_T57' for i in treatments]
    adata.obs['ann_cluster'] = adata.obs['leiden'].map(new_cluster_names)
    h5sample = 'sample'
    h5cluster = 'ann_cluster'

    exclude_groups = ['Empty']
    adata = adata[list(map(lambda x: not x in exclude_groups,adata.obs[h5cluster])),:].copy()
    adata.obs[h5cluster] = adata.obs[h5cluster].astype('category')
    clusterList = sorted(adata.obs[h5cluster].unique())
    totalCells = adata.obs[h5sample].value_counts()
    print(sum(totalCells))

    colors = sns.color_palette('tab10')
    hue_colors = {}
    for i,drug in enumerate(['vehicle','FSH','FSH+hCG','4221','0353','2599','5690']):
        hue_colors.update({drug:colors[i]})
        
    sns.set(font_scale=.9)
    fig,axes = plt.subplots(len(clusterList),1,figsize=(8,3*len(clusterList)))
    for n,cluster in enumerate(clusterList):
        cells = adata.obs[adata.obs[h5cluster]==cluster][h5sample].value_counts()
        df = pd.concat([cells,totalCells,meta],join='outer',axis=1)
        df['frac_cells'] = df.apply(lambda x: x[0]/float(x[1]),axis=1)
        df['xpos'] = df.apply(lambda x: treatment_order.index(x['Treatment']),axis=1)
        df1 = df.groupby('Treatment').agg(np.mean)
        df1['Drug'] = df1.apply(lambda x: x.name.split('_')[0] if not re.search(r':',x.name) else x.name.split(':')[0],axis=1)
        '''
        b = sns.barplot(x=df1.index,y='frac_cells',hue='Drug',order=treatment_order,data=df1,dodge=False,lw=3,palette=hue_colors,ax=axes[n])
        for patch in b.patches:
            color = patch.get_facecolor()
            patch.set_edgecolor(color)
            patch.set_facecolor('white')
            patch.set_width(.65)
        '''
        b = sns.boxplot(x=df1.index,y='frac_cells',hue='Drug',order=treatment_order,data=df1,dodge=False,palette=hue_colors,ax=axes[n])
        b.legend_.remove()
        #plt.setp(b.get_legend().get_texts(),fontsize='8')
        sns.swarmplot(x='xpos',y='frac_cells',hue='Drug',data=df,legend=False,palette=hue_colors,ax=axes[n])
        axes[n].set_xticklabels(treatments*2,rotation=40,ha='right')
        for xtick,color in zip(axes[n].get_xticklabels(),['Blue']*10+['Red']*10):
            xtick.set_color(color)
        axes[n].set_facecolor('white')
        axes[n].set_title(cluster)        
    plt.tight_layout()
    plt.savefig(prefix+'.Bar.ClusterDistBySample3.pdf')
    plt.close()


def plotCellComposition(adata,prefix,meta,celltype):
    celltypes = {'0':'GC','1':'GC','2':'GC','3':'GC','4':'GC','5':'Empty','6':'TC','7':'GC','8':'Unkown','9':'Stromal','10':'GC','11':'Empty','12':'GC','13':'Empty','14':'TC','15':'EpiC','16':'BEC','17':'Immune','18':'RBC','19':'LEC'}
    adata.obs['cell_type'] = adata.obs['leiden'].map(celltypes)
    new_cluster_names = {'0':'Proliferating GC','1':'Antral mural GC','2':'Luteinizing GC','3':'Atretic GC','4':'Preantral GC','5':'Empty','6':'Theca','7':'Luteinizing GC','8':'Luteinizing GC 2','9':'Stroma','10':'Luteinizing GC','11':'Empty','12':'Proliferating GC','13':'Empty','14':'Luteinizing GC','15':'EpiC','16':'BEC','17':'Immune','18':'RBC','19':'LEC'}
    adata.obs['ann_cluster'] = adata.obs['leiden'].map(new_cluster_names)
    adata = adata[adata.obs['cell_type']==celltype,:].copy()
    removeSampleList = ['Ov12','Ov18','Ov52','Ov76','Ov92']
    adata = adata[list(map(lambda x: not x in removeSampleList,adata.obs['sample'])),:].copy()

    treatments = ['vehicle','FSH','FSH+hCG','4221:6','5690:0.3','5690:1','0353:10','0353:30','2599:0.5','2599:6']
    treatment_order = [i+'_T28' for i in treatments] + [i+'_T57' for i in treatments]
    tab10 = sns.color_palette('tab10')
    colors = {'Proliferating GC':tab10[0],'Antral mural GC':tab10[3],'Atretic GC':tab10[1],'Preantral GC':tab10[2],'Luteinizing GC':tab10[4]}
    abbrs = ['Atr','Pra','Pro','Ant','Lut']
    
    sns.set(context='talk',font_scale=.85)
    #sns.set(context='talk',font_scale=0.6)
    for time in ['T28','T57']:
        array = []
        adata1 = adata[adata.obs['Time']==time,:].copy()
        totalCells = adata1.obs['sample'].value_counts()
        for cluster in adata1.obs['ann_cluster'].unique():
            cells = adata1.obs[adata1.obs['ann_cluster']==cluster]['sample'].value_counts()
            df1 = pd.concat([cells,totalCells],join='outer',axis=1).fillna(0)
            df1[cluster] = df1.apply(lambda x: x[0]/float(x[1]),axis=1)
            array.append(df1[cluster])
            
        df = pd.concat(array,join='outer',axis=1).fillna(0)
        df = pd.concat([df,meta],join='inner',axis=1)
        df['xpos'] = df.apply(lambda x: treatment_order.index(x['Treatment']),axis=1)
        df['sample'] = df.index
        df = df.sort_values(['xpos','sample'],ascending=True)
        df = df[['Atretic GC','Preantral GC','Proliferating GC','Antral mural GC','Treatment']]
        df.iloc[:,:-1] = df.iloc[:,:-1].apply(lambda x: x*100,axis=0)
        
        #colors = ['orange','green','blue','red','purple']
        #df.iloc[:,:-1].transpose().plot(kind='pie',autopct=None,labeldistance=None,colors=colors,subplots=True,ylabel='',normalize=True,startangle=90,counterclock=False,figsize=(10,10 if time=='T57' else 6),title=list(map(lambda x: x.split('_')[0],df['Treatment'])),layout=(5 if time=='T57' else 3,4),legend=False)
        #fig,ax = plt.subplots(figsize=(16 if time=='T57' else 12,6))
        #df.iloc[:,:-1].plot(kind='bar',stacked=True,xlabel='Treatment',ylabel='Fraction of cells',edgecolor='white',color=colors,linewidth=1.5,width=.85,ax=ax)
        fig,ax = plt.subplots(figsize=(30 if time=='T57' else 20,6))
        df.iloc[:,:-1].plot(kind='bar',stacked=False,xlabel='Treatment',ylabel='Percentage of cells',edgecolor='white',color=colors,linewidth=1,width=.88,ax=ax)
        for i,container in enumerate(ax.containers):
            labels = list(map(lambda x: abbrs[i]+'\n%.0f' % x,df.iloc[:,i]))
            ax.bar_label(container,labels,fontsize=12)
            #ax.bar_label(container,fmt='%.0f',fontsize=10)
        for bkpt in range(df.shape[0]-1):
            ax.axvline(bkpt + 0.5,lw=.5,color='lightgrey')

        #ax.set_xticklabels(list(map(lambda x: x.split('_')[0],df['Treatment'])),rotation=40,ha='right')
        #xlabels = [i.split('_')[0].replace('1','1.0') if not re.search(r'5690:1',i)==None else i.split('_')[0] for i in df['Treatment']]
        xlabels = [i.split('_')[0] for i in df['Treatment']]
        xlabels = [i.split(':')[0]+' (%.1f)' % float(i.split(':')[1]) if not re.search(r':',i)==None else i for i in xlabels]
        ax.set_xticklabels(xlabels,rotation=0)
        ax.set_facecolor('white')
        #ax.set_title(celltype+' sub population composition')
        ax.set_title('')
        ax.legend(bbox_to_anchor=(1.01,1),loc='upper left')
        plt.tight_layout()
        plt.savefig(prefix+'.CompositionBySample.%s.pdf' % time)
        plt.close()

        
def plotCellComposition2(adata,prefix,meta,celltype):
    celltypes = {'0':'GC','1':'GC','2':'GC','3':'GC','4':'GC','5':'Empty','6':'TC','7':'GC','8':'Unkown','9':'Stromal','10':'GC','11':'Empty','12':'GC','13':'Empty','14':'TC','15':'EpiC','16':'BEC','17':'Immune','18':'RBC','19':'LEC'}
    adata.obs['cell_type'] = adata.obs['leiden'].map(celltypes)
    new_cluster_names = {'0':'Proliferating GC','1':'Antral mural GC','2':'Luteinizing GC','3':'Atretic GC','4':'Preantral GC','5':'Empty','6':'Theca','7':'Luteinizing GC','8':'Luteinizing GC 2','9':'Stroma','10':'Luteinizing GC','11':'Empty','12':'Proliferating GC','13':'Empty','14':'Luteinizing GC','15':'EpiC','16':'BEC','17':'Immune','18':'RBC','19':'LEC'}
    adata.obs['ann_cluster'] = adata.obs['leiden'].map(new_cluster_names)
    adata = adata[adata.obs['cell_type']==celltype,:].copy()
    SampleList = ['Ov1','Ov41','Ov12','Ov18','Ov52','Ov76','Ov92']
    adata = adata[list(map(lambda x: x in SampleList,adata.obs['sample'])),:].copy()

    treatments = ['vehicle','FSH','FSH+hCG','4221:6','5690:0.3','5690:1','0353:10','0353:30','2599:0.5','2599:6']
    treatment_order = [i+'_T28' for i in treatments] + [i+'_T57' for i in treatments]
    tab10 = sns.color_palette('tab10')
    colors = {'Proliferating GC':tab10[0],'Antral mural GC':tab10[3],'Atretic GC':tab10[1],'Preantral GC':tab10[2],'Luteinizing GC':tab10[4]}
    #colors = {'Proliferating GC':'blue','Antral mural GC':'red','Atretic GC':'orange','Preantral GC':'green','Luteinizing GC':'purple'}
    abbrs = ['Atr','Pra','Pro','Ant','Lut']
    
    sns.set(context='talk',font_scale=.8)
    #sns.set(context='talk',font_scale=0.6)
    array = []
    totalCells = adata.obs['sample'].value_counts()
    for cluster in adata.obs['ann_cluster'].unique():
        cells = adata.obs[adata.obs['ann_cluster']==cluster]['sample'].value_counts()
        df1 = pd.concat([cells,totalCells],join='outer',axis=1).fillna(0)
        df1[cluster] = df1.apply(lambda x: x[0]/float(x[1]),axis=1)
        array.append(df1[cluster])
            
    df = pd.concat(array,join='outer',axis=1).fillna(0)
    df = pd.concat([df,meta],join='inner',axis=1)
    df['xpos'] = df.apply(lambda x: treatment_order.index(x['Treatment']),axis=1)
    df['sample'] = df.index
    df = df.sort_values(['xpos','sample'],ascending=True)
    df = df[['Atretic GC','Preantral GC','Proliferating GC','Antral mural GC',\
             'Luteinizing GC','Treatment']]
    df.iloc[:,:-1] = df.iloc[:,:-1].apply(lambda x: x*100,axis=0)
        
    #colors = ['orange','green','blue','red','purple']
    #df.iloc[:,:-1].transpose().plot(kind='pie',autopct=None,labeldistance=None,colors=colors,subplots=True,ylabel='',normalize=True,startangle=90,counterclock=False,figsize=(10,10 if time=='T57' else 6),title=list(map(lambda x: x.split('_')[0],df['Treatment'])),layout=(5 if time=='T57' else 3,4),legend=False)
    fig,ax = plt.subplots(figsize=(15,6))
    df.iloc[:,:-1].plot(kind='bar',stacked=False,xlabel='Treatment',ylabel='Percentage of cells',edgecolor='white',color=colors,linewidth=1,width=.88,ax=ax)
    for i,container in enumerate(ax.containers):
        labels = list(map(lambda x: abbrs[i]+'\n%.0f' % x,df.iloc[:,i]))
        ax.bar_label(container,labels,fontsize=10)
    for bkpt in range(df.shape[0]-1):
        ax.axvline(bkpt + 0.5,lw=.5,color='lightgrey')

    xlabels = [i.split('_')[0] for i in df['Treatment']]
    xlabels = [i.split(':')[0]+' (%.1f)' % float(i.split(':')[1]) if not re.search(r':',i)==None else i for i in xlabels]
    ax.set_xticklabels(xlabels,rotation=0)
    ax.set_facecolor('white')
    ax.set_title('')
    ax.legend(bbox_to_anchor=(1.01,1),loc='upper left')
    plt.tight_layout()
    plt.savefig(prefix+'.CompositionBySample.luteinized.pdf')
    plt.close()
    

def plotTreatmentDEgenes(adata,meta,celltype,control,treatment):
    new_cluster_names = {'0':'GC','1':'GC','2':'GC','3':'GC','4':'GC','5':'Empty','6':'TC','7':'GC','8':'unknown','9':'Stromal','10':'GC','11':'Empty','12':'GC','13':'Empty','14':'unknown','15':'EpiC','16':'BEC','17':'Immune','18':'Erythroid','19':'LEC'}
    adata.obs['new_cluster'] = adata.obs['leiden'].map(new_cluster_names)
    adata = adata[adata.obs['new_cluster']==celltype,:].copy()
    if not 'base' in adata.uns['log1p']:
        adata.uns['log1p']['base'] = None

    comparison = treatment + '_vs_' + control
    sc.tl.rank_genes_groups(adata,'Treatment',groups=[treatment],reference=control,method='wilcoxon',pts=True,key_added=celltype)

    index = list(map(lambda x: x[0],adata.uns[celltype]['names']))
    logfc = pd.Series(map(lambda x: x[0],adata.uns[celltype]['logfoldchanges']),index=index,name='logFC')
    pv = pd.Series(map(lambda x: x[0],adata.uns[celltype]['pvals_adj']),index=index,name='pvals_adj')
    markers = pd.concat([logfc,pv],join='inner',axis=1)
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
    sig = markers[(markers['-logP']>10) & ((markers['logFC']>1) | (markers['logFC']<-1))] \
        #if not re.search(r'5690',treatment) else \
        #   markers[(markers['-logP']>150) & ((markers['logFC']>1.5) | (markers['logFC']<-1.5))]
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
    #plotClustDist(adata,h5adFile.replace('.h5ad',''),meta)
    plotCellComposition(adata,h5adFile.replace('.h5ad',''),meta,'GC')
    #plotCellComposition2(adata,h5adFile.replace('.h5ad',''),meta,'GC')
    #plotViolin2(adata,h5adFile.replace('.h5ad',''),['Ccnd2','Top2a','Cdk6','Mki67','Inhbb','Lhcgr','Cyp19a1','Pappa','Nppc','Mitf','Pawr','Pik3ip1','Ctgf','Gls','Amh','Amhr2','Pcsk6','Igfbp5','Ereg','Areg','Ptgs2','Pgr','Mrap','Parm1','Prss35','Cemip','Runx2'],'GC',meta)
    #plotViolin(adata,h5adFile.replace('.h5ad',''),['Lhcgr','Cyp17a1','Cyp11a1','Star','Fdx1','Fdxr','Por','Ldlr','Scarb1','Ldb3'],'TC',meta)

    #plotTreatmentDEgenes(adata,meta,'GC','vehicle_T28','FSH_T28')
    #plotTreatmentDEgenes(adata,meta,'GC','vehicle_T28','4221:6_T28')
    #plotTreatmentDEgenes(adata,meta,'GC','vehicle_T28','2599:0.5_T28')
    #plotTreatmentDEgenes(adata,meta,'GC','vehicle_T28','2599:6_T28')
    #plotTreatmentDEgenes(adata,meta,'GC','vehicle_T28','0353:10_T28')
    #plotTreatmentDEgenes(adata,meta,'GC','vehicle_T28','0353:30_T28')
    #plotTreatmentDEgenes(adata,meta,'GC','vehicle_T28','5690:1_T28')
    #plotTreatmentDEgenes(adata,meta,'GC','vehicle_T28','5690:0.3_T28')
    

if __name__=='__main__':
    main()
