#! /usr/bin/python

import re
import os
import sys
import numpy as np
import pandas as pd
import scipy.stats as st
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import NMF
from sklearn.decomposition import non_negative_factorization
from sklearn import preprocessing
from scipy.cluster.hierarchy import linkage,fcluster

def plotGeneHeatmap(Hmat):
    df = pd.read_csv('expr.masterTbl.csv',low_memory=False,index_col=0).astype('float')
    sn = pd.read_csv('snRNAseq.celltype.avgExpr.csv',index_col=0)
    sn = sn.loc[sn.apply(lambda x: any(x > .005),axis=1),:]
    info = pd.read_csv('metadata.csv',index_col=0,dtype={'FSH_dose':'int'})
    H = pd.read_csv(Hmat,sep='\t',index_col=0).astype('float')
    H.columns = map(lambda x: 'comp'+x,H.columns)
    H = H.apply(lambda x: x/sum(x),axis=1)
    markerGenes = ['Bmp15','Gdf9']
    markers = pd.read_csv('/data/jzhuang/Optigon_snRNAseq_RatOvaries20_June2022/scanpy.marker.genes.csv',index_col=0)
    for cluster,group in markers.groupby('cluster'):
        if not cluster in [0,4,6,11]:
            markerGenes += list(group.index[:50])

    from pybiomart import Dataset
    dataset = Dataset(name='rnorvegicus_gene_ensembl',host='http://www.ensembl.org')
    names = dataset.query(attributes=['ensembl_gene_id','external_gene_name'])
    names = names.dropna(subset=['Gene name'],axis=0)
    names.index = names['Gene stable ID']

    sns.set(context='talk')
    from matplotlib.backends.backend_pdf import PdfPages
    pdf = PdfPages(Hmat.replace('.mat.tsv','.genesHeatmap2.pdf'))
    for comp in H.columns:
        geneList = list(H[H[comp]>.375].index)
        df1 = df.loc[map(lambda x: x in geneList and x in sn.index,df.index),:].copy()
        df1 = pd.concat([df1.transpose(),info[['FSH_dose','hCG_dose','time','treatment']]],join='inner',axis=1)
        timepoint = 'T57' if comp=='comp1' else 'T28'
        df1 = df1[df1['time']==timepoint].drop('time',axis=1)
        df1.loc[:,'treatment'] = df1.apply(lambda x: re.sub(r'_T\d\d$','',x['treatment']),axis=1)
        df1 = df1.groupby('treatment').agg(np.median)
        df1 = df1.sort_values(['FSH_dose','hCG_dose'],ascending=True)
        df2 = df1.iloc[:,:-2].copy().transpose()
        #print(df2)
        df3 = df1.iloc[:,-2:].copy()
        FSH_colors = dict(zip(df3['FSH_dose'].unique(),sns.color_palette('Blues',4)))
        hCG_colors = dict(zip(df3['hCG_dose'].unique(),sns.color_palette('Greens',4)))
        df3['FSH_dose'] = df3['FSH_dose'].map(FSH_colors)
        df3['hCG_dose'] = df3['hCG_dose'].map(hCG_colors)

        g = sns.clustermap(df2,cmap="Reds",vmax=1.0,vmin=0,metric='correlation',z_score=None,standard_scale=0,row_cluster=True,col_cluster=False,cbar_kws={'ticks':[0,1]},linecolor='black',linewidths=0,xticklabels=False,yticklabels=False,col_colors=df3,colors_ratio=.01,figsize=(8,18))
        g.ax_heatmap.vlines([4,8,12],*g.ax_heatmap.get_ylim(),colors='k')
        #g.ax_row_dendrogram.set_visible(False)
        g.ax_col_colors.yaxis.set_ticks_position('left')
        plt.setp(g.ax_col_colors.yaxis.get_majorticklabels(),fontsize=12)
        plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(),rotation='horizontal')
        plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(),fontsize=8)
        plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(),style='italic')
        g.fig.suptitle(comp + '  ' + timepoint)
        plt.savefig(pdf,format='pdf')
        plt.close()

        sn1 = sn.loc[list(df2.iloc[g.dendrogram_row.reordered_ind,:].index),:].copy()
        sn1 = sn1.apply(lambda x: x-min(x),axis=1)
        sn1 = sn1.apply(lambda x: x/max(x),axis=1)
        sn1 = pd.concat([sn1,names['Gene name']],join='inner',axis=1)
        sn1.index = sn1['Gene name']
        sn1 = sn1.drop('Gene name',axis=1)
        fig,axes = plt.subplots(1,2,figsize=(12,18),gridspec_kw={'wspace':0})
        #sns.heatmap(sn1,cmap='mako',vmax=1,vmin=0,linewidths=0,cbar_kws={'ticks':[0,1]},xticklabels=True,yticklabels=False,annot=False,cbar=False,ax=axes[0])
        sns.heatmap(sn1,cmap='Purples',vmax=1,vmin=0,linewidths=0,cbar_kws={'ticks':[0,1]},xticklabels=True,yticklabels=False,annot=False,cbar=False,ax=axes[0])
        plt.setp(axes[0].yaxis.get_majorticklabels(),rotation='horizontal')
        plt.setp(axes[0].yaxis.get_majorticklabels(),fontsize=8)
        plt.setp(axes[0].yaxis.get_majorticklabels(),style='italic')
        axes[0].xaxis.set_ticks_position('top')
        plt.setp(axes[0].xaxis.get_majorticklabels(),rotation='45')
        plt.setp(axes[0].xaxis.get_majorticklabels(),ha='left')
        plt.setp(axes[0].xaxis.get_majorticklabels(),fontsize=15)
        mask = np.zeros_like(sn1)
        mask[:,:] = True
        sns.heatmap(sn1,cbar=False,mask=mask,xticklabels=False,yticklabels=False,ax=axes[1])
        axes[1].set_facecolor('white')
        axes[1].set_ylabel('')

        #sn2 = sn1.loc[sn1.apply(lambda x: len([i for i in x if i > .33]) <= 1,axis=1),:].copy()
        sn2 = sn1.loc[map(lambda x: x in markerGenes,sn1.index),:].copy()
        specificGeneIndex = [i for i,gene in enumerate(sn1.index) if gene in list(sn2.index)]
        print(specificGeneIndex)
        from adjustText import adjust_text
        rightmost = axes[0].get_xlim()[1]
        texts = []
        for label,y in zip(sn2.index,specificGeneIndex):
            texts.append(axes[1].text(0,y+.5,label,size=10,style='italic'))
            #g.ax_heatmap.annotate(label,xy=(rightmost,y+.5),xytext=(1,0),textcoords='offset points',fontsize=8,style='italic',ha='left',va='center',annotation_clip=False)
        #print(adjust_text(texts,ax=axes[1],autoalign=False,expand_text=(1.2,1.2),force_text=(0.1,0.5),ha='left',va='center',arrowprops=dict(arrowstyle='-',color='k',lw=1),only_move={'text':'xy','points':'y','objects':'x'},annotation_clip=False))
        print(adjust_text(texts,ax=axes[1],autoalign=False,expand_text=(1.05,1.05),ha='left',va='center',arrowprops=dict(arrowstyle='-',color='k',lw=1),only_move={'text':'xy','points':'y','objects':'x'},lim=1000,annotation_clip=False))
        #plt.tight_layout()
        plt.savefig(pdf,format='pdf')
        plt.close()
    pdf.close()

    
def main():
    alpha = 0
    p = int(sys.argv[1])
    prefix = 'tpm'

    plotGeneHeatmap('H_%s_p%d_a%d.mat.tsv' % (prefix,p,alpha))
        

if __name__=='__main__':
    main()

