#! /usr/bin/python

import re
import os
import sys
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as st

def getWBmat():
    df = pd.read_csv('/mnt/nfs/analysis/180418_JP022NS083/180418_JP022NS083.congregated.tsv',index_col=0,sep='\t')
    df = df.loc[:,map(lambda x: not re.search(r'NS83',x)==None or x=='Description',df.columns)]
    df = df.drop('NTC-NS83-25',axis=1)
    df.index = map(lambda x: re.sub(r'\.\d+$','',x),df.index)
    names = df['Description']
    df = df.iloc[:,:-1].transpose()
    df.index = map(lambda x: x.split('-')[0],df.index)
    info = pd.read_csv('/home/jzhuang@ms.local/BloodFrac/83_samples.csv',usecols=[1,3,7],dtype=str).iloc[:-1,:]
    info.index = info['Aliquot ID']
    info = info[(info['Type']=='peripheral blood') | (info['Type']=='plasma')]
    df = pd.concat([df,info],join='inner',axis=1)
    df['condition'] = df.apply(lambda x: '-'.join([x['Patient'],x['Type']]),axis=1)
    df = df.groupby('condition').agg(np.median)
    return df

def getBMlist():
    df1 = getWBmat()
    df1['frac'] = map(lambda x: 'WB' if not re.search(r'peripheral blood',x)==None else 'Plasma',df1.index)
    df1 = df1.groupby('frac').agg(np.median).transpose()
    df2 = pd.read_csv('/mnt/shares/Users/jzhuang/blood_cell_types/BloodCellTypes.exprPabst.tsv',sep='\t',index_col=0)
    df = pd.concat([df1,df2[['BM','Description']]],join='inner',axis=1)
    df1 = df[df.apply(lambda x: (x['BM']+5)/(x['WB']+5) > 5,axis=1)]
    df1['FC1'] = df1.apply(lambda x: (x['BM']+5)/(x['WB']+5),axis=1)
    df1['FC2'] = df1.apply(lambda x: (x['Plasma']+1)/(x['WB']+1),axis=1)
    df2 = df[df.apply(lambda x: (x['WB']+5)/(x['BM']+5) > 5,axis=1)]
    df2['FC1'] = df2.apply(lambda x: (x['WB']+5)/(x['BM']+5),axis=1)
    df2['FC2'] = df2.apply(lambda x: (x['WB']+1)/(x['Plasma']+1),axis=1)
    print df1
    print df2
    return df1,df2


def exprDist(geneList,kind):
    #df = pd.read_csv('/home/jzhuang@ms.local/CoExpression/GTEx.tissues.median.tsv',sep='\t',index_col=0)
    #blood = pd.read_csv('/mnt/shares/Users/jzhuang/blood_cell_types/BloodCellTypes.expr4.tsv',sep='\t',index_col=0).iloc[:,:-2]
    #blood.index = map(lambda x: re.sub(r'\.\d+$','',x),blood.index)
    #df = pd.concat([blood,df],join='inner',axis=1)
    df = pd.read_csv('/mnt/shares2/annotations/hg38/Blueprint_cellType_medians.tsv',sep='\t',index_col=0).iloc[:-1,:]
    names = df['Description']
    df = df[df.iloc[:,:-1].apply(lambda x: any(x>25),axis=1)]
    df = df.iloc[:,:-1].apply(lambda x: x/max(x),axis=1)
    df = df[df.apply(lambda x: len([i for i in x if i > 0.2]) <= 10,axis=1)]

    sns.set(font_scale=2)
    from matplotlib.backends.backend_pdf import PdfPages
    pdf = PdfPages('BM%s.bluprint.exprMap.pdf' % kind)
    #pdf = PdfPages('BM.GTEx.exprMap2.pdf' % kind)
    df2 = df[map(lambda x: x in set(geneList.index),df.index)]
    df2 = pd.concat([df2,names],join='inner',axis=1)
    df2.index = df2['Description']
    df2 = df2.drop('Description',axis=1).transpose()
    print df2.shape

    g = sns.clustermap(df2,cmap="RdBu_r",vmax=1.0,vmin=0,metric='correlation',z_score=None,standard_scale=None,row_cluster=False,col_cluster=True,linecolor='black',linewidths=1,figsize=(24,18))
    plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(),rotation='horizontal')
    plt.setp(g.ax_heatmap.xaxis.get_majorticklabels(),rotation=45)
    plt.setp(g.ax_heatmap.xaxis.get_majorticklabels(),ha='right')
    g.fig.suptitle('BM genes')
    g.savefig(pdf,format='pdf')
    plt.close()
    pdf.close()

def exprDist2(geneList,kind):
    df1 = pd.read_csv('/home/jzhuang@ms.local/CoExpression/GTEx.tissues.median.tsv',sep='\t',index_col=0)
    blood = pd.read_csv('/mnt/shares/Users/jzhuang/blood_cell_types/BloodCellTypes.expr4.tsv',sep='\t',index_col=0).iloc[:,:-2]
    blood.index = map(lambda x: re.sub(r'\.\d+$','',x),blood.index)
    df1 = pd.concat([blood,df1],join='inner',axis=1)
    df2 = pd.read_csv('/mnt/shares2/annotations/hg38/Blueprint_cellType_medians.tsv',sep='\t',index_col=0).iloc[:-1,:]
    names = df2['Description']
    df3 = df1[df1.iloc[:,:-1].apply(lambda x: any(x>25),axis=1)]
    df3 = df3.iloc[:,:-1].apply(lambda x: x/max(x),axis=1)
    df3 = df3[df3.apply(lambda x: len([i for i in x if i > 0.2]) <= 10,axis=1)]

    sns.set(font_scale=1.4)
    from matplotlib.backends.backend_pdf import PdfPages
    pdf = PdfPages('BM%s.exprBar.pdf' % kind)
    for gene in geneList.index:
        if not gene in names or not gene in df3.index:
            continue
        fig,axes = plt.subplots(1,2,figsize=(30,12))
        if gene in df1.index:
            df1.loc[gene][:-1].sort_values(ascending=False)[:15].plot(kind='barh',color='blue',ax=axes[0],title='GTEx',fontsize=18)
            axes[0].set_xlabel('TPM')
        if gene in df2.index:
            df2.loc[gene][:-1].sort_values(ascending=False)[:15].plot(kind='barh',color='blue',ax=axes[1],title='Blueprint',fontsize=18)
            axes[1].set_xlabel('TPM')

        plt.suptitle('%s\n%.4f\n%.4f' % (names[gene],geneList.loc[gene,'FC1'],geneList.loc[gene,'FC2']),fontsize=20)
        plt.tight_layout(pad=5)
        plt.savefig(pdf,format='pdf')
        plt.close()

    pdf.close()


def plotScat(df1,df2):
    df1['color'] = ['red'] * df1.shape[0]
    df2['color'] = ['blue'] * df2.shape[0]
    df = pd.concat([df1,df2],join='inner',axis=0)
    comb = df[(df['Plasma']>5)&(df['WB']>5)][['WB','Plasma','color']]
    #comb = df[df['Plasma']>5][['WB','Plasma','color']]
    #comb.iloc[:,:-1] = comb.iloc[:,:-1].apply(lambda x: x+1,axis=0)
    pv = st.ranksums(np.log2(df1['FC2']),-np.log2(df2['FC2']))[1]

    sns.set(font_scale=2,context='talk')
    comb.plot('WB','Plasma',kind='scatter',loglog=True,c=comb['color'],s=100,fontsize=30,title='Plasma cf-RNA vs Whole blood\nP-value = %.3f' % pv,figsize=(15,12),alpha=0.8)
    plt.plot([0.1,100000],[0.1,100000],lw=3,color='red',ls='--')
    #comb1 = pd.concat([comb,df['Description']],join='inner',axis=1)
    #comb1 = comb1[(comb1['BM']>500) | (comb1['serum']>500)]
    #comb1 = comb1.loc[comb1.apply(lambda x: x['serum'] > x['BM']*2 or x['BM'] > x['serum']*10,axis=1),:]
    #pack = zip(comb1['Description'],comb1['serum'],comb1['BM'])
    #for label,x,y in pack:
    #    plt.annotate(label,xy=(x,y),xytext=(8,0),textcoords='offset points',ha='left',va='center')

    plt.xlim(1,10000)
    plt.ylim(1,10000)
    plt.xlabel('Periphral blood (TPM)')
    plt.ylabel('Plasma cf-RNA (TPM)')
    plt.savefig('BM.WBvsPlasma.scatter.pdf')
    plt.close()


def plotKDE(df1,df2):
    df1['color'] = ['red'] * df1.shape[0]
    df2['color'] = ['blue'] * df2.shape[0]
    df = pd.concat([df1,df2],join='inner',axis=0)
    comb = df[(df['Plasma']>5)&(df['WB']>5)][['WB','Plasma','color']]
    comb.iloc[:,:-1] = comb.iloc[:,:-1].apply(lambda x: np.log10(x),axis=0)
    pv = st.ranksums(np.log2(df1['FC2']),-np.log2(df2['FC2']))[1]
    comb1 = comb[comb['color']=='red']
    comb2 = comb[comb['color']=='blue']

    sns.set(font_scale=2,context='talk')
    '''
    g = sns.jointplot(x='WB',y='Plasma',data=comb1,kind='kde',color='r')
    g.plot_joint(plt.scatter,c='darkred',marker='o')
    g.x = comb2['WB']
    g.y = comb2['Plasma']
    g.plot_joint(sns.kdeplot,color='b')
    g.plot_joint(plt.scatter,c='darkblue',marker='s')
    '''
    fig,ax = plt.subplots(figsize=(15,12))
    comb.plot('WB','Plasma',kind='scatter',c=comb['color'],s=75,fontsize=30,title='Plasma cf-RNA vs Whole blood\nP-value = %.3f' % pv,alpha=0.65,ax=ax)
    plt.plot([-1,5],[-1,5],lw=3,color='red',ls='--')
    sns.kdeplot(comb1['WB'],comb1['Plasma'],cmap='Reds',shade=False,shade_lowest=False,n_levels=6,ax=ax)
    sns.kdeplot(comb2['WB'],comb2['Plasma'],cmap='Blues',shade=False,shade_lowest=False,n_levels=7,ax=ax)
    plt.xlim(0,4)
    plt.ylim(0,4)
    plt.xlabel('Periphral blood (Log10 TPM)')
    plt.ylabel('Plasma cf-RNA (Log10 TPM)')
    plt.savefig('BM.WBvsPlasma.kde.pdf')
    plt.close()


def plotViolin(df1,df2):
    df1['type'] = ['BM Enriched'] * df1.shape[0]
    df2['type'] = ['BM Depleted'] * df2.shape[0]
    pv = st.ranksums(np.log2(df1['FC2']),-np.log2(df2['FC2']))[1]
    df = pd.concat([df1,df2],join='inner',axis=0)
    df['log2FC'] = df.apply(lambda x: np.log2(x['FC2']) if x['type']=='BM Enriched' else -np.log2(x['FC2']),axis=1)
    sns.set(font_scale=2,context='talk',style='whitegrid')
    fig,ax = plt.subplots(figsize=(15,12))
    sns.violinplot(x='type',y='log2FC',data=df,inner='quartile',saturation=0.8,linewidth=0.6,width=0.6,ax=ax)
    df['x'] = df.apply(lambda x: np.random.normal(0,0.03) if x['type']=='BM Enriched' else np.random.normal(1,0.03),axis=1)
    df.plot('x','log2FC',kind='scatter',s=40,color='black',alpha=0.8,ax=ax)
    plt.title('Plasma cf-RNA vs Whole blood\nP-value = %.3f' % pv)
    plt.xlabel('Gene type')
    plt.ylabel('Log2 fold enrichment')
    plt.ylim(-10,10)
    plt.savefig('BM.WBvsPlasma.violin.pdf')
    plt.close()


def main():
    BMlist,WBlist = getBMlist()
    '''
    BMlist.to_csv('BMenrich.genes.csv')
    WBlist.to_csv('BMdeplete.genes.csv')
    #plotScat(BMlist,WBlist)
    #plotKDE(BMlist,WBlist)
    plotViolin(BMlist,WBlist)
    '''
    exprDist(BMlist,'enrich')
    exprDist(WBlist,'deplete')
    #BMlist = BMlist[BMlist.apply(lambda x: (x['Plasma']+1)/(x['WB']+1) > 3,axis=1)]
    #WBlist = WBlist[WBlist.apply(lambda x: (x['WB']+1)/(x['Plasma']+1) > 3,axis=1)]
    #exprDist2(BMlist,'enrich')
    #exprDist2(WBlist,'deplete')



if __name__=='__main__':
    main()
