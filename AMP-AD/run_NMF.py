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

def getMatrix(fn,norm=True):
    df = pd.read_csv(fn,low_memory=False,index_col=0)
    df.columns = map(lambda x: re.sub(r'^X','',x),df.columns)
    df.index = map(lambda x: re.sub(r'\.\d+$','',x),df.index)
    df = df.astype('float')

    #for values > 5 s.d. above the mean, cap them at 5 s.d.
    df = df.apply(lambda x: pd.Series(np.where((x - np.mean(x))/np.std(x) > 5, np.mean(x)+5*np.std(x), x),index=df.columns),axis=1)
    #upper quartile normalization, make genes comparable to each other
    upperQuartile = df.apply(lambda x: np.quantile(x,.75),axis=1)
    if norm:
        df = df.apply(lambda x: x/np.quantile(x,.75)*1000.0,axis=1)
    return df,upperQuartile


def applyNMF(df,alpha,p,prefix):
    df = df.transpose()
    model = NMF(n_components=p,random_state=142,solver='cd',max_iter=1000000,alpha=alpha,l1_ratio=0,init='nndsvd')
    model.fit(df)
    W = pd.DataFrame(data=model.transform(df),index=df.index)
    H = pd.DataFrame(data=model.components_,columns=df.columns)
    #w, h, n_iter = non_negative_factorization(mat,n_components=p,init='nndsvd',random_state=142,alpha=alpha,l1_ratio=0,solver='cd',max_iter=1000000,regularization='components')
    #W = pd.DataFrame(data=w,index=mat.index)
    #H = pd.DataFrame(data=h,columns=mat.columns)
    W.to_csv('W_%s_p%d_a%d.mat.tsv' % (prefix,p,alpha),sep='\t')
    H.to_csv('H_%s_p%d_a%d.mat.tsv' % (prefix,p,alpha),sep='\t')
    print(W)
    print(H)


def plotComp(Wmat,cat):
    W = pd.read_csv(Wmat,sep='\t',index_col=0).astype('float')
    W.columns = list(map(lambda x: 'comp'+str(int(x)+1),W.columns))
    info = pd.read_csv('/home/jzhuang/AMP/Output/demo_clinical.csv',index_col=0)
    info = info[list(map(lambda x: x in W.index,info.index))]
    info.loc[:,'Diagnosis'] = list(map(lambda x: re.sub(r'\+$','',x),info['diagnosis']))
    assay = pd.read_csv('/home/jzhuang/AMP/metadata/ROSMAP_assay_rnaSeq_metadata.csv',index_col=0)
    info = pd.concat([info,assay['libraryBatch']],join='inner',axis=1)
    W = pd.concat([W,info],join='inner',axis=1)
    #W = W[W[cat]!='None']

    sns.set(font_scale=1.8)
    from matplotlib.backends.backend_pdf import PdfPages
    pdf = PdfPages(Wmat.replace('.mat.tsv','.%s.box.pdf' % cat))
    for col in W.columns[:-info.shape[1]]:
        fig,ax = plt.subplots(figsize=(15,12))
        sns.boxplot(x=cat,y=col,data=W,color='lightgrey',ax=ax)
        sns.swarmplot(x=cat,y=col,data=W,dodge=True,ax=ax)
        plt.title(col + ' Coefficient')
        plt.savefig(pdf,format='pdf')
        plt.close()        
    pdf.close()


def plotEntrHist(Hmat):
    H = pd.read_csv(Hmat,sep='\t',index_col=0).astype('float').transpose()
    H.columns = map(lambda x: 'comp'+str(x+1),H.columns)
    entropy = H.apply(lambda x: st.entropy(x,base=H.shape[1]),axis=1)

    sns.set(font_scale=1.2,context='talk')
    fig,ax = plt.subplots(figsize=(15,12))
    sns.histplot(entropy,binwidth=.05,ax=ax)
    plt.savefig(Hmat.replace('.mat.tsv','.entropy.hist.pdf'))
    plt.close()
    

def plotCompLoadHist(Hmat):
    H = pd.read_csv(Hmat,sep='\t',index_col=0).astype('float').transpose()
    H.columns = map(lambda x: 'comp'+str(x+1),H.columns)
    H = H.apply(lambda x: x/sum(x),axis=1)
    pcc = pd.read_csv('/home/jzhuang/AMP/Output/ROSMAP_geneLevel_limma_normed_counts.genewisePCC.csv',low_memory=False,index_col=0)

    sns.set(font_scale=1.2,context='talk')
    from matplotlib.backends.backend_pdf import PdfPages
    pdf = PdfPages(Hmat.replace('.mat.tsv','.compLoading.hist.pdf'))
    for comp in H.columns:
        fig,axes = plt.subplots(1,3,figsize=(45,12))
        sns.histplot(H[comp],binrange=(0,1),binwidth=.05,ax=axes[0])
        sns.histplot(H.sort_values(comp,ascending=False).iloc[:500,:][comp],binwidth=.05,ax=axes[1])
        #p = 5,6
        #geneList = set(H[H[comp]>.4].index).intersection(set(pcc.index))
        #p = 7,8
        geneList = set(H[H[comp]>.25].index).intersection(set(pcc.index))
        il = np.tril_indices(len(geneList),-1)
        sns.histplot(pcc.loc[geneList,geneList].values[il],binwidth=.05,binrange=(0,1),ax=axes[2])
        plt.suptitle(comp+'\n%d' % len(geneList))
        plt.savefig(pdf,format='pdf')
        plt.close()

    fig,ax = plt.subplots(figsize=(15,12))
    geneList = set(H[H.apply(lambda x: all(x<=.25),axis=1)].index).intersection(set(pcc.index))
    il = np.tril_indices(len(geneList),-1)
    sns.histplot(pcc.loc[geneList,geneList].values[il],binwidth=.05,binrange=(0,1),ax=ax)
    plt.title('Non enriched\n%d' % len(geneList))
    plt.savefig(pdf,format='pdf')
    plt.close()
    pdf.close()


def exprDist(Hmat,db):
    H = pd.read_csv(Hmat,sep='\t',index_col=0).transpose().astype('float')
    H.columns = map(lambda x: 'comp'+str(x+1),H.columns)

    df = pd.read_csv('/home/jzhuang/annotations/GTEx/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz',sep='\t',index_col=0,skiprows=2)
    df.index = map(lambda x: re.sub(r'\.\d+$','',x),df.index)
    names = df['Description']
    if db == 'GTEx':
        df = df[df.iloc[:,1:].apply(lambda x: any(x>10),axis=1)]
        df = df.iloc[:,1:].apply(lambda x: x/max(x),axis=1)
        df = df[df.apply(lambda x: len([i for i in x if i > 0.2]) <= 10,axis=1)]
    elif db == 'blueprint':
        df = pd.read_csv('/home/jzhuang/annotations/GTEx/Blueprint_cellType_medians.tsv',sep='\t',index_col=0)
        df = df[df.iloc[:,:-1].apply(lambda x: any(x>10),axis=1)]
        df = df.iloc[:,:-1].apply(lambda x: x/max(x),axis=1)
        df = df[df.apply(lambda x: len([i for i in x if i > .2]) <= 10,axis=1)]
    elif db == 'singleCell':
        sc = pd.read_csv('/home/jzhuang/annotations/GTEx/brain_single_cell_trimmed_means.csv',index_col=0)
        df = df[['Description']].merge(sc,left_on='Description',right_index=True,how='inner')
        meta = pd.read_csv('/home/jzhuang/annotations/GTEx/metadata.csv',index_col=5)
        meta = meta[~meta.index.duplicated()]
        df = pd.concat([df.iloc[:,1:].transpose(),meta['subclass_label']],join='inner',axis=1)
        #df = df.groupby('subclass_label').agg(max).transpose()
        df = df.groupby('subclass_label').agg(np.median).transpose()
        df = pd.concat([names,df],join='inner',axis=1)
        df = df[df.iloc[:,1:].apply(lambda x: any(x>2),axis=1)]
        df = df.iloc[:,1:].apply(lambda x: x/max(x),axis=1)
        df = df[df.apply(lambda x: len([i for i in x if i > 0.2]) <= 8,axis=1)]
    elif db == 'singleCell2':
        sc = pd.read_csv('/home/jzhuang/annotations/GTEx/brain_single_cell_PFC_Seurat_means.csv',index_col=0)
        df = df[['Description']].merge(sc,left_on='Description',right_index=True,how='inner')
        #df = df[df.iloc[:,1:].apply(lambda x: any(x>5),axis=1)]
        df = df.iloc[:,1:].apply(lambda x: x - max(x),axis=1)
        df = df.apply(lambda x: np.exp(x), axis=1)
        df = df[df.apply(lambda x: len([i for i in x if i > 0.2]) <= 3,axis=1)]
    elif db == 'singleCell3':
        sc = pd.read_csv('/home/jzhuang/annotations/GTEx/brain_single_cell_PFC_means.csv',index_col=0)
        df = df[['Description']].merge(sc,left_on='Description',right_index=True,how='inner')
        df = df[df.iloc[:,1:].apply(lambda x: any(x>1),axis=1)]
        df = df.iloc[:,1:].apply(lambda x: x/max(x),axis=1)
        df = df[df.apply(lambda x: len([i for i in x if i > 0.2]) <= 3,axis=1)]

    #sns.set(font_scale=1.2,context='talk')
    from matplotlib.backends.backend_pdf import PdfPages
    pdf = PdfPages(Hmat.replace('.mat.tsv','.%s.exprMap.pdf' % db))
    for comp in H.columns:
        df1 = H
        df1['perc'] = df1.apply(lambda x: x[comp]/sum(x),axis=1)
        geneList = list(df1[df1['perc']>0.35].index) #cutoff 0.4 for p=8
        #df1 = df1.apply(lambda x: x/max(x),axis=1)
        #df1 = df1[df1.apply(lambda x: len([i for i in x if i>0.66])==1,axis=1)]
        #geneList = list(df1[df1[comp]==1].index)
        df2 = df[list(map(lambda x: x in geneList,df.index))]
        df2 = pd.concat([df2,names],join='inner',axis=1)
        df2.index = df2['Description']
        df2 = df2.drop('Description',axis=1).transpose()
        print(df2.shape)
        if df2.shape[1] < 10:
            continue
        #fig,ax = plt.subplots(figsize=(24,18))
        #sns.heatmap(df2,cmap="RdBu_r",vmax=1.0,vmin=-1.0,linecolor='black',linewidths=1)
        #plt.xticks(rotation=45,ha='right')
        #plt.title(comp)
        #plt.tight_layout()
        g = sns.clustermap(df2,cmap="RdBu_r",vmax=1.0,vmin=-1.0,metric='correlation',z_score=None,standard_scale=None,row_cluster=False,col_cluster=True,linecolor='black',linewidths=0,yticklabels=True,figsize=(24,7))
        g.ax_heatmap.hlines(range(df2.shape[0]),*g.ax_heatmap.get_xlim())
        g.ax_col_dendrogram.set_visible(False)
        plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(),rotation='horizontal')
        plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(),fontsize=25)
        plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(),fontweight='bold')
        plt.setp(g.ax_heatmap.xaxis.get_majorticklabels(),rotation=45)
        plt.setp(g.ax_heatmap.xaxis.get_majorticklabels(),ha='right')
        plt.setp(g.ax_heatmap.xaxis.get_majorticklabels(),style='italic')
        g.fig.suptitle(comp)
        g.savefig(pdf,format='pdf')
        plt.close()
    pdf.close()


def cellTypeSpecificGenes(Hmat):
    H = pd.read_csv(Hmat,sep='\t',index_col=0).transpose().astype('float')
    H.columns = map(lambda x: 'comp'+str(x+1),H.columns)
    H = H.apply(lambda x: x/sum(x),axis=1)

    df = pd.read_csv('/home/jzhuang/annotations/GTEx/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz',sep='\t',index_col=0,skiprows=2)
    df.index = map(lambda x: re.sub(r'\.\d+$','',x),df.index)
    names = df['Description']    
    sc = pd.read_csv('/home/jzhuang/annotations/GTEx/brain_single_cell_trimmed_means.csv',index_col=0)
    df = df[['Description']].merge(sc,left_on='Description',right_index=True,how='inner')
    meta = pd.read_csv('/home/jzhuang/annotations/GTEx/metadata.csv',index_col=5)
    meta = meta[~meta.index.duplicated()]
    df = pd.concat([df.iloc[:,1:].transpose(),meta['subclass_label']],join='inner',axis=1)
    df = df.groupby('subclass_label').agg(max).transpose()
    df = pd.concat([names,df],join='inner',axis=1)
    df = df[df.iloc[:,1:].apply(lambda x: any(x>5),axis=1)]
    df = df.iloc[:,1:].apply(lambda x: x/max(x),axis=1)

    colors = ['red','blue','purple','orange','green']
    genes = {'Astrocyte':['SOX9','AQP4','BMPR1B','GJA1','TGFB2'],
             'Endothelial':['VWF','CLDN5','SLC2A1','CLEC14A','ESAM'],
             'Microglia':['LAPTM5','C3','FYB','C1QA','AQP1'],
             'Neuron':['NRGN','KCNQ5','GAD1','GABRA3','NELL2'],
             'Oligodendrocyte':['PLP1','MYRF','KLK6','MOG','MBP']}
    array = []
    for celltype in ('Astrocyte','Endothelial','Microglia','Neuron','Oligodendrocyte'):
        #df1 = df[df[celltype]==1]
        #df1['2nd'] = df1.apply(lambda x: max(x[x<1]),axis=1)
        #print(df1.sort_values('2nd'))
        #H1 = pd.concat([H[list(map(lambda x: x in df1[df1['2nd']<.1].index,H.index))],names],join='inner',axis=1)
        names1 = names[list(map(lambda x: x in genes[celltype],names))]
        H1 = pd.concat([H,names1],join='inner',axis=1)
        H1['colors'] = colors.pop()
        g = sns.clustermap(H1.iloc[:,:-2],metric='euclidean',z_score=None,standard_scale=None,row_cluster=True,col_cluster=False)
        array.append(H1.iloc[g.dendrogram_row.reordered_ind,:])
        
    H2 = pd.concat(array,join='inner',axis=0)
    H2.index = H2['Description']

    sns.set(font_scale=1.2,context='talk')
    fig,ax = plt.subplots(figsize=(10,20))
    sns.heatmap(H2.drop(['Description','colors'],axis=1),cmap="RdBu_r",vmax=.8,vmin=-.8,square=True)
    for color,tick in zip(H2['colors'],ax.yaxis.get_major_ticks()):
        tick.label1.set_color(color)
    plt.xlabel('Component')
    plt.ylabel('Gene')
    plt.tight_layout()
    plt.savefig(Hmat.replace('.mat.tsv','.genes_loading.heatmap.pdf'))
    plt.close()


def geneCompCorr(df,Wmat):
    W = pd.read_csv(Wmat,sep='\t',index_col=0)
    W.columns = list(map(lambda x: 'comp'+str(int(x)+1),W.columns))
    df2 = df.transpose()
    df2 = df2.loc[list(W.index),:]
    array = []
    for gene in df2.columns:
        S = pd.Series(W.apply(lambda x: st.pearsonr(x,df2[gene])[0],axis=0),name=gene)
        array.append(S)

    corr = pd.concat(array,join='inner',axis=1).transpose()
    corr.to_csv(Wmat.replace('.mat.tsv','.compGeneCorr.csv'))

    sns.set(font_scale=1.2,context='talk')
    fig,axes = plt.subplots(2,2,figsize=(30,24))
    sns.histplot(data=corr.values.flatten(),binwidth=.05,binrange=(-1,1),kde=True,ax=axes[0][0])
    axes[0][0].set_title('Gene vs Pathway PCC distribution')
    counts = corr.apply(lambda x: len([i for i in x if i>.8]),axis=1)
    sns.histplot(data=counts,binwidth=.99,kde=False,ax=axes[0][1])
    axes[0][1].set_title('# of pathway with PCC > 0.8 for each gene')
    counts = corr.apply(lambda x: len([i for i in x if i>.75]),axis=1)
    sns.histplot(data=counts,binwidth=.99,kde=False,ax=axes[1][0])
    axes[1][0].set_title('# of pathway with PCC > 0.75 for each gene')
    counts = corr.apply(lambda x: len([i for i in x if i>.7]),axis=1)
    sns.histplot(data=counts,binwidth=.99,kde=False,ax=axes[1][1])
    axes[1][1].set_title('# of pathway with PCC > 0.7 for each gene')
    plt.tight_layout()
    plt.savefig(Wmat.replace('.mat.tsv','.compGeneCorr.dist.pdf'))
    plt.close()


def main():
    alpha = 0
    p = 8

    #df,upperQuartile = getMatrix('/home/jzhuang/AMP/Output/ROSMAP_geneLevel_limma_normed_counts.csv')
    prefix = 'limma'
    #applyNMF(df,alpha,p,prefix)

    #plotComp('W_%s_p%d_a%d.mat.tsv' % (prefix,p,alpha),'Braak')
    #plotComp('W_%s_p%d_a%d.mat.tsv' % (prefix,p,alpha),'CERAD')
    #plotComp('W_%s_p%d_a%d.mat.tsv' % (prefix,p,alpha),'diagnosis')
    #plotComp('W_%s_p%d_a%d.mat.tsv' % (prefix,p,alpha),'Diagnosis')
    #plotComp('W_%s_p%d_a%d.mat.tsv' % (prefix,p,alpha),'libraryBatch')
    #exprDist('H_%s_p%d_a%d.mat.tsv' % (prefix,p,alpha),'GTEx')
    exprDist('H_%s_p%d_a%d.mat.tsv' % (prefix,p,alpha),'singleCell3')
    #exprDist('H_%s_p%d_a%d.mat.tsv' % (prefix,p,alpha),'blueprint')
    #plotCompLoadHist('H_%s_p%d_a%d.mat.tsv' % (prefix,p,alpha))
    #cellTypeSpecificGenes('H_%s_p%d_a%d.mat.tsv' % (prefix,p,alpha))
    '''
    df,upperQuartile = getMatrix('/home/jzhuang/AMP/Output/ROSMAP_geneLevel_limma_normed_counts.csv',False)
    df = df.apply(lambda x: np.log2(x),axis=0)
    geneCompCorr(df,'W_limma_p%d_a%d.mat.tsv' % (p,alpha))
    '''

if __name__=='__main__':
    main()

