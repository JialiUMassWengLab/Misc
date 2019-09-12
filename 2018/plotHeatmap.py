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

def getMatrix():
    df = pd.read_csv('../GatesADboth_r1_0frag_DESeq2rlog.csv',index_col=0,low_memory=False).iloc[:-3,:].astype(float)
    df.index = map(lambda x: re.sub(r'\.\d+$','',x),df.index)
    #df = df.drop(['17944','17934'],axis=1)
    sig = pd.read_csv('../DESeq2_analysis/DESeq2.center.sig',index_col=0,header=None)
    sig.columns = ['logFC','FDR','Description','tissue']
    geneList = list(sig[sig['FDR']<0.1].index)
    df = df[map(lambda x: x in geneList,df.index)]
    maxval = df.apply(max,axis=1)
    df = df.apply(lambda x: x/max(x)*1000,axis=1)
    return df


def getPCCmat(df):
    array = []
    header = []
    df = df.transpose()
    print df.shape
    df = df.apply(lambda x: np.log2(x+0.5),axis=0)
    for gene in df.columns:
        pcc = df.apply(lambda x: st.pearsonr(x,df[gene])[0],axis=0)
        pcc.index = df.columns
        array.append(pcc)
        header.append(gene)

    PCC = pd.concat(array,axis=1)
    PCC.columns = header
    PCC.to_csv('pcc.matrix.csv')
    print PCC


def plotHeatmap():
    df = pd.read_csv('pcc.matrix.csv0',index_col=0)
    #H = pd.read_csv('H_LDA_p7.mat.tsv',sep='\t',index_col=0).transpose()
    H = pd.read_csv('H_rlog_p6_a0.mat.tsv',sep='\t',index_col=0).transpose()
    #H = H.apply(lambda x: x/sum(x),axis=1)
    #H = H[H.apply(lambda x: any(x>0.45),axis=1)]
    H = H.apply(lambda x: x/max(x),axis=1)
    H = H[H.apply(lambda x: len([i for i in x if i>0.4])==1,axis=1)]
    df = df[map(lambda x: x in H.index,df.index)]
    df = df.loc[:,map(lambda x: x in H.index,df.columns)]
    print df.shape

    colors = {0:'orange',1:'blue',2:'red',3:'black',4:'purple',5:'green'}
    sns.set(context='talk',font_scale=1.8,style='white')
    row_colors = H.apply(lambda x: [i for i,j in enumerate(x) if j==max(x)][0],axis=1).map(colors)
    g0 = sns.clustermap(df,method='average',metric='euclidean',row_cluster=True,col_cluster=True,xticklabels=False,yticklabels=False)
    s1 = pd.DataFrame({'0':range(df.shape[0]),'1':g0.dendrogram_row.reordered_ind})
    s2 = pd.DataFrame({'0':range(df.shape[0]),'1':g0.dendrogram_col.reordered_ind})
    row_order = list(s1.sort_values('1',ascending=True)['0'])
    col_order = list(s2.sort_values('1',ascending=True)['0'])
    mask = np.zeros_like(df)
    for i in range(mask.shape[0]):
        for j in range(mask.shape[1]):
            if col_order[j] >= row_order[i]:
                mask[i,j] = True
    #mask[np.triu_indices_from(pd.DataFrame(mask).iloc[row_order,col_order])] = True
    #print pd.DataFrame(mask).iloc[g0.dendrogram_row.reordered_ind,g0.dendrogram_col.reordered_ind]
    #g = sns.clustermap(df,method='average',metric='euclidean',row_cluster=True,col_cluster=True,mask=np.array(mask),figsize=(30,28),row_colors=row_colors,xticklabels=False,yticklabels=False)
    g = sns.clustermap(df,method='average',metric='euclidean',row_cluster=True,col_cluster=True,figsize=(30,28),col_colors=row_colors,xticklabels=False,yticklabels=False)
    plt.savefig('heatmap0.pdf')
    plt.close()


def plotLines(df0,cat):
    H0 = pd.read_csv('H_rlog_p6_a0.mat.tsv',sep='\t',index_col=0).transpose()
    #H0 = H0.apply(lambda x: x/sum(x),axis=1)
    H0 = H0.apply(lambda x: x/max(x),axis=1)
    H0 = H0[H0.apply(lambda x: len([i for i in x if i>0.5])==1,axis=1)]
    info = pd.read_csv('../sampleInfo.csv',index_col=0)
    info.index = info.index.astype(str)
    info = info[info[cat]!='None']
    if cat == 'MMSE':
        info['MMSE'] = pd.cut(info['MMSE'].astype(int),5)
    df = pd.concat([df0.transpose(),info[cat]],join='inner',axis=1)
    df = df.groupby(cat).agg(np.median)
    df = df.iloc[:,:-1].apply(lambda x: x/max(x),axis=0)
    for comp in H0.columns:
        #H = H0[H0[comp]>0.35]
        H = H0[H0[comp]==1]
        print H.shape
        if H.shape[0] < 50:
            continue
        df1 = df.loc[:,map(lambda x: x in H.index,df.columns)]
        avg = df1.apply(np.median,axis=1)
        sns.set(context='talk',font_scale=1.8,style='whitegrid')
        fig,ax = plt.subplots(figsize=(24,12))
        #df1.plot(kind='line',lw=1,color='lightgrey',xticks=None,ax=ax,legend=None)
        avg.plot(kind='line',lw=5,style=['ro-'],xticks=None,ax=ax,markersize=25,fontsize=50)
        #ax.plot(range(1,6),list(avg),'ro-',lw=4,markersize=15)
        plt.xticks(range(0,5),avg.index)
        plt.xlim(-0.25,4.25)
        plt.ylim(0.15,1.05)
        plt.ylabel('Normalized expression',fontsize=50)
        plt.xlabel(cat,fontsize=50)
        plt.box(on=None)
        plt.savefig('%s.comp%s.lines.png' % (cat,comp))
        plt.close()


def main():
    df = getMatrix()
    #getPCCmat(df)
    plotHeatmap()
    #plotLines(df,sys.argv[1])


if __name__=='__main__':
    main()
