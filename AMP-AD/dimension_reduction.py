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
from sklearn import preprocessing
from sklearn import decomposition
from sklearn.metrics import silhouette_score
from sklearn.metrics import silhouette_samples

def readMatrix(fn,sum_gene_types):
    analysis_types = ('protein_coding','lincRNA','processed_transcript','retained_intron','antisense')
    df = pd.read_csv(fn,sep='\t',index_col=0,low_memory=False)
    #remove samples without metadata
    df = df.drop(['800_130701','150_120419_0_merged', '764_130520'],axis=1)
    #remove outliers
    df = df.drop(['380_120503','500_120515','629_120524','367_120502','507_120515'],axis=1)
    df.loc[:,'gene_id'] = list(map(lambda x: x.split('|')[1],df.index))
    df.loc[:,'type'] = list(map(lambda x: x.split('|')[-2],df.index))

    if sum_gene_types:
        df2 = df.groupby('type').agg(sum)
        df2.to_csv(os.path.basename(fn).replace('.txt.gz','.geneTypes.sum.csv'))
        print(df2.apply(np.median,axis=1))

    df = df[list(map(lambda x: x in analysis_types,df['type']))]
    df = df.groupby('gene_id').agg(sum)
    counts = df
    df = df.apply(lambda x: x/sum(x)*1000000,axis=0)
    df = df[df.apply(lambda x: len([i for i in x if i>=5])/float(df.shape[1]) > .75,axis=1)]
    #print(df)
    #write count matrix for analysis_types genes; no outliers
    counts = counts[list(map(lambda x: x in list(df.index),counts.index))]
    counts.to_csv(fn.replace('_all_quant_counts-matrix.txt.gz','_geneLevel_counts.csv'))
    return df


def performPCA(df,prefix,normalization='std'):
    df = df.transpose()
    sampleList = df.index
    if normalization == 'std':
        df = preprocessing.StandardScaler().fit_transform(df)
    elif normalization == 'qNorm':
        df = df.apply(lambda x: np.log2(x/np.quantile(x,.75)+1),axis=0)
    pca = decomposition.PCA(n_components=4)
    comp = pd.DataFrame(pca.fit_transform(df),index=sampleList)
    comp.columns = map(lambda x: 'PC'+str(x+1),comp.columns)

    demo = pd.read_csv('/home/jzhuang/AMP/Output/demo_clinical.csv',index_col=0)
    assay = pd.read_csv('/home/jzhuang/AMP/metadata/ROSMAP_assay_rnaSeq_metadata.csv',index_col=0)
    comp = pd.concat([comp.iloc[:,:2],demo,assay[['libraryBatch','sequencingBatch','libraryPrep']]],join='inner',axis=1)
    print(comp)
    #print(silhouette_score(comp.iloc[:,:2],comp['libraryBatch']))
    comp.loc[:,'ss'] = pd.Series(silhouette_samples(comp.iloc[:,:2],comp['libraryBatch']),index=comp.index)
    print(np.mean(comp['ss']), np.mean(comp[comp['libraryBatch']=='7']['ss']))

    sns.set(context='talk',font_scale=1.2)
    from matplotlib.backends.backend_pdf import PdfPages
    pdf = PdfPages('%s.%s.PCA.pdf' % (prefix,normalization))
    for category in ('sex','race','diagnosis','apoeGenotype','libraryPrep','libraryBatch','sequencingBatch','RIN','pmi','ageDeath','Braak','CERAD'):
        tmp = comp
        if category == 'ageDeath':
            tmp = tmp.dropna(axis=0,subset=['ageDeath'])
            tmp.loc[:,'ageDeath'] = list(map(lambda x: 90.1 if x=='90+' else float(x),tmp['ageDeath']))
        fig,ax = plt.subplots(figsize=(15,12))
        sns.scatterplot(x='PC1',y='PC2',hue=category,data=tmp,ax=ax)
        #for label,x,y in zip(low.index,low['PC1'],low['PC2']):
        #    plt.annotate(label,xy=(x,y),xytext=(5,0),textcoords='offset points',ha='left',va='center')
        plt.tight_layout()
        plt.savefig(pdf,format='pdf')
        plt.close()
    pdf.close()


def main():
    #rawDf = readMatrix('/home/jzhuang/AMP/Source/ROSMAP_all_quant_counts-matrix.txt.gz',False)
    #performPCA(rawDf,'qNorm')
    rawDf = pd.read_csv('/home/jzhuang/AMP/Output/ROSMAP_geneLevel_limma_normed_counts.csv',index_col=0)
    rawDf.columns = list(map(lambda x: re.sub(r'^X','',x),rawDf.columns))
    performPCA(rawDf,'limmaNormCounts')
    print(rawDf.shape[0])


if __name__=='__main__':
    main()
    
