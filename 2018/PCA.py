#! /usr/bin/python

import re
import os
import sys
import glob
import numpy as np
import pandas as pd
import scipy.stats as st
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from sklearn import preprocessing
from sklearn import decomposition
import seaborn as sns
from adjustText import adjust_text

def getMatrix():
    df = pd.read_csv('1804Revalidation.congregated.tsv',sep='\t',index_col=0,low_memory=False).iloc[:,:-2].transpose()
    df = df[map(lambda x: not re.search(r'^4275',x) and not re.search(r'^10800',x),df.index)]
    df = df.loc[:,df.apply(lambda x: len([i for i in x if i > 10]) > df.shape[0] * 0.8, axis=0)]
    #df = df.apply(lambda x: x/sum(x)*1000000,axis=1)
    return df
    

def performPCA(df):
    colors = {'A1': 'red', 'A2':'blue'}
    df['rep'] = map(lambda x: x.split('-')[1],df.index)

    df.iloc[:,:-1] = preprocessing.StandardScaler().fit_transform(df.iloc[:,:-1])
    pca = decomposition.PCA()
    comp = pd.DataFrame(pca.fit_transform(df.iloc[:,:-1]),index=df.index)
    comp.columns = map(lambda x: 'PC'+str(x+1),comp.columns)
    comp = pd.concat([comp,df['rep']],join='inner',axis=1)
    comp['color'] = map(lambda x: colors[x], comp['rep'])
    print comp

    comp['alq'] = map(lambda x: x.split('-')[0],comp.index)
    withRep = []
    for aliquot,group in comp.groupby('alq'):
        if group.shape[0]>1:
            withRep.append(aliquot)
    comp = comp[map(lambda x: x in withRep,comp['alq'])]
    #comp.plot.scatter(x=0,y=1,alpha=0.8,c=comp['color'])
    sns.set(font_scale=2)
    fig, ax = plt.subplots(figsize=(18,15))
    grouped = comp.groupby('rep')
    for key,group in grouped:
        group.plot.scatter(x=0,y=1,alpha=0.8,color=colors[key],label=key,s=180,ax=ax,fontsize=20)

    #comp = comp[(comp['PC1']>-10) & (comp['PC2']>0) & ((comp['disease']=='valid1')|(comp['disease']=='valid2'))]
    #texts = []
    #for label,x,y in zip(comp.index,comp['PC1'],comp['PC2']):
    #    texts.append(plt.text(x,y,label,ha='center',size=25))
    #adjust_text(texts,autoalign='xy',arrowprops=dict(arrowstyle="simple, head_width=0.25, tail_width=0.05", color='r', lw=0.5, alpha=0.8))
    for label,x,y in zip(comp['alq'],comp['PC1'],comp['PC2']):
        plt.annotate(label,xy=(x,y),xytext=(5,0),textcoords='offset points',ha='left',va='center')
    plt.legend(fontsize=35)
    plt.xlabel('PC1',fontsize=20)
    plt.ylabel('PC2',fontsize=20)
    plt.savefig('pca.png')


def main():
    df = getMatrix()
    print df
    performPCA(df)


if __name__=='__main__':
    main()

