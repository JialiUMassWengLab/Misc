#! /usr/bin/python

import os
import re
import sys
import csv
import glob
import pandas as pd
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from adjustText import adjust_text

def plotV(fn):
    #F01vsF34,NAFLDvsNASH
    fc_cutoff = 0.5
    pv_cutoff = 1.3
    fc_cutoff2 = 1
    pv_cutoff2 = 3
    df = pd.read_csv(fn,header=None)
    df.columns = ['geneID','logFC','pvalue','geneName','tissue']
    df['logP'] = -np.log10(df['pvalue'])
    df['logFC'] = -df['logFC']
    df = df.loc[map(lambda x: not re.search(r'^ERCC-',x),df['geneID']),:]
    sns.set(font_scale=2,context='talk')
    df.plot('logFC','logP',kind='scatter',alpha=0.7,c=df[['logFC','logP']].apply(lambda x: 'blue' if x[1]<pv_cutoff or abs(x[0]) <= fc_cutoff else 'red',axis=1),figsize=(18,15),fontsize=30,s=120)
    sig1 = df[(df['logP']>pv_cutoff) & (abs(df['logFC'])>fc_cutoff)]
    sig = df[(df['logP']>pv_cutoff2) & (abs(df['logFC'])>fc_cutoff2)]
    print sig1.shape

    #for label,x,y in zip(sig['geneName'],sig['logFC'],sig['logP']):
    #    plt.annotate(label,xy=(x,y),xytext=(2,0),textcoords='offset points',ha='left',va='center')
    texts = []
    #for label,x,y in zip(sig['geneName'],sig['logFC'],sig['logP']):
    #    texts.append(plt.text(x,y,label,ha='center',size=18))
    #adjust_text(texts,autoalign='xy',arrowprops=dict(arrowstyle="simple, head_width=0.25, tail_width=0.05", color='r', lw=0.5, alpha=0.8))
    plt.xlabel('log2 FoldChange')
    plt.ylabel('-log10 adjusted P-value')
    #plt.xlim(-1,1)
    plt.savefig(fn.replace('.sig','.volcano.pdf'))
    plt.close()


def plotS(fn):
    #F01vsF34,NAFLDvsNASH
    fc_cutoff = 0.5
    pv_cutoff = 2
    #fc_cutoff = 1
    #pv_cutoff = 3
    cat = 'Disease'
    cat1 = 'Normal Control'
    cat2 = 'NASH'
    pv = pd.read_csv(fn,header=None,index_col=0)
    pv.columns = ['logFC','pvalue','geneName','tissue']
    pv['logP'] = -np.log10(pv['pvalue'])
    pv['logFC'] = -pv['logFC']
    pv = pv.loc[map(lambda x: not re.search(r'^ERCC-',x),pv.index),:]
    df = pd.read_csv('../NASHcomb_avg_1frag_TPM.csv',index_col=0).transpose()
    df.columns = map(lambda x: re.sub(r'\.\d+$','',x),df.columns)
    info = pd.read_csv('../sampleInfo.csv',index_col=0)
    info.index = info.index.astype(str)
    info = info[(info[cat]!='None') & (info[cat]!='') & (info[cat]!='F2')]
    #info['PathFibrosis'] = map(lambda x: 'F34' if x=='F3' or x=='F4' else 'F01',info['PathFibrosis'])
    df = pd.concat([df,info[cat]],join='inner',axis=1)
    df = df.groupby(cat).agg(np.mean)
    df = pd.concat([df.transpose(),pv],join='inner',axis=1)
    df = df[df.apply(lambda x: x[cat1]+x[cat2] > 10,axis=1)]
    df1 = df[df.apply(lambda x: x['logP'] <= pv_cutoff or abs(x['logFC']) <= fc_cutoff,axis=1)]
    df2 = df[df.apply(lambda x: x['logP'] > pv_cutoff and abs(x['logFC']) > fc_cutoff,axis=1)]

    sns.set(font_scale=2,context='talk')
    fig,ax = plt.subplots(figsize=(18,15))
    #df.plot(cat1,cat2,kind='scatter',alpha=0.7,c=df[['logFC','logP']].apply(lambda x: 'blue' if x[1]<pv_cutoff or abs(x[0]) <= fc_cutoff else 'red',axis=1),fontsize=30,s=120,loglog=True,ax=ax)
    df1.plot(cat1,cat2,kind='scatter',alpha=0.7,color='blue',fontsize=30,s=120,loglog=True,ax=ax)
    df2.plot(cat1,cat2,kind='scatter',alpha=0.7,color='red',fontsize=30,s=120,loglog=True,ax=ax)
    ax.plot([0.01,1000000],[0.01,1000000],lw=2,color='red',ls='--')
    sig = df[(df['logP']>pv_cutoff) & (abs(df['logFC'])>fc_cutoff)]
    #for label,x,y in zip(sig['geneName'],sig['logFC'],sig['logP']):
    #    plt.annotate(label,xy=(x,y),xytext=(2,0),textcoords='offset points',ha='left',va='center')
    texts = []
    #for label,x,y in zip(sig['geneName'],sig['logFC'],sig['logP']):
    #    texts.append(plt.text(x,y,label,ha='center',size=18))
    #adjust_text(texts,autoalign='xy',arrowprops=dict(arrowstyle="simple, head_width=0.25, tail_width=0.05", color='r', lw=0.5, alpha=0.8))
    plt.xlabel(cat1+' TPM')
    plt.ylabel(cat2+' TPM')
    plt.xlim(1,100000)
    plt.ylim(1,100000)
    plt.savefig(fn.replace('.sig','.scatter.pdf'))
    plt.close()


def main():
    plotV(sys.argv[1])
    #plotS(sys.argv[1])


if __name__=='__main__':
    main()
