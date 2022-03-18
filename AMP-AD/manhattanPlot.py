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
import seaborn as sns
import statsmodels.api as sm

def plotManhattan(fn):
    cs = pd.read_csv('/home/jzhuang/annotations/hg19.chrom.sizes',sep='\t',header=None,index_col=0)
    cs = cs[1]
    chromSizes = list(map(lambda x: sum([cs['chr'+str(i)] for i in range(1,x+1)]),range(23)))
    df = pd.read_csv(fn,sep='\t',low_memory=False)
    df['CHR'] = df['CHR'].astype(int)
    df = df[df['CHR'] <= 22]
    df.loc[:,'-logP'] = list(map(lambda x: 0-np.log10(x),df['P']))
    df.loc[:,'pos'] = df.apply(lambda x: x['BP'] + chromSizes[x['CHR']-1],axis=1)
    df.loc[:,'color'] = list(map(lambda x: 'darkgrey' if x % 2 == 0 else 'black',df['CHR']))

    mids = []
    for i in range(len(chromSizes)-1):
        mids.append( (chromSizes[i]+chromSizes[i+1])/2 )
    sns.set(font_scale=2,context='talk')
    fig,ax = plt.subplots(figsize=(25,12))
    df.plot(x='pos',y='-logP',color=df['color'],kind='scatter',alpha=0.8,ax=ax)
    ax.axhline(-np.log10(5e-8),color='red',linestyle='--',linewidth=2)
    ax.set_xticks(mids)
    ax.set_xticklabels(list(map(str,range(1,23,1))))
    ax.set_facecolor('white')
    plt.xlabel('Chromosome',fontsize=35,fontweight='bold')
    plt.ylabel('-log10(P-value)',fontsize=35,fontweight='bold')
    plt.tight_layout()
    plt.savefig(fn.replace('.tsv','.manhattan.png'))
    plt.close()
    

def plotQQ(fn):
    df = pd.read_csv(fn,sep='\t',low_memory=False)
    df['CHR'] = df['CHR'].astype(int)
    df = df[df['CHR'] <= 22]
    df.loc[:,'-logP'] = list(map(lambda x: 0-np.log10(x),df['P']))

    a = np.sort(df['-logP'])
    #loc = np.mean(df['P'])
    #scale = np.std(df['P'])
    #b = np.sort(list(map(lambda x: 0-np.log10(x),st.norm.rvs(loc=loc,scale=scale,size=df.shape[0]))))
    b = np.sort(list(map(lambda x: 0-np.log10(x),st.norm.cdf(st.norm.rvs(size=df.shape[0])))))
    df2 = pd.DataFrame({'obs':a,'exp':b})

    sns.set(font_scale=2.5,style='white')
    fig,ax = plt.subplots(figsize=(15,12))
    #sm.qqplot(df['-logP'],dist=st.norm,line='45',alpha=0.8,ax=ax)
    df2.plot(y='obs',x='exp',kind='scatter',alpha=0.8,color='k',ax=ax)
    ax.plot([0,15],[0,15],color='red',ls='--',lw=2)
    sns.despine()
    plt.xlim(0,8)
    plt.ylim(0,8)
    plt.xlabel('Expected[-log10(P-value)]',fontsize=25,fontweight='bold')
    plt.ylabel('Observed[-log10(P-value)]',fontsize=25,fontweight='bold')
    plt.tight_layout()
    plt.savefig(fn.replace('.tsv','.qq.png'))
    plt.close()


def main():
    plotManhattan(sys.argv[1])
    #plotQQ(sys.argv[1])


if __name__=='__main__':
    main()

