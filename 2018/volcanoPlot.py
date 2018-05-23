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
    fc_cutoff = 2
    pv_cutoff = 4
    df = pd.read_csv(fn,header=None)
    df.columns = ['geneID','logFC','pvalue','geneName','tissue']
    df['logP'] = -np.log10(df['pvalue'])
    df['logFC'] = -df['logFC']
    df = df.loc[map(lambda x: not re.search(r'^ERCC-',x),df['geneID']),:]
    sns.set(font_scale=2)
    df.plot('logFC','logP',kind='scatter',alpha=0.7,c=df[['logFC','logP']].apply(lambda x: 'blue' if x[1]<pv_cutoff or abs(x[0]) <= fc_cutoff else 'red',axis=1),figsize=(18,15),fontsize=20,s=100)
    sig = df[(df['logP']>pv_cutoff) & (abs(df['logFC'])>fc_cutoff)]

    #for label,x,y in zip(sig['geneName'],sig['logFC'],sig['logP']):
    #    plt.annotate(label,xy=(x,y),xytext=(2,0),textcoords='offset points',ha='left',va='center')
    texts = []
    for label,x,y in zip(sig['geneName'],sig['logFC'],sig['logP']):
        texts.append(plt.text(x,y,label,ha='center',size=18))
    adjust_text(texts,autoalign='xy',arrowprops=dict(arrowstyle="simple, head_width=0.25, tail_width=0.05", color='r', lw=0.5, alpha=0.8))
    plt.xlabel('log2 FoldChange')
    plt.ylabel('-log10 adjusted P-value')
    plt.xlim(-6,6)
    plt.savefig(fn.replace('.sig','.volcano.png'))
    plt.close()


def main():
    plotV(sys.argv[1])


if __name__=='__main__':
    main()
