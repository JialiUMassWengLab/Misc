#! /usr/bin/python

import re
import os
import sys
import glob
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns


def plotHeatmap(Hmat):
    colors = ['red','blue','purple','orange','green']
    H = pd.read_csv(Hmat,sep='\t',index_col=0).astype('float').transpose()
    H = H.apply(lambda x: x/sum(x),axis=1)
    array = []
    files = glob.glob('/home/jzhuang@ms.local/BloodFrac/*.specificGenes.csv')
    for fn in files:
        if os.path.basename(fn) == 'Bcell.specificGenes.csv' \
           or os.path.basename(fn) == 'monocyte_MF_dendritic.specificGenes.csv':
            continue
        df0 = pd.read_csv(fn,header=None,index_col=0)
        df0.columns = ['Description']
        H1 = pd.concat([H[map(lambda x: x in df0.index,H.index)],df0],join='inner',axis=1)
        H1['colors'] = colors.pop()
        g = sns.clustermap(H1.iloc[:,:-2],metric='euclidean',z_score=None,standard_scale=None,row_cluster=True,col_cluster=False,cmap="YlGnBu",vmax=1.0,vmin=0.0,square=False,yticklabels=False,figsize=(10,12))
        array.append(H1.iloc[g.dendrogram_row.reordered_ind,:])

    df = pd.concat(array,join='inner',axis=0)
    df = df[map(lambda x: not re.search(r'^IG[HLK]V',x),df['Description'])]
    df.index = df['Description']
    print df

    sns.set(font_scale=1.8)
    fig,ax = plt.subplots(figsize=(10,23))
    sns.heatmap(df.drop(['Description','colors'],axis=1),cmap="RdBu_r",vmax=1.0,vmin=-1.0,square=True)
    #ax.tick_params(axis='x',colors=df['colors'])
    for color,tick in zip(df['colors'][::-1],ax.yaxis.get_major_ticks()):
        tick.label1.set_color(color)
    plt.xlabel('Component')
    plt.ylabel('Gene')
    plt.savefig(Hmat.replace('.mat.tsv','.genes_loading.heatmap.pdf'))
    plt.close()



def main():
    plotHeatmap(sys.argv[1])


if __name__=='__main__':
    main()
