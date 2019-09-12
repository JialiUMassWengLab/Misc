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

def plotBar():
    geneList = map(lambda x: re.sub(r'\.\d+$','',x),getGeneList())
    highList = map(lambda x: re.sub(r'\.\d+$','',x),getMatrix())
    array = []
    for tissue in ['Brain','Kidney','Liver','Lung']:
        de = pd.read_csv('DESeq2.%s.4hLPS.sig' % tissue,index_col=0,header=None)
        de = de[de[1]<0]
        de = de[map(lambda x: x in geneList and x in highList,de.index)]
        de = de.iloc[:,:-2]
        de.columns = [tissue+' logFC',tissue,]
        array.append(de[tissue])

    df = pd.concat(array,join='outer',axis=1).fillna(1)
    names = pd.read_csv('../../mouseLPSv2/180502_APK102MouseLPSPlasma.congregated.tsv',sep='\t',index_col=0)['Description']
    names.index = map(lambda x: re.sub(r'\.\d+$','',x),names.index)
    df = pd.concat([df,names],join='inner',axis=1)
    print df
    df1 = df
    df1.index = df1['Description']
    df1 = df1.iloc[:,:-1].apply(lambda x: -np.log10(x),axis=0)
    sns.set(font_scale=1.5,context='talk',style='white')
    fig,axes = plt.subplots(2,1,sharex=True,gridspec_kw={'height_ratios':[1,3]},figsize=(24,12))
    plt.subplots_adjust(hspace=0.05)
    df1.plot(kind='bar',ax=axes[0],color=['blue','red','green','black'])
    df1.plot(kind='bar',ax=axes[1],color=['blue','red','green','black'])
    axes[1].hlines(-np.log10(0.05),*axes[0].get_xlim(),linestyles='dashed',color='red')
    axes[1].set_ylim(0,30)
    axes[1].set_xlabel('Genes')
    axes[0].set_ylim(45,120)
    axes[1].spines['top'].set_visible(False)
    axes[0].spines['bottom'].set_visible(False)
    axes[0].spines['top'].set_visible(False)
    axes[1].spines['right'].set_visible(False)
    axes[0].spines['right'].set_visible(False)
    axes[1].yaxis.set_ticks_position('left')
    axes[0].yaxis.set_ticks_position('left')
    axes[0].xaxis.set_ticks_position('none')
    axes[1].xaxis.tick_bottom()
    axes[1].set_xticklabels(df1.index,rotation=45,ha='right')
    #axes[0].tick_params(labelright='off')
    d = .01 # how big to make the diagonal lines in axes coordinates
    # arguments to pass plot, just so we don't keep repeating them
    kwargs = dict(transform=axes[0].transAxes, color='k', clip_on=False)
    axes[0].plot((-d,+d), (-d*3,+d*3), **kwargs)
    kwargs.update(transform=axes[1].transAxes)  # switch to the bottom axes
    axes[1].plot((-d,+d), (1-d,1+d), **kwargs)
    axes[1].get_legend().remove()
    #plt.tight_layout()
    plt.savefig('logP.tissue.4hLPS.pdf')
    plt.close()



def main():
    plotBar()


if __name__=='__main__':
    main()
