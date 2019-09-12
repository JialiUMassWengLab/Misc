#! /usr/bin/python

import re
import os
import sys
import glob
import psycopg2
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from adjustText import adjust_text

def plotLines():
    df = pd.read_csv('para_select_files_DESeq2topGenesFDR005_571_r1/RC.topGenes.Disease.0frag_TPM.label.prob.csv',index_col=0,header=None)
    df.columns = ['label','prob']
    info = pd.read_csv('sampleInfo.csv',index_col=0)
    df = pd.concat([df,info[['Center','Disease']]],join='inner',axis=1)
    df = df.sort_values('prob')
    df['color'] = map(lambda x: 'red' if x=='AD' else 'blue',df['Disease'])
    df['rank'] = range(0,df.shape[0])
    
    sns.set(font_scale=2,context='talk')
    fig,ax = plt.subplots(figsize=(18,12))
    #(df['prob']-0.44).plot(kind='bar',color=df['color'],ax=ax,xticks=None)
    (df['prob']-0.44).plot(kind='bar',stacked=True,bottom=0.44,color=df['color'],ax=ax,xticks=None)
    '''
    texts = []
    for label,x,y in zip(df2.index,df2['rank'],df2['prob']):
        texts.append(plt.text(x,y,label,ha='center',size=18))
    adjust_text(texts,autoalign='xy',arrowprops=dict(arrowstyle="simple, head_width=0.25, tail_width=0.05", color='r', lw=0.5, alpha=0.8))
    '''
    ax.set_xticklabels(['']*df.shape[0])
    #plt.ylim(0,1.1)
    #plt.ylim(-0.49,0.61)
    #ax.set_yticklabels(['0','0.1','0.3','0.5','0.7','0.9'])
    #ax.set_yticklabels(['-.06','0.04','0.24','0.44','0.64','0.84','1.04'])
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    plt.ylabel('Probability of AD')
    plt.savefig('waterfall.571.RC.pdf')
    plt.close()


def main():
    plotLines()


if __name__=='__main__':
    main()
