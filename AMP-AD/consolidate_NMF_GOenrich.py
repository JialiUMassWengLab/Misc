#! /usr/bin/python

import os
import re
import sys
import glob
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

def NMF():
    files = glob.glob('/home/jzhuang/AMP/Output/NMF_files_p8_a0/H_limma_p8_a0.comp*.GOenrich.txt')
    array = []
    for fn in files:
        df = pd.read_csv(fn,index_col=0)
        df = df[(df['Significant']>5) & (df['Annotated']<1000)]
        array.append(df[['Term','Significant','adjustP1']].head(n=5))

    consolid = pd.concat(array,join='inner',axis=0)
    consolid.columns = ['Term','Hit','FDR']
    print(consolid)
    consolid.to_csv('/home/jzhuang/AMP/Output/NMF_files_p8_a0/H_limma_p8_a0.consolidated.GOenrich.csv')


def WGCNA():
    colors = pd.read_csv('/home/jzhuang/AMP/Output/WGCNA_files_p8_a0/publish_glia_cutoff0.25/colorCode.csv',header=None,index_col=0)

    files = glob.glob('/home/jzhuang/AMP/Output/WGCNA_files_p8_a0/publish_glia_cutoff0.25/modules.module.module*.GOenrichment.txt')
    array = []
    for fn in files:
        module = int(os.path.basename(fn).split('.')[2].replace('module',''))
        df = pd.read_csv(fn,index_col=0)
        df = df[(df['Significant']>5) & (df['Annotated']<1000)]
        #df = df[['Term','Significant','adjustP1']].head(n=10)
        df = df[['Term','Significant','adjustP1']].head(n=5)
        df['module'] = [colors.loc[module,1]] * df.shape[0]
        array.append(df)

    consolid = pd.concat(array,join='inner',axis=0).fillna(1e-25)
    consolid.columns = ['Term','Hit','FDR','Module']
    print(consolid)
    #consolid.to_csv('/home/jzhuang/AMP/Output/WGCNA_files_p8_a0/publish_glia_cutoff0.25/consolidated.GOenrich.csv')

    consolid1 = consolid.loc[map(lambda x: x in ['cyan','brown','darkgrey'],consolid['Module']),:]
    consolid1['-logP'] = list(map(lambda x: -np.log10(x),consolid1['FDR']))
    sns.set(font_scale=1.5,style='white')
    fig,ax = plt.subplots(figsize=(12,10))
    #sns.barplot(x='-logP',y='Term',data=consolid1,ax=ax)
    consolid1.index = consolid1['Term']
    consolid1['-logP'].iloc[::-1].plot(kind='barh',color=['darkgrey']*5+['brown']*5+['cyan']*5,edgecolor='black',linewidth=2,ax=ax)
    ax.axvline(1.3,color='red',linestyle='--')
    ax.axhline(4.5,linestyle='--',color='black',linewidth=2.5)
    ax.axhline(9.5,linestyle='--',color='black',linewidth=2.5)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.set_yticklabels(consolid1.index[::-1],fontweight='bold')
    ax.set_xlabel('-log10 FDR',style='italic',fontweight='bold')
    plt.tight_layout()
    plt.savefig('barplot_WGCNA_GOenrich.pdf')
    plt.close()
    

def main():
    #NMF()
    WGCNA()


if __name__=='__main__':
    main()

        
