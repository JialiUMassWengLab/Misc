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
from sklearn.decomposition import NMF
from sklearn.decomposition import non_negative_factorization

def plotComp(Wmat,cat):
    W = pd.read_csv(Wmat,sep='\t',index_col=0).astype('float')
    info1 = pd.read_csv('/home/jzhuang@ms.local/NASH/2018FebValid/sampleInfo.csv',index_col=0)
    info2 = pd.read_csv('/home/jzhuang@ms.local/NASH/2018AprRevalMerged/sampleInfoAll.csv',index_col=0)
    info3 = pd.read_csv('/home/jzhuang@ms.local/NASH/PSC_1804/sampleInfo.csv',index_col=0)
    info = pd.concat([info1,info2,info3],join='inner',axis=0)
    W.index = W.index.astype(int)
    W = pd.concat([W,info],join='inner',axis=1)
    W = W[W[cat]!='None']
    print W

    plotList = ['Normal Control','NAFLD','NASH','PSC']
    grouped = W.groupby('Disease')
    sns.set(font_scale=1.8)
    from matplotlib.backends.backend_pdf import PdfPages
    pdf = PdfPages(Wmat.replace('.mat.tsv','.%s.box.pdf' % cat))
    for col in W.columns[:-4]:
        sArray = []
        fig,ax = plt.subplots(figsize=(12,10))
        for i,treat in enumerate(plotList):
            group = grouped.get_group(treat)
            sArray.append(pd.DataFrame({'x': np.random.normal(i,0.06,group.shape[0]),'y':group[col]}))
            ax.boxplot(group[col].values,positions=[i],boxprops={'linewidth':2},medianprops={'linewidth':2},whiskerprops={'linewidth':2},capprops={'linewidth':2},sym='',widths=0.5)
        ax.set_xticks(range(len(plotList)))
        ax.set_xticklabels(plotList)
        ax.set_xlim(xmin=-0.5)

        dfScat = pd.concat(sArray,join='inner',axis=0)
        dfScat.plot('x','y',kind='scatter',s=80,ax=ax,color='blue',fontsize=20,title='comp'+str(col))
        plt.xlabel(cat)
        plt.ylabel('Coefficient')
        plt.savefig(pdf,format='pdf')
        plt.close()
        
    pdf.close()


def main():
    files = glob.glob('W_TPM_*.mat.tsv')
    for fn in files:
        plotComp(fn,sys.argv[1])



if __name__=='__main__':
    main()

