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

def getMatrix(tbl):
    df = pd.read_csv(tbl,sep='\t',index_col=0).iloc[:,:-1]
    df = df.loc[map(lambda x: not re.search(r'^ERCC-',x),df.index),:]
    df.iloc[:,:-1] = df.iloc[:,:-1].apply(lambda x: x/sum(x)*1000000,axis=0)
    df.index = map(lambda x: re.sub(r'\.\d+$','',x),df.index)
    df = df.loc[df.iloc[:,:-1].apply(lambda x: max(x) > 10,axis=1),:]
    #df = df.loc[df.iloc[:,:-1].apply(lambda x: (max(x)+1)/(min(x)+1) > 5,axis=1),:]
    df.iloc[:,:-1] = df.iloc[:,:-1].apply(lambda x: x/max(x),axis=1)
    return df

def getMatrix2(tbl):
    df = pd.read_csv(tbl,sep='\t',index_col=0)
    df = df.loc[:,map(lambda x: not re.search(r'-MM02-',x)==None,df.columns)]
    df = df.loc[map(lambda x: not re.search(r'^ERCC-',x),df.index),:]
    df = df.apply(lambda x: x/sum(x)*1000000,axis=0)
    df.index = map(lambda x: re.sub(r'\.\d+$','',x),df.index)
    df = df.loc[df.apply(lambda x: max(x) > 0,axis=1),:]
    #df = df.loc[df.iloc[:,:-1].apply(lambda x: (max(x)+1)/(min(x)+1) > 5,axis=1),:]
    df = df.apply(lambda x: x/max(x),axis=1)
    return df
    
def getCBC():
    df = pd.read_csv('/mnt/shares/Users/jzhuang/2017Jul_MM/MM02_output/mm02.cbc.csv').iloc[-18:,:]
    df['doc'] = pd.to_datetime(df['doc'])
    df['doc'] = map(lambda x: x - df['doc'].iloc[2],df['doc'])
    df.index = map(lambda x: (x / np.timedelta64(1,'D')).astype(int),df['doc'])
    return df[['neut_abs','RBC','PLT']]


def plotLine(df,buffy,cbc,cellType):
    cbcType = {'neutrophil': 'neut_abs',
               'erythroid': 'RBC',
               'megakaryocyte': 'PLT',
           }
    sampleName = 'MM02'
    geneList = list(pd.read_csv('/home/jzhuang@ms.local/BloodFrac/%s.specificGenes.csv' % cellType,index_col=0).index)
    names = df['Description']
    df = df[map(lambda x: x in geneList,df.index)].transpose()
    df = df.iloc[:-1,:].astype('float')
    df.index = range(-2,16)
    buffy = buffy.transpose()
    buffy.index = map(lambda x: int('-'.join(x.split('-')[3:]).replace('day','')),buffy.index)

    from matplotlib.backends.backend_pdf import PdfPages
    pdf = PdfPages(sampleName + '.%s.normTPM.lines.pdf' % cellType)
    for col in df.columns:
        sns.set(font_scale=1.8,context='talk')
        fig,ax = plt.subplots(figsize=(15,12))
        df[col].plot(kind='line',lw=3,style='bo-',title=names.loc[col],ax=ax,fontsize=18,label='Plasma')
        buffy[col].plot(kind='line',lw=3,style='go-',ax=ax,fontsize=18,label='Buffy Coat')
        ax.tick_params(axis='y',labelcolor='blue')
        ax.set_ylabel('Fraction to Max TPM',color='blue')
        ax.set_xlabel('Days')
        ax.set_ylim(0,1.05)
        ax.legend(loc='best')
        sns.set(font_scale=1.8,context='talk',style='dark')
        ax2 = ax.twinx()
        cbc[cbcType[cellType]].plot(kind='line',color='red',lw=3,ls='--',ax=ax2,label='CBC')
        ax2.tick_params(axis='y',labelcolor='red')
        ax2.set_ylabel('CBC',color='red')
        plt.savefig(pdf,format='pdf')
        plt.close()
    pdf.close()


def main():
    df1 = getMatrix('../MM02_output/MM02.congregated.tsv')
    df2 = getMatrix2('180628_APK136MM02BuffyLP.congregated.tsv')
    cbc = getCBC()
    plotLine(df1,df2,cbc,sys.argv[1])


if __name__=='__main__':
    main()
