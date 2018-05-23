#! /usr/bin/python

import re
import os
import sys
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import scipy.stats as st
import seaborn as sns

def getMatrix():
    df = pd.read_csv('../180305_NASHvalidation_203_avg_TPM.csv',index_col=0,low_memory=False).transpose()
    #df = df[(df['PathFibrosis']=='F0') | (df['Disease']=='Normal Control')]
    #df['Disease'] = map(lambda x: 'Control' if x=='Normal Control' else 'Indiana F0',df['Disease'])
    df = df[(df['PathFibrosis']!='F2') & (df['Disease']!='Normal Control')]
    df['Disease'] = map(lambda x: 'F01' if x=='F0' or x=='F1' else 'F34',df['PathFibrosis'])
    df.columns = map(lambda x: re.sub(r'\.\d+$','',x),df.columns)
    df.iloc[:,:-3] = df.iloc[:,:-3].astype('float')
    return df


def plotBox(df,field,pdf,pv,name):
    sArray = []
    i = 1
    for key,group in df.groupby('Disease'):
        sArray.append(pd.DataFrame({'x': np.random.normal(i,0.03,group.shape[0]),'y':group[field]}))
        i += 1

    #pv = st.ranksums(df[df['Disease']=='Control'][field],df[df['Disease']=='Indiana F0'][field])[1]
    dfScat = pd.concat(sArray,join='inner',axis=0)
    sns.set(font_scale=1.6)
    fig,ax = plt.subplots(figsize=(12,10))
    bp = df.boxplot(column=field,by='Disease',boxprops={'linewidth':1.8},medianprops={'linewidth':1.8},whiskerprops={'linewidth':1.8},capprops={'linewidth':1.8},sym='',ax=ax,return_type='dict')
    [[item.set_color('blue') for item in bp[key]['boxes']] for key in bp.keys()]
    [[item.set_color('blue') for item in bp[key]['whiskers']] for key in bp.keys()]
    [[item.set_color('red') for item in bp[key]['medians']] for key in bp.keys()]
    dfScat.plot('x','y',kind='scatter',s=60,ax=ax,color='blue',fontsize=18,title='%s\nP-value = %.4f' % (name,pv))
    plt.xlabel('Disease')
    plt.ylabel('TPM')
    plt.savefig(pdf,format='pdf')
    plt.close()


def plotGeneBox(df):
    from matplotlib.backends.backend_pdf import PdfPages

    pdf = PdfPages('DESeq2.F01vsF34.box.pdf')
    df0 = pd.read_csv('DESeq2.F01vsF34.sig',index_col=0,header=None)
    for gene,row in df0.iterrows():
        #print gene,row[2],row[3]
        if row[2] < 0.1 and gene in df.columns:
            plotBox(df,gene,pdf,row[2],row[3])
    '''
    geneList = set([])
    with open('/home/jzhuang@ms.local/GTEx/liverGenes','rU') as infile:
        for line in infile:
            geneList.add(line[:-1])

    df = df.loc[:,map(lambda x: x in geneList or x=='NA',df.loc['Description'])]
    pct = pd.read_csv('/mnt/shares/Users/jzhuang/Plasma_ctrl_fusion/170715_run2_HD_CKD_HCC_CLL_MM_tissue_pct.csv',index_col=0).transpose()
    pct.index = map(lambda x: '_'.join(x.split('_')[::-1]+['XT2']),pct.index)
    pct.loc['Description'] = map(lambda x: x+' frac',pct.columns)
    df = pd.concat([pct['liver'],df],join='inner',axis=1)
    print df
    pdf = PdfPages('liverGenes.HCCvsCTRL.pdf')
    for gene in df.columns[:-1]:
        pv = st.ranksums(df[df['disease']=='Control'][gene],df[df['disease']=='HCC'][gene])[1]
        plotBox(df,gene,pdf,pv)
    '''
    pdf.close()


def main():
    df = getMatrix()
    plotGeneBox(df)


if __name__=='__main__':
    main()
