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

def getMatrix():
    df1 = pd.read_csv('170902_JRP104LG31.congregated.tsv',sep='\t',index_col=0).iloc[:,:-2]
    df2 = pd.read_csv('180417_AI066Human.congregated.tsv',sep='\t',index_col=0).iloc[:,:-1]
    df = pd.concat([df1,df2],join='inner',axis=1)
    df = df.loc[:,map(lambda x: not re.search(r'GCSF',x)==None or x=='Description',df.columns)]
    df = df.loc[map(lambda x: not re.search(r'^ERCC-',x),df.index),:]
    df.iloc[:,:-1] = df.iloc[:,:-1].apply(lambda x: x/sum(x)*1000000,axis=0)
    df.index = map(lambda x: re.sub(r'\.\d+$','',x),df.index)
    names = df['Description']
    df = df.iloc[:,:-1].transpose()
    df.index = map(lambda x: x.split('-')[0],df.index)

    return df,names
    

def writeMergedTbl(df):
    info = pd.read_csv('sampleInfo.csv',index_col=0)
    info.index = info.index.astype('str')
    info['DrawDate'] = map(lambda x: x.split()[0],info['DrawDate'])
    info['DrawDate'] = pd.to_datetime(info['DrawDate'])
    df = pd.concat([df,info[['DrawDate','PatientID']]],join='inner',axis=1)

    array = []
    for pid,group in df.groupby('PatientID'):
        group['DrawDate'] = map(lambda x: x - group['DrawDate'][0],group['DrawDate'])
        group['DrawDate'] = map(lambda x: (x / np.timedelta64(1,'D')).astype(int),group['DrawDate'])
        df1 = group.iloc[:,:-2]
        df1 = df1.loc[:,df1.apply(lambda x: any(x>10),axis=0)]
        df1 = pd.concat([df1.apply(lambda x: x/max(x),axis=0),group.iloc[:,-2:]],join='inner',axis=1)
        array.append(df1)

    norm = pd.concat(array,join='inner',axis=0)
    norm.transpose().to_csv('GCSF.normed.table.csv')
    

def plotLine(df,names,pid):
    remove = ['ARG2','ANKRD9','TRIM58','OSBP2','RHD','RNF208','OR2W3','PKLR','KCNH2','LCN2']
    df0 = pd.read_csv('/home/jzhuang@ms.local/Adipose/BoneMarrow/neutrophil_DE_genes.csv',index_col=0) 
    df0 = df0[map(lambda x: not x in remove,df0['Description'])]
    geneList = df0.index
    df = df.loc[:,map(lambda x: x in geneList,df.columns)]

    info = pd.read_csv('sampleInfo.csv',index_col=0)
    info.index = info.index.astype('str')
    info['DrawDate'] = map(lambda x: x.split()[0],info['DrawDate'])
    info = info[info['PatientID']==int(pid)]
    info['DrawDate'] = pd.to_datetime(info['DrawDate'])
    info['DrawDate'] = map(lambda x: x - info['DrawDate'][0],info['DrawDate'])
    df = df[map(lambda x: x in info.index,df.index)]
    df = df.loc[:,df.apply(lambda x: max(x) > 10,axis=0)]
    df = df.apply(lambda x: x/max(x),axis=0)
    df = pd.concat([df,info['DrawDate']],join='inner',axis=1)
    df.index = map(lambda x: (x / np.timedelta64(1,'D')).astype(int),df['DrawDate'])
    df = df.iloc[:,:-1].astype('float')
    print df

    sns.set(font_scale=1.8)
    from matplotlib.backends.backend_pdf import PdfPages
    pdf = PdfPages(pid + '.neutrophil.lines.pdf')
    for col in df.columns:
        df[col].plot(kind='line',color='blue',lw=3,title=names.loc[col],figsize=(15,12),fontsize=18)
        plt.xlabel('days')
        plt.ylabel('Fraction of Max')
        plt.savefig(pdf,format='pdf')
        plt.close()
    pdf.close()


def plotLine2(df,names):
    remove = ['ARG2','ANKRD9','TRIM58','OSBP2','RHD','RNF208','OR2W3','PKLR','KCNH2','LCN2']
    df0 = pd.read_csv('/home/jzhuang@ms.local/Adipose/BoneMarrow/neutrophil_DE_genes.csv',index_col=0) 
    df0 = df0[map(lambda x: not x in remove,df0['Description'])]
    geneList = df0.index
    df = df.loc[:,map(lambda x: x in geneList,df.columns)]

    info = pd.read_csv('sampleInfo.csv',index_col=0)
    info.index = info.index.astype('str')
    info['DrawDate'] = map(lambda x: x.split()[0],info['DrawDate'])
    info['DrawDate'] = pd.to_datetime(info['DrawDate'])
    df = pd.concat([df,info[['DrawDate','PatientID']]],join='inner',axis=1)    

    array = []
    for pid,group in df.groupby('PatientID'):
        group['DrawDate'] = map(lambda x: x - group['DrawDate'][0],group['DrawDate'])
        group.index = map(lambda x: (x / np.timedelta64(1,'D')).astype(int),group['DrawDate'])
        group1 = group.iloc[:,:-2]
        group1 = group1.loc[:,group1.apply(lambda x: max(x) > 10,axis=0)]
        #group1 = group1.apply(lambda x: x/max(x),axis=0)
        group1 = pd.concat([group1,group['PatientID']],join='inner',axis=1)
        array.append(group1)

    df = pd.concat(array,join='inner',axis=0)
    print df

    sns.set(font_scale=1.8)
    from matplotlib.backends.backend_pdf import PdfPages
    pdf = PdfPages('neutrophil.all.lines.pdf')
    for col in df.columns[:-1]:
        df1 = df.pivot(columns='PatientID',values=col)
        df1.plot(kind='line',lw=3,title=names.loc[col],figsize=(15,12),fontsize=18)
        plt.xlim(-1,11)
        plt.xlabel('days')
        #plt.ylabel('Fraction of Max')
        plt.ylabel('TPM')
        plt.savefig(pdf,format='pdf')
        plt.close()
    pdf.close()


def main():
    df,names = getMatrix()
    #print df
    writeMergedTbl(df)
    #plotLine(df,names,sys.argv[1])
    #plotLine2(df,names)


if __name__=='__main__':
    main()
