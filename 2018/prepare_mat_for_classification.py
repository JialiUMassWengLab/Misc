#! /usr/bin/python

import re
import os
import sys
import glob
import numpy as np
import pandas as pd
sys.path.append('/home/jzhuang@ms.local/NASH')
from norm_with_ERCC import normERCCIDT

def fragFilter(thred,rep_avg):
    table1 = '../2018DecAD/1812ADcohort.congregated.RC.tsv'
    table2 = '../2018NovAD/1811ADcohort.congregated.RC.tsv'
    df1 = pd.read_csv(table1,sep='\t',index_col=0,low_memory=False)
    df2 = pd.read_csv(table2,sep='\t',index_col=0,low_memory=False)
    df = pd.concat([df1,df2],join='inner',axis=1)
    df = df.loc[map(lambda x: not re.search(r'^ERCC-',x),df.index),:]
    df = df.transpose()

    if rep_avg == 'avg':
        df.index = map(lambda x: x.split('-')[0], df.index)
        withRepList = []
        for name,group in df.groupby(df.index):
            if group.shape[0] >= 2:
                withRepList.append(name)
        df = df[map(lambda x: x in withRepList,df.index)].astype('float')
        df = df.groupby(df.index).agg(sum)
    elif rep_avg == 'r1':
        df = df.loc[map(lambda x: not re.search(r'-r1$',x)==None,df.index),:]
        df.index = map(lambda x: x.split('-')[0], df.index)

    df = df.loc[:,df.apply(lambda x: any(x>0),axis=0)]
    #df = df.loc[:,df.apply(lambda x: len([i for i in x if i > thred]) >= df.shape[0]*0.8,axis=0)]
    df = df.loc[:,df.apply(lambda x: len([i for i in x if i >= thred]) >= df.shape[0]*0.5,axis=0)]
    print df.shape
    return list(df.columns)


def getTPMmat(table1,table2,geneList,rep_avg):
    df1 = pd.read_csv(table1,sep='\t',index_col=0,low_memory=False).iloc[:,:-2]
    df2 = pd.read_csv(table2,sep='\t',index_col=0,low_memory=False)
    df = pd.concat([df1,df2],join='inner',axis=1)
    names = df['Description']
    df = df.loc[map(lambda x: not re.search(r'^ERCC-',x),df.index),:]
    df.iloc[:,:-2] = df.iloc[:,:-2].apply(lambda x: x/sum(x)*1000000,axis=0)
    df = df.transpose()

    if rep_avg == 'avg':
        df.index = map(lambda x: x.split('-')[0], df.index)
        withRepList = []
        for name,group in df.groupby(df.index):
            if group.shape[0] >= 2:
                withRepList.append(name)
        df = df[map(lambda x: x in withRepList,df.index)].astype('float')
        df = df.groupby(df.index).agg(np.median)
    elif rep_avg == 'r1':
        df = df.loc[map(lambda x: not re.search(r'-r1$',x)==None,df.index),:]
        df.index = map(lambda x: x.split('-')[0], df.index)

    df = df.loc[:,map(lambda x: x in geneList,df.columns)]
    info = pd.read_csv('sampleInfo.csv',index_col=0)
    info.index = info.index.astype('str')
    table = pd.concat([ df,info[['PatientID','Center','Disease']] ],axis=1,join='inner')
    #table = table.append(names)
    return table


def getERCCmat(table1,table2,geneList,rep_avg):
    df1 = pd.read_csv(table1,sep='\t',index_col=0,low_memory=False).iloc[:,:-2]
    df2 = pd.read_csv(table2,sep='\t',index_col=0,low_memory=False)
    df = pd.concat([df1,df2],join='inner',axis=1)
    df = df.transpose()
    df = normERCCIDT(df.iloc[:-2,:].astype(float))

    if rep_avg == 'avg':
        df.index = map(lambda x: x.split('-')[0], df.index)
        withRepList = []
        for name,group in df.groupby(df.index):
            if group.shape[0] >= 2:
                withRepList.append(name)
        df = df[map(lambda x: x in withRepList,df.index)].astype('float')
        df = df.groupby(df.index).agg(np.median)
    elif rep_avg == 'r1':
        df = df.loc[map(lambda x: not re.search(r'-r1$',x)==None,df.index),:]
        df.index = map(lambda x: x.split('-')[0], df.index)

    df = df.loc[:,map(lambda x: x in geneList,df.columns)]
    info = pd.read_csv('sampleInfo.csv',index_col=0)
    info.index = info.index.astype('str')
    table = pd.concat([ df,info[['PatientID','Center','Disease']] ],axis=1,join='inner')
    return table


def getDESeq2mat(table1,table2,geneList,rep_avg):
    df1 = pd.read_csv(table1,sep='\t',index_col=0,low_memory=False)
    df1.columns = map(lambda x: '-'.join([x.split('-')[0],x.split('-')[-1]]),df1.columns)
    df2 = pd.read_csv(table2,sep='\t',index_col=0,low_memory=False)
    df2.columns = map(lambda x: '-'.join([x.split('-')[0],x.split('-')[-2]]),df2.columns)
    df = pd.concat([df1,df2],join='inner',axis=1)
    df = df.apply(lambda x: 2**x,axis=0)
    df = df.transpose()

    if rep_avg == 'avg':
        df.index = map(lambda x: x.split('-')[0], df.index)
        withRepList = []
        for name,group in df.groupby(df.index):
            if group.shape[0] >= 2:
                withRepList.append(name)
        df = df[map(lambda x: x in withRepList,df.index)].astype('float')
        df = df.groupby(df.index).agg(np.median)
    elif rep_avg == 'r1':
        df = df.loc[map(lambda x: not re.search(r'-r1$',x)==None,df.index),:]
        df.index = map(lambda x: x.split('-')[0], df.index)

    df = df.loc[:,map(lambda x: x in geneList,df.columns)]
    info = pd.read_csv('sampleInfo.csv',index_col=0)
    info.index = info.index.astype('str')
    table = pd.concat([ df,info[['PatientID','Center','Disease']] ],axis=1,join='inner')
    return table


def main():
    thred = int(sys.argv[1])
    rep_avg = sys.argv[2]
    geneList = fragFilter(thred,rep_avg)

    table2 = '../2018DecAD/1812ADcohort.congregated.tsv'
    table1 = '../2018NovAD/1811ADcohort.congregated.tsv'
    #df = getTPMmat(table1,table2,geneList,rep_avg)
    #df.transpose().to_csv('GatesADboth_%s_%dfrag_TPM.csv' % (rep_avg,thred))
    #df = getERCCmat(table1,table2,geneList,rep_avg)
    #df.transpose().to_csv('GatesADboth_%s_%dfrag_ERCCnorm.csv' % (rep_avg,thred))
    #df = getDESeq2mat('../2018NovAD/1811ADcohort.DESeq2VST.tsv','../2018DecAD/1812ADcohort.DESeq2VST.tsv',geneList,rep_avg)
    #df.transpose().to_csv('GatesADboth_%s_%dfrag_DESeq2VST.csv' % (rep_avg,thred))
    df = getDESeq2mat('../2018NovAD/1811ADcohort.DESeq2rlog.tsv','../2018DecAD/1812ADcohort.DESeq2rlog.tsv',geneList,rep_avg)
    df.transpose().to_csv('GatesADboth_%s_%dfrag_DESeq2rlog.csv' % (rep_avg,thred))
    

if __name__=='__main__':
    main()
