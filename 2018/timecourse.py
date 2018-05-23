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
import numpy.polynomial.polynomial as poly
from check_const_genes1 import getConstGenes
from check_const_genes1 import getMatrix
from check_const_genes1 import getMatrix2
import seaborn as sns

def getVmat(chainType):
    files = glob.glob('*.IGgenes.readCounts')
    data = {'1st':{}, '2nd':{}, 'total':{}}
    for fn in files:
        with open(fn,'rU') as infile:
            total = 0
            m = 0
            n = 0
            sample = fn.replace('.IGgenes.readCounts','')
            if re.search(r'-IP0',sample):
                sample = 'MM-'+sample.split('-')[0]

            for line in infile:
                a = line.split()
                if re.search(r'^%s' % chainType,a[0]):
                    a[2] = int(a[2])
                    total += a[2]
                    if a[2] > m:
                        m = a[2]
                    elif a[2] > n:
                        n = a[2]
            data['1st'].update({sample: m})
            data['2nd'].update({sample: n})
            data['total'].update({sample: total})

    df = pd.DataFrame.from_dict(data)
    df['ratio1'] = df['1st']/df['total']
    df['ratio2'] = df['1st']/df['2nd']
    return df

def getVmat1(Vgene):
    files = glob.glob('*.IGgenes.readCounts')
    data = {'Vgene':{}, 'total':{}}
    for fn in files:
        with open(fn,'rU') as infile:
            total = 0
            m = 0
            sample = fn.replace('.IGgenes.readCounts','')
            if re.search(r'-IP0',sample):
                sample = 'MM-'+sample.split('-')[0]

            for line in infile:
                a = line.split()
                if re.search(r'^%s' % Vgene[:4],a[0]):
                    a[2] = int(a[2])
                    total += a[2]
                    if Vgene == a[0]:
                        m += a[2]
            data['Vgene'].update({sample: m})
            data['total'].update({sample: total})

    df = pd.DataFrame.from_dict(data)
    df['ratio'] = df['Vgene']/df['total']
    return df
    

def plotLineKLratio():
    IGTRgeneList = getConstGenes()
    df = pd.read_csv('MM02.congregated.tsv',sep='\t',index_col=0).iloc[:,:-1]
    df = df.loc[map(lambda x: x in IGTRgeneList,df['Description']),:]
    df.index = df['Description']
    df = df.drop('Description',axis=1).transpose()
    df = df.astype('float')
    df = df.loc[:,map(lambda x: not re.search(r'^IG',x)==None,df.columns)]
    df['ratio1'] = df.apply(lambda x: x['IGKC']/sum(x[['IGLC1','IGLC2','IGLC3','IGLC7']]),axis=1)
    df['days'] = range(-2,16)
    #return df
    #df[['ratio1','days']].plot('days','ratio1',kind='line',title='Kappa/Lambda ratio time course',legend=False,lw=2.5)
    #plt.savefig('MM02.KLratio.timecourse.png')
    #plt.close()

    sns.set(font_scale=2)
    polyF = poly.Polynomial(poly.polyfit(df['days'],df['ratio1'],8))
    intvals = map(lambda x: x/10.0,range(-20,151))
    fit = pd.DataFrame({'days': intvals, 'Kappa/Lambda ratio': polyF(intvals)})
    fig,ax = plt.subplots(figsize=(12,10))
    df.plot('days','ratio1',kind='scatter',title='Kappa/Lambda ratio time course',legend=False,ax=ax,s=100,fontsize=18)
    fit.plot('days','Kappa/Lambda ratio',kind='line',lw=4,legend=False,ax=ax,alpha=0.5,fontsize=18)
    plt.ylabel('Ratio')
    plt.xlabel('Days')
    plt.savefig('MM02.KLratio.timecourse.fitted.png')
    plt.close()

def plotLine(chainType):
    df = getVmat(chainType)
    df['days'] = range(-2,16)
    return df
    #df[['ratio1','days']].plot('days','ratio1',kind='line',title='Dominant %s gene frac time course' % chainType,legend=False,lw=2.5)
    #plt.savefig('MM02.%s.ratio1.timecourse.png' % chainType)

def plotLine1(Vgene):
    df = getVmat1(Vgene)
    df['days'] = range(-2,16)
    return df
    #df[['ratio','days']].plot('days','ratio',kind='line',title='Dominant gene %s frac time course' % Vgene,legend=False,lw=2.5)
    #plt.savefig('MM02.%s.ratio.timecourse.png' % Vgene)

def plotLineAll():
    df1 = plotLine('IGHV')
    df2 = plotLine('IGKV')
    df3 = plotLineKLratio()
    df = pd.concat([df1['ratio1'],df2['ratio1'],df3['ratio1']],join='inner',axis=1)
    df.columns = ['IGHV','IGKV','K/L ratio']
    df['days'] = range(-2,16)

    sns.set(font_scale=2)
    fig,ax = plt.subplots(figsize=(12,10))
    df[['IGHV','days']].plot('days','IGHV',kind='line',lw=2.5,color='blue',ax=ax,fontsize=18)
    df[['IGKV','days']].plot('days','IGKV',kind='line',lw=2.5,color='red',ax=ax,fontsize=18)
    #df[['K/L ratio','days']].plot('days','K/L ratio',kind='line',lw=2.5,color='black',ax=ax)
    plt.title('MM02 dominant gene fraction time course')
    plt.ylim(0,1)
    plt.ylabel('Fraction')
    plt.savefig('MM02.ratio1.timecourse.png')
    plt.close()


def plotLineAll1():
    df1 = plotLine1('IGHV3-15')
    df2 = plotLine1('IGKV2-24')
    df = pd.concat([df1['ratio'],df2['ratio']],join='inner',axis=1)
    df.columns = ['IGHV3-15','IGKV2-24']
    df['days'] = range(-2,16)

    IGHVwt = ['IGHV1-18','IGHV1-46','IGHV1-3']
    IGKVwt = ['IGKV1-5','IGKV3-15']
    dfs = []
    for gene in IGHVwt + IGKVwt:
        tmp = plotLine1(gene)
        dfs.append(tmp['ratio'])
    df3 = pd.concat(dfs,join='outer',axis=1).fillna(0)
    df3.columns = IGHVwt + IGKVwt
    df3.index = range(-2,16)

    poly1 = poly.Polynomial(poly.polyfit(df[~(df['days']==3)]['days'],df[~(df['days']==3)]['IGHV3-15'],8))
    poly2 = poly.Polynomial(poly.polyfit(df['days'],df['IGKV2-24'],10))
    intvals = map(lambda x: x/10.0,range(-20,151))
    fit = pd.DataFrame({'days': intvals, 'IGHV3-15': poly1(intvals), 'IGKV2-24': poly2(intvals)})
    #df['HVfit'] = df.apply(lambda x: poly1(x['days']), axis=1)
    #df['KVfit'] = df.apply(lambda x: poly2(x['days']), axis=1)
    sns.set(font_scale=2)
    fig,ax = plt.subplots(figsize=(12,10))
    #df[['IGHV3-15','days']].plot('days','IGHV3-15',kind='line',lw=2.5,color='blue',ax=ax,fontsize=18)
    #df[['IGKV2-24','days']].plot('days','IGKV2-24',kind='line',lw=2.5,color='red',ax=ax,fontsize=18)
    df.plot('days','IGHV3-15',kind='scatter',color='blue',ax=ax,fontsize=18,s=100)
    df.plot('days','IGKV2-24',kind='scatter',color='red',ax=ax,fontsize=18,s=100)
    fit.plot('days','IGHV3-15',kind='line',lw=4.5,color='blue',ax=ax,fontsize=18,alpha=0.5)
    fit.plot('days','IGKV2-24',kind='line',lw=4.5,color='red',ax=ax,fontsize=18,alpha=0.5)
    df3[IGHVwt].plot(kind='line',lw=2,ls='--',ax=ax,fontsize=18)
    df3[IGKVwt].plot(kind='line',lw=2,ls='-.',ax=ax,fontsize=18)
    plt.title('MM02 Immunoglobulin V gene fraction time course')
    plt.ylabel('Fraction')
    plt.xlabel('Days')
    plt.ylim(-0.1,1.0)
    plt.savefig('MM02.ratio.timecourse.fitted.png')
    plt.close()

    
def main():
    pd.options.display.max_rows = 100
    #df = getVmat('IGHV')
    #df = df.loc[map(lambda x: not re.search(r'-dep',x),df.index),:]
    #print df
    #plotLine()
    #plotLine1('IGHV1-69')
    plotLineKLratio()
    #plotLineAll1()


if __name__=='__main__':
    main()
