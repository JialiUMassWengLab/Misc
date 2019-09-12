#! /usr/bin/python

import re
import os
import sys
import numpy as np
import pandas as pd
import scipy.stats as st
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels.stats.multitest as multi
import statsmodels.api as sm
from sklearn.preprocessing import LabelEncoder

def getPvalue(df,names,dtype):
    info = pd.read_csv('sampleInfo.csv',index_col=0).iloc[:,[1,2,3,5,6,7]]
    info.index = info.index.astype(str)
    df0 = pd.concat([df,info],join='inner',axis=1)

    df = df0[df0['MMSE']!='None'].iloc[:,:-6]
    df = pd.concat([df.loc[:,df.apply(lambda x: len([i for i in x if i > 5]) >= 0.8*df.shape[0],axis=0)],df0.iloc[:,-6:]],join='inner',axis=1)
    demo = pd.get_dummies(df[['Gender','Center']])
    print df.shape
    X = pd.concat([demo,df['MMSE'],df['Age']],join='inner',axis=1).astype(float)
    results = df.iloc[:,:-6].apply(lambda y: sm.OLS(y,X).fit(),axis=0)
    results = df.iloc[:,:-6].apply(lambda y: sm.OLS(y,X).fit(),axis=0)
    pv1 = results.apply(lambda x: x.pvalues['MMSE'])
    rs1 = results.apply(lambda x: x.rsquared)
    coef1 = results.apply(lambda x: x.params['MMSE'])
    fdr1 = pd.Series(multi.multipletests(pv1,method='fdr_bh')[1],index=pv1.index)
    out = pd.concat([coef1,rs1,pv1,fdr1,names],join='inner',axis=1)
    out.columns = ['Coef','R^2','Pvalue','FDR','Description']
    out.index = map(lambda x: re.sub(r'\.\d+$','',x),out.index)
    out.sort_values(['FDR','Pvalue']).to_csv('MMSE_cor.%s.centerAgeCorrected.csv' % dtype)

    df = df0[df0['CDR']!='None'].iloc[:,:-6]
    df = pd.concat([df.loc[:,df.apply(lambda x: len([i for i in x if i > 5]) >= 0.8*df.shape[0],axis=0)],df0.iloc[:,-6:]],join='inner',axis=1)
    demo = pd.get_dummies(df[['Gender','Center']])
    df['CDR'] = LabelEncoder().fit_transform(df['CDR'].astype(float))
    print df.shape
    X = pd.concat([demo,df['CDR'],df['Age']],join='inner',axis=1).astype(float)
    results = df.iloc[:,:-6].apply(lambda y: sm.OLS(y,X).fit(),axis=0)
    pv2 = results.apply(lambda x: x.pvalues['CDR'])
    rs2 = results.apply(lambda x: x.rsquared)
    coef2 = results.apply(lambda x: x.params['CDR'])
    fdr2 = pd.Series(multi.multipletests(pv2,method='fdr_bh')[1],index=pv2.index)
    out = pd.concat([coef2,rs2,pv2,fdr2,names],join='inner',axis=1)
    out.columns = ['Coef','R^2','Pvalue','FDR','Description']
    out.index = map(lambda x: re.sub(r'\.\d+$','',x),out.index)
    out.sort_values(['FDR','Pvalue']).to_csv('CDR_cor.%s.centerAgeCorrected.csv' % dtype)

    df = df0[df0['APOE']!='None'].iloc[:,:-6]
    df = pd.concat([df.loc[:,df.apply(lambda x: len([i for i in x if i > 5]) >= 0.8*df.shape[0],axis=0)],df0.iloc[:,-6:]],join='inner',axis=1)
    demo = pd.get_dummies(df[['Gender','Center']])
    df['APOE'] = LabelEncoder().fit_transform(df['APOE'])
    print df.shape
    X = pd.concat([demo,df['APOE'],df['Age']],join='inner',axis=1)
    results = df.iloc[:,:-6].apply(lambda y: sm.OLS(y,X).fit(),axis=0)
    pv3 = results.apply(lambda x: x.pvalues['APOE'])
    rs3 = results.apply(lambda x: x.rsquared)
    coef3 = results.apply(lambda x: x.params['APOE'])
    fdr3 = pd.Series(multi.multipletests(pv3,method='fdr_bh')[1],index=pv3.index)
    out = pd.concat([coef3,rs3,pv3,fdr3,names],join='inner',axis=1)
    out.columns = ['Coef','R^2','Pvalue','FDR','Description']
    out.index = map(lambda x: re.sub(r'\.\d+$','',x),out.index)
    out.sort_values(['FDR','Pvalue']).to_csv('APOE_cor.%s.centerAgeCorrected.csv' % dtype)

    pvs = pd.concat([pv1,pv2,pv3],join='inner',axis=1)
    pvs.columns = ['MMSE','CDR','APOE']
    pvs = pd.concat([pvs,names],join='inner',axis=1)
    #pvs.to_csv('clinc_cor.%s.centerAgeCorrected.pvalue.csv' % dtype)
    return pvs


def main():
    dtype = sys.argv[1]
    df = pd.read_csv('../GatesADboth_r1_%s.csv' % dtype,index_col=0,low_memory=False).transpose()
    #df = df[df['Center']=='University of Kentucky']
    df = df.iloc[:,:-3].astype('float')
    names = pd.read_csv('../../2018DecAD/1812ADcohort.congregated.tsv',sep='\t',index_col=0)['Description']

    pvs = getPvalue(df,names,dtype)
    #pvs = pd.read_csv('clinc_cor.%s.centerCorrected.pvalue.csv' % dtype,index_col=0)
    #pvs = pd.read_csv('clinc_cor.%s.centerAgeCorrected.pvalue.csv' % dtype,index_col=0)
    #pvs = pd.read_csv('clinc_cor.%s.KentuckyAgeCorrected.pvalue.csv' % dtype,index_col=0)
    #plotBox(df,pvs,'CDR',dtype)
    #plotBox(df,pvs,'MMSE',dtype)
    #plotBox(df,pvs,'APOE',dtype)


if __name__=='__main__':
    main()

