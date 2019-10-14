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
import statsmodels.stats.multitest as multi
import statsmodels.api as sm
from sklearn.preprocessing import LabelEncoder

def getPvalue(df,names,dtype,cat):
    info = pd.read_csv('../sampleInfo.csv',index_col=0)
    info = info[info['Center']=='Indiana University']
    info['NAS'] = info.apply(lambda x: 'None' if x['PathSteatosis']=='None' or x['PathInflammation']=='None' or x['PathBallooning']=='None' \
                             else int(x['PathSteatosis']) + int(x['PathInflammation']) + int(x['PathBallooning']),axis=1)
    info.index = info.index.astype(str)
    df0 = pd.concat([df,info.iloc[:,[6,7,8,9,12]]],join='inner',axis=1)
    #demo = pd.get_dummies(df0['Gender'])
    df = df0[(df0[cat]!='None') & (df0[cat]!='F01')].iloc[:,:-5]
    df = pd.concat([df.loc[:,df.apply(lambda x: len([i for i in x if i > 5]) >= 0.8*df.shape[0],axis=0)],df0.iloc[:,-5:]],join='inner',axis=1)
    print df.shape
    if cat == 'PathFibrosis':
        df[cat] = map(lambda x: int(x[1:]),df[cat])
    #X = pd.concat([demo,df['MMSE'],df['Age']],join='inner',axis=1).astype(float)
    X = LabelEncoder().fit_transform(df[cat])
    X = sm.add_constant(X)
    #print pd.concat([df0[cat],pd.Series(X,index=df.index)],join='inner',axis=1)
    data = {'Coef':{},'R^2':{},'Pvalue':{}}
    for gene in df.columns[:-5]:
        result = sm.OLS(df[gene],X).fit()
        data['Coef'].update({gene: result.params['x1']})
        data['R^2'].update({gene: result.rsquared})
        data['Pvalue'].update({gene: result.pvalues['x1']})
    out = pd.DataFrame.from_dict(data)
    out['FDR'] = pd.Series(multi.multipletests(out['Pvalue'],method='fdr_bh')[1],index=out.index)
    out = pd.concat([out,names],join='inner',axis=1)
    out.index = map(lambda x: re.sub(r'\.\d+$','',x),out.index)
    #print out.loc['ENSG00000075618']
    out.sort_values(['FDR','Pvalue']).to_csv('%s.cor.%s.IndianaOnly.csv' % (cat,dtype))
    return out


def main():
    dtype = 'TPM'
    cat = sys.argv[1]
    df = pd.read_csv('../NASHcomb_avg_0frag_%s.csv' % dtype,index_col=0,low_memory=False).transpose()
    df = df.astype('float')
    names = pd.read_csv('../18NASHcombined.congregated.tsv',sep='\t',index_col=0,low_memory=False)['Description']
    pvs = getPvalue(df,names,dtype,cat)


if __name__=='__main__':
    main()

