#! /usr/bin/python

import re
import os
import sys
import math
import numpy as np
import pandas as pd
import scipy.stats as st
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn import preprocessing
from sklearn import decomposition
from sklearn import metrics
import statsmodels.api as sm
import statsmodels.formula.api as smf
import statsmodels.stats.weightstats as ws
import statsmodels.stats.multitest as mul
from scipy.cluster.hierarchy import linkage,fcluster

def regressOnCluster(moduleGenes,clusterDf,prefix):
    names = pd.read_csv('/home/jzhuang/annotations/GTEx/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz',sep='\t',index_col=0,skiprows=2)['Description']
    names.index = map(lambda x: re.sub(r'\.\d+$','',x),names.index)

    info = pd.read_csv('/home/jzhuang/AMP/Output/demo_clinical.csv',index_col=0)
    info = info[list(map(lambda x: x in clusterDf.index,info.index))]
    info.loc[:,'Diagnosis'] = list(map(lambda x: 'AD' if re.sub(r'\+$','',x)=='AD' else 'ctrl',info['diagnosis']))
    #info.loc[:,'Diagnosis'] = list(map(lambda x: 'AD' if x<=2 else 'MCI' if x==3 else 'NCI',info['CERAD']))
    S = pd.cut(info[info['ageDeath']!='90+']['ageDeath'].astype(float),3).astype(str)
    info.loc[:,'deathAge'] = info.apply(lambda x: '90+' if x['ageDeath']=='90+' else S[x.name],axis=1)
    #info['pmi'] = pd.cut(info['pmi'],4)
    info = info[['RIN','sex','apoeGenotype','pmi','Braak','CERAD','Diagnosis']]

    for module in ['glia14']:
        names1 = names[list(map(lambda x: x in moduleGenes[module],names.index))].copy()
        labels = clusterDf[module].value_counts().keys()
        df = pd.concat([info,clusterDf[module]],join='inner',axis=1)
        df = df[(df[module]==labels[0]) | (df[module]==labels[1])]
        #df[module] = df[module]-min(df[module])
        print(df.dtypes)

        stats = {'statistic':{},'pvalue':{}}
        for variable in df.columns[:-1]:
            df1 = df.copy().dropna(subset=[variable])
            if variable == 'apoeGenotype':
                df1.loc[:,'apoeGenotype'] = list(map(lambda x: len([i for i in x.split('/') if i=='4']),df1['apoeGenotype']))
            if df1.dtypes[variable] == object:
                tbl = pd.crosstab(df1[module],df1[variable])
                print(tbl)
                oddsratio,pv = st.fisher_exact(np.array(tbl))
                stats['statistic'].update({variable: oddsratio})
                stats['pvalue'].update({variable: pv})
            else:
                print(pd.crosstab(df1[module],df1[variable]))
                U,pv = st.mannwhitneyu(df1[df1[module]==labels[0]][variable],df1[df1[module]==labels[1]][variable])
                stats['statistic'].update({variable: U})
                stats['pvalue'].update({variable: pv})

        statsDf = pd.DataFrame.from_dict(stats)
        statsDf.loc[:,'FDR'] = mul.multipletests(statsDf['pvalue'],method='fdr_bh')[1]
        print(statsDf)
        return statsDf


def plotAssoc(statsDf):
    statsDf['info'] = statsDf.index
    statsDf['-logP'] = list(map(lambda x: -np.log10(x),statsDf['FDR']))

    sns.set(font_scale=1.6,style='white')
    fig,ax = plt.subplots(figsize=(18,8))
    t_ax = sns.barplot(x='info',y='-logP',data=statsDf,color='lightgrey',linewidth=2.5,ax=ax)
    for patch in t_ax.patches:
        patch.set_edgecolor('black')
    plt.ylabel('-log10 corrected P-value',fontweight='bold',style='italic')
    plt.xlabel('')
    ax.axhline(1.3,color='red',linewidth=2,linestyle='--')
    ax.set_xticklabels(['RIN','Sex',r'APOE$\epsilon$4 copy','PMI','Braak','CERAD','Diagnosis'],fontweight='bold')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.savefig('ROSMAP_glia14_cluster_association.barplot.pdf')
    plt.close()    


def main():
    fn_prefix = sys.argv[1]
    '''
    df = pd.read_csv(fn_prefix+'.csv',low_memory=False,index_col=0)
    df.columns = map(lambda x: re.sub(r'^X','',x),df.columns)
    df.index = map(lambda x: re.sub(r'\.\d+$','',x),df.index)
    df = df.astype('float')
    #for values > 5 s.d. above the mean, cap them at 5 s.d.
    df = df.apply(lambda x: pd.Series(np.where((x - np.mean(x))/np.std(x) > 5, np.mean(x)+5*np.std(x), x),index=df.columns),axis=1)
    df = df.drop('492_120515',axis=1)
    '''
    module_gene_dir = '/home/jzhuang/AMP/Output/WGCNA_files_p8_a0'
    clusterDf = pd.read_csv(fn_prefix+'.patientCorrCluster.csv',index_col=0)
    geneList = pd.read_csv('../FP_MSBB_Code/glia14.genes.txt',index_col=0,header=None).index
    moduleGenes = {'glia14':list(geneList)}
    statsDf = regressOnCluster(moduleGenes,clusterDf,fn_prefix)
    plotAssoc(statsDf)


if __name__=='__main__':
    main()
