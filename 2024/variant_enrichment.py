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
import scipy.stats as st
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
from statsmodels.tools.tools import add_constant

def pQTLEnrich(prefix,chrom):
    import polars as pl
    from statsmodels.stats.weightstats import ztest
    with open('impactful_variant.list.txt','r') as file:
        variants = [i[:-1] for i in file.readlines()]
        
    modules = pd.read_csv(os.path.join('../VAE_protein_embedding',prefix+'_modules.txt'),sep=' ',index_col=0).astype('float')
    #df = pd.read_csv('impactful_variant.sumStats.csv.gz',low_memory=False)
    #df.columns = ['Variant','ASSAY','BETA','SE','LOG10P']
    df = pl.read_csv('tmp.%s.csv' % chrom).select(Variant=pl.col('CHROM').cast(str)+'_'+pl.col('GENPOS').cast(str)+'_'+pl.col('ALLELE0')+'_'+pl.col('ALLELE1'),ASSAY=pl.col('ASSAY'),BETA=pl.col('BETA'),SE=pl.col('SE'),LOG10P=pl.col('LOG10P')).filter(pl.col('Variant').is_in(variants)).to_pandas()
    data = {'statistic':{},'p-value':{}}
    for variant in df['Variant'].unique():
        print(variant)
        df1 = df.query('Variant == @variant')[['ASSAY','BETA']].dropna(axis=0)
        m = df1['BETA'].mean()
        s = df1['BETA'].std()
        df1['BETA'] = df1['BETA'].map(lambda x: np.where(abs(x - m)/s > 5, m + 5*np.sign(x-m)*s, x))
        df1.index = df1['ASSAY'].map(lambda x: x.split(':')[0]+'_'+x.split(':')[2])
        for module in modules.columns:
            pair = variant+'-'+module
            comb = pd.concat([df1,modules[module].astype(bool)],join='inner',axis=1)
            if comb[comb[module]].shape[0] > 0:
                diff,pv = ztest(comb[comb[module]]['BETA'],comb[~comb[module]]['BETA'],value=0)
            else:
                diff = 0
                pv = 1
            data['statistic'].update({pair: diff})
            data['p-value'].update({pair: pv})

    outDf = pd.DataFrame.from_dict(data)
    outDf['variant'] = outDf.index.map(lambda x: x.split('-')[0])
    outDf['module'] = outDf.index.map(lambda x: x.split('-')[1].replace('VAE_',''))
    outDf['FDR'] = multipletests(outDf['p-value'],method='fdr_bh')[1]
    #outDf.sort_values('FDR',ascending=True).to_csv('impactful_variant.pQTLEnrich.csv')
    outDf.sort_values('FDR',ascending=True).to_csv(chrom+'.pQTLEnrich.ztest.csv')
    print(outDf)

    
def pQTLEnrichNorm(prefix,chrom):
    import polars as pl
    from statsmodels.stats.weightstats import ztest
    with open('impactful_variant.list.txt','r') as file:
        variants = [i[:-1] for i in file.readlines()]
        
    modules = pd.read_csv(os.path.join('../VAE_protein_embedding',prefix+'_modules.txt'),sep=' ',index_col=0).astype('float')
    #df = pd.read_csv('impactful_variant.sumStats.csv.gz',low_memory=False)
    #df.columns = ['Variant','ASSAY','BETA','SE','LOG10P']
    df = pl.read_csv('tmp.%s.csv' % chrom).select(Variant=pl.col('CHROM').cast(str)+'_'+pl.col('GENPOS').cast(str)+'_'+pl.col('ALLELE0')+'_'+pl.col('ALLELE1'),ASSAY=pl.col('ASSAY'),BETA=pl.col('BETA'),SE=pl.col('SE'),LOG10P=pl.col('LOG10P')).filter(pl.col('Variant').is_in(variants)).to_pandas()
    data = {'statistic':{},'p-value':{},'mean-diff':{}}
    for variant in df['Variant'].unique():
        print(variant)
        df1 = df.query('Variant == @variant')[['ASSAY','BETA','SE']].dropna(axis=0)
        df1['Z'] = df1['BETA']/df1['SE']
        #df1['Z'] = df1['Z'].map(lambda x: np.where(abs(x) > 10, 10 * np.sign(x), x))
        m = df1['Z'].mean()
        s = df1['Z'].std()
        df1['Z'] = df1['Z'].map(lambda x: np.where(abs(x - m)/s > 5, m + 5*np.sign(x-m)*s, x))
        df1.index = df1['ASSAY'].map(lambda x: x.split(':')[0]+'_'+x.split(':')[2])
        for module in modules.columns:
            pair = variant+'-'+module
            comb = pd.concat([df1,modules[module].astype(bool)],join='inner',axis=1)
            if comb[comb[module]].shape[0] > 0:
                diff,pv = ztest(comb[comb[module]]['Z'],comb[~comb[module]]['Z'],value=0)
                diff2 = np.mean(comb[comb[module]]['Z']) - np.mean(comb[~comb[module]]['Z'])
            else:
                diff = 0
                pv = 1
            data['statistic'].update({pair: diff})
            data['p-value'].update({pair: pv})
            data['mean-diff'].update({pair: diff2})

    outDf = pd.DataFrame.from_dict(data)
    outDf['variant'] = outDf.index.map(lambda x: x.split('-')[0])
    outDf['module'] = outDf.index.map(lambda x: x.split('-')[1].replace('VAE_',''))
    outDf['FDR'] = multipletests(outDf['p-value'],method='fdr_bh')[1]
    outDf.sort_values('FDR',ascending=True).to_csv(chrom+'.pQTLEnrich.ztest2.csv')
    print(outDf)

    
def pQTLEnrichMW(prefix,chrom):
    import polars as pl
    from scipy.stats import mannwhitneyu
    with open('impactful_variant.list.txt','r') as file:
        variants = [i[:-1] for i in file.readlines()]
        
    modules = pd.read_csv(os.path.join('../VAE_protein_embedding',prefix+'_modules.txt'),sep=' ',index_col=0).astype('float')
    #df = pd.read_csv('impactful_variant.sumStats.csv.gz',low_memory=False)
    #df.columns = ['Variant','ASSAY','BETA','SE','LOG10P']
    df = pl.read_csv('tmp.%s.csv' % chrom).select(Variant=pl.col('CHROM').cast(str)+'_'+pl.col('GENPOS').cast(str)+'_'+pl.col('ALLELE0')+'_'+pl.col('ALLELE1'),ASSAY=pl.col('ASSAY'),BETA=pl.col('BETA'),SE=pl.col('SE'),LOG10P=pl.col('LOG10P')).filter(pl.col('Variant').is_in(variants)).to_pandas()
    data = {'statistic':{},'p-value':{}}
    for variant in df['Variant'].unique():
        print(variant)
        df1 = df.query('Variant == @variant')[['ASSAY','BETA','SE']].dropna(axis=0)
        df1['Z'] = df1['BETA']/df1['SE']
        df1.index = df1['ASSAY'].map(lambda x: x.split(':')[0]+'_'+x.split(':')[2])
        for module in modules.columns:
            pair = variant+'-'+module
            comb = pd.concat([df1,modules[module].astype(bool)],join='inner',axis=1)
            if comb[comb[module]].shape[0] > 0:
                diff,pv = mannwhitneyu(comb[comb[module]]['Z'],comb[~comb[module]]['Z'],method='asymptotic')
            else:
                diff = -1
                pv = 1
            data['statistic'].update({pair: diff})
            data['p-value'].update({pair: pv})

    outDf = pd.DataFrame.from_dict(data)
    outDf['variant'] = outDf.index.map(lambda x: x.split('-')[0])
    outDf['module'] = outDf.index.map(lambda x: x.split('-')[1].replace('VAE_',''))
    outDf['FDR'] = multipletests(outDf['p-value'],method='fdr_bh')[1]
    outDf.sort_values('FDR',ascending=True).to_csv(chrom+'.pQTLEnrich.MW.csv')
    print(outDf)

    
def pQTLEnrichNormARD(prefix,chrom):
    import polars as pl
    from sklearn.linear_model import ARDRegression
    with open('impactful_variant.list.txt','r') as file:
        variants = [i[:-1] for i in file.readlines()]
        
    modules = pd.read_csv(os.path.join('../VAE_protein_embedding',prefix+'_modules.txt'),sep=' ',index_col=0).astype('float')
    df = pl.read_csv('tmp.%s.csv' % chrom).select(Variant=pl.col('CHROM').cast(str)+'_'+pl.col('GENPOS').cast(str)+'_'+pl.col('ALLELE0')+'_'+pl.col('ALLELE1'),ASSAY=pl.col('ASSAY'),BETA=pl.col('BETA'),SE=pl.col('SE'),LOG10P=pl.col('LOG10P')).filter(pl.col('Variant').is_in(variants)).to_pandas()
    data = {'statistic':{},'p-value':{}}
    for variant in df['Variant'].unique():
        print(variant)
        df1 = df.query('Variant == @variant')[['ASSAY','BETA','SE']].dropna(axis=0)
        df1['Z'] = df1['BETA']/df1['SE']
        m = df1['Z'].mean()
        s = df1['Z'].std()
        df1['Z'] = df1['Z'].map(lambda x: np.where(abs(x - m)/s > 5, m + 5*np.sign(x-m)*s, x))
        df1.index = df1['ASSAY'].map(lambda x: x.split(':')[0]+'_'+x.split(':')[2])
        df1 = df1.loc[df1.index.isin(modules.index)]
        #y = df1.loc[modules.index]['Z']
        X = modules.copy().loc[df1.index]
        X = X.loc[:,X.apply(lambda x: any(x>0),axis=0)]
        clf = ARDRegression(threshold_lambda = 400)
        clf.fit(X,df1['Z'])
        w = pd.Series(clf.coef_, index=X.columns, name='posterior_w')
        w = w[w!=0]
        var = pd.Series(np.diagonal(clf.sigma_), index=w.index, name='var')
        dist = pd.concat([w,var],join='inner',axis=1)
        for module,row in dist.iterrows():
            pair = variant+'-'+module
            #cdf = st.norm.cdf(x=0,loc=row['posterior_w'],scale=np.sqrt(row['var']))
            pv = st.norm.cdf(x=0,loc=row['posterior_w'],scale=np.sqrt(row['var'])) if row['posterior_w'] > 0 else st.norm.sf(x=0,loc=row['posterior_w'],scale=np.sqrt(row['var']))
            data['statistic'].update({pair: row['posterior_w']})
            data['p-value'].update({pair: pv})

    outDf = pd.DataFrame.from_dict(data)
    outDf['variant'] = outDf.index.map(lambda x: x.split('-')[0])
    outDf['module'] = outDf.index.map(lambda x: x.split('-')[1].replace('VAE_',''))
    outDf.sort_values('p-value',ascending=True).to_csv(chrom+'.pQTLEnrich.ARD.csv')
    print(outDf)

    
def main():
    pQTLEnrichNormARD('protein-3layer-512-64-10-run2.leiden2.5-seed92-merged',sys.argv[1])
    #pQTLEnrichNorm('protein-3layer-512-64-10-run2.leiden2.5-seed92-merged',sys.argv[1])
    
    
if __name__=='__main__':
    main()
