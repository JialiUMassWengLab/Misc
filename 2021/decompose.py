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
from scipy.cluster.hierarchy import linkage,fcluster

def fitComponent(df,geneList):
    df = df[list(map(lambda x: x in geneList,df.index))]
    df = pd.melt(df,ignore_index=False,var_name='Subject',value_name='expr')
    df.loc[:,'gene'] = df.index
    #print(df.shape)
    
    model = smf.mixedlm("expr ~ C(gene)",data=df,groups='Subject',re_formula="1")
    fit = model.fit(method='lbfgs')
    #print(fit.summary())

    coef = pd.Series(fit.fe_params,name='coef')
    coef['Intercept'] = 0.
    coef = coef.rename({'Intercept': 'C(gene)[T.%s]' % df.iloc[0,2]})
    coef.index = list(map(lambda x: re.search(r'ENSG\w+',x).group(),coef.index))
    intercept = {}
    for key,value in fit.random_effects.items():
        intercept.update({key: value['Subject']})
    intercept = pd.Series(intercept)

    return coef,intercept


def plotFit(df0,module_gene_dir,prefix):
    gene1 = pd.read_csv(os.path.join(module_gene_dir,'modules.neuron.genes.txt'),index_col=0,sep=' ')
    gene2 = pd.read_csv(os.path.join(module_gene_dir,'modules.glia.genes.txt'),index_col=0,sep=' ')
    genes = pd.concat([gene1,gene2],join='outer',axis=1).fillna(0).astype(int)
    df1 = df0.copy().transpose()
    df1 = pd.DataFrame(preprocessing.StandardScaler().fit_transform(df1),index=df0.columns,columns=df0.index)
    df0 = df0.apply(lambda x: np.log2(x),axis=0)

    sns.set(context='talk',font_scale=1.6)
    from matplotlib.backends.backend_pdf import PdfPages
    pdf = PdfPages(prefix + '.fitScatter.pdf')
    for module in genes.columns:
        geneList = selectRepGenes(df1,list(genes[genes[module]==1].index))
        print(module,len(geneList))
        try:
            coef, intercept = fitComponent(df0,geneList)
        except:
            print('%s failed fitting' % module)
        df = df0[list(map(lambda x: x in geneList,df0.index))]
        df = pd.melt(df,ignore_index=False,var_name='Subject',value_name='expr')
        df.loc[:,'gene'] = df.index

        quantiles = intercept.quantile([.1,.5,.9],interpolation='nearest').values
        samples = list(map(lambda x: intercept[intercept == x].index[0],quantiles))
        #print(samples)
        df2 = df[list(map(lambda x: x in samples,df['Subject']))]
        df2 = df2.merge(coef,left_on='gene',right_index=True)
        colors = {}
        c = ['blue','orange','green']
        for sub in intercept[samples].sort_values().index:
            colors.update({sub: c.pop()})
            
        fig,ax = plt.subplots(figsize=(15,12))
        sns.scatterplot(x='coef',y='expr',hue='Subject',palette=colors,data=df2,alpha=.8,ax=ax)
        for q in quantiles:
            ax.axline((0,q),slope=1,color=colors[intercept[intercept==q].index[0]],linestyle='--',alpha=.8)
        plt.title(module)
        plt.tight_layout()
        plt.savefig(pdf,format='pdf')
        plt.close()
    pdf.close()
        

def plotCorr(df0,module_gene_dir,prefix):
    gene1 = pd.read_csv(os.path.join(module_gene_dir,'modules.neuron.genes.txt'),index_col=0,sep=' ')
    gene2 = pd.read_csv(os.path.join(module_gene_dir,'modules.glia.genes.txt'),index_col=0,sep=' ')
    genes = pd.concat([gene1,gene2],join='outer',axis=1).fillna(0).astype(int)
    df1 = df0.copy().transpose()
    df1 = pd.DataFrame(preprocessing.StandardScaler().fit_transform(df1),index=df0.columns,columns=df0.index)
    df0 = df0.apply(lambda x: np.log2(x),axis=0)

    sns.set(context='talk',font_scale=1.6)
    from matplotlib.backends.backend_pdf import PdfPages
    pdf = PdfPages(prefix + '.PCCvsIntercept.pdf')
    for module in genes.columns:
        geneList = selectRepGenes(df1,list(genes[genes[module]==1].index))
        print(module,len(geneList))
        try:
            coef, intercept = fitComponent(df0,geneList)
        except:
            print('%s failed fitting' % module)
            continue

        df = df0[list(map(lambda x: x in geneList,df0.index))].copy()
        df = pd.melt(df,ignore_index=False,var_name='Subject',value_name='expr')
        df.loc[:,'gene'] = df.index

        stats = {'pcc':{},'rmse':{}}
        for sub in df['Subject'].unique():
            comb = pd.concat([df[df['Subject']==sub]['expr'],coef],join='inner',axis=1)
            pcc, pv = st.pearsonr(comb['coef'],comb['expr'])
            stats['pcc'].update({sub: pcc})
            stats['rmse'].update({sub: math.sqrt(metrics.mean_squared_error(comb['expr'],comb['coef']+intercept[sub]))})

        #PCCs = pd.Series(corr,name='pcc')
        result = pd.DataFrame.from_dict(stats)
        result = pd.concat([result,intercept],join='inner',axis=1)
        print(result.sort_values('rmse',ascending=True))
        result.columns = ['pcc','rmse','intercept']

        fig,ax = plt.subplots(figsize=(15,12))
        sns.scatterplot(x='intercept',y='rmse',data=result,alpha=.8,ax=ax)
        plt.title(module)
        plt.tight_layout()
        plt.savefig(pdf,format='pdf')
        plt.close()
    pdf.close()


def patientCluster(df0,module_gene_dir,prefix):
    gene1 = pd.read_csv(os.path.join(module_gene_dir,'modules.neuron.genes.txt'),index_col=0,sep=' ')
    gene2 = pd.read_csv(os.path.join(module_gene_dir,'modules.glia.genes.txt'),index_col=0,sep=' ')
    genes = pd.concat([gene1,gene2],join='outer',axis=1).fillna(0).astype(int)
    df1 = df0.copy().transpose()
    df1 = pd.DataFrame(preprocessing.StandardScaler().fit_transform(df1),index=df0.columns,columns=df0.index)
    #df0 = df0[df0.apply(lambda x: np.mean(x) > 5,axis=1)]
    df0 = df0.apply(lambda x: np.log2(x),axis=0)
    print(df0.shape)

    sns.set(context='talk',font_scale=1.6)
    from matplotlib.backends.backend_pdf import PdfPages
    pdf = PdfPages(prefix + '.patientCorrCluster0.pdf')
    for module in genes.columns:
        geneList = selectRepGenes(df1,list(genes[genes[module]==1].index))
        #geneList = list(genes[genes[module]==1].index)
        df = df0[list(map(lambda x: x in geneList,df0.index))].copy()
        Z = linkage(df.transpose(),method='average',metric='correlation')
        print(df.shape)

        fig,ax = plt.subplots(figsize=(15,12))
        dist = df.corr(method='pearson')
        g = sns.clustermap(dist,row_linkage=Z,col_linkage=Z,vmin=.8,figsize=(30,25))
        g.ax_heatmap.xaxis.set_ticks_position('bottom')
        plt.title(module)
        g.savefig(pdf,format='pdf')
        plt.close()
    pdf.close()


def plotScatter(df,gene1,gene2):
    df = df.loc[[gene1,gene2],:].transpose()
    sns.set(context='talk',font_scale=1.5)
    fig,ax = plt.subplots(figsize=(15,12))
    sns.scatterplot(x=gene1,y=gene2,data=df,ax=ax)
    plt.savefig('%s_vs_%s.scatter.pdf' % (gene1,gene2))
    plt.close()


def calCompLevels(df,module_gene_dir,prefix):
    gene1 = pd.read_csv(os.path.join(module_gene_dir,'modules.neuron.genes.txt'),index_col=0,sep=' ')
    gene2 = pd.read_csv(os.path.join(module_gene_dir,'modules.glia.genes.txt'),index_col=0,sep=' ')
    genes = pd.concat([gene1,gene2],join='outer',axis=1).fillna(0).astype(int)
    df1 = df.copy().transpose()
    df1 = pd.DataFrame(preprocessing.StandardScaler().fit_transform(df1),index=df.columns,columns=df.index)
    df = df.apply(lambda x: np.log2(x),axis=0)
    '''
    compLevels = {}
    for module in genes.columns:
        geneList = selectRepGenes(df1,list(genes[genes[module]==1].index))
        try:
            coef,intercept = fitComponent(df,geneList)
            compLevels.update({module: intercept})
        except:
            print('%s fitting failed' % module)

    levels = pd.DataFrame.from_dict(compLevels)
    levels.to_csv(prefix+'.pathLevels2.csv')

    pcc = levels.corr(method='pearson')
    #pcc = levels.apply(lambda x: 2**x,axis=0).corr(method='pearson')
    sns.set(context='talk')
    sns.clustermap(pcc,method='average',metric='euclidean',cmap='RdBu_r',row_cluster=True,col_cluster=True,figsize=(30,28),xticklabels=True,yticklabels=False)
    plt.tight_layout()
    plt.savefig(fn.replace('.cliques.csv','.pathCorrHeatmap.pdf'))
    plt.close()
    '''
    levels = pd.read_csv(prefix+'.pathLevels2.csv',index_col=0)
    df2 = df.transpose()
    df2 = df2.loc[list(levels.index),:]
    array = []
    for gene in df2.columns:
        S = pd.Series(levels.apply(lambda x: st.pearsonr(x,df2[gene])[0],axis=0),name=gene)
        array.append(S)

    corr = pd.concat(array,join='inner',axis=1).transpose()
    corr.to_csv(prefix+'.pathGeneCorr.csv')
    print(corr.apply(lambda x: len([i for i in x if i > .7]),axis=0))


def selectRepGenes(df,moduleGenes):
    df = df.loc[:,moduleGenes]
    pca = decomposition.PCA(n_components=2)
    #coor = pd.DataFrame(pca.fit_transform(df),index=sampleList)
    pca.fit(df)
    print(pca.explained_variance_ratio_)
    loading = pd.DataFrame(pca.components_,columns=df.columns).transpose()
    loading['abs'] = loading.apply(lambda x: abs(x[0]),axis=1)
    loading = loading.sort_values('abs',ascending=False)
    thred = max(loading['abs']) * 0.67
    loadingHigh = loading[loading['abs']>thred]
    #print(loading)

    return list(loading.index[:min(100,loading.shape[0],loadingHigh.shape[0])])


def main():
    fn_prefix = sys.argv[1]
    df = pd.read_csv(fn_prefix+'.csv',low_memory=False,index_col=0)
    df.columns = map(lambda x: re.sub(r'^X','',x),df.columns)
    df.index = map(lambda x: re.sub(r'\.\d+$','',x),df.index)
    df = df.astype('float')
    #for values > 5 s.d. above the mean, cap them at 5 s.d.
    df = df.apply(lambda x: pd.Series(np.where((x - np.mean(x))/np.std(x) > 5, np.mean(x)+5*np.std(x), x),index=df.columns),axis=1)
    df = df.drop('492_120515',axis=1)

    module_gene_dir = '/home/jzhuang/AMP/Output/WGCNA_files_p8_a0'
    #calCompLevels(df,module_gene_dir,fn_prefix)
    #plotFit(df,module_gene_dir,fn_prefix)
    #plotCorr(df,module_gene_dir,fn_prefix)
    patientCluster(df,module_gene_dir,fn_prefix)


if __name__=='__main__':
    main()
