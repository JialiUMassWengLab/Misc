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
from scipy.cluster.hierarchy import linkage,fcluster

def fitComponent(df,geneList,cluster):
    df = df[list(map(lambda x: x in geneList,df.index))]
    df = pd.melt(df,ignore_index=False,var_name='Subject',value_name='expr')
    df.loc[:,'gene'] = df.index
    cluster.name = 'cluster'
    labels = cluster.value_counts().keys()
    df = df.merge(pd.DataFrame(cluster),left_on='Subject',right_index=True,how='inner')
    #print(df.shape)
    df = df[(df['cluster']==labels[0]) | (df['cluster']==labels[1])]
    model = smf.mixedlm("expr ~ C(gene)*C(cluster)",data=df,groups='Subject',re_formula="1")
    fit = model.fit(method='bfgs')
    #print(fit.summary())
    '''
    for cls in labels[:2]:
        model = smf.mixedlm("expr ~ C(gene)",data=df[df['cluster']==cls],groups='Subject',re_formula="1")
        #model = smf.mixedlm("expr ~ C(gene)",data=df,groups='Subject',re_formula="1")
        fit = model.fit(method='bfgs')
        print(fit.summary())
    '''
    coef = pd.Series(fit.fe_params,name='coef')
    coef['Intercept'] = 0.
    coef['C(cluster)[T.%d]' % max(labels[:2])] = 0.
    coef = coef.rename({'Intercept': 'C(gene)[T.%s]' % df.iloc[0,2]})
    #coef.index = list(map(lambda x: re.search(r'ENSG\w+',x).group() if not re.search(r'C\(cluster\)',x) else re.search(r'\[T.\d+\]',x).group()[3:-1] + '_' + re.search(r'ENSG\w+',x).group(),coef.index))
    coef.index = list(map(lambda x: '_'.join(list(map(lambda y: re.search(r'\[T\.\w+\]',y).group()[3:-1],x.split(':')))),coef.index))
    coef.index = list(map(lambda x: df.iloc[0,2]+'_'+x if not re.search(r'^ENSG',x) else x,coef.index))
    print(coef)
    intercept = {}
    for key,value in fit.random_effects.items():
        intercept.update({key: value['Subject']})
    intercept = pd.Series(intercept)

    return coef,intercept


def patientCluster(df0,module_gene_dir,prefix):
    genes = pd.read_csv(os.path.join(module_gene_dir,'modules.assigned.genes.txt'),index_col=0,sep=' ')
    df1 = df0.copy().transpose()
    df1 = pd.DataFrame(preprocessing.StandardScaler().fit_transform(df1),index=df0.columns,columns=df0.index)
    highGenes = set(df0[df0.apply(lambda x: np.mean(x) > 5,axis=1)])
    df0 = df0.copy().apply(lambda x: np.log2(x),axis=0)

    clustNum = {'glia8':9, 'glia14':9, 'glia18':23}
    array = []
    moduleGenes = {}
    sns.set(context='talk',font_scale=1.6)
    from matplotlib.backends.backend_pdf import PdfPages
    #pdf = PdfPages(prefix + '.patientCorrCluster.pdf')
    for module in ['glia8','glia14','glia18']:
        #geneList = selectRepGenes(df1,set(genes[genes[module]==1].index).intersection(highGenes))
        #geneList = selectRepGenes(df1,set(genes[genes[module]==1].index))
        geneList = list(pd.read_csv(module+'.genes.txt',header=None,index_col=0).index)
        df = df0[list(map(lambda x: x in geneList,df0.index))].copy()
        #print(df.shape)
        Z = linkage(df.transpose(),method='average',metric='correlation')
        S = pd.Series(fcluster(Z,t=clustNum[module],criterion='maxclust'),index=df0.columns,name=module)
        print(S.value_counts())
        array.append(S)
        moduleGenes.update({module: list(df.index)})
    '''
        fig,ax = plt.subplots(figsize=(15,12))
        dist = df.corr(method='pearson')
        g = sns.clustermap(dist,row_linkage=Z,col_linkage=Z,vmin=.85,figsize=(30,25))
        g.ax_heatmap.xaxis.set_ticks_position('bottom')
        plt.title('%s\nNumber of genes: %d' % (module,df.shape[0]))
        g.savefig(pdf,format='pdf')
        plt.close()
    pdf.close()
    '''
    return pd.concat(array,join='inner',axis=1),moduleGenes


def plotScatter(df0,moduleGenes,clusterDf,prefix):
    df0 = df0.copy().apply(lambda x: np.log2(x),axis=0)
    names = pd.read_csv('/home/jzhuang/annotations/GTEx/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz',sep='\t',index_col=0,skiprows=2)['Description']
    names.index = map(lambda x: re.sub(r'\.\d+$','',x),names.index)

    sns.set(context='talk',font_scale=1.6)
    from matplotlib.backends.backend_pdf import PdfPages
    pdf = PdfPages(prefix + '.patientClusterScatter.pdf')
    #for module in ['assigned8','assigned23','assigned13']:
    for module in ['glia8','glia14','glia18']:
        coef,intercept = fitComponent(df0,moduleGenes[module],clusterDf[module])
        df = pd.DataFrame(coef)
        df['gene'] = df.apply(lambda x: x.name.split('_')[0],axis=1)
        df['cluster'] = df.apply(lambda x: 'cluster1' if not re.search(r'_',x.name) else 'cluster2',axis=1)
        df = df.pivot(index='gene',columns='cluster',values='coef')
        df.iloc[1,1] = df.iloc[0,0]
        df['deviation'] = abs(df['cluster2'])
        df['cluster2'] = df['cluster1'] + df['cluster2']
        df = pd.concat([df,names],join='inner',axis=1)
        df = df.sort_values('deviation',ascending=False)

        fig,ax = plt.subplots(figsize=(15,12))
        sns.scatterplot(x='cluster1',y='cluster2',data=df,color='blue',alpha=.8,linewidth=2,s=100,ax=ax)
        ax.axline((0,0),slope=1,color='red',linestyle='--',alpha=.8)
        top = df.iloc[:10,:]
        for label,x,y in zip(top['Description'],top['cluster1'],top['cluster2']):
            plt.annotate(label,xy=(x,y),xytext=(5,0),textcoords='offset points',ha='left',va='center',fontsize=15)
        plt.title(module)
        plt.savefig(pdf,format='pdf')
        plt.close()
    pdf.close()


def pickRatio(df0,moduleGenes,clusterDf,prefix):
    df0 = df0.copy().apply(lambda x: np.log2(x),axis=0)
    names = pd.read_csv('/home/jzhuang/annotations/GTEx/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz',sep='\t',index_col=0,skiprows=2)['Description']
    names.index = map(lambda x: re.sub(r'\.\d+$','',x),names.index)

    sns.set(context='talk',font_scale=1.2)
    from matplotlib.backends.backend_pdf import PdfPages
    pdf = PdfPages(prefix + '.patientClusterRatioBox.pdf')
    for module in ['assigned8','assigned23','assigned13']:
        df = df0[list(map(lambda x: x in moduleGenes[module],df0.index))].copy()
        ratios = []
        for i in range(df.shape[0]):
            gene1 = df.index[i]
            for j in range(i+1,df.shape[0]):
                gene2 = df.index[j]
                ratios.append((df.loc[gene1,:]-df.loc[gene2,:]).rename(gene1+'_'+gene2))
                
        dRatio = pd.concat(ratios,join='inner',axis=1)
        labels = clusterDf[module].value_counts().keys()
        dRatio = pd.concat([dRatio,clusterDf[module]],join='inner',axis=1)
        pvalues = {}
        for pair in dRatio.columns[:-1]:
            d1 = ws.DescrStatsW(dRatio[dRatio[module]==labels[0]][pair])
            d2 = ws.DescrStatsW(dRatio[dRatio[module]==labels[1]][pair])
            cm_obj = ws.CompareMeans(d1,d2)
            pvalues.update({pair: cm_obj.ztest_ind()[1]})
        PVs = pd.Series(pvalues,name='pvalues').sort_values(ascending=True)
        print(PVs)

        dRatio = dRatio[list(map(lambda x: x in labels[:4],dRatio[module]))]
        for pair in PVs.index[:5]:
            gene1,gene2 = pair.split('_')
            fig,ax = plt.subplots(figsize=(15,12))
            sns.violinplot(x=module,y=pair,data=dRatio,color='lightgrey',ax=ax)
            sns.stripplot(x=module,y=pair,data=dRatio,color='blue',linewidth=1,edgecolor='white',ax=ax)
            plt.title('%s  %s\n%s over %s  P-value: %.3g' % (module,pair,names[gene1],names[gene2],PVs[pair]))
            plt.savefig(pdf,format='pdf')
            plt.close()
    pdf.close()


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

    return list(loading.index[:min(100,loading.shape[0])])


def main():
    fn_prefix = sys.argv[1]
    df = pd.read_csv(fn_prefix+'.csv',low_memory=False,index_col=0)
    df.columns = map(lambda x: re.sub(r'^X','',x),df.columns)
    df.index = map(lambda x: re.sub(r'\.\d+$','',x),df.index)
    df = df.astype('float')
    #for values > 5 s.d. above the mean, cap them at 5 s.d.
    df = df.apply(lambda x: pd.Series(np.where((x - np.mean(x))/np.std(x) > 5, np.mean(x)+5*np.std(x), x),index=df.columns),axis=1)

    module_gene_dir = '/home/jzhuang/AMP/Output/MSBB_PHG_WGCNA_files_p7_a0'
    clusterDf,moduleGenes = patientCluster(df,module_gene_dir,fn_prefix)
    #print(clusterDf)
    plotScatter(df,moduleGenes,clusterDf,fn_prefix)
    #pickRatio(df,moduleGenes,clusterDf,fn_prefix)


if __name__=='__main__':
    main()
