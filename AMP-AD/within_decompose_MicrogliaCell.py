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
    coef = coef.rename({'Intercept': 'C(gene)[T.%s]' % df.iloc[0,2]})
    #coef.index = list(map(lambda x: re.search(r'ENSG\w+',x).group() if not re.search(r'C\(cluster\)',x) else re.search(r'\[T.\d+\]',x).group()[3:-1] + '_' + re.search(r'ENSG\w+',x).group(),coef.index))
    coef.index = list(map(lambda x: '_'.join(list(map(lambda y: re.search(r'\[T\.\w+\]',y).group()[3:-1],x.split(':')))),coef.index))
    print(coef)
    intercept = {}
    for key,value in fit.random_effects.items():
        intercept.update({key: value['Subject']})
    intercept = pd.Series(intercept)

    return coef,intercept


def patientCluster(df0,prefix):
    df0 = df0.copy().apply(lambda x: np.log2(x+.5),axis=0)

    clustNum = {'glia8':3, 'glia14':5, 'glia18':3}
    array = []
    moduleGenes = {}
    sns.set(context='talk',font_scale=1.6)
    from matplotlib.backends.backend_pdf import PdfPages
    pdf = PdfPages(prefix + '.patientCorrCluster.pdf')
    for module in ['glia8','glia14','glia18']:
        geneList = list(pd.read_csv(module+'.genes.txt',header=None,index_col=0).index)
        df = df0[list(map(lambda x: x in geneList,df0.index))].copy()
        df = df[df.apply(lambda x: all(x>=0),axis=1)]
        print(len(geneList),df.shape)
        Z = linkage(df.transpose(),method='average',metric='correlation')
        S = pd.Series(fcluster(Z,t=clustNum[module],criterion='maxclust'),index=df0.columns,name=module)
        print(S.value_counts())
        array.append(S)
        moduleGenes.update({module: list(df.index)})

        colors = S.map({1:'black',3:'grey'}).rename('cluster')
        fig,ax = plt.subplots(figsize=(15,12))
        dist = df.corr(method='pearson')
        g = sns.clustermap(dist,row_linkage=Z,col_linkage=Z,vmin=.8,row_colors=colors,colors_ratio=0.015,figsize=(30,25),tree_kws=dict(linewidths=2))
        g.ax_heatmap.xaxis.set_visible(False)
        g.ax_row_dendrogram.set_visible(False)
        #plt.suptitle('%s\nNumber of genes: %d' % (module,df.shape[0]))
        g.savefig(pdf,format='pdf')
        plt.close()
    pdf.close()

    return pd.concat(array,join='inner',axis=1),moduleGenes


def corrWithAvg(df0,prefix):
    df0 = df0.copy().apply(lambda x: np.log2(x+.5),axis=0)
    names = pd.read_csv('/home/jzhuang/annotations/GTEx/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz',sep='\t',index_col=0,skiprows=2)['Description']
    names.index = map(lambda x: re.sub(r'\.\d+$','',x),names.index)

    sns.set(font_scale=1.2,context='talk')
    from matplotlib.backends.backend_pdf import PdfPages
    #pdf = PdfPages(prefix+'.patientClusterPCC.pdf')
    for module in ['glia8','glia14','glia18']:
        moduleGenes = list(pd.read_csv(module+'.genes.txt',header=None,index_col=0).index)
        avg = pd.read_csv(module+'.avgExpr.csv',index_col=0)
        df = df0[list(map(lambda x: x in moduleGenes,df0.index))].copy()
        avg = avg[list(map(lambda x: x in df.index,avg.index))]
        pcc1 = df.apply(lambda x: st.pearsonr(x,avg['cluster1'])[0],axis=0).transpose()
        pcc2 = df.apply(lambda x: st.pearsonr(x,avg['cluster2'])[0],axis=0).transpose()
        PCC = pd.concat([pcc1,pcc2],join='inner',axis=1)
        #labels = clusterDf[module].value_counts().keys()
        #PCC = PCC[(PCC[module]==labels[0]) | (PCC[module]==labels[1])].sort_values(module)
        print(PCC)


def plotScatter(df0,moduleGenes,clusterDf,prefix):
    df0 = df0.copy().apply(lambda x: np.log2(x+.5),axis=0)
    names = pd.read_csv('/home/jzhuang/annotations/GTEx/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz',sep='\t',index_col=0,skiprows=2)['Description']
    names.index = map(lambda x: re.sub(r'\.\d+$','',x),names.index)

    sns.set(context='talk',font_scale=1.6)
    from matplotlib.backends.backend_pdf import PdfPages
    pdf = PdfPages(prefix + '.patientClusterScatter.pdf')
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


def plotBox(df0,moduleGenes,clusterDf,prefix):
    df0 = df0.copy().apply(lambda x: np.log2(x+.5),axis=0)
    names = pd.read_csv('/home/jzhuang/annotations/GTEx/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz',sep='\t',index_col=0,skiprows=2)['Description']
    names.index = map(lambda x: re.sub(r'\.\d+$','',x),names.index)
    pvs = pd.read_csv('DESeq2.cluster.sig',index_col=0,header=None)

    module = 'glia14'
    labels = clusterDf[module].value_counts().keys()
    sns.set(context='talk',font_scale=1.2)
    from matplotlib.backends.backend_pdf import PdfPages
    pdf = PdfPages(prefix + '.patientClusterBox.pdf')
    df = df0[list(map(lambda x: x in moduleGenes[module],df0.index))].copy()
    df = pd.concat([df.transpose(),clusterDf[module]],join='inner',axis=1)
    df = df[(df[module]==labels[0]) | (df[module]==labels[1])]
    df['tissue'] = list(map(lambda x: x.split('_')[0],df.index))
    #for gene in ('C1QA','C1QB','C1QC','FCGBP','FCGR3A','CD14','FYB','C3','RUNX1','INPP5D','CSF3R','RCSD1','HAMP','TYMP'):
    fig,axes = plt.subplots(2,4,figsize=(32,15))
    for i,gene in enumerate(['C1QA','C1QB','C1QC','FCGBP','FCGR3A','CD14','HAMP','TYMP']):
        geneID = names[names==gene].index[0]
        #fig,ax = plt.subplots(figsize=(15,12))
        sns.boxplot(x=module,y=geneID,data=df,color='lightgrey',fliersize=0,width=.4,ax=axes[int(i/4)][i%4])
        sns.swarmplot(x=module,y=geneID,data=df,hue='tissue',linewidth=1,edgecolor='white',s=12,ax=axes[int(i/4)][i%4])
        axes[int(i/4)][i%4].set_title(gene+'\nFDR = %.3g' % pvs.loc[geneID,2],style='italic',fontweight='bold')
        axes[int(i/4)][i%4].set_xlabel('')
        axes[int(i/4)][i%4].set_xticklabels(['cluster B','cluster A'])
        axes[int(i/4)][i%4].set_ylabel('log2 TPM')
        axes[int(i/4)][i%4].legend().set_visible(False)

    plt.tight_layout()
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

    clusterDf,moduleGenes = patientCluster(df,fn_prefix)
    #print(clusterDf)
    #clusterDf['tissue'] = clusterDf.apply(lambda x: x.name.split('_')[0],axis=1)
    #clusterDf['donorID'] = clusterDf.apply(lambda x: x.name.split('_')[1],axis=1)
    #clusterDf.to_csv(fn_prefix+'.patientCorrCluster.csv')
    #corrWithAvg(df,fn_prefix)
    #plotScatter(df,moduleGenes,clusterDf,fn_prefix)
    #plotBox(df,moduleGenes,clusterDf,fn_prefix)


if __name__=='__main__':
    main()
