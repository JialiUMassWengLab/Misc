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
import seaborn as sns
from pandas_plink import read_plink

def getGT(varInfo,study):
    varID = varInfo[0]
    prefix = '/home/jzhuang/AMP/Derived/ROSMAP_maf0.01/chr%s_ROSMAP_maf0.01' % varInfo[1]
    (bim,fam,bed) = read_plink(prefix)
    target = bim[bim['snp']==varID]
    a0 = target['a0']
    a1 = target['a1']

    S = pd.Series(bed[target.i,:].compute()[0],index=fam['fid'],name='alleleCount')
    df0 = pd.read_csv('/home/jzhuang/AMP/Derived/rnaseq2wgs.specimenMap.%s.csv' % study,index_col=2)['specimenID_x']
    gt = pd.concat([S,df0],join='inner',axis=1)
    gt.loc[:,varID] = list(map(lambda x: a0+'/'+a0 if x==0 else a0+'/'+a1 if x==1 else a1+'/'+a1,gt['alleleCount']))
    gt.index = gt['specimenID_x']
    return gt[varID],gt['alleleCount']


def plotScatterByGenotype(prefix,gt,varID,pathway,gene):
    names = pd.read_csv('/home/jzhuang/annotations/GTEx/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz',sep='\t',index_col=0,skiprows=2)['Description']
    names.index = map(lambda x: re.sub(r'\.\d+$','',x),names.index)
    geneName = names[gene] if gene in names else 'NA'

    levels = pd.read_csv(prefix+'.pathLevels.csv',index_col=0)
    df = pd.read_csv(prefix+'.csv',low_memory=False,index_col=0).astype(float)
    df = df.apply(lambda x: np.log2(x),axis=0)
    df.columns = map(lambda x: re.sub(r'^X','',x),df.columns)
    df.index = map(lambda x: re.sub(r'\.\d+$','',x),df.index)

    from matplotlib.backends.backend_pdf import PdfPages
    pdf = PdfPages(prefix+'.byGenotype.scatter.%s_by_%s.%s.pdf' % (geneName,pathway,varID))
    sns.set(context='talk',font_scale=1.2)
    df1 = pd.concat([df.loc[gene],levels[pathway],gt],join='inner',axis=1)
    fig,ax = plt.subplots(figsize=(15,12))
    sns.scatterplot(x=pathway,y=gene,hue=varID,data=df1,alpha=.75,ax=ax)
    plt.title('%s  %s genotype' % (geneName,varID))
    plt.savefig(pdf,format='pdf')
    plt.close()
    pdf.close()


def plotBoxByGenotype(Wmat,gt,varID):
    W = pd.read_csv(Wmat,sep='\t',index_col=0).astype('float')
    W.columns = list(map(lambda x: 'comp'+str(int(x)+1),W.columns))
    df = pd.concat([W,gt],join='inner',axis=1)

    from matplotlib.backends.backend_pdf import PdfPages
    pdf = PdfPages(os.path.basename(Wmat).replace('.mat.tsv','')+'.byGenotype.box.%s.pdf' % varID)
    sns.set(context='talk',font_scale=1.2)
    for comp in W.columns:
        fig,ax = plt.subplots(figsize=(15,12))
        sns.boxplot(x=varID,y=comp,data=df,color='lightgrey',ax=ax)
        sns.swarmplot(x=varID,y=comp,data=df,dodge=True,ax=ax)
        plt.title('%s  %s genotype' % (comp,varID))
        plt.savefig(pdf,format='pdf')
        plt.close()
    pdf.close()


def plotBoxByGenotype2(levelFn,gt,varID):
    W = pd.read_csv(levelFn,index_col=0).astype('float')
    W.index = list(map(lambda x: re.sub(r'^X','',x),W.index))
    df = pd.concat([W,gt],join='inner',axis=1)
    df = df.sort_values(varID)
    print(df.shape)

    from matplotlib.backends.backend_pdf import PdfPages
    pdf = PdfPages(os.path.basename(levelFn).replace('.moduleLevels.csv','')+'.byGenotype.box.%s.pdf' % varID)
    sns.set(context='talk',font_scale=1.2)
    for comp in W.columns:
        fig,ax = plt.subplots(figsize=(15,12))
        sns.boxplot(x=varID,y=comp,data=df,color='lightgrey',ax=ax)
        sns.swarmplot(x=varID,y=comp,data=df,dodge=True,ax=ax)
        plt.title('%s  %s genotype' % (comp,varID))
        plt.savefig(pdf,format='pdf')
        plt.close()
    pdf.close()


def plotAFBarByGenotype(prefix,study,gt,varID):
    gt = gt.rename(varID)
    cluster = pd.read_csv(prefix+'.labels.csv',index_col=0)
    info = pd.read_csv('/home/jzhuang/AMP/Output/demo_clinical.csv',index_col=0)
    info = info[info['Study']==study]
    info.loc[:,'Diagnosis'] = list(map(lambda x: 'AD' if x<=2 else 'MCI' if x==3 else 'NCI',info['CERAD']))    
    
    cluster = pd.concat([cluster,info['Diagnosis'],gt],join='outer',axis=1).fillna(-1)
    cluster['cluster'] = cluster['cluster'].astype(int)
    cluster['Cluster'] = cluster.apply(lambda x: 'cluster'+str(x['cluster']+1) if x['cluster']>=0 else x['Diagnosis'],axis=1)
    cluster = cluster[cluster[varID]!=-1]
    df = cluster.groupby('Cluster').agg(np.mean)
    df = df.drop('AD',axis=0)
    df.loc[:,'AlleleFreq'] = df.apply(lambda x: x[varID]/2,axis=1)
    if df.loc['NCI','AlleleFreq'] > .5:
        df['AlleleFreq'] = 1 - df['AlleleFreq']
    df.loc[:,'cluster'] = df.index

    sns.set(font_scale=1.2,context='talk')
    fig,ax = plt.subplots(figsize=(15,12))
    sns.barplot(x='cluster',y='AlleleFreq',order=['NCI','MCI','cluster1','cluster2','cluster3','cluster4','cluster5'],data=df,ax=ax)
    plt.savefig(prefix+'.%s.AlleleFreqBar.pdf' % varID)
    plt.close()


def printTableByGenotype(prefix,module,gt,varID):
    cluster = pd.read_csv(prefix+'.exprRatio.%s.clusters.csv' % module,index_col=0)
    df = pd.concat([cluster,gt],join='inner',axis=1)
    print(pd.crosstab(df['cluster'],df[varID]))


def main():
    study = 'ROSMAP'
    #varInfo = ('rs1130474','8',100904277)
    #varInfo = ('rs1065691','12',6646320)
    varInfo = (sys.argv[1],sys.argv[2])
    gt,alleleCount = getGT(varInfo,study)
    print(gt)

    #plotScatterByGenotype('/home/jzhuang/AMP/Output/%s_geneLevel_limma_normed_counts' % study,gt,varInfo[0],'path4','ENSG00000164919')
    #plotScatterByGenotype('/home/jzhuang/AMP/Output/%s_geneLevel_limma_normed_counts' % study,gt,varInfo[0],'path4','ENSG00000269968')
    plotBoxByGenotype2(sys.argv[3],gt,varInfo[0])
    #printTableByGenotype('/home/jzhuang/AMP/Output/%s_geneLevel_limma_normed_counts' % study,'glia8',gt,varInfo[0])

    #plotAFBarByGenotype('/home/jzhuang/AMP/Output/%s_geneLevel_limma_normed_counts.hierarchy.AD.5clusters' % study,study,alleleCount,varInfo[0])


if __name__=='__main__':
    main()

