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
from sklearn.cluster import KMeans
from sklearn import preprocessing
from sklearn import decomposition

def normalize(tbl):
    df = pd.read_csv(tbl,sep='\t',index_col=0).iloc[:,:-1]
    df = df.loc[map(lambda x: not re.search(r'^ERCC-',x),df.index),:]
    df.index = map(lambda x: re.sub(r'\.\d+$','',x),df.index)
    print df.shape
    df = df.loc[df.iloc[:,:-1].apply(lambda x: max(x) > 50,axis=1),:]
    df = df.loc[df.iloc[:,:-1].apply(lambda x: (max(x)+1)/(min(x)+1) > 5,axis=1),:]
    df.iloc[:,:-1] = df.iloc[:,:-1].apply(lambda x: x/max(x),axis=1)
    #df.iloc[:,:-1] = preprocessing.StandardScaler().fit_transform(df.iloc[:,:-1])
    return df
    

def optNumClt():
    inertias = {}
    df = normalize('MM02_output/MM02.congregated.tsv')
    for clt in range(2,21):
        kmeans = KMeans(n_clusters=clt,random_state=132,algorithm='full',n_jobs=10,max_iter=50000,tol=0.00001).fit(df.iloc[:,:-1])
        inertias.update({clt: kmeans.inertia_})
    inert = pd.Series(inertias)
    sns.set(font_scale=1.8)
    inert.plot(lw=2,figsize=(12,10),fontsize=18,title='Elbow plot')
    plt.xlabel('Number of clusters (k)',fontsize=18)
    plt.ylabel('Inertia',fontsize=18)
    plt.savefig('elbow.png')
    plt.close()
    print inert


def printClusters(k):
    df = normalize('MM02_output/MM02.congregated.tsv')
    kmeans = KMeans(n_clusters=k,random_state=132,algorithm='full',n_jobs=10,max_iter=50000,tol=0.00001).fit(df.iloc[:,:-1])
    df['cluster'] = kmeans.labels_
    df['cluster'].to_csv('clusterNums.csv')


def hierCluster():
    df = normalize('MM02_output/MM02.congregated.tsv').iloc[:,:-1]
    df.columns = range(-2,16)
    print df.shape
    sns.set(font_scale=1.5)
    g = sns.clustermap(df,metric='correlation',z_score=None,standard_scale=None,row_cluster=True,col_cluster=False,cmap="YlGnBu",vmax=1.0,vmin=0.0,square=False,yticklabels=False,figsize=(10,12))
    plt.savefig('hierCluster.png')
    plt.close()
    

def plotClusters(k):
    df = normalize('MM02_output/MM02.congregated.tsv')
    kmeans = KMeans(n_clusters=k,random_state=132,algorithm='full',n_jobs=10,max_iter=50000,tol=0.00001).fit(df.iloc[:,:-1])
    centers = pd.DataFrame(kmeans.cluster_centers_,columns = df.columns[:-1])
    df['cluster'] = kmeans.labels_
    print df.shape
    '''
    from matplotlib.backends.backend_pdf import PdfPages
    pdf = PdfPages('clusters.pdf')
    sns.set(font_scale=1.8)
    for key,group in df.groupby('cluster'):
        group = group.transpose().iloc[:-2,:]
        group = pd.concat([group,centers.transpose().iloc[:,key]],join='inner',axis=1)
        group.index = range(-2,16)
        print group.shape
        fig,ax = plt.subplots(figsize=(12,10))
        group.iloc[:,:-1].plot(alpha=0.3,ax=ax,fontsize=18,color='grey',title='cluster'+str(key),legend=False)
        group[key].plot(ax=ax,color='blue',lw=2)
        plt.xlabel('time (days)',fontsize=18)
        plt.ylabel('Expression normalized to max',fontsize=18)
        plt.savefig(pdf,format='pdf')
        #plt.savefig('cluster'+str(key)+'.png')
        plt.close()
    pdf.close()
    '''
    from matplotlib.backends.backend_pdf import PdfPages
    pdf = PdfPages('clusters.heatmap.pdf')
    sns.set(font_scale=1.5)
    for key,group in df.groupby('cluster'):
        group = group.iloc[:,:-2]
        group = pd.concat([group,pd.DataFrame(centers.iloc[key,:]).transpose()],join='inner',axis=0)
        group.columns = range(-2,16)
        print group.shape
        dfSort = group.iloc[:-1,:]
        dfSort.loc[:,'dist'] = dfSort.apply(lambda x: np.sum((x - group.iloc[-1,:])**2),axis=1)
        dfSort = dfSort.sort_values(by='dist',ascending=True)
        fig,axes = plt.subplots(2,1,gridspec_kw={'height_ratios':[1,4]},figsize=(8,10))
        group.iloc[-1,:].transpose().plot(lw=2,ax=axes[0])
        axes[0].set_title('cluster'+str(key))
        #axes[0].tick_params(axis='x',bottom='off',top='off',which='both')
        axes[0].set_xticks([],[])
        axes[0].set_xlim(-2.5,15.5)
        sns.heatmap(dfSort.iloc[:,:-1],cmap="YlGnBu",vmax=1.0,vmin=0.0,square=False,yticklabels=False,cbar=False,ax=axes[1])
        plt.subplots_adjust(hspace=0.05)
        plt.xlabel('Time (days)',fontsize=18)
        plt.ylabel('Expression normalized to max',fontsize=18)
        plt.savefig(pdf,format='pdf')
        plt.close()
    pdf.close()

def plotDESeq2NormCounts():
    clusterNum = pd.read_csv('clusterNums.csv',header=None,index_col=0)
    df = pd.read_csv('DESeq2_normed_counts.MM02.tsv',index_col=0,sep='\t')
    df = df.loc[map(lambda x: not re.search(r'^ERCC-',x),df.index),:]
    df.index = map(lambda x: re.sub(r'\.\d+$','',x),df.index)
    df = df.apply(lambda x: 2**(x-max(x)),axis=1)
    df = pd.concat([df,clusterNum],join='inner',axis=1)
    print df

    from matplotlib.backends.backend_pdf import PdfPages
    pdf = PdfPages('clusters.normedCounts.heatmap.pdf')
    sns.set(font_scale=1.5)
    for key,group in df.groupby(1):
        group = group.iloc[:,:-1]
        group.loc['avg',:] = group.apply(np.mean,axis=0)
        group.columns = range(-2,16)
        print group.shape
        dfSort = group.iloc[:-1,:]
        dfSort.loc[:,'dist'] = dfSort.apply(lambda x: np.sum((x - group.iloc[-1,:])**2),axis=1)
        dfSort = dfSort.sort_values(by='dist',ascending=True)
        fig,axes = plt.subplots(2,1,gridspec_kw={'height_ratios':[1,4]},figsize=(8,10))
        group.iloc[-1,:].transpose().plot(lw=2,ax=axes[0])
        axes[0].set_title('cluster'+str(key))
        #axes[0].tick_params(axis='x',bottom='off',top='off',which='both')
        axes[0].set_xticks([],[])
        axes[0].set_xlim(-2.5,15.5)
        sns.heatmap(dfSort.iloc[:,:-1],cmap="YlGnBu",vmax=1.0,vmin=0.0,square=False,yticklabels=False,cbar=False,ax=axes[1])
        plt.subplots_adjust(hspace=0.05)
        plt.xlabel('Time (days)',fontsize=18)
        plt.ylabel('Expression normalized to max',fontsize=18)
        plt.savefig(pdf,format='pdf')
        plt.close()
    pdf.close()


def main():
    #plotClusters(10)
    #optNumClt()
    #printClusters(10)
    #hierCluster()
    plotDESeq2NormCounts()


if __name__=='__main__':
    main()
