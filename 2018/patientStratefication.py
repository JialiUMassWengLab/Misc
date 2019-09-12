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
from sklearn import preprocessing
from sklearn import decomposition
from sklearn.cluster import AgglomerativeClustering
from sklearn.cluster import KMeans
from scipy.cluster.hierarchy import linkage,fcluster

def getMatrix():
    df = pd.read_csv('../GatesADboth_r1_0frag_TPM.csv',index_col=0,low_memory=False).iloc[:-3,:].astype(float)
    df.index = map(lambda x: re.sub(r'\.\d+$','',x),df.index)
    df = df.drop(['17951','17949','17943','17944','17933','17934','17937','17938','1900','1914'],axis=1)
    #df = df.loc[df.apply(lambda x: np.median(x)>=1,axis=1),:]
    df = df.loc[df.apply(lambda x: np.mean(x)>=2,axis=1),:]
    #df = df.loc[df.apply(lambda x: len([i for i in x if i > 0]) >= df.shape[1]*0.5,axis=1),:]
    df['CV'] = df.apply(lambda x: np.std(x)/np.mean(x),axis=1)
    df = df[df['CV'] > 0.2].drop('CV',axis=1)
    maxval = df.apply(max,axis=1)
    df = df.apply(lambda x: x/max(x)*1000,axis=1)
    return df,maxval


def HierachyCluster():
    #df = pd.read_csv('W_TPM_p11_a0.mat.tsv',sep='\t',index_col=0)
    df = pd.read_csv('DESeq2rlog_DEgenesFDR0.10_files3/W_rlog_p6_a0.mat.tsv',sep='\t',index_col=0)
    info = pd.read_csv('../sampleInfo.csv',index_col=0)
    info = info[info['Disease']=='AD']
    df = df[map(lambda x: x in info.index,df.index)]
    #df = preprocessing.StandardScaler().fit_transform(df)
    df = df.apply(lambda x: x/x.quantile(0.75),axis=0)

    #cl = AgglomerativeClustering(n_clusters=5,affinity='l1',linkage='average')
    #cl = KMeans(n_clusters=5,random_state=132,algorithm='full',n_jobs=10,max_iter=50000,tol=0.00001)
    #y = cl.fit_predict(df)
    #methods = ['single','complete','weighted','median']
    #methods = ['centroid','ward']
    Z = linkage(df,method='average',metric='correlation')
    df['cluster'] = pd.Series(fcluster(Z,t=5,criterion='maxclust')-1,index=df.index)
    return df

    cols = {0:'red',1:'blue',2:'orange',3:'purple',4:'green',5:'black',6:'pink'}
    ccols = {0:'orange',1:'blue',2:'red',3:'black',4:'purple',5:'green'}
    rcols = dict(zip(range(5),sns.color_palette('Blues',5)))
    sns.set(context='talk',font_scale=1.8,style='white')
    g = sns.clustermap(df.iloc[:,:-1],row_linkage=Z,col_cluster=False,figsize=(30,28),vmax=3,row_colors=df['cluster'].map(rcols),col_colors=map(lambda x: ccols[x],range(6)))
    #g = sns.clustermap(df.iloc[:,:-1],method='average',metric='correlation',row_cluster=True,col_cluster=False,figsize=(30,28),vmax=4,row_colors=df['cluster'].map(cols))
    g.ax_heatmap.xaxis.set_ticks_position('bottom')
    g.savefig('cluster2.pdf')
    plt.close()


def plotBox(df):
    info = pd.read_csv('/home/jzhuang@ms.local/Neuro/GatesADcombined2/correlation/sampleInfo.csv',index_col=0)
    info = info[info['Disease']=='AD']
    info2 = pd.read_csv('/home/jzhuang@ms.local/Neuro/GatesADcombined2/NMF_files/MMSEInfo.csv',index_col=0)
    #df = pd.concat([info2[['diff','ratio']],df['cluster']],join='inner',axis=1)
    df = pd.concat([info[['Age','MMSE']].astype(int),df],join='inner',axis=1)
    #df = pd.concat([info['Center'],df],join='inner',axis=1)
    print df
    grouped = df.groupby('cluster')
    sns.set(font_scale=1.2,context='talk')
    from matplotlib.backends.backend_pdf import PdfPages
    #pdf = PdfPages('response_clust.pdf')
    pdf = PdfPages('dist_clust.pdf')
    for col in df.columns[:-1]:
        df1 = df[df[col]>0]
        pv = st.f_oneway(grouped.get_group(0)[col],grouped.get_group(1)[col],grouped.get_group(2)[col],grouped.get_group(3)[col],grouped.get_group(4)[col])[1]
        fig,ax = plt.subplots(figsize=(15,12))
        sns.boxplot(x='cluster',y=col,data=df1,fliersize=0,linewidth=2,width=0.4,color='black',ax=ax)
        sns.swarmplot(x='cluster',y=col,data=df1,size=9,edgecolor='black',linewidth=1,palette='Blues',ax=ax)
        #sns.boxplot(x='cluster',y=col,hue='Center',data=df,fliersize=0,linewidth=1.8,width=0.5,ax=ax)
        #sns.swarmplot(x='cluster',y=col,hue='Center',data=df,dodge=True,size=8,edgecolor='white',ax=ax)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.yaxis.set_ticks_position('left')
        ax.xaxis.set_ticks_position('bottom')
        for i,artist in enumerate(ax.artists):
            col = artist.get_facecolor()
            artist.set_edgecolor(col)
            artist.set_facecolor('None')
            for j in range(i*6,i*6+6):
                line = ax.lines[j]
                line.set_color(col)
                line.set_mfc(col)
                line.set_mec(col)

        plt.title('P-value = %.4f' % pv)
        plt.tight_layout()
        plt.savefig(pdf,format='pdf')
        plt.close()
    pdf.close()


def twowayANOVA(df):
    import statsmodels.formula.api as smf
    from statsmodels.stats.anova import anova_lm

    df.columns = map(lambda x: 'comp'+str(x) if not re.search(r'^cluster',x) else x,df.columns)
    info = pd.read_csv('/home/jzhuang@ms.local/Neuro/GatesADcombined2/correlation/sampleInfo.csv',index_col=0)
    info = info[info['Disease']=='AD']
    #df = pd.concat([info[['Age','MMSE']].astype(int),df],join='inner',axis=1)
    df = pd.concat([df,info['Center']],join='inner',axis=1)
    print df
    for comp in df.columns[:-2]:
        print comp
        formula = '%s ~ C(cluster) + C(Center) + C(cluster):C(Center)' % comp
        model = smf.ols(formula=formula,data=df).fit()
        aov_table = anova_lm(model, typ=2)
        print aov_table
        #print model.summary()


def main():
    df = HierachyCluster()
    #plotBox(df)
    twowayANOVA(df)


if __name__=='__main__':
    main()
