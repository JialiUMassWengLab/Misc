import re
import os
import sys
import glob
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

def plotTopicDistOnImage(adata,pixelPropFn):
    n = int(pixelPropFn.split('.')[1][1:])
    topics = pd.read_csv(pixelPropFn,index_col=0)
    topics.index = topics.apply(lambda x: x.name.replace('-','-1-'), axis=1)
    topics.columns = list(map(lambda x: 'Topic'+str(x),topics.columns))
    topicList = list(topics.columns)
    adata.obs = pd.concat([adata.obs,topics],join='inner',axis=1)

    sns.set(font_scale=1.2)
    from matplotlib.backends.backend_pdf import PdfPages
    pdf = PdfPages(pixelPropFn.replace('.csv','.dist.pdf'))
    fig,ax = plt.subplots(figsize=(15,15))
    sc.pl.umap(adata, color=topicList, ncols=3)
    plt.savefig(pdf,format='pdf')
    plt.close()

    for topic in topicList:
        fig,axes = plt.subplots(2,2,figsize=(15,15))
        for n,sample in enumerate(['JY97A','JY98A','JY101A','JY102A']):
            i = int(n/2)
            j = n % 2
            upperbound = max(adata.obs[topic])
            sc.pl.spatial(adata[adata.obs['library_id']==sample,:].copy(), color=topic, vmin=0,
                          vmax = upperbound, library_id=sample, ncols=3, size=1.2, title=sample,
                          ax=axes[i][j])
        plt.suptitle(topic)
        plt.savefig(pdf,format='pdf')
        plt.close()

    pdf.close()
    

def plotTopicCorrelation(pixelPropFn):
    from scipy import stats
    topics = pd.read_csv(pixelPropFn,index_col=0)
    topics.index = topics.apply(lambda x: x.name.replace('-','-1-'), axis=1)
    topics.columns = list(map(lambda x: 'Topic'+str(x),topics.columns))
    topicList = list(topics.columns)
    
    corr = topics.corr()
    print(corr)
    fig,ax = plt.subplots(figsize=(15,15))
    sns.heatmap(corr,square=True,ax=ax,annot=True)
    plt.savefig(pixelPropFn.replace('.csv','.corrHeatmap.pdf'))
    
    
def main():
    sc.logging.print_header()
    print(f"squidpy=={sq.__version__}")
    #adata = sc.read_h5ad('squidpy.h5ad')
    #plotTopicDistOnImage(adata,sys.argv[1])

    plotTopicCorrelation(sys.argv[1])


if __name__=='__main__':
    main()
