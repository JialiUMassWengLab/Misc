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

def getMatrix(geneList):
    info = pd.read_csv('../../mouseLPSv2/sampleInfo.csv',index_col=0)
    info['analyte'] = info.apply(lambda x: '|'.join([x['AnimalName'],x['Tissue'],x['Treatment']]),axis=1)
    info = info.drop_duplicates('analyte',keep='last').drop('analyte',axis=1)
    info = info[info['uniqFragNum']>50000]
    info = info[map(lambda x: not re.search(r'Tofa',x),info['Treatment'])]
    info = info[info['Tissue']=='Plasma+16000g']
    info.index = info.index.astype(str)
    info1 = pd.read_csv('../sampleInfo1.csv',index_col=0)
    info1 = info1[(info1['Tissue']=='Serum') & (info1['uniqFragNum']>50000)]
    info1.index = info1.index.astype(str)
    info = pd.concat([info,info1],join='inner',axis=0)

    df = pd.read_csv('/mnt/nfs/analysis/180614_APK129LPS2PBMCRBC/180614_APK129LPS2PBMCRBC.congregated.tsv',sep='\t',index_col=0)
    names = df['Description']
    df2 = pd.read_csv('/mnt/nfs/analysis/180502_APK102MouseLPSPlasma/180502_APK102MouseLPSPlasma.congregated.tsv',sep='\t',index_col=0).iloc[:,:-2]
    df3 = pd.read_csv('/mnt/nfs/analysis/180616_APK129LPS2PlasmaPlt/180616_APK129LPS2PlasmaPlt.congregated.tsv',sep='\t',index_col=0).iloc[:,:-2]
    df4 = pd.read_csv('/mnt/nfs/analysis/180528_AI72MouseLPS2Plt/180528_AI72MouseLPS2Plt.congregated.tsv',sep='\t',index_col=0).iloc[:,:-2]
    df5 = pd.read_csv('/mnt/nfs/analysis/180120_AI055MouseLPS/180120_AI055MouseLPS.congregated.tsv',sep='\t',index_col=0).iloc[:,:-2]
    #df = pd.concat([df.iloc[:,:-2],df2,df3,df4,df5],join='inner',axis=1)
    df = pd.concat([df.iloc[:,:-2],df2,df3,df4],join='inner',axis=1)
    df.columns = map(lambda x: x.replace('_','-').split('-')[0],df.columns)
    df = df.loc[:,map(lambda x: x in list(info.index),df.columns)]
    df = df[map(lambda x: not re.search(r'^ERCC-',x),df.index)]
    df = df.apply(lambda x: x/sum(x)*1000000,axis=0)
    df = df[map(lambda x: names[x] in geneList,df.index)]
    return df,info

def getMatrixLiver(geneList):
    info = pd.read_csv('../../mouseLPSv2/sampleInfo.csv',index_col=0)
    info['analyte'] = info.apply(lambda x: '|'.join([x['AnimalName'],x['Tissue'],x['Treatment']]),axis=1)
    info = info.drop_duplicates('analyte',keep='last').drop('analyte',axis=1)
    info = info[info['uniqFragNum']>50000]
    info = info[map(lambda x: not re.search(r'Tofa',x),info['Treatment'])]
    info = info[info['Tissue']=='Liver']
    info.index = info.index.astype(str)
    #print info.sort_values('AnimalName')
    df = pd.read_csv('/mnt/nfs/analysis/180503_APK104MouseLPSLiver/180503_APK104MouseLPSLiver.congregated.tsv',sep='\t',index_col=0)
    names = df['Description']
    df = df.iloc[:,:-2].transpose()
    df['Alq'] = map(lambda x: x.replace('_','-').split('-')[0],df.index)
    df = df.groupby('Alq').agg(np.mean).transpose()
    df = df.loc[:,map(lambda x: x in list(info.index),df.columns)]
    df = df[map(lambda x: not re.search(r'^ERCC-',x),df.index)]
    df = df.apply(lambda x: x/sum(x)*1000000,axis=0)
    #df = df[df.apply(lambda x: any(x>0),axis=1)]
    df = df[map(lambda x: names[x] in geneList,df.index)]
    return df,info


def plotComp(df,info,tissue):
    W = pd.concat([df.transpose(),info[['Treatment','Tissue']]],join='inner',axis=1)
    #W = W[W['Tissue']=='Plasma+16000g']
    W['Treatment'] = map(lambda x: re.sub('post_','',x),W['Treatment'])
    array = []
    for ti,group in W.groupby('Tissue'):
        group = group.groupby('Treatment').agg(np.mean)
        group['time'] = map(lambda x: 0 if not re.search(r'^\d',x) else int(x.split('_')[0]),group.index)
        group['Drug'] = map(lambda x: 'CTRL' if not re.search(r'^\d',x) else 'LPS' if not re.search(r'AZD$',x) else 'LPS+AZD',group.index)
        if ti == 'Serum' or ti == 'Liver':
            group['Drug'] = map(lambda x: ti+' '+x,group['Drug'])
        array.append(group)
    W = pd.concat(array,join='inner',axis=0)
    W = W.sort_values(['time','Drug'])
    W = W[W['time']<=8]
    W1 = W[map(lambda x: not re.search(r'^Liver',x),W['Drug'])].iloc[:,:-2]
    geneList = W1.loc[:,W1.apply(lambda x: any(x>30),axis=0)].columns
    W = W.loc[:,map(lambda x: x in geneList or x=='time' or x=='Drug',W.columns)]
    #print W

    colors = {'LPS':'blue','LPS+AZD':'orange','Serum LPS':'red','Liver LPS':'blue','Liver LPS+AZD':'orange'}
    sns.set(font_scale=1.8,context='talk',style='white')
    #from matplotlib.backends.backend_pdf import PdfPages
    #pdf = PdfPages(Wmat.replace('.mat.tsv','.lines.pdf'))
    fig,ax = plt.subplots(figsize=(15,12))
    #for drug in ['LPS','LPS+AZD','Liver LPS','Liver LPS+AZD']:
    #for drug in ['LPS','LPS+AZD']:
    for drug in ['Liver LPS']:
        if not re.search(r'^Liver',drug):
            df1 = W[(W['Drug']==drug) | (W['Drug']=='CTRL')]
            s = 'o-'
        else:
            df1 = W[(W['Drug']==drug) | (W['Drug']=='Liver CTRL')]            
            s = 'x-'
        df1.iloc[:,:-2] = df1.iloc[:,:-2].apply(lambda x: np.log2((x+1)/(x[0]+1)),axis=0)
        #df1.iloc[:,:-2] = df1.iloc[:,:-2].apply(lambda x: (x+1)/(x[0]+1),axis=0)
        df1['avg'] = df1.iloc[:,:-2].apply(np.mean,axis=1)
        for gene in df1.columns[:-3]:
            df1.plot('time',gene,style='-',lw=2,c='lightblue',alpha=0.6,legend=False,ax=ax)
        df1.plot('time','avg',style=s,lw=4,c=colors[drug],markersize=15,markeredgecolor=colors[drug],markeredgewidth=3,label=drug,legend=False,ax=ax)

    #plt.xlim(-1,73)
    plt.xlim(-1,9)
    #plt.ylim(-0.1,W[col].max()*1.1)
    plt.xlabel('Time after treatment (hours)')
    plt.ylabel('Log2 Fold change relative to baseline')
    #plt.legend(loc='best',numpoints=1)
    plt.suptitle('')
    plt.tight_layout()
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    plt.savefig('%s.avg.timecourse3.pdf' % tissue)
    plt.close()
    #pdf.close()


def main():
    tissue = sys.argv[1]
    gene = pd.read_csv('../candidate.specific.genes3.csv',index_col=0)
    gene = gene[gene['tissue']==tissue]
    print gene
    df,info = getMatrix(gene.index)
    df2,info2 = getMatrixLiver(gene.index)
    df = pd.concat([df,df2],join='inner',axis=1)
    info = pd.concat([info,info2],join='inner',axis=0)
    plotComp(df,info,tissue)



if __name__=='__main__':
    main()

