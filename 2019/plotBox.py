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

def getMatrix():
    info = pd.read_csv('../mouseLPSv2/sampleInfo.csv',index_col=0)
    info['analyte'] = info.apply(lambda x: '|'.join([x['AnimalName'],x['Tissue'],x['Treatment']]),axis=1)
    info = info.drop_duplicates('analyte',keep='last').drop('analyte',axis=1)
    info = info[info['uniqFragNum']>50000]
    info = info[map(lambda x: not re.search(r'Tofa',x),info['Treatment'])]
    info = info[(info['Tissue']=='Plasma+16000g') | (info['Tissue']=='PBMC')]
    info['Tissue'] = map(lambda x: x.split('+')[0],info['Tissue'])
    info.index = info.index.astype(str)
    #print info.sort_values('AnimalName')
    df = pd.read_csv('/mnt/nfs/analysis/180614_APK129LPS2PBMCRBC/180614_APK129LPS2PBMCRBC.congregated.tsv',sep='\t',index_col=0)
    names = df['Description']
    df2 = pd.read_csv('/mnt/nfs/analysis/180502_APK102MouseLPSPlasma/180502_APK102MouseLPSPlasma.congregated.tsv',sep='\t',index_col=0).iloc[:,:-2]
    df3 = pd.read_csv('/mnt/nfs/analysis/180616_APK129LPS2PlasmaPlt/180616_APK129LPS2PlasmaPlt.congregated.tsv',sep='\t',index_col=0).iloc[:,:-2]
    df4 = pd.read_csv('/mnt/nfs/analysis/180528_AI72MouseLPS2Plt/180528_AI72MouseLPS2Plt.congregated.tsv',sep='\t',index_col=0).iloc[:,:-2]
    df = pd.concat([df.iloc[:,:-2],df2,df3,df4],join='inner',axis=1)
    df.columns = map(lambda x: x.replace('_','-').split('-')[0],df.columns)
    df = df.loc[:,map(lambda x: x in list(info.index),df.columns)]
    df = df[map(lambda x: not re.search(r'^ERCC-',x),df.index)]
    df = df.apply(lambda x: x/sum(x)*1000000,axis=0)
    df = df[df.apply(lambda x: any(x>0),axis=1)]
    return df,info,names    
    
def getMatrixLiver():
    info = pd.read_csv('../mouseLPSv2/sampleInfo.csv',index_col=0)
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
    df = df[df.apply(lambda x: any(x>0),axis=1)]
    return df,info,names
    
    
def plotBox(df,info,names):
    df = pd.concat([df.transpose(),info[['AnimalName','Tissue','Treatment']]],join='inner',axis=1)
    df['time'] = map(lambda x: 0 if not re.search(r'^\d',x) else int(x.split('_')[0]),df['Treatment'])
    df = df.sort_values(['time','Treatment']).drop('AnimalName',axis=1)
    df['Treatment'] = map(lambda x: re.sub('post_','',x),df['Treatment'])
    df = df[df['Tissue']=='Plasma']
    df['drug'] = map(lambda x: 'LPS' if not re.search(r'^\d',x) else x.split('_')[-1],df['Treatment'])
    df['drug'] = map(lambda x: x+'+LPS' if x=='AZD' else x,df['drug'])

    df2 = df.groupby('Treatment').agg(np.median)
    df2['time'] = map(lambda x: 0 if not re.search(r'^\d',x) else int(x.split('_')[0]),df2.index)
    df2['drug'] = map(lambda x: 'CTRL' if not re.search(r'^\d',x) else x.split('_')[-1],df2.index)
    df2['drug'] = map(lambda x: x+'+LPS' if x=='AZD' else x,df2['drug'])
    df2 = df2.sort_values('time')

    df = df.drop(['Treatment','Tissue'],axis=1)
    df0 = df
    geneList = list(df0.columns[:-2])
    df = pd.melt(df,id_vars=['drug','time'])
    df['name'] = df['gene_id'].map(names)
    print df[['gene_id','name']].drop_duplicates()

    sns.set(font_scale=1.5,style='white',context='talk')
    g = sns.catplot(x='time',y='value',hue='drug',data=df,col='name',col_wrap=3,kind='swarm',sharey=False,palette=['blue','orange'],edgecolor='white',dodge=True,aspect=1.2,height=6,legend_out=True,s=8)
    #p = g._get_palette(g.data,'drug',None,None)
    p = ['blue','orange']
    for i,ax in enumerate(g.axes):
        sns.boxplot(x='time',y=geneList[i],hue='drug',data=df0,linewidth=1.8,width=0.8,fliersize=0,palette=p,ax=ax)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.yaxis.set_ticks_position('left')
        ax.xaxis.set_ticks_position('bottom')
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.get_legend().remove()
        #ax.set_xticklabels(list(df['time'].unique()),rotation=45,ha='right')
        for j,artist in enumerate(ax.artists):
            col = artist.get_facecolor()
            artist.set_edgecolor(col)
            artist.set_facecolor('None')
            for k in range(j*6,j*6+6):
                line = ax.lines[k]
                line.set_color(col)
                line.set_mfc(col)
                line.set_mec(col)

    g.savefig('Fig1E.genes.box.pdf')
    #g.savefig('lympho.genes.box.pdf')
    #g.savefig('immune.genes.box.pdf')
    plt.close()
    #pdf.close()


def main():
    #geneList = ['Ifit1','Cxcl10','Ifit3','Oasl1','Isg15','Ifit3b']
    #geneList = ['Ifit1','Cxcl10','Ifit3','Fcrl1','Mzb1','Gzma']
    #geneList = ['Iglc1','Fcrl1','Gzma','Ly9','Mzb1','Bcl11a']
    geneList = ['Il1b','Ccl5','Cxcl2','Tnf','Ifng','Csf3']
    df,info,names = getMatrix()
    #df,info,names = getMatrixLiver()
    df = df[map(lambda x: names[x] in geneList,df.index)]
    print df.shape
    plotBox(df,info,names)



if __name__=='__main__':
    main()
