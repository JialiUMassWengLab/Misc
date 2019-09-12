#! /usr/bin/python

import re
import os
import sys
import glob
import numpy as np
import pandas as pd
import scipy.stats as stat
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

def getMatrix():
    info = pd.read_csv('../sampleInfo1.csv',index_col=0)
    info['analyte'] = info.apply(lambda x: '|'.join([x['AnimalName'],x['Tissue'],x['Treatment']]),axis=1)
    info = info.drop_duplicates('analyte',keep='last').drop('analyte',axis=1)
    info = info[info['uniqFragNum']>50000]
    info.index = info.index.astype(str)
    #print info.sort_values('AnimalName')
    df = pd.read_csv('/mnt/nfs/analysis/180120_AI055MouseLPS/180120_AI055MouseLPS.congregated.tsv',sep='\t',index_col=0)
    names = df['Description']
    df2 = pd.read_csv('/mnt/nfs/analysis/180223_MouseTisAIFeb/180223_MouseTisAIFeb.congregated.tsv',sep='\t',index_col=0).iloc[:,:-2]
    df3 = pd.read_csv('/mnt/nfs/analysis/180310_APK76AI62mouseLPSmouseEPO/180310_APK76AI62mouseLPSmouseEPO.congregated.tsv',sep='\t',index_col=0).iloc[:,:-2]
    df = pd.concat([df.iloc[:,:-2],df2,df3],join='inner',axis=1)
    df = df.loc[:,map(lambda x: not re.search(r'-22FEB2018',x),df.columns)]
    df.columns = map(lambda x: x.replace('_','-').split('-')[0],df.columns)
    df = df.loc[:,map(lambda x: x in list(info.index),df.columns)]
    df = df[map(lambda x: not re.search(r'^ERCC-',x),df.index)]
    df = df.apply(lambda x: x/sum(x)*1000000,axis=0)
    df = df[df.apply(lambda x: any(x>0),axis=1)]
    return df,info,names
    
    
def getNBgenes(ratio):
    #ratio = ratio.loc[:,map(lambda x: not re.search(r'AZD$',x),ratio.columns)]
    #ratio = ratio[ratio.apply(lambda x: len([i for i in x if i > 3]) == len([i for i in x if i > 0]),axis=1)]
    #ratio = ratio[ratio.apply(lambda x: len([i for i in x if i > 3]) >= len([i for i in x if i > 0])-2,axis=1)]
    ratio = ratio[ratio.apply(lambda x: len([i for i in x if i > 3]) >= len([i for i in x if i > 0]) * 0.9,axis=1)]
    ratio = ratio[ratio.apply(lambda x: len([i for i in x if i > 0]) >= 6,axis=1)]
    return ratio
    

def getFC(df,info,tissue):
    #regulators = ['lipopolysaccharide','IFNG','Interferon alpha','IL1B','STAT1']
    regulators = ['Acute Phase Response Signaling','IL-10 Signaling','Interferon Signaling','IL-6 Signaling','Toll-like Receptor Signaling']
    geneInfo = pd.read_csv('mouseLPSv1.%s.4hLPS.geneInfo.txt' % tissue,sep='\t',index_col=0,skiprows=2)['Ensembl']
    #genesDf = pd.read_csv('mouseLPSv1.%s.4hLPS.regulator.txt' % tissue,sep='\t',index_col=0,skiprows=2)
    genesDf = pd.read_csv('mouseLPSv1.%s.4hLPS.pathway.txt' % tissue,sep='\t',index_col=0,skiprows=2)
    genesDf = genesDf[map(lambda x: x in regulators,genesDf.index)]
    genes = []
    genesDict = {}
    #for i,g in genesDf['Target molecules in dataset'].iteritems():
    for i,g in genesDf['Molecules'].iteritems():
        genes += map(lambda x: geneInfo[x],g.split(','))
        genesDict.update({i:map(lambda x: geneInfo[x],g.split(','))})
    df = df[map(lambda x: re.sub(r'\.\d+$','',x) in genes,df.index)]
    df = pd.concat([df.transpose(),info[['AnimalName','Tissue','Treatment']]],join='inner',axis=1)
    df = df[df['Tissue']==tissue]
    df = df.groupby('Treatment').agg(np.median)

    array = []
    for r in regulators:
        df1 = df.iloc[:,map(lambda x: re.sub(r'\.\d+$','',x) in genesDict[r],df.columns)]
        df1 = df1.apply(lambda x: (x+10),axis=0)
        df1 = df1.apply(lambda x: x/x['pre_dose'],axis=0)
        df1[r] = df1.apply(np.median,axis=1)
        array.append(df1[r])

    fc = pd.concat(array,join='inner',axis=1)
    fc['tissue'] = [tissue] * fc.shape[0]
    return fc


def plotBar(df,names):
    from mpl_toolkits.mplot3d import axes3d
    labelDict = {'lipopolysaccharide':'LPS','IFNG':'Ifng','Interferon alpha':'Ifna','IL1B':'Il1b','STAT1':'Stat1'}
    c = sns.color_palette('Blues_r',4) + sns.color_palette('Reds_r',4) + sns.color_palette('Greens_r',4) + sns.color_palette('Greys_r',4)
    L = []
    color = []
    for i in df.columns:
        L.append(df[i])
        color += c
    z = np.hstack(L).ravel()
    xlabels = df.index.unique()
    ylabels = df.columns.unique()
    x = np.arange(df.shape[0])
    y = np.arange(df.shape[1])
    print xlabels
    print ylabels
    x_M, y_M = np.meshgrid(x, y, copy=False)

    sns.set(font_scale=1.8,context='talk',style='white')
    fig = plt.figure(figsize=(26,18))
    ax = fig.add_subplot(111, projection='3d')
    ax.w_xaxis.set_ticks(x + 0.5/2.)
    ax.w_yaxis.set_ticks(y + 0.5/2.)
    ax.w_xaxis.set_ticklabels(map(lambda x: str(x[1])+'h',xlabels))
    ax.w_yaxis.set_ticklabels(ylabels,ha='left',fontsize=24)
    #ax.w_yaxis.set_ticklabels(map(lambda x: labelDict[x],ylabels),ha='left')
    ax.w_zaxis.set_ticklabels(map(lambda x: x/2.0,range(2,10)),ha='left')
    #ax.set_xlabel('Time')
    #ax.set_ylabel('Upstream regulator')
    #ax.set_zlabel('Median fold change')
    ax.bar3d(x_M.ravel(), y_M.ravel(), np.ones_like(z), dx=0.5, dy=0.2, dz=z-1, color=color, alpha=0.6)
    ax.view_init(elev=32., azim=282)
    plt.tight_layout()
    #plt.savefig('regulator.FC.bar3d.pdf')
    plt.savefig('pathway.FC.bar3d.pdf')
    plt.close()


def plotBar2():
    labelDict = {'lipopolysaccharide':'LPS','IFNG':'Ifng','Interferon alpha':'Ifna','IL1B':'Il1b','STAT1':'Stat1'}
    tissues = ('Brain','Kidney','Liver','Lung')
    colors = ['blue']*5+['red']*5+['green']*5+['black']*5
    array = []
    for t in tissues:
        df0 = pd.read_csv('mouseLPSv1.%s.4hLPS.regulator.txt' % t,sep='\t',index_col=0,skiprows=2)
        array.append(df0.iloc[:5,:]['p-value of overlap'])
        
    S = pd.concat(array,axis=0)
    S.index = map(lambda x: x if not x in labelDict.keys() else labelDict[x],S.index)
    S = -np.log10(S)
    #print S
    sns.set(font_scale=1.6,context='talk',style='white')
    fig,ax = plt.subplots(figsize=(12,18))
    S.iloc[::-1].plot(kind='barh',color=colors[::-1],ax=ax)
    ax.hlines(map(lambda x: x+0.5,range(4,20,5)),*ax.get_xlim(),linewidth=1.8)
    ax.vlines(-np.log10(0.05),*ax.get_ylim(),linewidth=2.5,linestyles='dashed',color='red')
    #ax.spines['right'].set_visible(False)
    #ax.spines['top'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    ax.set_xlabel('-log P-value')
    plt.tight_layout()
    plt.savefig('regulator.pvalue.bar.pdf')
    plt.close()


def main():

    df,info,names = getMatrix()
    tissues = ['Liver','Kidney','Lung','Brain']
    array = []
    for tissue in tissues:
        array.append(getFC(df,info,tissue))
    df = pd.concat(array,join='inner',axis=0)
    df['time'] = map(lambda x: x.split('_')[0],df.index)
    df = df[df['time']!='pre']
    df.index = pd.MultiIndex.from_tuples(list(zip(df['tissue'],df['time'].astype(int))),names=['tissue','time'])
    df = df.iloc[:,:-2].sort_index(level=('tissue','time'))
    print df 
    plotBar(df,names)
    '''
    plotBar2()
    '''

if __name__=='__main__':
    main()
