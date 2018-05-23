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

def getMatrix():
    df = pd.read_csv('180305_NASHvalidation_203_avg_TPM.csv',index_col=0,low_memory=False).transpose()
    df1 = df[df['Disease']=='Normal Control'].iloc[:,:-3].astype('float').transpose()
    #df2 = df[df['PathFibrosis']=='F0'].iloc[:,:-3].astype('float').transpose()
    #df1 = df[(df['PathFibrosis']=='F0') | (df['PathFibrosis']=='F1')].iloc[:,:-3].astype('float').transpose()
    #df2 = df[(df['PathFibrosis']=='F3') | (df['PathFibrosis']=='F4')].iloc[:,:-3].astype('float').transpose()
    df2 = df[df['PathFibrosis']=='F4'].iloc[:,:-3].astype('float').transpose()
    return df1,df2


def main():
    cellType = sys.argv[1]
    names = pd.read_csv('1802NASHserum.congregated.tsv',sep='\t',index_col=0,usecols=['gene_id','Description'])
    df0 = pd.read_csv('/mnt/shares/Users/jzhuang/blood_cell_types/BloodCellTypes.expr4.tsv',sep='\t',index_col=0)
    df0.index = map(lambda x: re.sub(r'\.\d+$','',x),df0.index)

    df1,df2 = getMatrix()
    df1['median'] = df1.iloc[:,:-1].apply(np.median,axis=1)
    df1['sd'] = df1.iloc[:,:-2].apply(st.sem,axis=1)
    df2['median'] = df2.iloc[:,:-1].apply(np.median,axis=1)
    df2['sd'] = df2.iloc[:,:-2].apply(st.sem,axis=1)
    comb = pd.concat([df1['median'],df2['median'],names],join='inner',axis=1)
    #lab1 = 'F01'
    #lab2 = 'F34'
    lab1 = 'CTRL'
    #lab2 = 'F0IND'
    lab2 = 'F4IND'
    comb.columns = [lab1,lab2,'Description']

    sns.set(font_scale=4)
    if not cellType == 'Liver':
        geneList = list(df0[df0.iloc[:,:-1].apply(lambda x: len([i for i in x if x[cellType] < i*50])==1 and x[cellType] > 50,axis=1)].index)
    else:
        geneList = list([])
        with open('/home/jzhuang@ms.local/CoExpression/NMF_serumNASHset2_l1ratio1_out/H1_serum_NB_TPM_p6_a2000.comp0.genes','rU') as infile:
            for line in infile:
                geneList.append(line[:-1])

    comb1 = comb[map(lambda x: re.sub(r'\.\d+$','',x) in geneList,comb.index)]            
    comb1.index = comb1['Description']
    comb1 = comb1.drop('Description',axis=1)
    #comb1.to_csv('%s.genes.NASHvsCTRL.csv' % cellType)
    pcc = st.pearsonr(comb1[lab1],comb1[lab2])[0]
    comb1.apply(lambda x: x+0.01,axis=0).plot(lab1,lab2,loglog=True,kind='scatter',figsize=(30,25),s=400,fontsize=50,color='blue')
    plt.plot([0.1,100000],[0.1,100000],ls='--',lw=2.5,color='red')
    comb1 = comb1[(comb1[lab1]>1000) & (comb1[lab2]>1000)]
    #comb1 = comb1.loc[comb1.apply(lambda x: (x[lab1]+1)/(x[lab2]+1) > 1.5,axis=1),:]
    pack = zip(comb1.index,comb1[lab1],comb1[lab2])
    for label,x,y in pack:
        plt.annotate(label,xy=(x,y),xytext=(8,0),textcoords='offset points',ha='left',va='center')
    plt.xlim(1,30000)
    plt.ylim(1,30000)
    plt.title(cellType+'\nPearson Correlation Coefficient: '+'{0:.3f}'.format(pcc))
    plt.xlabel('SDBB Control Median TPM')
    plt.ylabel('Indiana F4 Median TPM')
    plt.savefig(cellType+'.cor.F4INDvsCTRL.png')
    #plt.xlabel('F01 Median TPM')
    #plt.ylabel('F34 Median TPM')
    #plt.savefig(cellType+'.cor.F01vsF34.png')
    plt.close()


if __name__=='__main__':
    main()

