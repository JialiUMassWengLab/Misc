#! /usr/bin/python

import os
import re
import glob
import sys
import csv
import numpy as np
import scipy.stats as st
import pandas as pd

def getClasses():
    classes = {}
    names = {}
    allowed_types = ['protein_coding','lincRNA','processed_transcript','unprocessed_transcript','IG_C_gene','TR_C_gene']
    with open('/mnt/shares2/annotations/hg38/GeneInfo.csv','rU') as infile:
        reader = csv.DictReader(infile,fieldnames=['txnID','geneID','geneName','type','tissue','nonBlood'])
        for row in reader:
            if not row['tissue']=='' and row['type'] in allowed_types:
                if row['tissue']=='Whole':
                    row['tissue'] = 'Whole blood'
                classes[row['geneID']] = row['tissue']+'_ts'
            else:
                classes[row['geneID']] = row['type']

            names[row['geneID']] = row['geneName']

    return classes,names

def congregate(dirName,classes,names):
    title = os.path.basename(dirName)
    if title == '':
        title = os.path.basename(os.path.dirname(dirName))

    files = glob.glob(os.path.join(dirName,'*.dedup.genes.results'))
    matList = []
    matListRC = []
    for fn in files:
        name = os.path.basename(fn).replace('.dedup.genes.results','')
        mat = pd.read_csv(fn,delimiter='\t',index_col=0,usecols=[0,4,5])
        mat.columns = [name,name]
        matList.append(mat.iloc[:,1])
        matListRC.append(mat.iloc[:,0])
        
    df = pd.concat(matList,axis=1)
    df = df[map(lambda x: x in names, df.index)]
    df['Description'] = map(lambda x: names[x], df.index)
    df = df[map(lambda x: x in classes and not classes[x]=='rRNA' and not classes[x]=='misc_RNA' and \
                not classes[x]=='miRNA' and not classes[x]=='pseudogene' and not re.search(r'^s[\w]+RNA$',classes[x]) and \
                not re.search(r'^Mt_',classes[x]),df.index)]
    df = df[map(lambda x: not x=='Metazoa_SRP' and not x=='Y_RNA',df['Description'])]
    df.iloc[:,:-1] = df.iloc[:,:-1].apply(lambda x: x/sum(x)*1000000, axis=0)
    df['tissue'] = map(lambda x: re.sub(r'_ts$','',classes[x]) if x in classes and not re.search(r'_ts$',classes[x])==None else '',df.index)
    df.to_csv(os.path.join(dirName,title+'.congregated.tsv'),sep='\t')

    series = []
    df.iloc[:,:-2] = df.iloc[:,:-2].apply(lambda x: np.log10(x+0.01),axis=0)
    for col in df.columns[:-2]:
        vec = df[col]
        s = pd.Series(df.iloc[:,:-2].apply(lambda x: st.pearsonr(vec,x)[0], axis=0),name=col)
        series.append(s)
    cor = pd.concat(series,axis=1)
    cor.to_csv(os.path.join(dirName,title+'.exprPCC.tsv'),sep='\t')

    df = pd.concat(matListRC,axis=1)
    df = df[map(lambda x: x in names, df.index)]
    df['Description'] = map(lambda x: names[x], df.index)
    df = df[map(lambda x: x in classes and not classes[x]=='rRNA' and not classes[x]=='misc_RNA' and \
                not classes[x]=='miRNA' and not classes[x]=='pseudogene' and not re.search(r'^s[\w]+RNA$',classes[x]) and \
                not re.search(r'^Mt_',classes[x]),df.index)]
    df = df[map(lambda x: not x=='Metazoa_SRP' and not x=='Y_RNA',df['Description'])]
    df.iloc[:,:-1] = df.iloc[:,:-1].astype('int')
    df['tissue'] = map(lambda x: re.sub(r'_ts$','',classes[x]) if x in classes and not re.search(r'_ts$',classes[x])==None else '',df.index)
    df.to_csv(os.path.join(dirName,title+'.congregated.RC.tsv'),sep='\t')


def main():
    classes,names = getClasses()
    congregate(sys.argv[1],classes,names)


if __name__=='__main__':
    main()
