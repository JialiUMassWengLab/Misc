#! /usr/bin/python

import re
import os
import sys
import glob
import numpy as np
import pandas as pd

def getGeneType():
    #gtf = pd.read_csv('/home/jzhuang/annotations/homo_sapiens/Homo_sapiens.GRCh38.96.gtf',sep='\t',skiprows=5,header=None,low_memory=False)
    #gtf = gtf[gtf[2]=='gene']
    #print(gtf)
    types = {}
    with open('/home/jzhuang/annotations/homo_sapiens/Homo_sapiens.GRCh38.96.gtf','r') as infile:
        for line in infile.readlines():
            if not re.search(r'^#!',line)==None:
                continue
            a = line.split('\t')
            if a[2] == 'gene':
                b = a[8].split(';')
                gene_id = re.search(r'\"\w+\"',b[0]).group()[1:-1]
                biotype = re.search(r'\"\w+\"',b[-2]).group()[1:-1]
                types.update({gene_id: biotype})
                
    return types


def getMatrix():
    analysis_types = ('protein_coding','lincRNA','processed_transcript','retained_intron','antisense')
    biotypes = getGeneType()
    files = glob.glob('*/abundance.tsv')
    geneMap = pd.read_csv('/home/jzhuang/annotations/homo_sapiens/transcripts_to_genes.txt',sep='\t',header=None)
    geneMap.columns = ['txn_id','gene_id','gene_name']
    geneMap.index = geneMap['txn_id']

    array = []
    for fn in files:
        sample = os.path.basename(os.path.dirname(fn))
        df = pd.read_csv(fn,sep='\t',index_col=0)
        df = pd.concat([df,geneMap['gene_id']],join='inner',axis=1)
        df = df.groupby('gene_id').agg(sum)
        array.append(df['tpm'].rename(sample))
        
    tpm = pd.concat(array,join='inner',axis=1)
    tpm.index = list(map(lambda x: re.sub(r'.\d+$','',x),tpm.index))
    tpm['biotype'] = tpm.index.map(biotypes)
    tpm = tpm[list(map(lambda x: x in analysis_types,tpm['biotype']))]
    tpm = tpm.iloc[:,:-1].apply(lambda x: x/sum(x)*1000000,axis=0)
    return tpm


def getMatrix2():
    analysis_types = ('protein_coding','lincRNA','processed_transcript','retained_intron','antisense')
    biotypes = getGeneType()
    files = glob.glob('/home/jzhuang/AMP/Microglia/bulk/*.genes.txt.gz')
    meta = pd.read_csv('/home/jzhuang/AMP/Microglia/bulk/metadata.csv',index_col=0,header=None)
    meta.columns = ['donorID','disease','region']
    #meta['donorID'] = meta.apply(lambda x: x['donorID'].replace('-','_'),axis=1)
    meta['donorID'] = meta.apply(lambda x: x['donorID'].split('-')[1],axis=1)

    array = []
    for fn in files:
        sample = os.path.basename(fn).split('_')[0]
        df = pd.read_csv(fn,sep='\t',index_col=0,header=None)
        df.columns = [meta.loc[sample,'region']+'_'+meta.loc[sample,'donorID']]
        array.append(df)

    tpm = pd.concat(array,join='outer',axis=1).fillna(0)
    tpm['biotype'] = tpm.index.map(biotypes)
    print(tpm.groupby('biotype').agg(sum))
    tpm = tpm[list(map(lambda x: x in analysis_types,tpm['biotype']))].iloc[:,:-1]
    #tpm = tpm.apply(lambda x: x/sum(x)*1000000,axis=0)
    return tpm


def main():
    df = getMatrix2()
    #df.to_csv('RNAseq_tpm_matrix.csv')
    df.to_csv('RNAseq_counts_matrix.csv')


if __name__=='__main__':
    main()
