#! /usr/bin/python

import os
import re
import sys
import glob
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from tabulate import tabulate

def joinClinical():
    specimens = []
    files = glob.glob('/home/jzhuang/AMP/Source/*abundance-matrix.txt.gz')
    for fn in files:
        d = pd.read_csv(fn,sep='\t',nrows=1,index_col=0)
        if not re.search(r'^MSSM',os.path.basename(fn))==None:
            d = d.rename(columns={'hB_RNA_8475':'hB_RNA_8475_L43C014',
                                  'hB_RNA_9208_resequenced':'hB_RNA_9208_L43C014',
                                  'hB_RNA_10892':'hB_RNA_10892_K77C014',
                                  'hB_RNA_8255':'hB_RNA_8255_L43C014',
                                  'hB_RNA_4991':'hB_RNA_4991_B82C014'})

        specimens += list(d.columns)
        print(fn,d.shape[1])

    array = []
    files = glob.glob(os.path.join('/home/jzhuang/AMP/metadata','*_biospecimen_metadata.csv'))
    for fn in files:
        d = pd.read_csv(fn,index_col=1)
        d = d[list(map(lambda x: not x=='NA',d['individualID'].astype(str)))]
        d = d[~d.index.duplicated()]
        #print(fn)
        #print(d[d.apply(lambda x: x.name in specimens and x['exclude']==True,axis=1)].iloc[:,-2:])
        array.append(d[['individualID','tissue']])
    tissue = pd.concat(array,join='inner',axis=0)
    tissue = tissue[list(map(lambda x: x in specimens,tissue.index))]

    array = []
    files = glob.glob(os.path.join('/home/jzhuang/AMP/metadata','*_assay_[rR]*.csv'))
    for fn in files:
        d = pd.read_csv(fn,index_col=0)
        d = d[~d.index.duplicated()]
        array.append(d[['RIN']])
    assay = pd.concat(array,join='inner',axis=0)
    tissue = pd.concat([tissue,assay],join='inner',axis=1)

    mayo_dx_code = {'Alzheimer Disease':'AD','control':'NCI','pathological aging':'PA','progressive supranuclear palsy':'PSP'}
    rosmap_dx_code = {1:'NCI',2:'MCI',3:'MCI+',4:'AD',5:'AD+',6:'Other'}
    array = []
    d = pd.read_csv('/home/jzhuang/AMP/metadata/ROSMAP_clinical.csv',dtype={'apoe_genotype':str},index_col=17)
    d = d.rename(columns={'msex':'sex','apoe_genotype':'apoeGenotype','educ':'yearsEducation','age_death':'ageDeath','cogdx':'diagnosis','braaksc':'Braak','ceradsc':'CERAD'})
    d.loc[:,'sex'] = list(map(lambda x: 'male' if x==1 else 'female',d['sex']))
    d.loc[:,'diagnosis'] = d['diagnosis'].map(rosmap_dx_code)
    d.loc[:,'Study'] = ['ROSMAP'] * d.shape[0]
    array.append(d)
    d = pd.read_csv('/home/jzhuang/AMP/metadata/MSBB_individual_metadata.csv',dtype={'apoeGenotype':str},index_col=0)
    d.loc[:,'diagnosis'] = list(map(lambda x: 'NCI' if x==0 else 'MCI' if x<=1 else 'AD',d['CDR']))
    d.loc[:,'Study'] = ['MSBB'] * d.shape[0]
    array.append(d)
    d = pd.read_csv('/home/jzhuang/AMP/metadata/MayoRNAseq_individual_metadata.csv',dtype={'apoeGenotype':str},index_col=0)
    d.loc[:,'diagnosis'] = d['diagnosis'].map(mayo_dx_code)
    d.loc[:,'Study'] = ['Mayo'] * d.shape[0]
    array.append(d)
    clinical = pd.concat(array,join='inner',axis=0)
    clinical.loc[:,'apoeGenotype'] = list(map(lambda x: str(x)[0]+'/'+str(x)[1],clinical['apoeGenotype']))

    df = tissue.merge(clinical,left_on='individualID',right_index=True,how='inner')
    df.loc[:,'race'] = list(map(lambda x: 'W' if x==1 or x=='White' else 'B' if x==2 else 'U' if x==3 else x,df['race']))
    #print(df)
    print(tabulate(pd.crosstab(df['tissue'],df['Study']), headers='keys', tablefmt='psql'))
    print(df.drop_duplicates(subset='individualID')['Study'].value_counts())
    print(set(specimens).difference(set(df.index)))
    #df.to_csv('demo_clinical.csv')


def plotDist():
    df = pd.read_csv('/home/jzhuang/AMP/Output/demo_clinical.csv',dtype={'individualID':str},index_col=0)
    patient = df.drop(df.columns[1:3],axis=1).drop_duplicates()
    patient.index = patient['individualID']
    patient = patient.drop('individualID',axis=1)

    print(tabulate(pd.crosstab(patient['diagnosis'],patient['Study']), headers='keys', tablefmt='psql'))
    print(tabulate(pd.crosstab(patient['apoeGenotype'],patient['Study']), headers='keys', tablefmt='psql'))
    print(tabulate(pd.crosstab(patient['sex'],patient['Study']), headers='keys', tablefmt='psql'))
    print(tabulate(pd.crosstab(patient['race'],patient['Study']), headers='keys', tablefmt='psql'))
    print(tabulate(pd.crosstab(patient['CERAD'],patient['Study']), headers='keys', tablefmt='psql'))
    print(tabulate(pd.crosstab(patient['Braak'],patient['Study']), headers='keys', tablefmt='psql'))

    sns.set(context='talk',font_scale=1.2)
    from matplotlib.backends.backend_pdf import PdfPages
    pdf = PdfPages('demo_dist.pdf')
    fig,axes = plt.subplots(3,1,figsize=(15,15))
    i = 0
    for study,group in df.groupby('Study'):
        sns.histplot(data=group,x='RIN',ax=axes[i])
        axes[i].set_title(study)
        i += 1
    plt.suptitle('RIN')
    plt.tight_layout()
    plt.savefig(pdf,format='pdf')
    plt.close()

    for col in ('ageDeath','pmi'):
        fig,axes = plt.subplots(3,1,figsize=(15,15))
        i = 0
        for study,group in patient.groupby('Study'):
            if col=='ageDeath':
                group = group.dropna(axis=0,subset=[col])
                group.loc[:,col] = list(map(lambda x: 90.1 if x=='90+' else float(x),group[col]))
            sns.histplot(data=group,x=col,ax=axes[i])
            axes[i].set_title(study)
            i += 1
        plt.suptitle(col)
        plt.tight_layout()
        plt.savefig(pdf,format='pdf')
        plt.close()
    pdf.close()


def intersectWGS():
    df = pd.read_csv('/home/jzhuang/AMP/Output/demo_clinical.csv',dtype={'individualID':str},index_col=0)
    files = glob.glob('/home/jzhuang/AMP/metadata/*_biospecimen_metadata.csv')
    wgsPatients = []
    for fn in files:
        specimen = pd.read_csv(fn,index_col=0)
        specimen = specimen[specimen['assay']=='wholeGenomeSeq']        
        wgsPatients += list(specimen.index.astype(str))

    print(df[list(map(lambda x: not x in wgsPatients,df['individualID']))])
    df = df[list(map(lambda x: x in wgsPatients,df['individualID']))]
    print(tabulate(pd.crosstab(df['tissue'],df['Study']), headers='keys', tablefmt='psql'))
    print(df.drop_duplicates(subset='individualID')['Study'].value_counts())
        
def mapRNAseq2WGS(study):
    info = pd.read_csv('/home/jzhuang/AMP/Output/demo_clinical.csv',dtype={'individualID':str})
    samples = pd.read_csv('/home/jzhuang/AMP/Source/%s_all_quant_abundance-matrix.txt.gz' % study,nrows=1,sep='\t',index_col=0).columns
    info = info[list(map(lambda x: x in list(samples),info['specimenID']))]
    
    specimen = pd.read_csv('/home/jzhuang/AMP/metadata/%s_biospecimen_metadata.csv' % study)
    specimen = specimen[(specimen['assay']=='wholeGenomeSeq') & (specimen['exclude']==False)]
    df = info.iloc[:,:2].merge(specimen.iloc[:,:2],on='individualID')
    print(df.drop_duplicates(subset='individualID').shape)
    df.to_csv('/home/jzhuang/AMP/Derived/rnaseq2wgs.specimenMap.%s.csv' % study,index=False)
    

def ROSMAPbatch():
    demo = pd.read_csv('/home/jzhuang/AMP/Output/demo_clinical.csv',index_col=0)
    demo = demo[['RIN','pmi','Braak','CERAD']]
    assay = pd.read_csv('/home/jzhuang/AMP/metadata/ROSMAP_assay_rnaSeq_metadata.csv',index_col=0)
    picard = pd.read_csv('/home/jzhuang/AMP/metadata/ROSMAP_all_metrics_matrix.txt',sep='\t',index_col=0)
    picard = picard[['AlignmentSummaryMetrics__PCT_PF_READS','AlignmentSummaryMetrics__PCT_PF_READS_ALIGNED',\
                     'AlignmentSummaryMetrics__TOTAL_READS','RnaSeqMetrics__PCT_INTERGENIC_BASES',\
                     'RnaSeqMetrics__PCT_INTRONIC_BASES','RnaSeqMetrics__PCT_MRNA_BASES',\
                     'RnaSeqMetrics__PCT_RIBOSOMAL_BASES','RnaSeqMetrics__PCT_UTR_BASES']]
    df = pd.concat([assay['sequencingBatch'],demo,picard],join='inner',axis=1)
    df = df[df['sequencingBatch']!='0, 6, 7']
    df = df.sort_values('sequencingBatch',ascending=True)
    sns.set(context='talk',font_scale=1.2)
    from matplotlib.backends.backend_pdf import PdfPages
    pdf = PdfPages('ROSMAP_by_batch.pdf')
    for category in df.columns[1:]:
        fig,ax = plt.subplots(figsize=(15,12))
        sns.violinplot(x='sequencingBatch',y=category,data=df,color='lightgrey',ax=ax)
        sns.swarmplot(x='sequencingBatch',y=category,data=df,dodge=True,ax=ax)
        plt.title(category)
        plt.tight_layout()
        plt.savefig(pdf,format='pdf')
        plt.close()
    pdf.close()


def main():
    pd.options.display.max_columns = 100
    #joinClinical()
    #plotDist()

    #intersectWGS()
    mapRNAseq2WGS('ROSMAP')

    #ROSMAPbatch()


if __name__=='__main__':
    main()
    
        
