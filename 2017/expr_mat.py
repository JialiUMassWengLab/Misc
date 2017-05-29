#! /usr/bin/python

import re
import csv
import random
import pandas as pd
import numpy as np
from pyearth import Earth
import sklearn.feature_selection as fs
from sklearn.linear_model import LinearRegression
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

ALB_median = 27334.12169886384

def getSampleTissues():
    tissues = {}
    with open('/usr/local/bioinfo/MS/tissueExpnDB/GTEx_Data_V6_Annotations_SampleAttributesDS.txt','rU') as infile:
        reader = csv.DictReader(infile,delimiter='\t')
        for row in reader:
            tissues[row['SAMPID']] = row['SMTSD']

    return tissues

def getSamples(tissue,tissues,allSamples):
    samples = []
    for key,value in tissues.iteritems():
        if value == tissue:
            samples.append(key)

    samples = [sample for sample in samples if sample in list(allSamples)]
    return samples

def generateClfMatrix(tissue,nsamp,seed):
    tissues = getSampleTissues()
    filename = '/usr/local/bioinfo/MS/tissueExpnDB/GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct.gz'
    allSamples = pd.read_csv(filename,sep='\t',skiprows=[0,1],nrows=3).columns

    samples = getSamples(tissue,tissues,allSamples)
    data = pd.read_csv(filename,sep='\t',skiprows=[0,1],usecols=['Name','Description']+samples)
    data.iloc[:,2:] = data.iloc[:,2:].apply(lambda x: x/sum(x)*1000000, axis=0)
    data = data[data.iloc[:,2:].apply(lambda x:any(x>10), axis=1)]
    data = data[data.apply(lambda x: re.search(r'^MT-',x['Description'])==None, axis=1)]
    nfields = data.shape[1]

    random.seed(a=seed)
    ctrl_samps = []
    for t in set(tissues.values()):
        if not t == tissue:
            samps = getSamples(t,tissues,allSamples)
            if nsamp == 0:
                n = len(samps)
            else:
                n = min(nsamp,len(samps))
            ctrl_samps += random.sample(samps,n)
    ctrl = pd.read_csv(filename,sep='\t',skiprows=[0,1],usecols=['Name']+ctrl_samps)
    ctrl.iloc[:,1:] = ctrl.iloc[:,1:].apply(lambda x: x/sum(x)*1000000, axis=0)

    merge = data.merge(ctrl,how='inner',on='Name')
    merge.index = merge['Name']
    merge = merge.drop('Name',axis=1)
    merge = merge.transpose()
    merge['Target'] = [None]+map(lambda x: x in samples, merge.index[1:])

    return merge

def getSpecificGenes(tissue):
    for i in range(50):
        data = getMatrix(tissue,5,19+i*12)
        coef = LRclassify(data.iloc[:,:-1],list(data['Target']))
        coef.to_csv('LR.genes.%d.tsv' % i, sep='\t')
        coef = featureSelect(data.iloc[:,:-1],list(data['Target']))
        coef.to_csv('Ftest.genes.%d.tsv' % i, sep='\t')


def generateSynTxn(liverSamples,bloodSamples,seed):
    tissues = getSampleTissues()
    filename = '/usr/local/bioinfo/MS/tissueExpnDB/GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct.gz'
    allSamples = pd.read_csv(filename,sep='\t',skiprows=[0,1],nrows=3).columns

    selected_samples = []
    selected_fracs = {}
    tissue_fracs = {}
    allTissues = set(tissues.values())

    random.seed(a=seed)
    sample = random.choice(liverSamples)
    frac = random.randint(1,50)/1000.0
    left = 1.0 - frac
    selected_samples.append(sample)
    selected_fracs[sample] = frac
    tissue_fracs['Liver'] = frac

    random.seed(a=seed+23)
    sample = random.choice(bloodSamples)
    selected_samples.append(sample)
    selected_fracs[sample] = left
    tissue_fracs['Whole Blood'] = left
    '''
    allTissues.discard('Whole Blood')
    random.shuffle(list(allTissues))
    for tissue in allTissues:
        samples = getSamples(tissue,tissues,allSamples)
        if len(samples) < 10:
            continue
        seed += 2
        ceiling = min(50,left)
        if ceiling < 2:
            break
        random.seed(a=seed)
        sample = random.choice(samples)
        frac = random.randint(2,ceiling)/1000.0
        left -= int(frac*1000)
        selected_samples.append(sample)
        selected_fracs[sample] = frac
        tissue_fracs[tissue] = frac
    ''' 
    return selected_samples,selected_fracs,tissue_fracs
    

def featureSelect(df,cat):
    model = fs.SelectKBest(k=1000).fit(df.iloc[1:,:],cat[1:])
    df.loc['scores',:] = model.scores_
    df.loc['pvalues',:] = model.pvalues_

    coef_mat = df.loc[['Description','scores','pvalues'],:].transpose()
    coef_mat = coef_mat.sort_values('pvalues',ascending=True)
    return coef_mat.iloc[:50,:]

def LRclassify(data,cat):
    import sklearn.linear_model as skln
    lr = skln.LogisticRegression(penalty='l1',C=0.05).fit(data.iloc[1:,:],cat[1:])
    print lr.predict_proba(data.iloc[1:,:])

    data.loc['coef',:] = lr.coef_
    coef_mat = data.loc[['Description','coef'],].transpose()
    coef_mat['abs'] = abs(coef_mat['coef'])
    coef_mat = coef_mat[coef_mat['abs']>0].sort_values('abs',ascending=False)
    return coef_mat

def DecisionTreeClassify(data,cat):
    from sklearn.tree import DecisionTreeClassifier
    clf = DecisionTreeClassifier(random_state=19,max_depth=6)
    clf.fit(data.iloc[1:,:],cat[1:])
    data.loc['importance',:] = clf.feature_importances_

    coef_mat = data.loc[['Description','importance'],:].transpose()
    coef_mat = coef_mat.sort_values('importance',ascending=False)
    return coef_mat.iloc[:50,:]


def main():
    tissues = getSampleTissues()
    filename = '/usr/local/bioinfo/MS/tissueExpnDB/GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct.gz'
    allSamples = pd.read_csv(filename,sep='\t',skiprows=[0,1],nrows=3).columns

    selected_genes = []
    with open('/home/jzhuang@ms.local/GTEx/liverGenes','rU') as infile:
    #with open('/home/jzhuang@ms.local/GTEx/LR.geneList.txt','rU') as infile:
        for line in infile:
            selected_genes.append(line[:-1])
    #with open('/home/jzhuang@ms.local/GTEx/Ftest.geneList.txt','rU') as infile:
    #    for line in infile:
    #        selected_genes.append(line[:-1])

    print 'Generate random dataset...'
    liverSamples0 = getSamples('Liver',tissues,allSamples)
    bloodSamples0 = getSamples('Whole Blood',tissues,allSamples)
    random.seed(a=192)
    random.shuffle(liverSamples0)
    random.shuffle(bloodSamples0)
    sampleSet = set([])
    sampleList = []
    fracList = []
    liverFracList = []

    liverSamples = liverSamples0[:72]
    bloodSamples = bloodSamples0[:236]
    for i in range(120):
        out = generateSynTxn(liverSamples,bloodSamples,96+i*3)
        map(sampleSet.add,out[0])
        sampleList.append(out[0])
        fracList.append(out[1])
        liverFracList.append(out[2]['Liver'])

    liverSamples = liverSamples0[72:]
    bloodSamples = bloodSamples0[236:]
    for i in range(80):
        out = generateSynTxn(liverSamples,bloodSamples,319+i*3)
        map(sampleSet.add,out[0])
        sampleList.append(out[0])
        fracList.append(out[1])
        liverFracList.append(out[2]['Liver'])

    print 'Start reading table...'
    data = pd.read_csv(filename,sep='\t',skiprows=[0,1],usecols=['Name','Description']+list(sampleSet))
    data.iloc[:,2:] = data.iloc[:,2:].apply(lambda x: x/sum(x)*1000000, axis=0)
    data = data[map(lambda x: x in set(selected_genes),data['Description'])]
    mat = data[['Name','Description']]
    for i in range(200):
        colname = 'synthetic'+str(i)
        mat[colname] = 0.0
        for col in sampleList[i]:
            mat[colname] += data[col]*fracList[i][col]

    mat.index = mat['Name']
    mat = mat.drop('Name',axis=1)
    mat = mat.transpose()
    print mat.shape
    #mat['Real_frac'] = ['NA'] + liverFracList
    #mat = mat.transpose()
    #mat.to_csv('Simulated_liver_blood.tsv',delimiter='\t')

    print 'Start fitting...'
    from sklearn.metrics import mean_squared_error
    from math import sqrt
    #for i in range(2,31):
    #model = LinearRegression()
    model = Earth(max_degree=1,penalty=9)
    model.fit(mat.iloc[1:121,:],liverFracList[:120])
    out_mat = pd.DataFrame({'Real': liverFracList[120:],
                            'Predicted': model.predict(mat.iloc[121:,:]),
                            'ALB': mat['ENSG00000163631.12'][121:]/ALB_median
                        })
    rms1 = sqrt(mean_squared_error(out_mat['Real'],out_mat['Predicted']))
    rms2 = sqrt(mean_squared_error(out_mat['Real'],out_mat['ALB']))
    #print i,rms1

    fig, ax = plt.subplots()
    out_mat.plot.scatter(x=2,y=1,ax=ax)
    plt.xlabel('Real Perc')
    plt.ylabel('Predicted Perc')
    plt.xlim(-0.01,0.06)
    plt.ylim(-0.01,0.06)
    plt.title('MARS predictions on simulated dataset\n'+str(rms1))
    plt.savefig('pred_vs_real.pdf')
    plt.close()
    '''
    fig, ax = plt.subplots()
    out_mat.plot.scatter(x=2,y=0,color='red',ax=ax)
    plt.xlabel('Real Perc')
    plt.ylabel('Predicted Perc')
    plt.title('ALB predictions on simulated dataset\n'+str(rms2))
    plt.savefig('ALB_vs_real.pdf')
    plt.close()
    '''

if __name__=='__main__':
    main()
