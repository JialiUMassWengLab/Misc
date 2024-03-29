#! /usr/bin/python

import re
import os
import sys
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn import preprocessing
from sklearn import model_selection
from sklearn import linear_model


def getOrganSpecificGenes(organ):
    df = pd.read_csv('references/bulk-gex_v8_rna-seq_GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz',sep='\t',index_col=1,skiprows=2)
    organDf = pd.read_csv('references/organ_designation.tsv',sep='\t',index_col=0)
    df = df.transpose()
    df = df.loc[map(lambda x: re.subn(r'[\s\(\)-]','.',x)[0] in list(organDf.index),df.index),:]
    df.loc[:,'organ'] = list(map(lambda x: organDf.loc[re.subn(r'[\s\(\)-]','.',x)[0],'organ'], df.index))
    df = df.groupby('organ').agg(max).transpose()
    #df = df.loc[df.apply(lambda x: x[organ] > 5 and len([i for i in x if i/x[organ] >= .25])==1,axis=1),:]
    df = df.loc[df.apply(lambda x: len([i for i in x if (i+1)/(x[organ]+1) >= .25])==1,axis=1),:]
    df[organ].to_csv('references/%s_specific_genes.csv' % organ)
    print(df)
    

def UKB_organ_specific_genes(prefix,organ,plotCorr=False):
    import scipy.stats as st
    geneDf = pd.read_csv('references/%s_specific_genes.csv' % organ,index_col=0)
    '''
    #distribution among WGCNA modules
    df = pd.read_csv('WGCNA_results/%s.modules.genes.txt' % prefix,sep=' ',index_col=0)
    df = df.loc[map(lambda x: x.split('_')[0] in list(geneDf.index),df.index),:]
    print(df.melt().groupby('variable').sum())
    '''
    df = pd.read_csv('combined_pheno_forconsortium_v1_NPX.tsv',sep='\t',index_col=0,low_memory=False)
    df = df.loc[:,map(lambda x: x.split(':')[0] in list(geneDf.index),df.columns)]
    #df.to_csv('%s_specific_genes.GTEx.NPX_combined.tsv' % organ, sep='\t')

    if not plotCorr:
        return df

    moduleDf = pd.read_csv(os.path.join('WGCNA_results',prefix+'.moduleLevels.csv'),index_col=0)
    array = []
    for protein in df.columns:
        tmpDf = pd.concat([moduleDf,df[protein]],join='inner',axis=1)
        data = {'rho':{},'pv':{}}
        for module in tmpDf.columns[:-1]:
            tmpDf2 = tmpDf[[module,protein]].dropna(axis=0)
            rho,pv = st.spearmanr(tmpDf2[module],tmpDf2[protein])
            i = module.split('_')[-1]+' - '+protein
            data['rho'].update({i: rho})
            data['pv'].update({i: pv})
        resultDf = pd.DataFrame.from_dict(data)
        array.append(resultDf)

    df = pd.concat(array,join='inner',axis=0)
    df['protein'] = df.apply(lambda x: x.name.split(' - ')[1],axis=1)
    df['module'] = df.apply(lambda x: x.name.split(' - ')[0],axis=1)
    df = df.pivot(columns='module',index='protein',values='rho')
    print(df)

    #fig,ax = plt.subplots(figsize=(12,22))
    #sns.heatmap(df,cmap="RdBu_r",annot=False,linecolor='k',linewidths=.5,square=True,ax=ax)
    sns.set(font_scale=.8)
    #sns.clustermap(df,cmap="RdBu_r",metric='correlation',z_score=None,standard_scale=None,row_cluster=True,col_cluster=False,linecolor='k',linewidths=.8,yticklabels=True,figsize=(12,22))
    plt.rcParams['xtick.bottom'] = True
    df.index = df.apply(lambda x: x.name.split(':')[0],axis=1)
    g = sns.clustermap(df.transpose(),cmap="RdBu_r",metric='correlation',z_score=None,standard_scale=None,row_cluster=False,col_cluster=True,linecolor='k',linewidths=.8,yticklabels=True,figsize=(26,10))
    plt.setp(g.ax_heatmap.xaxis.get_majorticklabels(), rotation=40)
    plt.setp(g.ax_heatmap.xaxis.get_majorticklabels(), ha='right')
    plt.tight_layout()
    plt.savefig(prefix+'.%s_VS_module.rho.pdf' % organ)
    plt.close()
    return None
    

def UKB_module_specific_genes(prefix,module):
    geneDf = pd.read_csv('WGCNA_results/%s.modules.genes.txt' % prefix,index_col=0,sep=' ')
    geneDf.index = geneDf.index.map(lambda x: x.split('_')[0])
    geneDf.columns = geneDf.columns.map(lambda x: x.split('_')[-1])
    geneDf = geneDf[geneDf[module]==1]
    df = pd.read_csv('combined_pheno_forconsortium_v1_NPX.tsv',sep='\t',index_col=0,low_memory=False)
    df = df.loc[:,map(lambda x: x.split(':')[0] in list(geneDf.index),df.columns)]
    return df


def splitTrainTest(df,seed=523):
    df = df.dropna(axis=0)
    info = pd.read_csv('sampleInfo.csv',index_col=0)
    info = info[info.index.isin(df.index)]
    info = info.loc[df.index,:]
    X_train,X_test,y_train,y_test = model_selection.train_test_split(df,info['Age'],test_size=0.25,random_state=seed)
    return X_train,X_test,y_train,y_test
    

def crossValidation(X,y):
    lasso = linear_model.Lasso()
    params = {'alpha':[1e-6,1e-5,5e-5,1e-4,5e-4,1e-3,5e-3,0.01,0.02,0.05,0.1,0.2,0.5,1,2,10]}
    reg = model_selection.GridSearchCV(lasso,params,cv=5,scoring='r2')
    reg.fit(X,y)
    print(reg.cv_results_)
    
    param_series = pd.Series(reg.cv_results_['param_alpha'],name='alpha')
    mean_series = pd.Series(reg.cv_results_['mean_test_score'],name='mean_r2')
    df = pd.concat([param_series,mean_series],axis=1)
    fig,ax = plt.subplots(figsize=(12,10))
    #sns.scatterplot(x='alpha',y='mean_r2',data=df,ax=ax)
    df.plot(x='alpha',y='mean_r2',kind='scatter',logx=True,ax=ax)
    plt.savefig('crossValid.results.pdf')
    plt.close()
    

def predOnTest(X_train,X_test,y_train,y_test,alpha,geneSet):
    subjectIDs = X_test.index
    lasso = linear_model.Lasso(alpha=alpha)
    lasso.fit(X_train,y_train)
    y_pred = lasso.predict(X_test)
    age_array = np.stack((y_test,y_pred),axis=-1)
    
    import statsmodels.api as sm
    df = pd.DataFrame(age_array,columns=['age','pred_age'],index=subjectIDs)
    smoothed = sm.nonparametric.lowess(exog=df['age'], endog=df['pred_age'], frac=2./3, return_sorted=False)
    fig,ax = plt.subplots(figsize=(12,10))
    sns.scatterplot(x='age',y='pred_age',data=df,color='darkgrey',alpha=.7,s=10,ax=ax)
    #sns.regplot(x='age',y='pred_age',data=df,color='darkgrey',scatter_kws={'alpha':.7,'s':20},ax=ax)
    ax.plot(df['age'], smoothed, c="k")
    plt.savefig('chronological_age_pred.%s.png' % geneSet)
    plt.close()

    df.loc[:,'lowess_smoothed'] = smoothed
    df.loc[:,'age_gap'] = df.apply(lambda x: x['pred_age']-x['lowess_smoothed'],axis=1)
    gap_mean, gap_std = np.mean(df['age_gap']), np.std(df['age_gap'])
    print(gap_mean,gap_std)
    df.loc[:,'z_score'] = df.apply(lambda x: (x['age_gap']-gap_mean)/gap_std,axis=1)
    df.to_csv('age_gap.%s.csv' % geneSet)
    

def testEnrichment(geneSet):
    from scipy.stats import hypergeom
    df = pd.read_csv('age_gap.%s.csv' % geneSet,index_col=0)
    diag = pd.read_csv('sampleDiag.csv',index_col=0).dropna()
    diag = diag.loc[map(lambda x: x in list(df.index),diag.index),:]
    #diseaseList = ['G10','G20','G21','G23','G24','G25','G30','G31','G32','G35','G36','G37','G40','G43','G44','G46','G47','G93','G98']
    diseaseList = ['G10','G20','G21','G30','G31','G32','G35','G36','G37','G40','G43','G44','G47','G52','G53','G54','G55']
    patientList = []
    for eid,row in diag.iterrows():
        diseases = row['diagnosis'].split(',')
        i = 0
        matched = False
        while not matched and i < len(diseases):
            if diseases[i][:3] in diseaseList:
                patientList.append(eid)
                matched = True
            i += 1

    df.loc[:,'affected'] = df.apply(lambda x: True if x.name in patientList else False,axis=1)
    k = df.query('affected == True & z_score > 2').shape[0]
    n = df.query('z_score > 2').shape[0]
    N = df.query('affected == True').shape[0]
    pv = 1 - hypergeom.cdf(k,df.shape[0],n,N)
    print(k,df.shape[0],n,N)
    print(pv)
    

def main():
    prefix = sys.argv[1]
    #getOrganSpecificGenes('Kidney')
    #df = UKB_organ_specific_genes(prefix,'Brain',True)
    
    #############################
    ##GTEx brain specific genes##
    #############################
    #getOrganSpecificGenes('Brain')
    #df = UKB_organ_specific_genes(prefix,'Brain',False)
    #X_train,X_test,y_train,y_test = splitTrainTest(df)
    #crossValidation(X_train,y_train)
    #predOnTest(X_train,X_test,y_train,y_test,0.05,'Brain')
    #testEnrichment('Brain')
    
    ###############################
    ##WGCNA module specific genes##
    ###############################
    #df = UKB_module_specific_genes(prefix,'module19')
    #X_train,X_test,y_train,y_test = splitTrainTest(df)
    #print(X_train.shape, X_test.shape)
    #crossValidation(X_train,y_train)
    #predOnTest(X_train,X_test,y_train,y_test,0.05,'module19')
    #testEnrichment('module19')
    
    df = UKB_module_specific_genes(prefix,'module14')
    X_train,X_test,y_train,y_test = splitTrainTest(df)
    #print(X_train.shape, X_test.shape)
    crossValidation(X_train,y_train)
    #predOnTest(X_train,X_test,y_train,y_test,0.05,'module19')

    
if __name__=='__main__':
    main()
