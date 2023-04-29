#! /usr/bin/python

import re
import os
import sys
import glob
import functools
import scipy
import numpy as np
import pandas as pd
import scipy.stats as ss
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn import preprocessing
from sklearn import model_selection
from sklearn import ensemble
from sklearn import metrics
from sklearn import linear_model
from sklearn import svm
from sklearn import neighbors
from sklearn import pipeline


def PRAUC_score(y_true,y_pred):
    precision, recall, _ = metrics.precision_recall_curve(y_true,y_pred)
    pr_auc = metrics.auc(recall, precision)
    return pr_auc


def splitTrainTest(table,target,standardize,seed):
    X_train,X_test,y_train,y_test = model_selection.train_test_split(table.iloc[:,:-3].astype(float),table[target],\
                                                                     test_size=0.4,random_state=seed)
    if standardize:
        X_train = preprocessing.StandardScaler().fit_transform(X_train)
        X_test = preprocessing.StandardScaler().fit_transform(X_test)

    names = pd.read_csv('../2018DecAD/1812ADcohort.congregated.tsv',sep='\t',index_col=0,usecols=['gene_id','Description'])
    table = pd.concat([table,names.transpose()],join='inner',axis=0)
    return X_train,X_test,y_train,y_test,table.iloc[-1,:]
    

def readMatrix(fn,pheno):
    tbl = pd.read_csv('/home/jzhuang@ms.local/files_for_Jerome/AD/bootstrap/summary.csv',index_col=0)
    geneList = list(tbl[(tbl['padj']<.05) & ((tbl['Prob_boot_lt_null']<.05) | (tbl['Prob_boot_gt_null']<.05))].index)
    geneList = list(map(lambda x: re.sub(r'\.\d+$','',x),geneList))

    df = pd.read_csv(fn,index_col=0,low_memory=False).transpose()
    #selectList = list(df.columns[map(lambda x: x in geneList,df.columns)])
    df1 = df.iloc[:,:-3].astype(float)
    selectList = list(df1.loc[:,df1.apply(lambda x: len([i for i in x if i>0]) > df1.shape[0]*0.8,axis=0)].columns)
    df = df.loc[:,selectList + ['PatientID','Center','Disease']]
    df.iloc[:,:-3] = df.iloc[:,:-3].astype('float')
    df = df[df['Center']=='UCSD1']
    if pheno == 'Disease':
        df = df[(df['Disease']=='AD') | (df['Disease']=='NCI')]
        df['Disease'] = map(lambda x: 1 if x=='AD' else 0, df['Disease'])
        
    return df


def RandomForest(X_train,X_test,y_train,y_test,n0,m0):
    rf = ensemble.RandomForestClassifier(min_samples_split=m0,n_estimators=n0,random_state=119)
    rf.fit(X_train,y_train)
    results = pd.concat([y_test,pd.Series(rf.predict(X_test),index=y_test.index,name='Predicted')],axis=1)
    #print results
    y_score = rf.predict_proba(X_test)
    fpr, tpr, _ = metrics.roc_curve(y_test,y_score[:,1])
    roc_auc = metrics.auc(fpr, tpr)
    precision, recall, _ = metrics.precision_recall_curve(y_test,y_score[:,1])
    pr_auc = metrics.auc(recall, precision)
    return fpr,tpr,roc_auc,precision,recall,pr_auc,y_score[:,1]

def RandomForest2(X,y,n0,m0,standardize):
    rf = ensemble.RandomForestClassifier(min_samples_split=m0,n_estimators=n0,random_state=119)
    clf = pipeline.make_pipeline(preprocessing.StandardScaler(),rf) if standardize else rf
    #roc_auc = model_selection.cross_val_score(clf,X,y,cv=5,scoring='roc_auc')
    pr_auc = model_selection.cross_val_score(clf,X,y,cv=5,scoring=metrics.make_scorer(PRAUC_score,needs_proba=True))
    pred = model_selection.cross_val_predict(clf,X,y,cv=5,method='predict_proba')
    return pr_auc,pred[:,1]


def RidgeClassifier(X_train,X_test,y_train,y_test,a0):
    rc = linear_model.LogisticRegression(penalty='l2',C=a0,random_state=254,fit_intercept=False)
    rc.fit(X_train,y_train)
    results = pd.concat([y_test,pd.Series(rc.predict(X_test),index=y_test.index,name='Predicted')],axis=1)
    #print results
    y_score = rc.predict_proba(X_test)
    fpr, tpr, _ = metrics.roc_curve(y_test,y_score[:,1])
    roc_auc = metrics.auc(fpr, tpr)
    precision, recall, _ = metrics.precision_recall_curve(y_test,y_score[:,1])
    pr_auc = metrics.auc(recall, precision)
    return fpr,tpr,roc_auc,precision,recall,pr_auc,y_score[:,1]

def RidgeClassifier2(X,y,a0,standardize):
    rc = linear_model.LogisticRegression(penalty='l2',C=a0,random_state=254,fit_intercept=False)
    clf = pipeline.make_pipeline(preprocessing.StandardScaler(),rc) if standardize else rc
    #roc_auc = model_selection.cross_val_score(clf,X,y,cv=5,scoring='roc_auc')
    pr_auc = model_selection.cross_val_score(clf,X,y,cv=5,scoring=metrics.make_scorer(PRAUC_score,needs_proba=True))
    pred = model_selection.cross_val_predict(clf,X,y,cv=5,method='predict_proba')
    return pr_auc,pred[:,1]


def plotBox(resultFile,classifier):
    data = {}
    fig,ax = plt.subplots(figsize=(14,10))
    df = pd.read_csv(resultFile,index_col=0)
    df.plot(kind='box',legend=False,title=classifier,fontsize=16,ax=ax)
    plt.setp(ax.get_xticklabels(), ha="right", rotation=45)
    plt.ylim(0.1,1.05)
    plt.ylabel('AUC')
    plt.xlabel('Parameter')
    plt.savefig(resultFile.replace('csv','png'))
    plt.close()


def plotPR(y_test,y_pred,paras,pdf):
    precision, recall, _ = metrics.precision_recall_curve(y_test,y_pred)
    pr_auc = metrics.auc(recall, precision)
    no_skill = y_test[y_test==1].shape[0] / float(y_test.shape[0])

    fig,ax = plt.subplots(figsize=(15,12))
    ax.plot(recall,precision,lw=3,label='Precision-recall curve (area = %0.3f)' % pr_auc,color='black')
    ax.plot([0,1],[no_skill,no_skill],color='black',lw=2,ls='--')
    plt.xlim([-.02,1])
    plt.ylim([0,1.05])
    plt.title(str(paras))
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.legend(loc='best')
    plt.savefig(pdf,format='pdf')
    plt.close()


def main():
    pheno = sys.argv[1]
    dtype = sys.argv[2]
    classifier = sys.argv[3]
    standardize = False
    df = readMatrix('ADDF_run1.normalized.%s.csv' % dtype,pheno)
    print(df.shape)

    #geneset='high_conf_DEGs'
    geneset='allGenes'
    sns.set(font_scale=1.6)
    '''
    if classifier == 'RC':
        #Cvalues = [1,5,10,50,100,200,500,1000,2000,5000]
        Cvalues = [0.1,0.2,0.5,1,5,10,50,100,200,500,1000,2000]
        data = {}
        for i in range(len(Cvalues)):
            AUCs = []
            C0 = Cvalues[i]
            for j in range(15):
                X_train,X_test,y_train,y_test,description = splitTrainTest(df,pheno,standardize,129+j*3)
                AUCs.append(RidgeClassifier(X_train,X_test,y_train,y_test,C0)[2])
            data.update({C0: AUCs})
        df1 = pd.DataFrame.from_dict(data)
        print(df1)
        dtype = dtype+'z' if standardize else dtype
        df1.to_csv('RidgeClf.%s.%s.%s.csv' % (geneset,pheno,dtype))
        plotBox('RidgeClf.%s.%s.%s.csv' % (geneset,pheno,dtype),'Logistic Regression L2')
    '''
    if classifier == 'RC':
        #Cvalues = [0.1,0.2,0.5,1,5,10,50,100,200,500,1000,2000]
        Cvalues = [1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,0.0001,0.001,0.01,0.1,1,10,100]
        if standardize:
            dtype = dtype+'z'
        AUCs = {}
        from matplotlib.backends.backend_pdf import PdfPages
        pdf = PdfPages('RidgeClf.%s.%s.%s.CV_PRcurve.pdf' % (geneset,pheno,dtype))
        for i in range(len(Cvalues)):
            C0 = Cvalues[i]
            AUC,pred = RidgeClassifier2(df.iloc[:,:-3],df[pheno],C0,standardize)
            AUCs.update({C0: AUC})
            plotPR(df[pheno],pred,C0,pdf)
        df1 = pd.DataFrame.from_dict(AUCs)
        print(df1)
        pdf.close()
        df1.to_csv('RidgeClf.%s.%s.%s.csv' % (geneset,pheno,dtype))
        plotBox('RidgeClf.%s.%s.%s.csv' % (geneset,pheno,dtype),'Logistic Regression L2')
    '''
    if classifier =='RF':
        ntrees = [10,20,50,80,100,250,500]
        leafSampleSizes = [2,5,8,10]
        data = {}
        for n in ntrees:
            for m in leafSampleSizes:
                AUCs = []
                for j in range(15):
                    X_train,X_test,y_train,y_test,description = splitTrainTest(df,pheno,standardize,129+j*3)
                    AUCs.append(RandomForest(X_train,X_test,y_train,y_test,n,m)[2])
                    data.update({'%d_%d' % (n,m): AUCs})
        df1 = pd.DataFrame.from_dict(data)
        print(df1)
        dtype = dtype+'z' if standardize else dtype
        df1.to_csv('RandomForest.%s.%s.%s.csv' % (geneset,pheno,dtype))
        plotBox('RandomForest.%s.%s.%s.csv' % (geneset,pheno,dtype),'Random Forest')
    '''    
    if classifier == 'RF':
        ntrees = [10,20,50,80,100,250,500]
        leafSampleSizes = [2,5,8,10]
        if standardize:
            dtype = dtype+'z'
        AUCs = {}
        from matplotlib.backends.backend_pdf import PdfPages
        pdf = PdfPages('RandomForest.%s.%s.%s.CV_PRcurve.pdf' % (geneset,pheno,dtype))
        for n in ntrees:
            for m in leafSampleSizes:
                AUC,pred = RandomForest2(df.iloc[:,:-3],df[pheno],n,m,standardize)
                AUCs.update({'%d_%d' % (n,m): AUC})
                plotPR(df[pheno],pred,'%d_%d' % (n,m),pdf)
        df1 = pd.DataFrame.from_dict(AUCs)
        print(df1)
        pdf.close()
        df1.to_csv('RandomForest.%s.%s.%s.csv' % (geneset,pheno,dtype))
        plotBox('RandomForest.%s.%s.%s.csv' % (geneset,pheno,dtype),'Random Forest')



if __name__=='__main__':
    main()
