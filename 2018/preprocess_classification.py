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
from compareRep import normERCC
from sklearn import preprocessing
from sklearn import model_selection
from sklearn import ensemble
from sklearn import metrics
from sklearn import linear_model
from sklearn import svm
from sklearn import neighbors
from sklearn import pipeline

def splitTrainTest(table,target,seed):
    X_train,X_test,y_train,y_test = model_selection.train_test_split(table.iloc[:,:-3],table[target],\
                                                                     test_size=0.4,random_state=seed)
    names = pd.read_csv('1802NASHserum.congregated.tsv',sep='\t',index_col=0,usecols=['gene_id','Description'])
    table = pd.concat([table,names.transpose()],join='inner',axis=0)
    return X_train,X_test,y_train,y_test,table.iloc[-1,:]
    

def readMatrix(fn,pheno):
    df = pd.read_csv(fn,index_col=0,low_memory=False).transpose()
    df.iloc[:,:-3] = df.iloc[:,:-3].astype('float')
    df = df.loc[:,df.apply(lambda x: len([i for i in x if i > 10]) > df.shape[0] * 0.8, axis=0)]
    if pheno == 'PathFibrosis':
        df = df[(df['Disease']=='NASH') | (df['Disease']=='NAFLD')]
        df = df[(df['PathFibrosis']=='F0') | (df['PathFibrosis']=='F3') | (df['PathFibrosis']=='F4')]
        df['PathFibrosis'] = map(lambda x: 1 if x=='F3' or x=='F4' else 0, df['PathFibrosis'])
    elif pheno == 'Disease':
        df = df[(df['Disease']=='NASH') | (df['Disease']=='Normal Control')]
        df = df[(df['PathFibrosis']=='F0') | (df['PathFibrosis']=='F1') | (df['PathFibrosis']=='None')]
        df['Disease'] = map(lambda x: 1 if x=='NASH' else 0, df['Disease'])
    '''
    selected_genes = set([])
    with open('/home/jzhuang@ms.local/CoExpression/Serum_vs_Plasma/nonBlood.csv','rU') as infile:
        for line in infile:
            a = line.split(',')
            selected_genes.append(a[0])
    panel_genes = {}
    with open('../qPCR_genes.tsv','rU') as infile:
        for line in infile:
            a = line.split()
            panel_genes[a[0]] = a[1]
    pharma_genes = set([])
    with open('GenesisCombinedGeneList.csv','rU') as infile:
        next(infile)
        for line in infile:
            pharma_genes.add(line.split(',')[0])
    #df = df.loc[:,map(lambda x: x in selected_genes,df.loc['Description',:])]
    #df = df.loc[:,map(lambda x: x in panel_genes,df.loc['Description',:])]
    df = pd.concat([df.loc[:,map(lambda x: re.sub('\.\d+$','',x) in selected_genes,df.columns)],df.iloc[:,-3:]],join='inner',axis=1)
    df.iloc[:-1,:-3] = df.iloc[:-1,:-3].apply(lambda x: x/sum(x)*1000000,axis=1)
    print df.shape
    '''
    return df

def readMatrix2(fn,pheno):
    #df = pd.read_csv(fn,index_col=1,low_memory=False).drop(['Unnamed: 0'],axis=1).transpose()
    df = pd.read_csv(fn,index_col=0,low_memory=False).transpose()
    df = df.astype('float')
    info = pd.read_csv('sampleInfo.csv',index_col=0)
    info.index = info.index.astype('str')
    df = pd.concat([df,info],join='inner',axis=1)
    if pheno == 'PathFibrosis':
        df = df[(df['Disease']=='NASH') | (df['Disease']=='NAFLD')]
        df = df[(df['PathFibrosis']=='F0') | (df['PathFibrosis']=='F3') | (df['PathFibrosis']=='F4')]
        df['PathFibrosis'] = map(lambda x: 1 if x=='F3' or x=='F4' else 0, df['PathFibrosis'])
    elif pheno == 'Disease':
        df = df[(df['Disease']=='NASH') | (df['Disease']=='Normal Control')]
        df = df[(df['PathFibrosis']=='F0') | (df['PathFibrosis']=='F1') | (df['PathFibrosis']=='None')]
        df['Disease'] = map(lambda x: 1 if x=='NASH' else 0, df['Disease'])
    elif pheno == 'CTRLvsF34':
        df = df[(df['Disease']=='Normal Control') | (df['PathFibrosis']=='F3') | (df['PathFibrosis']=='F4')]
        df['CTRLvsF34'] = map(lambda x: 0 if x=='Normal Control' else 1, df['Disease'])
        df.drop('PatientID',axis=1,inplace=True)
    elif pheno == 'CTRLvsNAFLD':
        df = df[(df['Disease']=='Normal Control') | (df['Disease']=='NAFLD')]
        df['CTRLvsNAFLD'] = map(lambda x: 1 if x=='NAFLD' else 0, df['Disease'])
        df.drop('PatientID',axis=1,inplace=True)
    elif pheno == 'NAFLDvsNASH':
        df = df[(df['Disease']=='NASH') | (df['Disease']=='NAFLD')]
        df['NAFLDvsNASH'] = map(lambda x: 0 if x=='NAFLD' else 1, df['Disease'])
        df.drop('PatientID',axis=1,inplace=True)
    elif pheno == 'NAFLDvsLowNASH':
        df = df[((df['Disease']=='NASH')&((df['PathFibrosis']=='F0')|(df['PathFibrosis']=='F1'))) | (df['Disease']=='NAFLD')]
        df['NAFLDvsLowNASH'] = map(lambda x: 0 if x=='NAFLD' else 1, df['Disease'])
        df.drop('PatientID',axis=1,inplace=True)
    elif pheno == 'NAFLDvsHighNASH':
        df = df[((df['Disease']=='NASH')&((df['PathFibrosis']=='F3')|(df['PathFibrosis']=='F4'))) | (df['Disease']=='NAFLD')]
        df['NAFLDvsHighNASH'] = map(lambda x: 0 if x=='NAFLD' else 1, df['Disease'])
        df.drop('PatientID',axis=1,inplace=True)
    elif pheno == 'F0vsF4':
        df = df[(df['PathFibrosis']=='F0') | (df['PathFibrosis']=='F4')]
        df['F0vsF4'] = map(lambda x: 0 if x=='F0' else 1, df['PathFibrosis'])
        df.drop('PatientID',axis=1,inplace=True)
    elif pheno == 'F01vsF34':
        df = df[(df['PathFibrosis']!='F2') & (df['PathFibrosis']!='None')]
        df['F01vsF34'] = map(lambda x: 0 if x=='F0' or x=='F1' else 1, df['PathFibrosis'])
        df.drop('PatientID',axis=1,inplace=True)
    elif pheno == 'F012vsF34':
        df = df[df['PathFibrosis']!='None']
        df['F012vsF34'] = map(lambda x: 1 if x=='F3' or x=='F4' else 0, df['PathFibrosis'])
        df.drop('PatientID',axis=1,inplace=True)
    elif pheno == 'F01vsF23':
        df = df[(df['PathFibrosis']!='F4') & (df['PathFibrosis']!='None')]
        df['F01vsF23'] = map(lambda x: 0 if x=='F0' or x=='F1' else 1, df['PathFibrosis'])
        df.drop('PatientID',axis=1,inplace=True)
        
    return df


def getMatrix(table):
    df = pd.read_csv(table,sep='\t',index_col=0,low_memory=False)
    #names = df['Description']
    df = df.loc[map(lambda x: not re.search(r'^ERCC-',x),df.index),:]
    df.iloc[:,:-2] = df.iloc[:,:-2].apply(lambda x: x/sum(x)*1000000,axis=0)
    #df = df.apply(lambda x: 2**x,axis=0)
    df = df.transpose()
    #df = normERCC(df.iloc[:,:-2].transpose())

    df.index = map(lambda x: x.replace('_','-').split('-')[0], df.index)
    withRepList = []
    for name,group in df.groupby(df.index):
        if group.shape[0] >= 2:
            withRepList.append(name)
    df = df[map(lambda x: x in withRepList,df.index)].astype('float')
    df = df.groupby(df.index).agg(np.median)
    print df.shape

    info = pd.read_csv('sampleInfo.csv',index_col=0)
    info.index = info.index.astype('str')
    #info = info[~(info['Source'].astype('str')=='nan')]
    table = pd.concat([ df,info ],axis=1,join='inner')
    #table = table.append(names)
    return table


def RandomForest(X_train,X_test,y_train,y_test,n0,m0):
    rf = ensemble.RandomForestClassifier(min_samples_split=m0,n_estimators=n0,random_state=119)
    rf.fit(X_train,y_train)
    results = pd.concat([y_test,pd.Series(rf.predict(X_test),index=y_test.index,name='Predicted')],axis=1)
    #print results
    y_score = rf.predict_proba(X_test)
    fpr, tpr, _ = metrics.roc_curve(y_test,y_score[:,1])
    roc_auc = metrics.auc(fpr, tpr)
    return fpr,tpr,roc_auc,y_score[:,1]


def GBClassifier(X_train,X_test,y_train,y_test,n0,d0):
    gb = ensemble.GradientBoostingClassifier(random_state=91,max_depth=d0,n_estimators=n0,learning_rate=2/float(n0))
    gb.fit(X_train,y_train)
    results = pd.concat([y_test,pd.Series(gb.predict(X_test),index=y_test.index,name='Predicted')],axis=1)
    #print results
    y_score = gb.predict_proba(X_test)
    fpr, tpr, _ = metrics.roc_curve(y_test,y_score[:,1])
    roc_auc = metrics.auc(fpr, tpr)
    return fpr,tpr,roc_auc,y_score[:,1]


def KNNClassifier(X_train,X_test,y_train,y_test,K):
    knn = neighbors.KNeighborsClassifier(n_neighbors=K)
    knn.fit(X_train,y_train)
    results = pd.concat([y_test,pd.Series(knn.predict(X_test),index=y_test.index,name='Predicted')],axis=1)
    #print results
    y_score = knn.predict_proba(X_test)
    fpr, tpr, _ = metrics.roc_curve(y_test,y_score[:,1])
    roc_auc = metrics.auc(fpr, tpr)
    return fpr,tpr,roc_auc,y_score[:,1]


def LRClassifier(X_train,X_test,y_train,y_test,C0):
    #lr = pipeline.Pipeline([('scaling',preprocessing.StandardScaler()),\
    #                        ('classifier',linear_model.LogisticRegression(penalty='l1',C=C0,random_state=254))])
    X_train = preprocessing.StandardScaler().fit_transform(X_train)
    X_test = preprocessing.StandardScaler().fit_transform(X_test)
    lr = linear_model.LogisticRegression(penalty='l1',C=C0,random_state=254,fit_intercept=False)
    lr.fit(X_train,y_train)
    print len([i for i in lr.coef_[0,:] if not i==0])
    results = pd.concat([y_test,pd.Series(lr.predict(X_test),index=y_test.index,name='Predicted')],axis=1)
    #print results
    y_score = lr.predict_proba(X_test)
    fpr, tpr, _ = metrics.roc_curve(y_test,y_score[:,1])
    roc_auc = metrics.auc(fpr, tpr)
    return fpr,tpr,roc_auc,y_score[:,1]


def RidgeClassifier(X_train,X_test,y_train,y_test,a0):
    X_train = preprocessing.StandardScaler().fit_transform(X_train)
    X_test = preprocessing.StandardScaler().fit_transform(X_test)
    #rc = linear_model.RidgeClassifier(alpha=a0,fit_intercept=True,random_state=254)
    rc = linear_model.LogisticRegression(penalty='l2',C=a0,random_state=254,fit_intercept=False)
    rc.fit(X_train,y_train)
    results = pd.concat([y_test,pd.Series(rc.predict(X_test),index=y_test.index,name='Predicted')],axis=1)
    #print results
    y_score = rc.predict_proba(X_test)
    fpr, tpr, _ = metrics.roc_curve(y_test,y_score[:,1])
    roc_auc = metrics.auc(fpr, tpr)
    return fpr,tpr,roc_auc,y_score[:,1]


def SVMClassifier(X_train,X_test,y_train,y_test,C0):
    X_train = preprocessing.StandardScaler().fit_transform(X_train)
    X_test = preprocessing.StandardScaler().fit_transform(X_test)
    svc = svm.SVC(C=C0,random_state=254,tol=0.00001,probability=True)
    svc.fit(X_train,y_train)
    results = pd.concat([y_test,pd.Series(svc.predict(X_test),index=y_test.index,name='Predicted')],axis=1)
    #print results
    y_score = svc.predict_proba(X_test)
    fpr, tpr, _ = metrics.roc_curve(y_test,y_score[:,1])
    roc_auc = metrics.auc(fpr, tpr)
    return fpr,tpr,roc_auc,y_score[:,1]



def plotLR1(resultFile,classifier):
    data = {}
    with open(resultFile,'rU') as infile:
        for line in infile:
            a = line[:-1].split()
            data[a[0]] = map(float,a[1].split(','))
        df = pd.DataFrame(data)
        df.plot(kind='box',legend=False,title=classifier,figsize=(14,10),fontsize=16)
        plt.ylim(0.1,1)
        plt.ylabel('AUC')
        plt.xlabel('C')
        plt.savefig(resultFile.replace('tsv','png'))
        plt.close()

def plotRandomForest1(resultFile,classifier):
    data = {}
    with open(resultFile,'rU') as infile:
        for line in infile:
            a = line[:-1].split()
            data[a[0]+'_'+a[1]] = map(float,a[2].split(','))
        df = pd.DataFrame(data)
        df.plot(kind='box',legend=False,title=classifier,figsize=(14,10),fontsize=16)
        plt.ylim(0.1,1)
        plt.ylabel('AUC')
        plt.xlabel('Para')
        plt.xticks(fontsize=15,rotation=45,ha='right')
        plt.savefig(resultFile.replace('tsv','png'))
        plt.close()


def main():
    #df = getMatrix('1802NASHserum.congregated.tsv')
    #df.transpose().to_csv('180305_NASHvalidation_203_avg_TPM.csv')

    pheno = sys.argv[1]
    dtype = sys.argv[2]
    classifier = sys.argv[3]
    #df = readMatrix('180305_NASHvalidation_146_avg_%s.csv' % dtype,pheno)
    df = readMatrix2('/mnt/shares/R&D/projects/201803_NASH_Validation/Data/All_data/20180305_NASH_Validation_204_All_%s_cutoff.csv' % dtype,pheno)
    print df.shape

    sns.set(font_scale=1.6)
    if classifier == 'LR':
        #Cvalues = [0.01,0.05,0.1,0.2,0.5,0.75,1,2.5,5,10]
        Cvalues = [0.2,0.5,1,5,10,50,100,200,500,1000,2000]
        out = open('LogisticRegre.allGenes.%s.%s.std.tsv' % (pheno,dtype),'w')
        for i in range(len(Cvalues)):
            AUCs = []
            C0 = Cvalues[i]
            for j in range(15):
                X_train,X_test,y_train,y_test,description = splitTrainTest(df,pheno,129+j*3)
                AUCs.append(LRClassifier(X_train,X_test,y_train,y_test,C0)[2])
            print AUCs
            out.write(str(C0)+'\t'+','.join(map(str,AUCs))+'\n')
        out.close()
        plotLR1('LogisticRegre.allGenes.%s.%s.std.tsv' % (pheno,dtype),'Logistic Regression')

    if classifier == 'RC':
        #Cvalues = [1,5,10,50,100,200,500,1000,2000,5000]
        Cvalues = [0.2,0.5,1,5,10,50,100,200,500,1000,2000]
        out = open('RidgeClf.allGenes.%s.%s.tsv' % (pheno,dtype),'w')
        for i in range(len(Cvalues)):
            AUCs = []
            C0 = Cvalues[i]
            for j in range(15):
                X_train,X_test,y_train,y_test,description = splitTrainTest(df,pheno,129+j*3)
                AUCs.append(RidgeClassifier(X_train,X_test,y_train,y_test,C0)[2])
            print AUCs
            out.write(str(C0)+'\t'+','.join(map(str,AUCs))+'\n')
        out.close()
        plotLR1('RidgeClf.allGenes.%s.%s.tsv' % (pheno,dtype),'Logistic Regression L2')

    elif classifier =='RF':
        ntrees = [10,20,50,80,100,250,500]
        leafSampleSizes = [2,5,8,10]
        out = open('RandomForest.allGenes.%s.%s.tsv' % (pheno,dtype),'w')
        for n in ntrees:
            for m in leafSampleSizes:
                AUCs = []
                for j in range(15):
                    X_train,X_test,y_train,y_test,description = splitTrainTest(df,pheno,129+j*3)
                    AUCs.append(RandomForest(X_train,X_test,y_train,y_test,n,m)[2])
                print AUCs
                out.write(str(n)+'\t'+str(m)+'\t'+','.join(map(str,AUCs))+'\n')
        out.close()
        plotRandomForest1('RandomForest.allGenes.%s.%s.tsv' % (pheno,dtype),'Random Forest')

    elif classifier =='GB':
        nestimator = [100,200,500,1000]
        maxDepth = [2,3]
        out = open('GradientBoost.allGenes.%s.%s.tsv' % (pheno,dtype),'w')
        for n in nestimator:
            for m in maxDepth:
                AUCs = []
                for j in range(15):
                    X_train,X_test,y_train,y_test,description = splitTrainTest(df,pheno,129+j*3)
                    AUCs.append(GBClassifier(X_train,X_test,y_train,y_test,n,m)[2])
                print AUCs
                out.write(str(n)+'\t'+str(m)+'\t'+','.join(map(str,AUCs))+'\n')
        out.close()
        plotRandomForest1('GradientBoost.allGenes.%s.%s.tsv' % (pheno,dtype),'Gradient Boosting')

    elif classifier == 'SVM':
        Cvalues = [0.001,0.01,0.05,0.1,0.5,1,5,10,50,100,500]
        out = open('SVM.allGenes.%s.%s.std.tsv' % (pheno,dtype),'w')
        for i in range(len(Cvalues)):
            AUCs = []
            C0 = Cvalues[i]
            for j in range(15):
                X_train,X_test,y_train,y_test,description = splitTrainTest(df,pheno,129+j*3)
                AUCs.append(SVMClassifier(X_train,X_test,y_train,y_test,C0)[2])
            print AUCs
            out.write(str(C0)+'\t'+','.join(map(str,AUCs))+'\n')
        out.close()
        plotLR1('SVM.allGenes.%s.%s.std.tsv' % (pheno,dtype),'SVM')

    elif classifier == 'KNN':
        Kvalues = [3,5,10,25]
        out = open('KNN.allGenes.%s.%s.tsv' % (pheno,dtype),'w')
        for i in range(len(Kvalues)):
            AUCs = []
            K = Kvalues[i]
            for j in range(15):
                X_train,X_test,y_train,y_test,description = splitTrainTest(df,pheno,129+j*3)
                AUCs.append(KNNClassifier(X_train,X_test,y_train,y_test,K)[2])
            print AUCs
            out.write(str(K)+'\t'+','.join(map(str,AUCs))+'\n')
        out.close()
        plotLR1('KNN.allGenes.%s.%s.tsv' % (pheno,dtype),'K Nearest Neighbors')



if __name__=='__main__':
    main()
