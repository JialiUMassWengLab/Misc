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
#from process_sequencing import readMatrix
#from process_sequencing import readMatrix2
#from process_sequencing2 import readMatrix2
#from process_sequencing2 import readMatrix3
from process_sequencing3 import readMatrix2
from process_sequencing import splitTrainTest

def WriteCoef(df,para,clf,pheno,dtype):
    X_train,X_test,y_train,y_test,names = splitTrainTest(df,pheno,249)
    X = pd.concat([X_train,X_test],axis=0)
    y = pd.concat([y_train,y_test],axis=0)
    names = pd.read_csv('../2018DecAD/1812ADcohort.congregated.tsv',sep='\t',index_col=0)['Description']
    if clf == 'LR':
        X = preprocessing.StandardScaler().fit_transform(X)
        lr = linear_model.LogisticRegression(penalty='l1',C=para[0],random_state=254,fit_intercept=False)
        lr.fit(X,y)
        coef = pd.Series(lr.coef_[0],index=df.columns[:-3],name='Coef')
    elif clf == 'RF':
        rf = ensemble.RandomForestClassifier(min_samples_split=para[1],n_estimators=para[0],random_state=119)
        rf.fit(X,y)
        coef = pd.Series(rf.feature_importances_,index=df.columns[:-3],name='Coef')
    elif clf == 'RC':
        X = preprocessing.StandardScaler().fit_transform(X)
        rc = linear_model.LogisticRegression(penalty='l2',C=para[0],random_state=254,fit_intercept=False)
        rc.fit(X,y)
        coef = pd.Series(rc.coef_[0],index=df.columns[:-3],name='Coef')

    coef = pd.concat([coef,names],join='inner',axis=1)
    coef = coef[coef['Coef']!=0]
    coef.index = map(lambda x: re.sub(r'\.\d+$','',x),coef.index)
    #coef.sort_values('Coef',ascending=False).to_csv('.'.join([clf,'allGenes',pheno,dtype])+'.Coefficient.csv')
    coef.sort_values('Coef',ascending=False).to_csv('.'.join([clf,'topGenes',pheno,dtype])+'.Coefficient.csv')
    print clf,pheno,dtype


def main():
    pheno = sys.argv[1]
    dtype = sys.argv[2]
    #files = glob.glob(os.path.join('para_select_files_KentuckyTopGenes_r1','*.%s.%s.*tsv' % (pheno,dtype)))
    #files = glob.glob(os.path.join('para_select_files_testKentucky_r1','*.%s.%s.*tsv' % (pheno,dtype)))
    #files = glob.glob(os.path.join('para_select_files_split171xUCSD20_r1','*.%s.%s.*tsv' % (pheno,dtype)))
    files = glob.glob(os.path.join('para_select_files_DESeq2topGenesxUCSD20_r1','*.%s.%s.*tsv' % (pheno,dtype)))
    #files = glob.glob(os.path.join('para_select_files_DESeq2topGenes_r1','*.%s.%s.*tsv' % (pheno,dtype)))
    for fn in files:
        clf = os.path.basename(fn).split('.')[0]
        #dtype = os.path.basename(fn).split('.')[3]
        if clf == 'RandomForest':
            clf = 'RF'
        elif clf == 'LogisticRegre':
            clf = 'LR'
        elif clf == 'RidgeClf':
            clf = 'RC'
        else:
            clf = ''
        if not clf == '':
            df = readMatrix2('GatesADboth_r1_%s.csv' % dtype,pheno)
            print df.shape
            tbl = pd.read_csv(fn,sep='\t',header=None)
            tbl['median'] = tbl.apply(lambda x: np.median(map(float,x[tbl.shape[1]-1].split(','))),axis=1)
            para = tuple(tbl.sort_values('median',ascending=False).iloc[0,:-2])
            WriteCoef(df,para,clf,pheno,dtype)
        
    '''
    X_train,X_test,y_train,y_test,description = splitTrainTest(df,pheno,123)
    if classifier == 'LR':
        coef = LRCoef(X_train,X_test,y_train,y_test,paras[comb],description)
    elif classifier == 'RF':
        coef = RFCoef(X_train,X_test,y_train,y_test,paras[comb],description)
    coef.to_csv('.'.join([classifier,'allGenes',pheno,dtype,'coef']),sep='\t')
    '''


if __name__=='__main__':
    main()
