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
from process_sequencing import readMatrix2

def splitTrainTest(table,target,seed):
    X_train,X_test,y_train,y_test = model_selection.train_test_split(table.iloc[:,:-3],table[target],\
                                                                     test_size=0.4,random_state=seed)
    names = pd.read_csv('1802NASHserum.congregated.tsv',sep='\t',index_col=0,usecols=['gene_id','Description'])
    table = pd.concat([table,names.transpose()],join='inner',axis=0)
    return X_train,X_test,y_train,y_test,table.iloc[-1,:]


def WriteCoef(df,para,clf,pheno,dtype):
    X_train,X_test,y_train,y_test,names = splitTrainTest(df,pheno,249)
    X = pd.concat([X_train,X_test],axis=0)
    y = pd.concat([y_train,y_test],axis=0)
    if clf == 'LR':
        X = preprocessing.StandardScaler().fit_transform(X)
        lr = linear_model.LogisticRegression(penalty='l1',C=para[0],random_state=254,fit_intercept=False)
        lr.fit(X,y)
        coef = pd.Series(lr.coef_[0],index=names,name='Coef')
    elif clf == 'RF':
        rf = ensemble.RandomForestClassifier(min_samples_split=para[1],n_estimators=para[0],random_state=119)
        rf.fit(X_train,y_train)
        coef = pd.Series(rf.feature_importances_,index=names,name='Coef')
    elif clf == 'RC':
        X = preprocessing.StandardScaler().fit_transform(X)
        rc = linear_model.LogisticRegression(penalty='l2',C=para[0],random_state=254,fit_intercept=False)
        rc.fit(X,y)
        coef = pd.Series(rc.coef_[0],index=names,name='Coef')

    coef = coef[~(coef==0)]
    coef.sort_values(ascending=False,inplace=True)
    coef.to_csv('.'.join([clf,'allGenes',pheno,dtype])+'.Coefficient.csv')
    print clf,pheno,dtype


def main():
    pheno = sys.argv[1]
    files = glob.glob(os.path.join('para_select_files','*.%s.*.*tsv' % pheno))
    #files = glob.glob(os.path.join('para_select_files_Mar8th_more_phenos','*.%s.*.*tsv' % pheno))
    for fn in files:
        clf = os.path.basename(fn).split('.')[0]
        dtype = os.path.basename(fn).split('.')[3]
        if clf == 'RandomForest':
            clf = 'RF'
        elif clf == 'LogisticRegre':
            clf = 'LR'
        elif clf == 'RidgeClf':
            clf = 'RC'
        else:
            clf = ''
        if not clf == '':
            df = readMatrix2('/mnt/shares/R&D/projects/201803_NASH_Validation/Data/All_data/20180305_NASH_Validation_204_All_%s_cutoff.csv' % dtype,pheno)
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
