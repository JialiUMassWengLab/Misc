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
sys.path.append('/home/jzhuang@ms.local/NASH/2018FebValid/')
from process_sequencing import RandomForest
from process_sequencing import KNNClassifier
from process_sequencing import LRClassifier
from process_sequencing import SVMClassifier
from process_sequencing import RidgeClassifier

def splitTrainTest(table,target,seed):
    X_train,X_test,y_train,y_test = model_selection.train_test_split(table.iloc[:,:-3],table[target],\
                                                                     test_size=0.4,random_state=seed)
    names = pd.read_csv('1802NASHserum.congregated.tsv',sep='\t',index_col=0,usecols=['gene_id','Description'])
    table = pd.concat([table,names.transpose()],join='inner',axis=0)
    return X_train,X_test,y_train,y_test,table.iloc[-1,:]

def splitTrainTest2(df1,df2,target):
    X_train = df1.iloc[:,:-3]
    X_test = df2.iloc[:,:-3]
    y_train = df1[target]
    y_test = df2[target]
    return X_train,X_test,y_train,y_test


def getTestMat(fn2,pheno):
    '''
    array = []
    files = glob.glob('*_avg_TPM.csv')
    for fn in files:
        df0 = pd.read_csv(fn,index_col=0,low_memory=False)
        array.append(df0)
    df1 = pd.concat(array,join='inner',axis=1).transpose()
    '''
    df1 = pd.read_csv('/mnt/shares/R&D/projects/201803_NASH_Validation/2018AprRevalMerged/NASH_Re_Validation_All_RUVNormCounts.csv',index_col=0,low_memory=False).transpose()
    df1['aliquot'] = map(lambda x: x.split('_')[0],df1.index)
    df1 = df1.groupby('aliquot').agg(np.mean)

    df2 = pd.read_csv(fn2,index_col=0,low_memory=False).transpose()
    df = pd.concat([df1,df2],join='inner',axis=0)
    df = df.astype('float')

    info = pd.read_csv('../2018FebValid/sampleInfo.csv',index_col=0)
    info.index = info.index.astype('str')
    infoT = pd.read_csv('sampleInfoAll.csv',index_col=0)
    infoT.index = infoT.index.astype('str')
    info0 = pd.concat([info,infoT],join='inner',axis=0)
    df = pd.concat([df,info0],join='inner',axis=1)
    if pheno == 'PathFibrosis':
        df = df[(df['PathFibrosis']=='F0') | (df['PathFibrosis']=='F3') | (df['PathFibrosis']=='F4') | (df['PathFibrosis']=='None')]
        df['PathFibrosis'] = map(lambda x: 1 if x=='F3' or x=='F4' else 0, df['PathFibrosis'])
    elif pheno == 'NAFLDvsNASH':
        df = df[(df['Disease']=='NASH') | (df['Disease']=='NAFLD')]
        df['NAFLDvsNASH'] = map(lambda x: 0 if x=='NAFLD' else 1, df['Disease'])
        df.drop('PatientID',axis=1,inplace=True)
    elif pheno == 'F01vsF34':
        df = df[(df['PathFibrosis']!='F2') & (df['PathFibrosis']!='None')]
        df['F01vsF34'] = map(lambda x: 1 if x=='F3' or x=='F4' else 0, df['PathFibrosis'])
        df.drop('PatientID',axis=1,inplace=True)
    elif pheno == 'F012vsF34':
        #df = df[df['PathFibrosis']!='None']
        df['F012vsF34'] = map(lambda x: 1 if x=='F3' or x=='F4' else 0, df['PathFibrosis'])
        df.drop('PatientID',axis=1,inplace=True)
    elif pheno == 'F01vsF23':
        df = df[df['PathFibrosis']!='F4']
        df['F01vsF23'] = map(lambda x: 1 if x=='F2' or x=='F3' else 0, df['PathFibrosis'])
        df.drop('PatientID',axis=1,inplace=True)

    df1 = df[map(lambda x: not x in infoT.index,df.index)]
    df2 = df[map(lambda x: x in infoT.index,df.index)]
    return df1,df2



def plotROC(X_train,X_test,y_train,y_test,para,clf,prefix):
    if clf == 'LR':
        fpr,tpr,auc,score = LRClassifier(X_train,X_test,y_train,y_test,para[0])
    elif clf == 'RF':
        fpr,tpr,auc,score = RandomForest(X_train,X_test,y_train,y_test,para[0],para[1])
    elif clf == 'SVM':
        fpr,tpr,auc,score = SVMClassifier(X_train,X_test,y_train,y_test,para[0])
    elif clf == 'KNN':
        fpr,tpr,auc,score = KNNClassifier(X_train,X_test,y_train,y_test,para[0])
    elif clf == 'RC':
        fpr,tpr,auc,score = RidgeClassifier(X_train,X_test,y_train,y_test,para[0])
    else:
        fpr,tpr,auc,score = (0.,0.,0.,0.)

    out = open(prefix+'.label.prob.csv','w')
    for alq,true,pred in zip(X_test.index,y_test,score):
        out.write(','.join(map(str,[alq,true,pred]))+'\n')
    out.close()

    fig,ax = plt.subplots(figsize=(15,12))
    ax.plot(fpr,tpr,lw=2.5,label='ROC curve (area = %0.3f)' % auc,color='black')
    ax.plot([0,1],[0,1],color='red',lw=2,ls='--')
    plt.xlim([0,1])
    plt.ylim([0,1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.legend(loc='best')
    plt.savefig(prefix+'.test.ROC.png')
    plt.close()

def plotMeanROC(df,target,para,prefix):
    FPR = []
    TPR = []
    for i in range(15):
        X_train,X_test,y_train,y_test,description = splitTrainTest(df,target,129+i*3)
        fpr,tpr = LRClassifier(X_train,X_test,y_train,y_test,para[0])[:2]
        #fpr,tpr = RandomForest(X_train,X_test,y_train,y_test,para[0],para[1])[:2]
        FPR.append(fpr)
        TPR.append(tpr)
    all_fpr = np.unique(np.concatenate([fpr for fpr in FPR]))
    mean_tpr = np.zeros_like(all_fpr)
    for i in range(15):
        mean_tpr += scipy.interp(all_fpr,FPR[i],TPR[i])
        plt.plot(FPR[i],TPR[i],color='grey',alpha=0.4)
        #plt.plot(all_fpr,scipy.interp(all_fpr,FPR[i],TPR[i]),color='grey',alpha=0.4)
    mean_tpr /= 15
    AUC = metrics.auc(all_fpr,mean_tpr)
    plt.plot(all_fpr,mean_tpr,lw=2,label='Mean ROC curve (area = %0.3f)' % AUC)
    plt.plot([0,1],[0,1],color='red',ls='--')
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.legend(loc='best')
    plt.savefig(prefix+'.mean_ROC.png')
    plt.close()

def plotMeanROC2(df,target,para,clf,prefix):
    FPR = []
    TPR = []
    out = open(prefix+'.samp.pred.csv','w')
    for i in range(15):
        X_train,X_test,y_train,y_test,description = splitTrainTest(df,target,129+i*3)
        if clf == 'LR':
            fpr,tpr,auc,score = LRClassifier(X_train,X_test,y_train,y_test,para[0])
        elif clf == 'RF':
            fpr,tpr,auc,score = RandomForest(X_train,X_test,y_train,y_test,para[0],para[1])
        elif clf == 'SVM':
            fpr,tpr,auc,score = SVMClassifier(X_train,X_test,y_train,y_test,para[0])
        elif clf == 'KNN':
            fpr,tpr,auc,score = KNNClassifier(X_train,X_test,y_train,y_test,para[0])
        elif clf == 'RC':
            fpr,tpr,auc,score = RidgeClassifier(X_train,X_test,y_train,y_test,para[0])
        else:
            fpr,tpr,auc,score = (0.,0.,0.,0.)
        FPR.append(fpr)
        TPR.append(tpr)
        out.write('iter%d\n' % i)
        for true,pred in zip(y_test,score):
            out.write(','.join(map(str,[true,pred]))+'\n')
    out.close()

    all_fpr = np.unique(np.concatenate([fpr for fpr in FPR]))
    mean_tpr = np.zeros_like(all_fpr)
    segs = np.zeros((15,len(all_fpr)))
    for i in range(15):
        #for j in range(len(FPR[i])):
        #    segs[int(FPR[i][j]/0.1)].append(TPR[i][j])
        segs[i,] = scipy.interp(all_fpr,FPR[i],TPR[i])
        mean_tpr += scipy.interp(all_fpr,FPR[i],TPR[i])
        
    mean_tpr /= 15
    AUC = metrics.auc(all_fpr,mean_tpr)
    means = np.apply_along_axis(np.mean,0,segs)
    stds = np.apply_along_axis(np.std,0,segs)
    ses = np.apply_along_axis(ss.sem,0,segs)
    df = pd.DataFrame({'tpr': means, 'std': stds, 'se': ses, 'xpos': all_fpr})
    df.plot('xpos','tpr',kind='line',lw=2,label='Mean ROC curve (area = %0.3f)' % AUC)
    plt.fill_between(all_fpr,means-ses,means+ses,color='b',alpha=0.2)
    plt.plot([0,1],[0,1],color='red',ls='--')
    #plt.xlim(-0.02,1)
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.legend(loc='lower right')
    plt.savefig(prefix+'.ROC_with_err.png')
    plt.close()

def writeScores(df,target,para,clf,prefix):
    repeats = []
    for i in range(15):
        X_train,X_test,y_train,y_test,description = splitTrainTest(df,target,129+i*3)
        if clf == 'LR':
            scores = LRClassifier(X_train,X_test,y_train,y_test,para[0])[3]
        elif clf == 'RF':
            scores = RandomForest(X_train,X_test,y_train,y_test,para[0],para[1])[3]
        elif clf == 'SVM':
            scores = SVMClassifier(X_train,X_test,y_train,y_test,para[0])[3]
        elif clf == 'KNN':
            scores = KNNClassifier(X_train,X_test,y_train,y_test,para[0])[3]
        elif clf == 'RC':
            scores = RidgeClassifier(X_train,X_test,y_train,y_test,para[0])[3]            
        else:
            scores = [0.] * len(y_test)

        repeat = {}
        for j in range(len(scores)):
            repeat.update({'sample'+str(j)+'label': y_test[j]})
            repeat.update({'sample'+str(j)+'prob': scores[j]})
        repeats.append(repeat)

    mat = pd.DataFrame(repeats)
    mat.to_csv(prefix+'.scores.csv')


def main():
    pheno = sys.argv[1]
    dtype = sys.argv[2]
    #classifier = sys.argv[3]
    df1,df2 = getTestMat('/mnt/shares/R&D/projects/201803_NASH_Validation/Data/All_data/20180305_NASH_Validation_204_All_%s_cutoff.csv' % dtype,pheno)
    print df1.shape,df2.shape

    sns.set(font_scale=1.8)
    paras = {}
    files = glob.glob(os.path.join('/home/jzhuang@ms.local/NASH/2018FebValid','para_select_files_Apr_revalid','*.%s.%s.*tsv' % (pheno,dtype)))
    for fn in files:
        clf = os.path.basename(fn).split('.')[0]
        print fn
        if clf == 'RandomForest':
            clf = 'RF'
        elif clf == 'LogisticRegre':
            clf = 'LR'
        elif clf == 'RidgeClf':
            clf = 'RC'
        tbl = pd.read_csv(fn,sep='\t',header=None)
        tbl['median'] = tbl.apply(lambda x: np.median(map(float,x[tbl.shape[1]-1].split(','))),axis=1)
        paras[clf] = tuple(tbl.sort_values('median',ascending=False).iloc[0,:-2])
        #writeScores(df,pheno,paras[clf],clf,'.'.join([clf,'allGenes',pheno,dtype]))
        #plotMeanROC2(df,pheno,paras[clf],clf,'.'.join([clf,'allGenes',pheno,dtype]))
        X_train,X_test,y_train,y_test = splitTrainTest2(df1,df2,pheno)
        plotROC(X_train,X_test,y_train,y_test,paras[clf],clf,'.'.join([clf,'allGenes',pheno,dtype]))
        
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
