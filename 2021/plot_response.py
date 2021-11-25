#! /usr/bin/python

import re
import os
import sys
import glob
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

def getData(metric):
    files = glob.glob('/home/jzhuang/DB109/Source/*/eff_%s_nom.csv' % metric)
    array = []
    for fn in files:
        d = pd.read_csv(fn,index_col=0)
        d = d[d['POP_FASFN']==1]
        d = d.dropna(subset=['VISIT'])
        if 'PARAM' in d.columns:
            d = d[d['PARAM']=='ADAS-cog Total Score']
        else:            
            d = d.rename(columns={'ADAS_TOT_DL':'AVAL_DL'})

        d = d.loc[:,list(map(lambda x: x=='VISIT' or x=='STUDYID' or x=='ANALYSIS_WEEK' or \
                             not re.search(r'_DL$',x)==None,d.columns))]
        array.append(d)

    df = pd.concat(array,join='inner',axis=0)
    df = df[df['ANALYSIS_WEEK']>=0]
    #print(df)
    return df


def getResponseTbl(metric):
    demo = pd.read_csv(os.path.join('/home/jzhuang/DB109/Output','patient_info.csv'),index_col=0)
    #demo.loc[:,'Treatment2'] = demo.apply(lambda x: x['Study']+' '+x['Treatment'],axis=1)
    df = getData(metric)
    #only keep subjects who finished the study
    subject_finished = list(df[(df['ANALYSIS_WEEK']==12) | (df['ANALYSIS_WEEK']==24)].index)
    df = df.loc[list(map(lambda x: x in subject_finished,df.index)),:]

    data = {}
    for subject,group in df.groupby(level=0):
        data.update({subject: group.sort_values('ANALYSIS_WEEK')['AVAL_DL'][-1]})

    S = pd.Series(data,name='ADAS_TOT_DL')
    return S


def plotTimeResponse(metric):
    demo = pd.read_csv(os.path.join('/home/jzhuang/DB109/Output','patient_info.csv'),index_col=0)
    demo = demo[demo['TRBTYP']=='DONEPEZIL']
    #demo.loc[:,'Treatment2'] = demo.apply(lambda x: x['Study']+' '+x['Treatment'],axis=1)
    df = getData(metric)

    sns.set(font_scale=1.5,context='talk')
    from matplotlib.backends.backend_pdf import PdfPages
    pdf = PdfPages('%s_response.pdf' % metric)
    for study in df['STUDYID'].unique():
        tmp = df.dropna(axis=0,subset=['AVAL_DL'])
        tmp = tmp[tmp['STUDYID']==study]
        tmp['Subject'] = tmp.index
        #S = tmp['Subject'].value_counts()
        #print(tmp[list(map(lambda x: x in list(S[S==2].index),tmp['Subject']))])

        tmp.index = tmp.apply(lambda x: str(x['Subject'])+' - '+str(x['ANALYSIS_WEEK']),axis=1)
        #remove duplicated measurements
        tmp = tmp[~tmp.index.duplicated(keep='first')]
        tmp = tmp[['Subject','ANALYSIS_WEEK','AVAL_DL']].merge(demo,left_on='Subject',right_index=True,how='left')
        print(pd.crosstab(tmp['ANALYSIS_WEEK'],tmp['ARM']))

        fig,ax = plt.subplots(figsize=(20,12))
        #sns.boxplot(x='Time',y=m2,hue='Treatment2',data=tmp,width=0.6,ax=ax)
        sns.pointplot(x='ANALYSIS_WEEK',y='AVAL_DL',hue='ARM',data=tmp,dodge=True,ci=68.47,join=True,ax=ax)
        plt.title(study)
        plt.savefig(pdf,format='pdf')
        plt.close()
    pdf.close()


def getAUC(metric):
    df = getData(metric)
    #only keep subjects who finished the study
    subject_finished = list(df[(df['ANALYSIS_WEEK']==12) | (df['ANALYSIS_WEEK']==24)].index)
    df = df.loc[list(map(lambda x: x in subject_finished,df.index)),:]

    tmp = df.dropna(axis=0,subset=['AVAL_DL']).copy()
    tmp['Subject'] = tmp.index
    tmp.index = tmp.apply(lambda x: str(x['Subject'])+' - '+str(x['ANALYSIS_WEEK']),axis=1)
    #remove duplicated measurements
    tmp = tmp[~tmp.index.duplicated(keep='first')]

    AUC = {}
    for pid,group in tmp[['Subject','ANALYSIS_WEEK','AVAL_DL']].groupby('Subject'):
        group = group.sort_values('ANALYSIS_WEEK',ascending=True)
        changes = list(group['AVAL_DL'])
        if not 24 in list(group['ANALYSIS_WEEK']):
            changes.append(group.loc[str(pid)+' - 12.0','AVAL_DL'])
        auc = 0
        for i in range(group.shape[0]-1):
            auc += (changes[i]+changes[i+1])/2.0*(group.iloc[i+1,-2]-group.iloc[i,-2])
            AUC.update({str(pid):auc})

    S = pd.Series(AUC,name=metric+'_AUC')
    demo = pd.read_csv(os.path.join('/home/jzhuang/DB109/Output','patient_info.csv'),index_col=0)
    demo = pd.concat([demo,S],join='inner',axis=1)
    array = []
    for study,group in demo.groupby('STUDYID'):
        pbo_mean = np.mean(group[(group['ARM']=='PBO') | (group['ARM']=='Placebo')][metric+'_AUC'])
        pbo_std = np.std(group[(group['ARM']=='PBO') | (group['ARM']=='Placebo')][metric+'_AUC'])
        group2 = group[(group['ARM']!='PBO') & (group['ARM']!='Placebo')].copy()
        group2[metric+'_AUC_z'] = group2.apply(lambda x: (x[metric+'_AUC']-pbo_mean)/pbo_std,axis=1)
        group2['STUDYID'] = [study] * group2.shape[0]
        array.append(group2[['ARM','STUDYID',metric+'_AUC_z']])

    zDf = pd.concat(array,join='inner',axis=0)
    zDf.to_csv('/home/jzhuang/DB109/Output/%s_AUC_zscore.csv' % metric)
    print(zDf)

    demo[['ARM','STUDYID',metric+'_AUC']].to_csv('/home/jzhuang/DB109/Output/%s_AUC.csv' % metric)


def main():
    plotTimeResponse('adas')
    #S = getResponseTbl('adas')
    #S.to_csv('/home/jzhuang/DB109/Output/response.csv')
    #getAUC('adas')


if __name__=='__main__':
    main()

