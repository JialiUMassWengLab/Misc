#! /usr/bin/python

import re
import os
import sys
import numpy as np
import pandas as pd
import scipy.stats as st
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import NMF
from sklearn.decomposition import non_negative_factorization

def getMatrix():
    df = pd.read_csv('/mnt/nfs/analysis/180408_APK88PSC/180408_APK88PSC.congregated.tsv',sep='\t',index_col=0)
    df = df.loc[:,map(lambda x: not re.search(r'-PSC-',x)==None,df.columns)]
    df.columns = map(lambda x: x.split('-')[0],df.columns)
    df2 = pd.read_csv('/home/jzhuang@ms.local/NASH/2018FebValid/180305_NASHvalidation_203_avg_TPM.csv',index_col=0,low_memory=False).iloc[:-3,:]
    df3 = pd.read_csv('/home/jzhuang@ms.local/NASH/2018AprRevalMerged/1804Revalidation_avg_TPM.csv',index_col=0,low_memory=False).iloc[:-3,:]
    df4 = pd.read_csv('/home/jzhuang@ms.local/NASH/2018AprRevalMerged/1804Revalidation2_avg_TPM.csv',index_col=0,low_memory=False).iloc[:-3,:]
    #df = pd.concat([df,df2,df3,df4],join='inner',axis=1)
    df = pd.concat([df2,df3,df4],join='inner',axis=1)
    df.index = map(lambda x: re.sub(r'\.\d+$','',x),df.index)
    df = df.astype('float')
    #df = df.apply(lambda x: 2**x, axis=0)
    df = df.loc[map(lambda x: not re.search(r'^ERCC-',x),df.index),:]
    df = df.apply(lambda x: x/sum(x)*1000000,axis=0)
    df = df.loc[df.apply(lambda x: any(x>50),axis=1),:]
    maxval = df.apply(max,axis=1)
    df = df.apply(lambda x: x/max(x)*100,axis=1)
    return df,maxval


def getMatrix2():
    df = pd.read_csv('/mnt/shares/R&D/projects/201803_NASH_Validation/Data/All_data/20180305_NASH_Validation_204_All_RUVNormCounts_2frag_cutoff.csv',index_col=0,low_memory=False)
    df2 = pd.read_csv('/mnt/shares/R&D/projects/201803_NASH_Validation/2018AprRevalMerged/NASH_Re_Validation_All_RUVNormCounts.csv',index_col=0,low_memory=False).transpose()
    df2['aliquot'] = map(lambda x: x.split('_')[0],df2.index)
    df2 = df2.groupby('aliquot').agg(np.mean).transpose()
    df = pd.concat([df,df2],join='inner',axis=1)
    df.index = map(lambda x: re.sub(r'\.\d+$','',x),df.index)
    df = df.astype('float')
    df = df.loc[map(lambda x: not re.search(r'^ERCC-',x),df.index),:]
    df = df.loc[df.apply(lambda x: any(x>20),axis=1),:]
    maxval = df.apply(max,axis=1)
    df = df.apply(lambda x: x/max(x)*100,axis=1)
    return df,maxval


def applyNMF(df,alpha,p):
    mat = df.transpose()
    #model = NMF(n_components=15,random_state=142,solver='cd',max_iter=1000000,alpha=alpha,l1_ratio=1,init='nndsvd')
    #W = pd.DataFrame(data=model.fit_transform(mat),index=mat.index)
    #H = pd.DataFrame(data=model.components_,columns=mat.columns)
    w, h, n_iter = non_negative_factorization(mat,n_components=p,init='nndsvd',random_state=142,alpha=alpha,l1_ratio=0,solver='cd',max_iter=1000000,regularization='components')
    W = pd.DataFrame(data=w,index=mat.index)
    H = pd.DataFrame(data=h,columns=mat.columns)
    #W.to_csv('W_TPM_p%d_a%d.mat.tsv' % (p,alpha),sep='\t')
    #H.to_csv('H_TPM_p%d_a%d.mat.tsv' % (p,alpha),sep='\t')
    W.to_csv('W_RUV_p%d_a%d.mat.tsv' % (p,alpha),sep='\t')
    H.to_csv('H_RUV_p%d_a%d.mat.tsv' % (p,alpha),sep='\t')
    print W
    print H
    print n_iter


def plotComp(Wmat,cat):
    W = pd.read_csv(Wmat,sep='\t',index_col=0).astype('float')
    info1 = pd.read_csv('/home/jzhuang@ms.local/NASH/2018FebValid/sampleInfo.csv',index_col=0)
    info2 = pd.read_csv('/home/jzhuang@ms.local/NASH/2018AprRevalMerged/sampleInfoAll.csv',index_col=0)
    info3 = pd.read_csv('/home/jzhuang@ms.local/NASH/PSC_1804/sampleInfo.csv',index_col=0)
    info = pd.concat([info1,info2,info3],join='inner',axis=0)
    W.index = W.index.astype(int)
    W = pd.concat([W,info],join='inner',axis=1)
    W = W[W[cat]!='None']
    print W

    sns.set(font_scale=1.8)
    from matplotlib.backends.backend_pdf import PdfPages
    pdf = PdfPages(Wmat.replace('.mat.tsv','.%s.box.pdf' % cat))
    for col in W.columns[:-4]:
        sArray = []
        i = 1
        for key,group in W.groupby(cat):
            sArray.append(pd.DataFrame({'x': np.random.normal(i,0.05,group.shape[0]), 'y': group[col]}))
            i += 1

        dfScat = pd.concat(sArray,join='inner',axis=0)
        fig,ax = plt.subplots(figsize=(12,10))
        bp = W.boxplot(column=col,by=cat,boxprops={'linewidth':1.8},medianprops={'linewidth':1.8},whiskerprops={'linewidth':1.8},capprops={'linewidth':1.8},sym='',ax=ax,return_type='dict')
        [[item.set_color('blue') for item in bp[key]['boxes']] for key in bp.keys()]
        [[item.set_color('blue') for item in bp[key]['whiskers']] for key in bp.keys()]
        [[item.set_color('red') for item in bp[key]['medians']] for key in bp.keys()]
        dfScat.plot('x','y',kind='scatter',s=80,ax=ax,color='blue',fontsize=20,title='comp'+str(col))
        plt.xlabel(cat)
        plt.ylabel('Coefficient')
        plt.savefig(pdf,format='pdf')
        plt.close()
        
    pdf.close()


def plotMapping(Hmat):
    H = pd.read_csv(Hmat,sep='\t',index_col=0).astype('float').transpose()
    H = H.apply(lambda x: x/sum(x),axis=1)
    df0 = pd.read_csv('/mnt/nfs/analysis/180408_APK88PSC/180408_APK88PSC.congregated.tsv',index_col=0,sep='\t').iloc[:,-2:]
    df0.index = map(lambda x: re.sub(r'\.\d+$','',x),df0.index)
    H = pd.concat([H,df0['Description']],join='inner',axis=1)
    H = H.sort_values(7)

    coef = pd.read_csv('/mnt/shares/R&D/projects/201803_NASH_Validation/Coefficients/RC.allGenes.F01vsF34.RUVNormCounts_10frag.Coefficient.csv',header=None).fillna('NA')
    coef = coef[coef[0]!='NA']
    coef.columns = ['Description','coef']
    #geneList = list(coef.iloc[:25,0]) + list(coef.iloc[-25:,0])
    #H = H.loc[map(lambda x: x in geneList,H.index)]
    H = H.merge(coef.iloc[:50,],on='Description').sort_values('coef',ascending=False)
    H.index = H['Description']

    sns.set(font_scale=2.0)
    fig,ax = plt.subplots(figsize=(26,8))
    sns.heatmap(H.iloc[:,:-2].transpose(),cmap="RdBu_r",vmax=0.7,vmin=-0.7,square=True)
    plt.yticks(rotation=0)
    plt.xticks(rotation=45,ha='right')
    plt.savefig(Hmat.replace('.mat.tsv','.RC_genes_loading.heatmap.png'))
    plt.close()


def corWithGTEx(Hmat,maxval):
    H = pd.read_csv(Hmat,sep='\t',index_col=0).astype('float').transpose()
    H.columns = map(lambda x: 'comp'+str(x),H.columns)
    H = H.apply(lambda x: x*maxval,axis=0)
    df = pd.read_csv('/home/jzhuang@ms.local/CoExpression/GTEx.tissues.median.tsv',sep='\t',index_col=0)
    blood = pd.read_csv('/mnt/shares/Users/jzhuang/blood_cell_types/BloodCellTypes.expr4.tsv',sep='\t',index_col=0).iloc[:,:-2]
    blood.index = map(lambda x: re.sub(r'\.\d+$','',x),blood.index)
    df = pd.concat([blood,df],join='inner',axis=1)
    #df = pd.read_csv('/mnt/shares2/annotations/hg38/Blueprint_cellType_medians.tsv',sep='\t',index_col=0)
    series = []
    for col in H.columns:
        df1 = H
        df1['perc'] = df1.apply(lambda x: x[col]/sum(x),axis=1)
        df1 = df1[df1['perc']>0.45]
        df1 = pd.concat([df1[col],df.iloc[:,:-1]],join='inner',axis=1)
        print col,df1.shape
        pccs = {}
        for tissue in df.columns[:-1]:
            comb = df1[[col,tissue]]
            comb = comb[comb[tissue]>50]
            comb = comb.apply(lambda x: np.log10(x+0.01),axis=0)
            if comb.shape[0] > 15:
                pcc = st.pearsonr(comb[tissue],comb[col])[0]
            else:
                pcc = 0.0
            pccs.update({tissue: pcc})

        s = pd.Series(pccs,name='comp'+str(col))
        series.append(s)

    df = pd.concat(series,axis=1)
    sns.set(font_scale=2.2)
    fig,ax = plt.subplots(figsize=(16,23))
    sns.heatmap(df,cmap="RdBu_r",vmax=1.0,vmin=-1.0,square=True)
    #plt.savefig(Hmat.replace('.mat.tsv','.bluprint.heatmap.png'))
    plt.savefig(Hmat.replace('.mat.tsv','.GTEx.heatmap.png'))
    plt.close()

    return df


def exprDist(Hmat):
    H = pd.read_csv(Hmat,sep='\t',index_col=0).astype('float').transpose()
    H.columns = map(lambda x: 'comp'+str(x),H.columns)
    #df = pd.read_csv('/home/jzhuang@ms.local/CoExpression/GTEx.tissues.median.tsv',sep='\t',index_col=0)
    #blood = pd.read_csv('/mnt/shares/Users/jzhuang/blood_cell_types/BloodCellTypes.expr4.tsv',sep='\t',index_col=0).iloc[:,:-2]
    #blood.index = map(lambda x: re.sub(r'\.\d+$','',x),blood.index)
    #df = pd.concat([blood,df],join='inner',axis=1)
    df = pd.read_csv('/mnt/shares2/annotations/hg38/Blueprint_cellType_medians.tsv',sep='\t',index_col=0).iloc[:-1,:]
    df = df.loc[:,map(lambda x: not re.search(r'endothelial',x) or x=='Description',df.columns)]
    names = df['Description']
    df = df[df.iloc[:,:-1].apply(lambda x: any(x>50),axis=1)]
    df = df.iloc[:,:-1].apply(lambda x: x/max(x),axis=1)
    df = df[df.apply(lambda x: len([i for i in x if i > 0.2]) <= 8,axis=1)]

    sns.set(font_scale=2)
    from matplotlib.backends.backend_pdf import PdfPages
    pdf = PdfPages(Hmat.replace('.mat.tsv','.bluprint.exprMap.pdf'))
    #pdf = PdfPages(Hmat.replace('.mat.tsv','.GTEx.exprMap.pdf'))
    for comp in H.columns:
        df1 = H
        df1['perc'] = df1.apply(lambda x: x[comp]/sum(x),axis=1)
        geneList = df1[df1['perc']>0.35].index
        df2 = df[map(lambda x: x in geneList,df.index)]
        df2 = pd.concat([df2,names],join='inner',axis=1)
        df2.index = df2['Description']
        df2 = df2.drop('Description',axis=1).transpose()
        print df2.shape
        if df2.shape[1] < 10:
            continue
        #fig,ax = plt.subplots(figsize=(24,18))
        #sns.heatmap(df2,cmap="RdBu_r",vmax=1.0,vmin=-1.0,linecolor='black',linewidths=1)
        #plt.xticks(rotation=45,ha='right')
        #plt.title(comp)
        #plt.tight_layout()
        g = sns.clustermap(df2,cmap="RdBu_r",vmax=1.0,vmin=-1.0,metric='correlation',z_score=None,standard_scale=None,row_cluster=False,col_cluster=True,linecolor='black',linewidths=1,figsize=(24,18))
        plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(),rotation='horizontal')
        plt.setp(g.ax_heatmap.xaxis.get_majorticklabels(),rotation=45)
        plt.setp(g.ax_heatmap.xaxis.get_majorticklabels(),ha='right')
        g.fig.suptitle(comp)
        g.savefig(pdf,format='pdf')
        plt.close()
    pdf.close()


def main():
    alpha = 0
    p = 12
    #df,maxval = getMatrix2()
    #print df
    #applyNMF(df,alpha,p)
    plotComp('W_RUV_p%d_a%d.mat.tsv' % (p,alpha),sys.argv[1])
    #corWithGTEx('H_TPM_p%d_a%d.mat.tsv' % (p,alpha),maxval)
    #exprDist('H_RUV_p%d_a%d.mat.tsv' % (p,alpha))
    #plotMapping('H_RUV_p%d_a%d.mat.tsv' % (p,alpha))




if __name__=='__main__':
    main()

