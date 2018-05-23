#! /usr/bin/python

import re
import os
import sys
import numpy as np
import pandas as pd
import scipy.stats as stat
from math import pi
from sklearn.svm import SVR
from sklearn.svm import LinearSVR
from scipy.optimize import minimize
from cvxopt import solvers,matrix,log,spdiag
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

def getMatrix():
    df = pd.read_csv('180417_AI066Human.congregated.tsv',sep='\t',index_col=0).iloc[:,:-1]
    df.index = map(lambda x: re.sub(r'\.\d+$','',x),df.index)
    name = df['Description']
    df = df.loc[:,map(lambda x: not re.search(r'-BMSF-',x)==None,df.columns)]
    #df = df.transpose()
    #df['disease'] = map(lambda x: 'BM' if not re.search(r'^9',x) else 'serum',df.index)
    '''
    df2 = pd.read_csv('/home/jzhuang@ms.local/NASH/2018FebValid/180305_NASHvalidation_203_avg_TPM.csv',index_col=0,low_memory=False)
    df2 = df2.loc[:,map(lambda x: x=='Normal Control',df2.loc['Disease'])].iloc[:-3,:].transpose()
    df2['disease'] = ['SDBB2']*df2.shape[0]
    df = pd.concat([df,df2],join='inner',axis=0)
    '''
    df = df.loc[map(lambda x: not re.search(r'^ERCC-',x),df.index)]
    df = df.apply(lambda x: x/sum(x)*1000000,axis=0).fillna(0)
    #df = df.loc[:,df.iloc[:-1,:].apply(lambda x: len([i for i in x if i > 1]) > (df.shape[0]-1) * 0.8, axis=0)]
    return df,name

def getBase():
    bp = pd.read_csv('/mnt/shares/Users/jzhuang/Blueprint/Blueprint_cellType_medians.tsv',sep='\t',index_col=0).iloc[:-1,:]
    bp = bp.drop(['myeloid_cell','peripheral_blood_mononuclear_cell','mononuclear_cell_of_bone_marrow'],axis=1)
    rbc = pd.read_csv('/mnt/shares/Users/jzhuang/blood_cell_types/GSE53983.erythroid.median.tsv',sep='\t',index_col=0)
    bp = bp.merge(rbc,left_on='Description',right_index=True,how='left').fillna(0).drop('Description',axis=1).astype('float')

    celltypes = {}
    celltypes.update({'RBC': list(rbc.columns)})
    celltypes.update({'megakaryocyte' : ['CD34-negative_CD41-positive_CD42-positive_megakaryocyte_cell','common_myeloid_progenitor','megakaryocyte-erythroid_progenitor_cell']})
    celltypes.update({'neutrophil' : list(bp.columns[map(lambda x: not re.search(r'neutrophi',x)==None,bp.columns)]) + ['common_myeloid_progenitor','granulocyte_monocyte_progenitor_cell']})
    #celltypes.update({'erythroid' : list(bp.columns[map(lambda x: not re.search(r'erythro',x)==None,bp.columns)]) + ['common_myeloid_progenitor']})
    #celltypes.update({'NK' : ['cytotoxic_CD56-dim_natural_killer_cell']})
    #celltypes.update({'monocyte' : list(bp.columns[map(lambda x: not re.search(r'monoc',x)==None or not re.search(r'macrophage',x)==None,bp.columns)]) + ['common_myeloid_progenitor']})
    #celltypes.update({'Tcell' : list(bp.columns[map(lambda x: not re.search(r'T_cell',x)==None or not re.search(r'thymocyte',x)==None,bp.columns)])})
    #celltypes.update({'Bcell' : list(bp.columns[map(lambda x: not re.search(r'B_cell',x)==None,bp.columns)]) + ['plasma_cell']})
    #celltypes.update({'MSC' : list(bp.columns[map(lambda x: not re.search(r'endothelial',x)==None,bp.columns)]) + ['mesenchymal_stem_cell_of_the_bone_marrow']})

    geneList = []
    cellTypes = []
    for lineage in celltypes.keys():
        bp1 = bp[celltypes[lineage]]
        bp2 = bp.drop(celltypes[lineage],axis=1)
        bp0 = pd.concat([bp2,bp1.apply(max,axis=1)],join='inner',axis=1)
        geneList += list(bp0[bp0.apply(lambda x: len([i for i in x if x[-1] < i * 20])==1 and x[-1]>50,axis=1)].index)
        cellTypes += celltypes[lineage]

    bp = bp[list(set(cellTypes))].drop(['common_myeloid_progenitor','granulocyte_monocyte_progenitor_cell','megakaryocyte-erythroid_progenitor_cell'],axis=1)
    bp = bp.loc[map(lambda x: x in geneList,bp.index)]
    #bp1 = bp.apply(lambda x: x/max(x),axis=1)
    #bp1 = bp1[bp1.apply(lambda x: len([i for i in x if i > 0.2])==1,axis=1)]
    #bp = bp.loc[map(lambda x: x in bp1.index,bp.index)]
    return bp


def SVR(B,expr):
    svr = LinearSVR(C=500.0,epsilon=0.1,max_iter=10000,fit_intercept=False)
    svr.fit(B,expr)
    frac = pd.Series(svr.coef_,index = B.columns)
    y = pd.concat([expr,pd.Series(svr.predict(B),index=B.index)],join='inner',axis=1)
    print frac
    print y
    y = y.apply(lambda x: np.log10(x+0.01),axis=1)
    print stat.pearsonr(y.iloc[:,0],y.iloc[:,1])[0]


def QP0(B,expr):
    solvers.options['reltol'] = 0.0000001
    P = matrix(B.transpose().dot(B).values,tc='d')
    q = matrix(-1*expr.transpose().dot(B).values,tc='d')
    G = matrix(np.concatenate((-1*np.identity(B.shape[1]),np.expand_dims(np.ones(B.shape[1]),axis=0)),axis=0),tc='d')
    #h = matrix(0.0000000001*np.ones(B.shape[1]),tc='d')
    h = matrix([matrix(np.zeros(B.shape[1])),1],tc='d')
    sol = solvers.qp(P,q,G,h)
    frac = pd.Series(sol['x'],index = B.columns)
    y = pd.concat([expr,pd.Series(B.dot(frac),index=B.index)],join='inner',axis=1)
    print frac
    y = y.apply(lambda x: np.log2(x+1),axis=0)
    pcc = stat.pearsonr(y.iloc[:,0],y.iloc[:,1])[0]
    print y
    return frac


def QP(B,expr,samp,name,pdf):
    solvers.options['reltol'] = 0.0000001
    P = matrix(B.transpose().dot(B).values,tc='d')
    q = matrix(-1*expr.transpose().dot(B).values,tc='d')
    G = matrix(np.concatenate((-1*np.identity(B.shape[1]),np.expand_dims(np.ones(B.shape[1]),axis=0)),axis=0),tc='d')
    #h = matrix(0.0000000001*np.ones(B.shape[1]),tc='d')
    h = matrix([matrix(np.zeros(B.shape[1])),1],tc='d')
    sol = solvers.qp(P,q,G,h)
    frac = pd.Series(sol['x'],index = B.columns)
    y = pd.concat([expr,pd.Series(B.dot(frac),index=B.index)],join='inner',axis=1)
    print frac
    y = y.apply(lambda x: np.log2(x+1),axis=0)
    pcc = stat.pearsonr(y.iloc[:,0],y.iloc[:,1])[0]
    print y
    y.plot(samp,0,kind='scatter',color='blue',figsize=(15,12),s=100,alpha=0.8,title=samp+'\nPCC = %.4f' % pcc,fontsize=18)
    y1 = pd.concat([y,name],join='inner',axis=1)
    y1 = y1[(y1[samp]>7) | (y1[0]>7)]
    #for label,x,y in zip(y1['Description'],y1[samp],y1[0]):
    #    plt.annotate(label,xy=(x,y),xytext=(8,0),textcoords='offset points',ha='left',va='center')
    plt.xlabel('real')
    plt.ylabel('predicted')
    plt.savefig(pdf,format='pdf')
    plt.close()
    return frac

    
def KLdiverg(B,expr):
    n = B.shape[1]
    p = matrix(expr.values/sum(expr.values),tc='d')
    #p = matrix(expr.values,tc='d')
    B1 = matrix(B.values,tc='d')
    solvers.options['maxiters'] = 10000
    #solvers.options['reltol'] = 0.00001
    def F(x=None, z=None):
        if x is None: return 0, matrix(1.0, (n,1))
        if min(x) <= 0: return None
        prod = B1 * x
        norm = matrix(np.ones(B.shape[0]),tc='d').T * prod
        f = p.T * log(prod) - log(norm)
        grad = p.T * spdiag(prod**-1) * B1 - norm**-1 * matrix(np.ones(B.shape[0]),tc='d').T * B1
        if z is None: return f, grad
        rows = []
        for i in B.columns:
            rows.append(-p.T * spdiag(prod**-2) * spdiag(matrix(B[i].values,tc='d')) * B1)
        H = z[0] * (matrix(rows,tc='d') + norm**-2 * B1.T * matrix(np.ones((B.shape[0],B.shape[0])),tc='d') * B1)
        return f, grad, H
        
    G = matrix(np.concatenate((-1*np.identity(n),np.identity(n),np.expand_dims(np.ones(n),axis=0)),axis=0),tc='d')
    h = matrix([matrix(np.zeros(n)),matrix(np.ones(n)),1],tc='d')
    sol = solvers.cp(F,G,h)
    frac = pd.Series(sol['x'],index = B.columns)
    pred = pd.Series(B.dot(frac),index=B.index)
    y = pd.concat([expr/sum(expr)*100,pred/sum(pred)*100],join='inner',axis=1)
    #print frac
    print y
    y = y.apply(lambda x: np.log10(x+0.01),axis=1)
    print stat.pearsonr(y.iloc[:,0],y.iloc[:,1])[0]
    return frac

def KLdiverg2(B,expr):
    n = B.shape[1]
    expr1 = expr/sum(expr)
    fun = lambda x: np.dot(expr1,map(np.log,B.dot(x))) - np.log(sum(B.dot(x)))
    #fun = lambda x: -np.dot(expr,map(np.log,B.dot(x)))
    cons = ({'type':'ineq', 'fun': lambda x: 1-sum(x)})
    res = minimize(fun,np.ones(n)/float(n),method='SLSQP',bounds=((0,1),)*n,constraints=cons)

    frac = pd.Series(res.x,index = B.columns)
    pred = B.dot(frac)
    y = pd.concat([expr1,pd.Series(pred/sum(pred),index=B.index)],join='inner',axis=1)
    print frac
    print y
    y = y.apply(lambda x: np.log10(x+0.01),axis=1)
    print stat.pearsonr(y.iloc[:,0],y.iloc[:,1])[0]
    return frac/sum(frac)


def main():
    pd.options.display.max_rows=200
    B = getBase()
    #print B
    df,name = getMatrix()
    df = df[map(lambda x: x in B.index,df.index)]
    #df = df.apply(lambda x: np.log2(x+0.125),axis=0)
    #B = B.apply(lambda x: np.log2(x+0.125),axis=0)
    print df
    #SVR(B,df['13506-Mouse-2019-3-d-LPS-AI055'])
    #QP(B,df['13506-Mouse-2019-3-d-LPS-AI055'])
    #KLdiverg(B,df['13506-Mouse-2019-3-d-LPS-AI055'])
    fracs = df.apply(lambda x: QP0(B,x),axis=0)
    #fracs = df.apply(lambda x: KLdiverg(B,x),axis=0)
    print fracs
    fracs.iloc[:,:3].to_csv('BM.cellType.perc.csv')
    #sns.set(font_scale=2)
    #fracs.columns = map(lambda x: '-'.join(x.split('-')[-2:]),fracs.columns)
    #fracs.plot.pie(subplots=True, figsize=(18,12), layout=(2,3))
    #plt.savefig('cellType.perc.BM.pie.png')
    #plt.close()
    '''
    sns.set(font_scale=2)
    colors = ['black','blue','red','orange','purple']
    angles = [n / float(fracs.shape[1]) * 2 * pi for n in range(fracs.shape[1])]
    angles += angles[:1]
    fig,ax = plt.subplots(figsize=(15,15),subplot_kw=dict(projection='polar'))
    ax.set_theta_offset(pi/2)
    ax.set_theta_direction(-1)
    plt.xticks(angles[:-1],fracs.columns)
    ax.set_rlabel_position(0)
    for index,row in fracs.iterrows():
        color = colors.pop()
        vec = row.values.flatten().tolist()
        vec += vec[:1]
        ax.plot(angles, vec, color=color, lw=2, ls='-', label=index)
        ax.fill(angles, vec, color=color, alpha=0.1)

    plt.legend()
    plt.savefig('cellType.perc.radar.png')
    plt.close()
    '''
    sns.set(font_scale=1.8)
    from matplotlib.backends.backend_pdf import PdfPages
    pdf = PdfPages('predict_frac_cor.pdf')
    for col in df.columns:
        QP(B,df[col],col,name,pdf)
    pdf.close()


if __name__=='__main__':
    main()
