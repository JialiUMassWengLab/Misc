#! /usr/bin/python

import re
import os
import sys
import glob
import numpy as np
import pandas as pd
from sklearn import preprocessing

def makeFam4PRsice(study):
    info = pd.read_csv('/home/jzhuang/AMP/Output/demo_clinical.csv',index_col=0)
    info = info[info['race']=='W']
    info = info[list(map(lambda x: not x==3,info['CERAD']))]
    #info = info[list(map(lambda x: x=='AD' or x=='AD+' or x=='NCI',info['diagnosis']))]
    #info = info[list(map(lambda x: x=='AD+' or x=='NCI',info['diagnosis']))]
    #info = info[list(map(lambda x: x<=2 or x>=5,info['Braak']))]
    info['diagnosis'] = list(map(lambda x: 'AD' if x<=2 else 'NCI',info['CERAD']))
    mapping = pd.read_csv('/home/jzhuang/AMP/Derived/rnaseq2wgs.specimenMap.%s.csv' % study,index_col=0)
    mapping = mapping.drop_duplicates(['individualID'],keep='last')
    df = pd.concat([info[['sex','diagnosis']],mapping['specimenID_y']],join='inner',axis=1)
    df.index = df['specimenID_y']
    df = df.drop(['specimenID_y'],axis=1)
    print(df.shape)
    ids = pd.read_csv('/home/jzhuang/AMP/Derived/ROSMAP_PRSice/chr21_ROSMAP_PRSice.fam',header=None,sep=' ',index_col=0).index
    df = df[list(map(lambda x: x in ids,df.index))]

    df.insert(0,'FID',df.index)
    df.insert(1,'IID',df.index)
    df.insert(2,'2',[0]*df.shape[0])
    df.insert(3,'3',[0]*df.shape[0])
    df.loc[:,'sex'] = list(map(lambda x: 1 if x=='male' else 2 if x=='female' else 0,df['sex']))
    df.loc[:,'diagnosis'] = list(map(lambda x: 0 if x=='NCI' else 1,df['diagnosis']))
    df.to_csv(os.path.join('/home/jzhuang/AMP/Derived','%s_PRSice.fam' % study),index=None,header=None,sep=' ')


def makeCov4PLINK(study):
    info = pd.read_csv('/home/jzhuang/AMP/Output/demo_clinical.csv',index_col=0)
    mapping = pd.read_csv('/home/jzhuang/AMP/Derived/rnaseq2wgs.specimenMap.%s.csv' % study,index_col=0)
    mapping = mapping.drop_duplicates(['individualID'],keep='last')
    df = pd.concat([info,mapping['specimenID_y']],join='inner',axis=1)
    df.index = df['specimenID_y']
    S = pd.cut(df[df['ageDeath']!='90+']['ageDeath'].astype(float),3).astype(str)
    df.loc[:,'deathAge'] = df.apply(lambda x: '90+' if x['ageDeath']=='90+' else S[x.name],axis=1)
    df = df.drop(['individualID','tissue','Study','yearsEducation','ageDeath','specimenID_y'],axis=1)
    df = df.round(2).fillna(-9)

    le = preprocessing.LabelEncoder()
    for obj_col in df.select_dtypes(include=['object']).columns:
        le.fit(df[obj_col])
        pd.Series(le.classes_,name=obj_col).to_csv(os.path.join('/home/jzhuang/AMP/Derived','%s.%s.coding.csv' % (study,obj_col)))
        df.loc[:,obj_col] = pd.Series(le.transform(df[obj_col])+1,name=obj_col,index=df.index)

    df.insert(0,'FID',df.index)
    df.insert(1,'IID',df.index)
    df.loc[:,'sex'] = list(map(lambda x: 1 if x==2 else 2 if x==1 else 0,df['sex']))
    df.to_csv(os.path.join('/home/jzhuang/AMP/Output','%s_covariate4PLINK.txt' % study),index=None,sep=' ')


def makeCluster4PLINK(clusterFn,study):
    cluster = pd.read_csv(clusterFn,index_col=0)
    mapping = pd.read_csv('/home/jzhuang/AMP/Derived/rnaseq2wgs.specimenMap.%s.csv' % study,index_col=0)
    mapping = mapping.drop_duplicates(['individualID'],keep='last')

    oh = preprocessing.OneHotEncoder()
    df = pd.DataFrame(oh.fit_transform(cluster).toarray(),index=cluster.index)
    df.columns = list(map(lambda x: 'cluster'+str(x+1),df.columns))
    df = pd.concat([df,mapping['specimenID_y']],join='outer',axis=1).fillna(-1)
    df = df[df['specimenID_y']!=-1]

    info = pd.read_csv('/home/jzhuang/AMP/Output/demo_clinical.csv',index_col=0)
    info = info[list(map(lambda x: x in df.index,info.index))]
    info.loc[:,'Diagnosis'] = list(map(lambda x: re.sub(r'\+$','',x),info['diagnosis']))
    df1 = df.drop(['specimenID_y'],axis=1).copy()
    df2 = df.drop(['specimenID_y'],axis=1).copy()
    #d1 compares cluster to all other AD patients, NCI & MCI samples are ignored
    df1 = df1.apply(lambda x: x-1 if all(x==0) else x,axis=1)
    #d2 compares cluster to only NCIs, MCI & other AD patients are ignored
    df2 = df2.apply(lambda x: pd.Series(map(lambda y: y-1 if not y==1 else y,x),index=df2.columns) if not info.loc[x.name,'Diagnosis']=='NCI' else x,axis=1)
    df1.columns = list(map(lambda x: x+'d',df1.columns))
    df2.columns = list(map(lambda x: x+'c',df2.columns))
    df = pd.concat([df,df1,df2],join='inner',axis=1)

    df.index = df['specimenID_y']
    df = df.drop(['specimenID_y'],axis=1).astype(int).apply(lambda x: x+1,axis=0)

    df.insert(0,'FID',df.index)
    df.insert(1,'IID',df.index)
    prefix = os.path.basename(clusterFn).replace('.labels.csv','')
    df.to_csv(os.path.join('/home/jzhuang/AMP/Output',prefix+'.%s.cluster4PLINK.txt' % study),index=None,sep=' ')


def makeCluster4PLINK2(clusterFn,study):
    cluster = pd.read_csv(clusterFn,index_col=0)
    mapping = pd.read_csv('/home/jzhuang/AMP/Derived/rnaseq2wgs.specimenMap.%s.csv' % study,index_col=0)
    mapping = mapping.drop_duplicates(['individualID'],keep='last')

    oh = preprocessing.OneHotEncoder()
    df = pd.DataFrame(oh.fit_transform(cluster).toarray(),index=cluster.index)
    df.columns = list(map(lambda x: 'cluster'+str(x+1),df.columns))
    df = pd.concat([df,mapping['specimenID_y']],join='outer',axis=1).fillna(-1)
    df = df[df['specimenID_y']!=-1]
    df.index = df['specimenID_y']
    df = df.drop(['specimenID_y'],axis=1).astype(int).apply(lambda x: x+1,axis=0)

    df.insert(0,'FID',df.index)
    df.insert(1,'IID',df.index)
    prefix = os.path.basename(clusterFn).replace('.labels.csv','')
    df.to_csv(os.path.join('/home/jzhuang/AMP/Output',prefix+'.labels4PLINK.txt'),index=None,sep=' ')


def makePathLevel4PLINK(levelFn,study):
    #study = os.path.basename(levelFn).split('_')[0]
    mapping = pd.read_csv('/home/jzhuang/AMP/Derived/rnaseq2wgs.specimenMap.%s.csv' % study,index_col=0)
    mapping = mapping.drop_duplicates(['individualID'],keep='last')
    df = pd.read_csv(levelFn,index_col=0)
    df.index = list(map(lambda x: re.sub(r'^X','',x),df.index))
    #df = df.apply(lambda x: pd.Series(np.where(abs((x-np.mean(x)))/np.std(x) > 3, np.mean(x)+np.sign(x-np.mean(x))*3*np.std(x), x), index=df.index), axis=0)
    #df = df.apply(lambda x: 2**x,axis=0)

    df = pd.concat([df,mapping['specimenID_y']],join='outer',axis=1).fillna(-9)
    df = df[df['specimenID_y']!=-9]
    df.index = df['specimenID_y']
    df = df.drop(['specimenID_y'],axis=1)
    df.insert(0,'FID',df.index)
    df.insert(1,'IID',df.index)
    df.to_csv(os.path.join('/home/jzhuang/AMP/Output','%s_moduleLevels4PLINK.txt' % study),index=None,sep=' ')


def makeCompLevel4PLINK(Wmat,study):
    mapping = pd.read_csv('/home/jzhuang/AMP/Derived/rnaseq2wgs.specimenMap.%s.csv' % study,index_col=0)
    mapping = mapping.drop_duplicates(['individualID'],keep='last')
    W = pd.read_csv(Wmat,sep='\t',index_col=0).astype('float')
    W.columns = list(map(lambda x: 'comp'+str(int(x)+1),W.columns))

    df = pd.concat([W,mapping['specimenID_y']],join='outer',axis=1).fillna(-9)
    df = df[df['specimenID_y']!=-9]
    df.index = df['specimenID_y']
    df = df.drop(['specimenID_y'],axis=1)
    df.insert(0,'FID',df.index)
    df.insert(1,'IID',df.index)
    df.to_csv(os.path.join('/home/jzhuang/AMP/Output','%s_compLevels4PLINK.txt' % study),index=None,sep=' ')


def makeWithinCluster4PLINK(study):
    mapping = pd.read_csv('/home/jzhuang/AMP/Derived/rnaseq2wgs.specimenMap.%s.csv' % study,index_col=0)
    mapping = mapping.drop_duplicates(['individualID'],keep='last')
    files = glob.glob('/home/jzhuang/AMP/Output/%s_geneLevel_limma_normed_counts.exprRatio.*.clusters.csv' % study)

    array = []
    for fn in files:
        df = pd.read_csv(fn,index_col=0)
        path = os.path.basename(fn).split('.')[-3]
        df.columns = list(map(lambda x: path+'_%s' % x,df.columns))
        array.append(df)

    df = pd.concat(array,join='outer',axis=1).fillna(-1)
    df = pd.concat([df,mapping['specimenID_y']],join='outer',axis=1).fillna(-1)
    df = df[df['specimenID_y']!=-1]
    df.index = df['specimenID_y']
    df = df.drop(['specimenID_y'],axis=1).astype(int).apply(lambda x: x+1,axis=0)

    df.insert(0,'FID',df.index)
    df.insert(1,'IID',df.index)
    df.to_csv(os.path.join('/home/jzhuang/AMP/Output','%s_withinPathClusters4PLINK.txt' % study),index=None,sep=' ')


def makeWithinCluster4PLINK2(study):
    mapping = pd.read_csv('/home/jzhuang/AMP/Derived/rnaseq2wgs.specimenMap.%s.csv' % study,index_col=0)
    mapping = mapping.drop_duplicates(['individualID'],keep='last')
    df = pd.read_csv('/home/jzhuang/AMP/Output/%s_geneLevel_limma_normed_counts.patientCorrCluster.csv' % study,index_col=0)
    for module in df.columns:
        labels = df[module].value_counts().keys()
        df[module] = list(map(lambda x: 0 if x==labels[0] else 1 if x==labels[1] else -1,df[module]))

    df = pd.concat([df,mapping['specimenID_y']],join='outer',axis=1).fillna(-1)
    df = df[df['specimenID_y']!=-1]
    df.index = df['specimenID_y']
    df = df.drop(['specimenID_y'],axis=1).astype(int).apply(lambda x: x+1,axis=0)

    df.insert(0,'FID',df.index)
    df.insert(1,'IID',df.index)
    df.to_csv(os.path.join('/home/jzhuang/AMP/Output','%s_withinModuleClusters4PLINK.txt' % study),index=None,sep=' ')


def main():
    #makeFam4PRsice('ROSMAP')
    #makeCov4PLINK('ROSMAP')
    #makeCluster4PLINK2(sys.argv[1],'ROSMAP')
    #makePathLevel4PLINK(sys.argv[1],'ROSMAP')
    #makeCompLevel4PLINK(sys.argv[1],'ROSMAP')
    #makeWithinCluster4PLINK('ROSMAP')
    makeWithinCluster4PLINK2('ROSMAP')
    

if __name__=='__main__':
    main()

    
