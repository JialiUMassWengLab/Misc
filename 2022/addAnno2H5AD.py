#! /usr/bin/python

import re
import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

def loadMetaData(fn):
    df = pd.read_csv(fn,sep='\t',header=None,index_col=0)
    df.index = list(map(lambda x: 'CB'+str(x),df.index))
    df.columns = ['FSH_dose','FSH','hCG_dose','hCG']
    df['FSH_dose'] = df['FSH_dose'].astype(int)
    df['treatment'] = df.apply(lambda x: '049:'+str(x['FSH_dose'])+'_'+'302:'+str(x['hCG_dose']),axis=1)
    mapping = {'049:0_302:0.0':'Vehicle','049:0_302:0.0375':'hCG_m','049:0_302:1.2':'hCG_hi',
               '049:1_302:0.0':'FSH_1','049:3_302:0.0':'FSH_3',
               '049:3_302:0.0375':'FSH_3 + hCG_m','049:3_302:1.2':'FSH_3 + hCG_hi',
               '049:10_302:0.0':'FSH_10','049:10_302:0.0375':'FSH_10 + hCG_m',
               '049:10_302:1.2':'FSH_10 + hCG_hi'}
    df['treatment1'] = df['treatment'].map(mapping)
    return df
    

def annotateH5AD(adata,meta):
    adata.obs['treatment'] = adata.obs['sample'].map(meta['treatment1'])
    dict1 = {'0':'Antral','1':'Proliferating','2':'Preantral','3':'Atretic'}
    dict2 = {'0':'Native','1':'Lhcgr-'}

    adata1 = sc.read_h5ad('scanpy.GC.recluster.h5ad')
    adata1.obs['GC_cluster'] = adata1.obs.apply(lambda x: dict1[x['leiden']] if x['leiden'] in dict1 else x['leiden'],axis=1)
    
    new_cluster_names = {'0':'GC','1':'GC','2':'GC','3':'GC','4':'Theca','5':'Empty','6':'Empty','7':'Stroma','8':'Epithelial','9':'BEC','10':'Immune','11':'Oocyte','12':'LEC','13':'Luteinizing GC','14':'Doublets'}
    adata.obs['ann_cluster'] = adata.obs['leiden'].map(new_cluster_names)
    shown_groups = ['GC','Theca','Stroma','Epithelial','BEC','Oocyte','Immune','LEC']
    #adata = adata[list(map(lambda x: x in shown_groups,adata.obs['ann_cluster'])),:]

    adata.obs = pd.concat([adata.obs,adata1.obs['GC_cluster']],join='outer',axis=1)
    adata.obs.loc[:,'new_cluster'] = adata.obs.apply(lambda x: 'GC-'+x['GC_cluster'] if x['ann_cluster']=='GC' else x['ann_cluster'],axis=1)
    print(adata)
    adata.write('scanpy.DIEM.annotated.h5ad')
    

def main():
    h5adFile = sys.argv[1]
    sc.settings.verbosity = 3
    sc.logging.print_header()
    adata = sc.read_h5ad(h5adFile)
    meta = loadMetaData('../metadata.tsv')
    
    annotateH5AD(adata,meta)
    
    
if __name__=='__main__':
    main()
