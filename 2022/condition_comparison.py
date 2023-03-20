#! /usr/bin/python

import re
import os
import sys
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

def plotVenn(celltype,comparison):
    df1 = pd.read_csv('../Optigon_snRNAseq_RatOvaries20_June2022/scanpy.%s.DE.%s.csv' % (celltype,comparison),index_col=3)
    conditions = comparison.split('_vs_')
    conditions[1] = re.sub(r'\.0$','',conditions[1])
    df2 = pd.read_csv('DESeq_results/%s_T28_vs_%s_T28/gene.results.merged.count.DESeq2.Wald.FC1.BH.P0.05.txt' % (conditions[0],conditions[1]),sep='\t',index_col=0)
    df2['significant'] = df2.apply(lambda x: abs(x[0])>.8 and x[-2]<.05,axis=1)
    df = pd.concat([df1,df2],join='inner',axis=1)
    df1 = df1.loc[map(lambda x: x in df.index,df1.index),:]
    df2 = df2.loc[map(lambda x: x in df.index,df2.index),:]

    up1 = set(df1[(df1['logFC']>0) & (df1['significant'])].index)
    up2 = set(df2[(df2['Log2FC %s_T28 vs %s_T28' % (conditions[0],conditions[1])]>0) & (df2['significant'])].index)
    down1 = set(df1[(df1['logFC']<0) & (df1['significant'])].index)
    down2 = set(df2[(df2['Log2FC %s_T28 vs %s_T28' % (conditions[0],conditions[1])]<0) & (df2['significant'])].index)
    sns.set(font_scale=1.5)
    from matplotlib_venn import venn2
    fig,axes = plt.subplots(2,1,figsize=(10,18))
    venn2([up1,up2],('snRNAseq_%s' % celltype,'bulk-RNAseq'),ax=axes[0])
    axes[0].set_title('Up-regulated genes')
    venn2([down1,down2],('snRNAseq_%s' % celltype,'bulk-RNAseq'),ax=axes[1])
    axes[1].set_title('Down-regulated genes')
    plt.savefig('sn_bulk_comparison.%s.%s.venn.pdf' % (comparison,celltype))
    plt.close()    
    

def compareFoldChange(celltype,compound1,compound2):
    time = compound1.split('_')[1]
    df1 = pd.read_csv('scanpy.%s.DE.%s_vs_vehicle_%s.csv' % (celltype,compound1,time),index_col=0)
    df2 = pd.read_csv('scanpy.%s.DE.%s_vs_vehicle_%s.csv' % (celltype,compound2,time),index_col=0)
    genes = set(list(df1[df1['significant']].index) + list(df2[df2['significant']].index))
    df = pd.concat([df1['logFC'],df2['logFC']],join='inner',axis=1)
    df.columns = [compound1,compound2]
    df = df.loc[map(lambda x: x in genes,df.index),:]
    df['deviated'] = df.apply(lambda x: True if abs(x[0]-x[1])>2 else False,axis=1)
    #df['concordant'] = df.apply(lambda x: True if x.name in df1[df1['significant']].index and x.name in df2[df2['significant']].index else False,axis=1)
    df['concordant'] = list(map(lambda x: 'concordant' if df1.loc[x,'significant']==df2.loc[x,'significant'] else 'discordant',df.index))
    print(df.shape)

    #sig = df[df['deviated']==True]
    #sig = df.loc[['Lhcgr','Cyp19a1','Cyp11a1','Hsd3b2','Hsd17b1','Inhba','Inhbb','Npas2','Sema3d','Tmtc2','Prlr','Gja1','Tox','Nr5a2','Foxo1','Foxp2','Tgfb2','Lox','Adamts5','Gls','Tmcc3'],:]
    sig = df.loc[df.apply(lambda x: abs(x[0])>3 or abs(x[1])>3,axis=1),:]
    if sig.shape[0] > 35:
        sig = sig.loc[sig.apply(lambda x: x[0]>3 or x[1]>3 or x[0]<-4 or x[1]<-4,axis=1),:]

    from adjustText import adjust_text
    sns.set(font_scale=1.2,context='talk')
    fig,ax = plt.subplots(figsize=(15,12))
    sns.scatterplot(x=compound1,y=compound2,hue='concordant',data=df,s=150,linewidth=1.5,palette={'discordant':'darkgrey','concordant':'red'},ax=ax)
    ax.axline((0,0),slope=1,color='red',ls='--',lw=2)
    texts = []
    for label,x,y in zip(sig.index,sig[compound1],sig[compound2]):
        texts.append(plt.text(x,y,label,ha='center',size=18))
    adjust_text(texts,autoalign='xy',arrowprops=dict(arrowstyle="simple, head_width=0.25, tail_width=0.05", color='r', lw=0.5, alpha=0.8))
    plt.legend(loc='upper left')
    plt.xlabel('Log2 Fold change %s_vs_vehicle' % compound1)
    plt.ylabel('Log2 Fold change %s_vs_vehicle' % compound2)
    plt.savefig('compound_comparison.%s.%s_VS_%s.logFC.scatter.png' % (celltype,compound1,compound2))
    plt.close()    
    

def main():
    #plotVenn('GC',sys.argv[1])
    #compareFoldChange('GC','FSH_T28','4221:6_T28')
    #compareFoldChange('GC','FSH_T28','2599:0.5_T28')
    #compareFoldChange('GC','FSH_T28','2599:6_T28')
    #compareFoldChange('GC','FSH_T28','0353:10_T28')
    #compareFoldChange('GC','FSH_T28','0353:30_T28')
    #compareFoldChange('GC','FSH_T28','5690:1_T28')
    #compareFoldChange('GC','FSH_T28','5690:0.3_T28')

    compareFoldChange('GC','FSH_T57','4221:6_T57')
    compareFoldChange('GC','FSH_T57','2599:0.5_T57')
    compareFoldChange('GC','FSH_T57','2599:6_T57')
    compareFoldChange('GC','FSH_T57','0353:10_T57')
    compareFoldChange('GC','FSH_T57','0353:30_T57')
    compareFoldChange('GC','FSH_T57','5690:1_T57')
    compareFoldChange('GC','FSH_T57','5690:0.3_T57')
    

if __name__=='__main__':
    main()
