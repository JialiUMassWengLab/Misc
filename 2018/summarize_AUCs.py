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

def summarize(pheno):
    mean_aucs = {}
    files = glob.glob(os.path.join('para_select_files_Apr_revalid','*.%s.*tsv' % pheno))
    #files = glob.glob('*.%s.*tsv' % pheno)
    for fn in files:
        clf = os.path.basename(fn).split('.')[0]
        dtype = os.path.basename(fn).split('.')[3]
        tbl = pd.read_csv(fn,sep='\t',header=None)
        tbl['median'] = tbl.apply(lambda x: np.median(map(float,x[tbl.shape[1]-1].split(','))),axis=1)
        tbl = tbl.sort_values('median',ascending=False)
        if not clf in mean_aucs:
            mean_aucs.update({clf: {}})
        mean_aucs[clf].update({dtype: tbl.iloc[0,-1]})

    df = pd.DataFrame.from_dict(mean_aucs)
    df.to_csv('%s.AUCs.summary.csv' % pheno)
    print df

    sns.set(font_scale=2.2)
    fig,ax = plt.subplots(figsize=(15,18))
    sns.heatmap(df,cmap="YlGnBu",vmax=1.0,vmin=0.4,square=True,ax=ax,annot=True)
    plt.title('Best AUC = %.3f' % max(df.max()))
    plt.xticks(rotation=45,ha='right')
    plt.savefig('%s.AUCs.heatmap.png' % pheno)
    plt.close()


def main():
    summarize(sys.argv[1])



if __name__=='__main__':
    main()
