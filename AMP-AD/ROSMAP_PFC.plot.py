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

def main():
    adata = sc.read_h5ad('ROSMAP_PFC.h5ad')
    sc.settings.set_figure_params(dpi=240, facecolor='white')

    sns.set(font_scale=1.2,context='talk')
    '''
    fig,ax = plt.subplots(figsize=(15,15))
    sc.pl.umap(adata, color='new_cluster', legend_loc='on data', legend_fontsize='x-large', size=12, ax=ax)
    ax.set_facecolor('white')
    plt.savefig('ROSMAP_PFC.umap.pdf',format='pdf')

    marker_genes = ['AQP4','GAD1','MOG','FYB','CLDN5','PDGFRA','SLC17A7','CUX2']
    fig,axes = plt.subplots(4,2,figsize=(20,32))
    for i,gene in enumerate(marker_genes):
        sc.pl.violin(adata, gene, groupby='new_cluster', ax=axes[int(i/2)][i%2], rotation=40, size=2)
        labels = axes[int(i/2)][i%2].get_xticklabels()
        axes[int(i/2)][i%2].set_xticklabels(labels,ha='right',fontweight='bold')
        axes[int(i/2)][i%2].set_ylabel(gene,style='italic',fontweight='bold')
        axes[int(i/2)][i%2].set_facecolor('white')
    plt.tight_layout()
    #plt.savefig('ROSMAP_PFC.violin.pdf',format='pdf')
    plt.savefig('ROSMAP_PFC.violin.png')
    '''
    genes = pd.read_csv('../../Source/brain_scRNAseq/glia14.genes.csv',index_col=0,header=None)
    genes = genes[genes[1]!='GIMAP1']
    #sns.set(rc={'figure.facecolor':'white','axes.facecolor':'white'})
    fig,ax = plt.subplots(figsize=(45,8))
    dp = sc.pl.dotplot(adata, genes[1], groupby='new_cluster', standard_scale='var', return_fig=True, colorbar_title='Normalized expression', size_title='Percent expressed', figsize=(45,10),ax=ax)
    dp.style(largest_dot=600, cmap='Blues', dot_max=.8, dot_min=0)
    dp.make_figure()
    #dp.show()
    #sc.pl.dotplot(adata, genes[1], groupby='new_cluster', cmap='Blues', dot_max=.8, standard_scale='var', ax=ax)

    ax = dp.get_axes()['mainplot_ax']
    ax2 = dp.get_axes()['size_legend_ax']
    ax3 = dp.get_axes()['color_legend_ax']
    ax2.set_facecolor('white')
    ax2.set_aspect(.6)
    for item in ([ax2.title,ax3.title] + ax2.get_xticklabels() + ax3.get_xticklabels()):
        item.set_fontsize(10)
    labels = ax.get_xticklabels()
    celltypes = ax.get_yticklabels()
    ax.set_xticklabels(labels,ha='right',rotation=45,style='italic')
    ax.set_yticklabels(celltypes,fontweight='bold')
    ax.tick_params(axis='y',which='major',labelsize=32)
    ax.tick_params(axis='x',which='major',labelsize=18)
    ax.spines['bottom'].set_visible(True)
    ax.spines['left'].set_visible(True)
    ax.xaxis.set_ticks_position('bottom')
    ax.set_facecolor('white')
    plt.tight_layout()
    plt.savefig('ROSMAP_PFC.glia14_dot.pdf',format='pdf')


if __name__=='__main__':
    main()

