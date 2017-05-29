#! /usr/bin/python

import re
import os
import sys
import glob
import pandas as pd

def consolidate():
    Ftest = {}
    files = glob.glob(os.path.join('specific_genes_simulation','Ftest.genes.*'))
    for fn in files:
        df = pd.read_csv(fn,delimiter='\t')
        for index,row in df.iterrows():
            if row['Description'] in Ftest:
                Ftest[row['Description']] += row['scores']
            else:
                Ftest[row['Description']] = row['scores']

    LR = {}
    files = glob.glob(os.path.join('specific_genes_simulation','LR.genes.*'))
    for fn in files:
        df = pd.read_csv(fn,delimiter='\t')
        for index,row in df.iterrows():
            if row['Description'] in LR:
                LR[row['Description']] += row['coef']
            else:
                LR[row['Description']] = row['coef']

    with open('Ftest.geneList.txt','w') as out:
        for gene in sorted(Ftest,key=Ftest.get,reverse=True)[:50]:
            print gene,Ftest[gene]
            out.write(gene+'\n')

    print '\n----------\n'

    with open('LR.geneList.txt','w') as out:
        for gene in sorted(LR,key=LR.get,reverse=True)[:50]:
            print gene,LR[gene]
            if LR[gene] > 0:
                out.write(gene+'\n')


def main():
    consolidate()


if __name__=='__main__':
    main()
