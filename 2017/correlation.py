#! /usr/bin/python

import sys
import os
import re
import glob
import numpy
import scipy.stats
import subprocess
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import seaborn

def normalize(file_dir):
    files = glob.glob(os.path.join(file_dir,'*.proteinCoding.*.results'))
    for fn in files:
        tpmSum = 0.0
        with open(fn,'rU') as infile:
            next(infile)
            for line in infile:
                a = line.split('\t')
                tpmSum += float(a[5])

        with open(fn,'rU') as infile:
            out = open(fn+'.norm','w')
            for line in infile:
                if re.search(r'gene_id',line):
                    out.write(line)
                else:
                    a = line[:-1].split('\t')
                    a[5] = '{0:.2f}'.format( float(a[5])/tpmSum*1000000 )
                    out.write('\t'.join(a)+'\n')


def getMatrix(file_dir, log=False):
    samples = []
    files = glob.glob(os.path.join(file_dir,'*.proteinCoding.genes.results.norm'))
    matrix = numpy.zeros([len(files),len(files)])
    for fn in files:
        sample = os.path.basename(fn).replace('.proteinCoding.genes.results.norm','')
        samples.append(sample)

    i = 0
    for fn1 in files:
        subprocess.check_output('cut -f1,6 %s | tail -n+2 | sort -k1,1 > tmp1' % fn1, shell=True, stderr=subprocess.STDOUT)
        row = []
        for fn2 in files:
            subprocess.check_output('cut -f1,6 %s | tail -n+2 | sort -k1,1 > tmp2' % fn2, shell=True, stderr=subprocess.STDOUT)
            subprocess.check_output("join -j 1 tmp1 tmp2 | awk '{OFS=\"\t\"; if (($2>0)||($3>0)) print}' > tmp3", shell=True, stderr=subprocess.STDOUT)
            e1 = []
            e2 = []
            with open('tmp3','rU') as infile:
                for line in infile:
                    a = line.split()
                    if log:
                        e1.append(numpy.log10(float(a[1])+0.5))
                        e2.append(numpy.log10(float(a[2])+0.5))
                    else:
                        e1.append(float(a[1]))
                        e2.append(float(a[2]))
                    
            print scipy.stats.pearsonr(e1,e2)
            row.append(scipy.stats.pearsonr(e1,e2)[0])

        print row
        matrix[i,] = row
        i += 1
        
    subprocess.check_output('rm tmp*', shell=True, stderr=subprocess.STDOUT)
    df = pd.DataFrame(matrix,index=samples,columns=samples)
    return df


def main():
    normalize(sys.argv[1])

    fig, ax = plt.subplots()
    df = getMatrix(sys.argv[1],True)
    seaborn.heatmap(df)
    plt.title('Pearson Correlation Coefficients log10 TPM')
    plt.savefig(os.path.join(sys.argv[1],sys.argv[2]+'.pdf'))
    df.to_excel(os.path.join(sys.argv[1],sys.argv[2]+'.xlsx'))


if __name__=='__main__':
    main()
