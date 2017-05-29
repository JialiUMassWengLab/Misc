#! /usr/bin/python

import re
import os
import sys
import csv
import glob
import pysam
import subprocess
from operator import add
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

def getClasses():
    classes = {}
    txns = {}
    allowed_types = ['protein_coding','lincRNA','processed_transcript','unprocessed_transcript']
    with open('/home/jzhuang@ms.local/TCGA/GeneInfo.csv','rU') as infile:
        reader = csv.DictReader(infile,fieldnames=['txnID','geneID','geneName','type','tissue','nonBlood'])
        for row in reader:
            if not row['tissue']=='' and row['type'] in allowed_types:
                classes[row['txnID']] = row['tissue']+'_ts'
            else:
                classes[row['txnID']] = row['type']

            txns[row['txnID']] = row['geneName']

    return classes,txns

def getGeneTypes(gtfFile):
    geneTypes = {}
    geneNames = {}
    with open(gtfFile,'rU') as infile:
        for line in infile:
            a = line.split()
            if not re.search(r'^#',line):
                match1 = re.search(r'gene_type "[\w]+";',line)
                match2 = re.search(r'transcript_id "[\w\.-]+";',line)
                match3 = re.search(r'gene_name "[\w]+";',line)
                if match2:
                    transcriptID = match2.group().split('"')[1]                    
                    if match1:
                        geneType = match1.group().split('"')[1]
                        if re.search(r'_pseudogene',geneType):
                            geneType = 'pseudogene'
                        if a[0]=='chrM' and not re.search(r'^Mt_',geneType):
                            geneType = 'Mt_' + geneType
                        if re.search(r'ERCC-',a[0]):
                            geneType = 'ERCC'

                        geneTypes[transcriptID] = geneType

                    if match3:
                        geneNames[transcriptID] = match3.group().split('"')[1]

    return geneTypes,geneNames

def getGeneList(dirName,sampleName,thre):
    genes = set([])
    with open(os.path.join(dirName,sampleName)+'.isoforms.results','rU') as infile:
        next(infile)
        for line in infile:
            a = line.split()
            if float(a[4]) >= thre:
                genes.add(a[0])

    return genes

def getPairs(bamFile):
    dirName = os.path.dirname(bamFile)
    sampleName = os.path.basename(bamFile).replace('Aligned.toTranscriptome.sorted.dedup.bam','')
    genes = getGeneList(dirName,sampleName,1)
    pairs = {}

    print len(genes)
    bamfile = pysam.AlignmentFile(bamFile,'rb')
    refs = bamfile.references
    lens = bamfile.lengths
    for i in range(len(refs)):
        if refs[i] in genes:
            pairs[refs[i]] = []

    i = 0
    for read in bamfile.fetch(until_eof=True):
        if read.is_proper_pair and read.is_read2 and not read.is_reverse and read.reference_name in genes:
            pairs[read.reference_name].append((read.reference_start+1,read.reference_start+read.template_length))
            i += 1
            if i % 1000000 == 0:
                print 'processed %s pairs' % str(i)
            
    bamfile.close()
    return pairs

def analyzePairs(pairs,geneTypes,geneNames,sampleName):
    with open(sampleName+'.pairs','w') as out1:
        from matplotlib.backends.backend_pdf import PdfPages
        #pdf = PdfPages(sampleName+'.pairHexbin.pdf')
        for key,value in pairs.iteritems():
            if not key in geneNames:
                geneNames[key] = 'NA'
            if not key in geneTypes:
                geneTypes[key] = 'NA'

            s = []
            e = []
            length = len(value)
            if length <= 100:
                continue
            for i in range(length):
                s.append(value[i][0])
                e.append(value[i][1])            
            out1.write('\t'.join([ key,','.join(map(str,s)),','.join(map(str,e)),str(length),geneTypes[key],geneNames[key] ])+'\n')
            '''
            plt.hexbin(s,e,extent=(1,length+1,1,length+1),gridsize=(100,100),cmap='inferno')
            plt.xlabel('Start position')
            plt.ylabel('End position')
            plt.title('%s (%s): %s fragment distribution' % (key,geneNames[key],geneTypes[key]))
            plt.savefig(pdf,format='pdf')
            plt.close()
            '''
        #pdf.close()
    subprocess.check_output('Rscript /home/jzhuang@ms.local/fragment/plot_hexbin.R %s' % sampleName,shell=True,stderr=subprocess.STDOUT)
    subprocess.check_output('rm %s.pairs' % sampleName,shell=True,stderr=subprocess.STDOUT)


def pairCov(pairs,geneTypes,geneNames,sampleName):
    with open(sampleName+'.bed','w') as out1:
        for key,value in pairs.iteritems():
            if not key in geneNames:
                geneNames[key] = 'NA'
            if not key in geneTypes:
                geneTypes[key] = 'NA'

            s = []
            e = []
            length = len(value)
            if length <= 100:
                continue
            for i in range(length):
                s.append(value[i][0])
                e.append(value[i][1])
                out1.write('\t'.join([key,str(s[i]-1),str(e[i])])+'\n')

    subprocess.check_output('sort -k1,1 %s.bed > tmp.bed' % sampleName,shell=True,stderr=subprocess.STDOUT)
    subprocess.check_output('mv tmp.bed %s.bed' % sampleName,shell=True,stderr=subprocess.STDOUT)
    subprocess.check_output('bedtools genomecov -i %s.bed -g /home/jzhuang@ms.local/fragment/txn.info -d | awk \'{if ($3>0) print}\' > %s.txn.cov.depth' % (sampleName,sampleName),shell=True,stderr=subprocess.STDOUT)

            
def main():
    #geneTypes,geneNames = getGeneTypes('/mnt/shares2/annotations/hg38/STAR/GRCh38_Gencode24_ERCC/gencode.v24.primary_assembly.annotation.with_ERCC.gtf')
    geneTypes,geneNames = getClasses()
    files = glob.glob(os.path.join(sys.argv[1],'*Aligned.toTranscriptome.sorted.dedup.bam'))
    all_pairs = {}

    bamfile = pysam.AlignmentFile(files[0],'rb')
    refs = bamfile.references
    lens = bamfile.lengths
    for i in range(len(refs)):
        all_pairs[refs[i]] = []

    for bamFile in files:
        pairs = getPairs(bamFile)
        for key,value in pairs.iteritems():
            all_pairs[key] += value

    analyzePairs(all_pairs,geneTypes,geneNames,'Swift-Prince')
    #pairCov(all_pairs,geneTypes,geneNames,'Genesis-97')


if __name__=='__main__':
    main()
