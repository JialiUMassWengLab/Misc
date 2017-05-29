#! /usr/bin/python

import re
import os
import sys
import csv
import glob
import pysam
import subprocess
from operator import add

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

def getEnds(bamFile):
    dirName = os.path.dirname(bamFile)
    sampleName = os.path.basename(bamFile).replace('Aligned.toTranscriptome.sorted.dedup.bam','')
    genes = getGeneList(dirName,sampleName,1)
    starts = {}
    ends = {}

    print len(genes)
    bamfile = pysam.AlignmentFile(bamFile,'rb')
    refs = bamfile.references
    lens = bamfile.lengths
    for i in range(len(refs)):
        if refs[i] in genes:
            starts[refs[i]] = [0] * lens[i]
            ends[refs[i]] = [0] * lens[i]

    i = 0
    for read in bamfile.fetch(until_eof=True):
        if read.is_proper_pair and read.is_read2 and not read.is_reverse and read.reference_name in genes:
            starts[read.reference_name][read.reference_start] += 1
            i += 1
            if i % 1000000 == 0:
                print 'processed %s pairs' % str(i)
        if read.is_proper_pair and read.is_read1 and read.is_reverse and read.reference_name in genes:
            ends[read.reference_name][read.reference_end-1] += 1
            
    bamfile.close()
    return starts,ends


def analyzeEnds(starts,ends,geneTypes,geneNames,sampleName):
    with open(sampleName+'.starts.counts','w') as out1:
        for key,value in starts.iteritems():
            if key in ends:
                if not key in geneNames:
                    geneNames[key] = 'NA'
                if not key in geneTypes:
                    geneTypes[key] = 'NA'
                out1.write('\t'.join([ key,','.join(map(str,starts[key])),','.join(map(str,ends[key])),str(sum(starts[key])),str(sum(ends[key])),geneTypes[key],geneNames[key] ])+'\n')

    subprocess.check_output('Rscript /home/jzhuang@ms.local/fragment/plot_spikes.R %s' % sampleName, shell=True, stderr=subprocess.STDOUT)
    #subprocess.check_output('rm %s' % sampleName+'.starts.counts', shell=True, stderr=subprocess.STDOUT)
    '''
    phase = {}
    for key,value in starts.iteritems():
        if not key in geneNames: geneNames[key] = 'NA'
        if not key in geneTypes: geneTypes[key] = 'NA'
        phase[key] = [0,0,0]
        for i in range(len(value)):
            phase[key][i % 3] += value[i]

        print '\t'.join([key]+map(str,phase[key]))
    '''
            
def main():
    #geneTypes,geneNames = getGeneTypes('/mnt/shares2/annotations/hg38/STAR/GRCh38_Gencode24_ERCC/gencode.v24.primary_assembly.annotation.with_ERCC.gtf')
    geneTypes,geneNames = getClasses()
    files = glob.glob(os.path.join(sys.argv[1],'*Aligned.toTranscriptome.sorted.dedup.bam'))
    all_starts = {}
    all_ends = {}

    bamfile = pysam.AlignmentFile(files[0],'rb')
    refs = bamfile.references
    lens = bamfile.lengths
    for i in range(len(refs)):
        all_starts[refs[i]] = [0] * lens[i]
        all_ends[refs[i]] = [0] * lens[i]

    for bamFile in files:
        starts,ends = getEnds(bamFile)
        #analyzeEnds(starts,ends,geneTypes,geneNames,sampleName)
        for key,value in starts.iteritems():
            if key in ends:
                all_starts[key] = map(add,starts[key],all_starts[key])
                all_ends[key] = map(add,ends[key],all_ends[key])

    analyzeEnds(all_starts,all_ends,geneTypes,geneNames,'Swift-Prince')


if __name__=='__main__':
    main()
