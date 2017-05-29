#! /usr/bin/python

import re
import sys
import os
import csv
import json
import glob
import pysam
import numpy
import pandas as pd
import statsmodels.api as sm

def getClasses():
    classes = {}
    txns = {}
    allowed_types = ['protein_coding','lincRNA','processed_transcript','unprocessed_transcript']
    with open('/mnt/shares2/annotations/hg38/GeneInfo.csv','rU') as infile:
        reader = csv.DictReader(infile,fieldnames=['txnID','geneID','geneName','type','tissue','nonBlood'])
        for row in reader:
            if not row['tissue']=='' and row['type'] in allowed_types:
                classes[row['geneName']] = row['tissue']+'_ts'
            else:
                classes[row['geneName']] = row['type']

            #txns[row['txnID']] = row['geneName']
            txns[row['geneID']] = row['geneName']

    return classes,txns


def countFrags(classes,names,BamFile):
    counts = {}
    bamfile = pysam.AlignmentFile(BamFile,'rb')
    for read in bamfile.fetch(until_eof=True):
        if read.is_read1 and read.is_proper_pair and read.reference_id >= 0 and read.reference_id == read.next_reference_id:
            if read.reference_name in counts:
                counts[read.reference_name] += 1
            else:
                counts[read.reference_name] = 1

    bamfile.close()

    geneFrags = {}
    for transcriptID,count in counts.iteritems():
        name = names[transcriptID]
        if not name in classes:
            continue
        if (name in geneFrags and count > geneFrags[name]) or not name in geneFrags:
            geneFrags[name] = count

    return geneFrags

def countFrags1(classes,names,exprFile):
    geneFrags = {}
    with open(exprFile,'rU') as infile:
        next(infile)
        for line in infile:
            a = line.split()
            name = names[a[0]]
            if not name in classes:
                continue
            elif not name in geneFrags or float(a[4]) > geneFrags[name]:
                geneFrags[name] = float(a[4])

    return geneFrags


def summarizeFrags(classes,names,dirName):
    runName = os.path.basename(dirName)
    if runName == '':
        runName = os.path.basename(os.path.dirname(dirName))
    
    liverGenes = set([])
    with open('/mnt/shares2/annotations/hg38/liverGenes','rU') as infile:
        for line in infile:
            liverGenes.add(line[:-1])
    '''
    nonBlood = set([])
    with open('/home/jzhuang@ms.local/TCGA/NonBlood.txt','rU') as infile:
        for line in infile:
            a = line.split()
            nonBlood.add(a[1])
    '''
    bamfiles = glob.glob(os.path.join(dirName,'*Aligned.toTranscriptome.sorted.dedup.bam'))
    summary = { 'total': {}, 'proteinCoding': {}, 'liverSpecific': {}, 'tissueSpecific': {} }
    for bamfile in bamfiles:
        exprFile = bamfile.replace('Aligned.toTranscriptome.sorted.dedup.bam','.dedup.genes.results')
        #fragCounts = countFrags(classes,names,bamfile)
        fragCounts = countFrags1(classes,names,exprFile)
        ts_f = []
        nb_f = []
        lv_f = []
        tissue_f = {}
        for gene,count in fragCounts.iteritems():
            #if gene in nonBlood:
            #    nb_f.append(count)
            if gene in liverGenes:
                lv_f.append(count)
            if re.search(r'_ts$',classes[gene]):
                ts_f.append(count)
                tissue = re.sub(r'_ts$','',classes[gene])
                if tissue in tissue_f:
                    tissue_f[tissue].append(count)
                else:
                    tissue_f[tissue] = [count]

        filename = os.path.basename(bamfile).replace('Aligned.toTranscriptome.sorted.dedup.bam','')
        summary['total'].update({ filename: str(int(sum(fragCounts.values()))) })
        summary['proteinCoding'].update({filename: ','.join([ str(sum(i>=2 for i in fragCounts.values())), str(sum(i>=5 for i in fragCounts.values())), str(sum(i>=10 for i in fragCounts.values())) ])})
        #summary['nonBlood'].update({filename: ','.join([ str(sum(i>0 for i in nb_f)), str(sum(i>5 for i in nb_f)), str(sum(i>10 for i in nb_f)), str(sum(i>20 for i in nb_f)), str(sum(i>50 for i in nb_f)) ])})
        summary['liverSpecific'].update({filename: ','.join([ str(sum(i>=2 for i in lv_f)), str(sum(i>=5 for i in lv_f)), str(sum(i>=10 for i in lv_f)) ])})
        summary['tissueSpecific'].update({filename: ','.join([ str(sum(i>=2 for i in ts_f)), str(sum(i>=5 for i in ts_f)), str(sum(i>=10 for i in ts_f)) ])})
        '''
        for tissue,vec in tissue_f.iteritems():
            if not tissue in summary: summary.update({ tissue: {} })
            if len(vec) > 0:
                summary[tissue].update({filename: ','.join([ str(sum(i>0 for i in vec)), str(sum(i>5 for i in vec)), str(sum(i>10 for i in vec)), str(sum(i>20 for i in vec)), str(sum(i>50 for i in vec)) ])})

        for tissue,geneNum in tsNum.iteritems():
            if not tissue in tissue_f:
                if not tissue in summary: summary.update({ tissue: {} })
                summary[tissue].update({filename: '0,0,0,0,0'})
        '''
    df = pd.DataFrame.from_dict(summary)
    '''
    newcol = []
    for tissue in df.columns:
        if tissue in tsNum:
            newcol.append(tissue+' (%s)' % str(tsNum[tissue]))
        else:
            newcol.append(tissue)

    df.columns = newcol
    '''
    df.to_csv(os.path.join(dirName,runName+'_genesDetected.summary.tsv'),sep=',')
    return summary        


def summarizeERCC(dirName):
    runName = os.path.basename(dirName)
    if runName == '':
        runName = os.path.basename(os.path.dirname(dirName))

    summary = { 'Sensitivity':{}, 'LOD':{}, 'ERCC multi':{}, 'ERCC detected':{}, 'Chrom multi':{} }
    files = glob.glob(os.path.join(dirName,'*.ercc.stat.txt'))
    for fn in files:
        filename = os.path.basename(fn).replace('.filtered.bam.ercc.counts.ercc.stat.txt','')
        try:
            infile = open(fn,'rU')
            a = infile.readlines()[1].split()
            summary['Sensitivity'].update({ filename: '{0:.3f}'.format(float(a[1])) })
            summary['LOD'].update({ filename: '{0:.3f}'.format(float(a[2])) })
            summary['ERCC multi'].update({ filename: '{0:.1f}'.format(float(a[3])) })
            summary['ERCC detected'].update({filename: int(a[4])})
            infile.close()
        except IOError:
            print 'Cannot open file %s!' % fn
            
        with open(fn.replace('filtered.bam.ercc.counts.ercc.stat.txt','stats.json'),'rU') as json_file:
            data = json.load(json_file)
            summary['Chrom multi'].update({ filename: '{0:.1f}'.format(float(data['multStats']['median'])) })
            
    df = pd.DataFrame.from_dict(summary)
    df.to_csv(os.path.join(dirName,runName+'_ercc.summary.tsv'),sep='\t')
    return summary


def main():
    classes,names = getClasses()
    summarizeFrags(classes,names,sys.argv[1])
    summarizeERCC(sys.argv[1])


if __name__=='__main__':
    main()
