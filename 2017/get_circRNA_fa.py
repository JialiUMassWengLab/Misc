#! /usr/bin/python

import csv
import sys
import os
import glob
import re
import math
import subprocess
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

def revcomp(s): 
    basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 
    letters = list(s) 
    letters = [basecomplement[base] for base in letters] 
    complement = ''.join(letters)
    return complement[::-1]

def bedOverlap(info,threshold,bedFile,options):
    with open('tmp.bed','w') as out:
        for key,value in info.iteritems():
            a = key.split('|')
            if value['juncClass']=='circular' and value['counts'] >= threshold:
                lower = a[1]
                upper = a[4]
                if int(a[1]) > int(a[4]):
                    lower = a[4]
                    upper = a[1]

                out.write('\t'.join([a[0],lower,upper,key,'.',a[5]])+'\n')

    subprocess.check_output('bedtools intersect -a tmp.bed -b %s %s | uniq > tmp' % (bedFile,options),shell=True,stderr=subprocess.STDOUT)

def getInfo(summaryFile):
    info = {}
    with open(summaryFile,'rU') as infile:
        for line in infile:
            a = line.split()
            name = '|'.join(a[1:7])
            juncClass = 'circular'
            if not a[1]==a[4] or not a[3]==a[6]:
                juncClass = 'fusion'
            if a[3]=='+' and a[6]=='+' and a[1]==a[4] and int(a[2])<int(a[5]):
                juncClass = 'intra-fusion'
            if a[3]=='-' and a[6]=='-' and a[1]==a[4] and int(a[2])>int(a[5]):
                juncClass = 'intra-fusion'
            chimericInfo = {}
            chimericInfo.update({'counts':int(a[0]),'juncType':a[7],'juncClass':juncClass})
            info.update({name:chimericInfo})

    return info

def getCircRNASeq(info,threshold,destPrefix,mappingFile,refSeqFile,twoBitFile):
    from Bio.Seq import Seq

    genes = {}
    with open(mappingFile,'rU') as infile:
        for line in infile:
            a = line.split(',')
            genes[a[1][:-1]] = a[0]

    circSeqs = {}
    bedOverlap(info,threshold,refSeqFile,'-s -f 0.95 -wo')
    with open('tmp','rU') as infile:
        for line in infile:
            a = line.split()
            starts = a[17].split(',')[:-1]
            ends = a[16].split(',')[:-1]
            start_match = -1
            end_match = -1
            for i in range(len(starts)):
                starts[i] = int(starts[i])+int(a[7])
                if abs(starts[i]-int(a[1])) <= 5:
                    start_match = i
            for i in range(len(ends)):
                ends[i] = starts[i]+int(ends[i])
                if abs(ends[i]-int(a[2])) <= 5:
                    end_match = i
            
            circKey = '|'.join([ '>circ',a[3],genes[a[9]] ])
            if start_match > -1 and end_match > -1 and not circKey in circSeqs:
                j = start_match
                #Generate exon fasta sequences
                while j <= end_match:
                    start = str(starts[j]+1)
                    end = str(ends[j]+1)
                    subprocess.check_output('twoBitToFa -noMask %s:%s:%s-%s stdout >> tmp.fa' % (twoBitFile,a[6],start,end),
                                            shell=True,stderr=subprocess.STDOUT)
                    j += 1
                #Read exon seqs and generate circRNA seq
                with open('tmp.fa','rU') as faFile:
                    lines = faFile.readlines()
                    exons = []
                    seq = ''
                    k = 0
                    while k < len(lines):
                        if re.search(r'^>',lines[k]):
                            if not seq == '':
                                exons.append(seq)
                            seq = ''
                        else:
                            seq += lines[k][:-1]
                        k += 1

                    exons.append(seq)
                    midseq = ''.join(exons[:-1])
                    lastseq = exons[-1]
                    if len(exons[-1]) > 200:
                        lastseq = exons[-1][:200]
                    circSeqs[circKey] = ''.join([exons[-1],midseq,lastseq])
                    if a[5] == '-':
                        #circSeqs[circKey] = Seq(circSeqs[circKey]).reverse_complement()
                        circSeqs[circKey] = revcomp(circSeqs[circKey])
                
                #print circKey
                #print circSeqs[circKey]
                os.remove('tmp.fa')
                    
    subprocess.check_output('rm tmp.bed tmp',shell=True,stderr=subprocess.STDOUT)
    with open(destPrefix+'.circRNA.fa','w') as out:
        for key,value in circSeqs.iteritems():
            out.write(key+'\n')
            out.write(value+'\n')

    with open(destPrefix+'.transcript2gene.mapping','w') as out:
        for key,value in genes.iteritems():
            out.write('\t'.join([value,key])+'\n')
        for key,value in circSeqs.iteritems():
            geneID = key.split('|')[-1]
            out.write('\t'.join([geneID,key[1:]])+'\n')


def main():
    '''
    files = glob.glob(os.path.join(sys.argv[1],'*Chimeric.out.junction'))
    for filename0 in files:
        filename = os.path.basename(re.sub(r'junction','summary',filename0))
        subprocess.check_output('cat %s | awk \'$1!="chrM" && $4!="chrM" && $7>0 && $8+$9<=5 {print $1,$2,$3,$4,$5,$6,$7,$8,$9}\' | sort | uniq -c | sort -k1,1rn > %s' % (filename0,filename),shell=True,stderr=subprocess.STDOUT)

        info=getInfo(filename)
        destPrefix = re.sub(r'Chimeric.out.junction','',filename0)
        print filename,circRatio(info,1,destPrefix,True)
    '''
    info = getInfo('TH08Chimeric.out.summary')
    getCircRNASeq(info,10,'TH08','/mnt/shares2/annotations/hg38/gencode_hg38.transcript2gene.mapping.csv','/mnt/shares2/annotations/hg38/gencode.v24.basic.annotation.gtf.gz.bed','/mnt/shares2/annotations/hg38/hg38.2bit')

if __name__=='__main__':
    main()
            
