#! /usr/bin/python

import re
import os
import sys
import glob
import math
import subprocess
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

def preprocess(sample):
    subprocess.check_output('export IGDATA=/mnt/shares2/blast_db/', shell=True, stderr=subprocess.STDOUT)
    subprocess.check_output('export BLASTDB=/mnt/shares2/blast_db:/mnt/shares2/blast_db/database', shell=True, stderr=subprocess.STDOUT)
    if not os.path.exists('%s.igBLAST.results' % sample):
        subprocess.check_output('~/bin/ncbi-igblast-1.6.1/bin/igblastn -germline_db_V imgt_v1.fa -germline_db_J imgt_j1.fa -germline_db_D imgt_d1.fa -organism human -domain_system imgt -query %s.candidate.reads.fa -out %s.igBLAST.results -outfmt 7 -auxiliary_data /mnt/shares2/blast_db/optional_file/human_gl.aux -num_threads 8 -evalue 0.000001' % (sample,sample), shell=True, stderr=subprocess.STDOUT)

    if not os.path.exists('%s.BLAST.results' % sample):
        subprocess.check_output('blastn -query %s.candidate.reads.fa -task blastn -db imgtrefseq -out %s.BLAST.results -outfmt "6 qseqid sseqid scomnames qstart qend length pident evalue score" -evalue 0.000001 -num_threads 8 -max_hsps 2 -max_target_seqs 10' % (sample,sample), shell=True, stderr=subprocess.STDOUT)


def parseBlastOutput0(BlastOut):
    anno = {}
    with open('/mnt/shares/Users/jzhuang/LIGM-DB/imgtrefseq.fasta','rU') as infile:
        for line in infile:
            if re.search(r'^>',line):
                a = line.split('|')
                if a[2] == 'Homo sapiens' and not a[4] == 'V-REGION' and not a[4] == 'D-REGION'\
                   and not a[4] == 'J-REGION':
                    b = a[0].split()
                    anno[b[0][1:]] = a[1]

    Cgene = {}
    CE = {}
    with open (BlastOut,'rU') as infile:
        for line in infile:
            a = line.split()
            a[7] = float(a[7])
            if not a[1] in anno: continue
            anno[a[1]] = anno[a[1]].split('*')[0]
            if not a[0] in CE or a[7] < CE[a[0]]: 
                CE[a[0]] = a[7]
                Cgene[a[0]] = set([anno[a[1]]])
            elif a[7] == CE[a[0]]: 
                Cgene[a[0]].add(anno[a[1]])

    return Cgene
    '''    
    IG_contig = set([])
    TR_contig = set([])
    with open(BlastOut,'rU') as infile:
        for line in infile:
            a = line.split()
            if not a[0] in IG_contig and not a[0] in TR_contig and a[1] in anno:
                if re.search(r'^IG',anno[a[1]]):
                    IG_contig.add(a[0])
                elif re.search(r'^TR',anno[a[1]]):
                    TR_contig.add(a[0])

    IG = subprocess.check_output('awk -F "\t" \'{OFS="\t"; if (($1~/^[VDJ]/)&&($3~/^IG/)) print}\' %s | cut -f2 | sort -u | wc -l' % igBlastOut, shell=True, stderr=subprocess.STDOUT)[:-1]
    TCR = subprocess.check_output('awk -F "\t" \'{OFS="\t"; if (($1~/^[VDJ]/)&&($3~/^TR/)) print}\' %s | cut -f2 | sort -u | wc -l' % igBlastOut, shell=True, stderr=subprocess.STDOUT)[:-1]
    with open(BlastOut.replace('results','summary'),'w') as out:
        out.write('\t'.join(['igBLAST IG:',IG])+'\n')
        out.write('\t'.join(['igBLAST TCR:',TCR])+'\n')
        out.write('\t'.join(['BLAST IG:',str(len(IG_contig))])+'\n')
        out.write('\t'.join(['BLAST TCR:',str(len(TR_contig))])+'\n')
    '''


def parseBlastOutput(igBlastOut):
    Vtop = {}
    Dtop = {}
    Jtop = {}
    VE = {}
    DE = {}
    JE = {}
    with open(igBlastOut,'rU') as infile:
        for line in infile:
            if not re.search(r'^[VDJ]',line): continue

            a = line.split()
            a[1] = a[1].replace('reversed|','')
            b = a[2].split('*')
            a[2] = b[0]
            if a[0] == 'V':
                if not a[1] in VE or float(a[12]) < VE[a[1]]:
                    Vtop[a[1]] = set([a[2]])
                    VE[a[1]] = float(a[12])
                elif float(a[12]) == VE[a[1]]:
                    Vtop[a[1]].add(a[2])
            if a[0] == 'D':
                if not a[1] in DE or float(a[12]) < DE[a[1]]:
                    Dtop[a[1]] = set([a[2]])
                    DE[a[1]] = float(a[12])
                elif float(a[12]) == DE[a[1]]:
                    Dtop[a[1]].add(a[2])
            if a[0] == 'J':
                if not a[1] in JE or float(a[12]) < JE[a[1]]:
                    Jtop[a[1]] = set([a[2]])
                    JE[a[1]] = float(a[12])
                elif float(a[12]) == JE[a[1]]:
                    Jtop[a[1]].add(a[2])
                
    return [Vtop,Dtop,Jtop]


def parseCombination(sampleName,IGgenes,igBlastOut):
    comb = {}
    anno = {}
    CDR3 = {}
    uniq = {}
    last_contig = ''
    follow = 0
    with open(igBlastOut,'rU') as infile:
        for line in infile:
            if follow > 0:
                a = line[:-1].split('\t')
                if a[0]=='N/A' or a[1]=='N/A' or (follow==2 and a[2]=='N/A') or last_contig == '': 
                    follow = 0
                    continue

                comb[last_contig] = '::'.join(a[:follow+1])
                anno[last_contig] = ','.join(a[follow+1:])
                if follow == 1 and len(IGgenes[0][last_contig]) == 1 and len(IGgenes[2][last_contig]) == 1:
                    uniq[last_contig] = comb[last_contig]
                elif follow == 2 and len(IGgenes[0][last_contig]) == 1 and len(IGgenes[1][last_contig]) == 1\
                     and len(IGgenes[2][last_contig]) == 1:
                    uniq[last_contig] = comb[last_contig]

                follow = 0
            elif follow == -5:
                a = line[:-1].split('\t')
                CDR3[last_contig] = a[2]
                follow = 0
            elif re.search(r'^# Query: ',line):
                a = line.split()
                last_contig = a[2]
            elif re.search(r'^# V-\(D\)-J rearrangement summary',line):
                follow = 1
                if re.search(r'Top D gene match',line): 
                    follow = 2
            elif re.search(r'^# Sub-region sequence details',line):
                follow = -5

    completeType = {}
    for contig,clonetype in comb.iteritems():
        if clonetype in completeType: 
            completeType[clonetype] += 1
        else: 
            completeType[clonetype] = 1

    cdr3expr = {}
    cdr3chain = {}
    for contig,cdr3 in CDR3.iteritems():
        if contig in anno and anno[contig].split(',')[-2] == 'Yes':
            if cdr3 in cdr3expr: cdr3expr[cdr3] += 1
            else: 
                cdr3expr[cdr3] = 1
                cdr3chain[cdr3] = anno[contig].split(',')[0]

    with open(sampleName+'.IGgenes.clonetypes','w') as out:
        for idiotype in sorted(completeType,key=completeType.get,reverse=True):
            expr = int(completeType[idiotype])
            out.write('\t'.join( [idiotype,str(expr)] )+'\n')

    with open(sampleName+'.IGgenes.CDR3.abundance','w') as out:
        for cdr3 in sorted(cdr3expr,key=cdr3expr.get,reverse=True):
            if cdr3expr[cdr3] > 10:
                out.write('\t'.join( [cdr3,str(int(cdr3expr[cdr3])),cdr3chain[cdr3]] )+'\n')
                

def quantify(sampleName,IGgenes):
    '''
    counts = {}
    with open(sampleName+'.contigs.genes.results','rU') as infile:
        next(infile)
        for line in infile:
            a = line.split()
            counts[a[0]] = float(a[4])
    '''
    IGexpr = {}
    uniq = {}
    for top in IGgenes:
        for contig,matches in top.iteritems():
            for match in matches:
                #value = counts[contig]/float(len(matches))
                value = 1/float(len(matches))
                if match in IGexpr:
                    IGexpr[match] += value
                else:
                    IGexpr[match] = value

                if len(matches) == 1:
                    if match in uniq:
                        uniq[match] += 1
                    else:
                        uniq[match] = 1


    with open(sampleName+'.IGgenes.readCounts','w') as out:
        for key in sorted(IGexpr,key=IGexpr.get,reverse=True):
            value = int(IGexpr[key])
            if value > 0:
                if not key in uniq: uniq[key] = 0
                out.write('\t'.join([key,str(value),str(uniq[key])])+'\n')
    ''' 
    #return IGexpr,uniq
    '''
    consolid = {}
    for key,value in IGexpr.iteritems():
        if re.search(r'^IG[HKL][VDJ]',key):
            gene = key.split('-')[0]
            if gene in consolid:
                consolid[gene] += value
            else:
                consolid[gene] = value
    print consolid
    with open(sampleName+'.IGgenes.consolid.counts','w') as out:
        for key,value in consolid.iteritems():
            out.write('%s\t%d\n' % (key,value))


def plotLenDist(sampleName,IGgenes,faFile):
    contigs = set([])
    for top in IGgenes:
        for contig in top.keys():
            contigs.add(contig)

    lens = []
    seq = ''
    contig = ''
    length = 0
    out = open(faFile.replace('fasta','IGgenes.fa'),'w')
    with open(faFile,'rU') as infile:
        for line in infile:
            if re.search(r'^>',line):
                a = line.split()
                if not contig == '' and not contig in contigs:
                    out.write('>'+contig+'\n')
                    out.write(seq+'\n')
                    lens.append(math.log10(length))
                seq = ''
                contig = a[0].strip('>')
                length = int(a[1].replace('len=',''))
            else:
                seq += line[:-1]
    out.close()

    fig, ax = plt.subplots()
    plt.hist(lens,bins=50)
    plt.ylabel('Count')
    plt.xlabel('Log10 Contig length')
    plt.title('Length Distribution of contigs covering IG/TCR genes\n(%s contigs)' % str(len(lens)))
    plt.savefig(sampleName + '.IGgenes.contigs.lendist.pdf')
    plt.close()


def main():
    files = glob.glob('[78]*-IP*.candidate.reads.fa')
    for fn in files:
        sample = fn.replace('.candidate.reads.fa','')
        preprocess(sample)
        IGgenes = parseBlastOutput(sample+'.igBLAST.results')
        IGgenes.append(parseBlastOutput0(sample+'.BLAST.results'))
        #counts = quantify(sys.argv[1],IGgenes)
        quantify(sample,IGgenes)
        #parseCombination(sample,IGgenes,sample+'.igBLAST.results')
        #plotLenDist(sys.argv[1],IGgenes,'Trinity.fasta')


if __name__=='__main__':
    main()
