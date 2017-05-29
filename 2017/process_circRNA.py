#! /usr/bin/python

import csv
import sys
import os
import glob
import re
import math
import subprocess
import operator
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

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

def plotCountDist(title,Info):
    logCountList = []
    for key,value in Info.iteritems():
        logCountList.append(math.log10(value['counts']))
        
    fig, ax = plt.subplots()
    plt.hist(logCountList,bins=20)
    plt.savefig(title + '.pdf')
    plt.close()

def plotLenDist(title,Info,threshold):
    lenList = []
    for key,value in Info.iteritems():
        a = key.split('|')
        if value['juncClass']=='circular' and value['counts'] >= threshold:
            lenList.append(math.log10(abs(int(a[1])-int(a[4]))))
        
    fig, ax = plt.subplots()
    plt.hist(lenList,bins=20)
    plt.xlabel('log10 circRNA length')
    plt.ylabel('Number of circRNAs')
    plt.title(title)
    plt.savefig(title + '.pdf')
    plt.close()

def plotExonLenDist(title,info,threshold):
    bedOverlap(info,threshold,'/mnt/shares2/annotations/hg38/refFlat_hg38.bed','-s -f 0.95 -wo')
    exon_nums = {}
    transcript_lens = {}
    with open('tmp','rU') as infile:
        for line in infile:
            a = line.split()
            starts = a[17].split(',')[:-1]
            ends = a[16].split(',')[:-1]
            lens = a[16].split(',')[:-1]
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
            if start_match > -1 and end_match > -1:
                exon_nums[a[3]] = end_match-start_match+1
                transcript_lens[a[3]] = 0
                for i in range(start_match,end_match+1):
                    transcript_lens[a[3]] += int(lens[i])

    subprocess.check_output('rm tmp.bed tmp',shell=True,stderr=subprocess.STDOUT)

    exonNumList = []
    for value in exon_nums.values():
        if value > 20:
            value = 20
        exonNumList.append(value)
    fig, ax = plt.subplots()
    plt.hist(exonNumList,bins=20)
    plt.xlabel('Number of exons')
    plt.ylabel('Number of circRNAs')
    plt.title(title)
    plt.savefig(title + 'Num.pdf')

    exonLenList = []
    for value in transcript_lens.values():
        exonLenList.append(math.log10(value))
    fig, ax = plt.subplots()
    plt.hist(exonLenList,bins=20)
    plt.xlabel('Log 10 exon lengths')
    plt.ylabel('Number of circRNAs')
    plt.title(title)
    plt.savefig(title + 'Len.pdf')    

def plotRecurrence(occ):
    count = occ.values()
    fig, ax = plt.subplots()
    plt.hist(count,bins=20)
    plt.xlabel('Number of samples')
    plt.ylabel('Number of circRNAs')
    plt.title('circRNA recurrence across samples')
    plt.savefig('circRNA_recurrence.pdf')

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

def circBaseValidate(info,threshold):
    bedOverlap(info,threshold,'hsa_hg38_circRNA.bed','-s -f 0.9 -wo')
    validated_set = {}
    with open('tmp','rU') as infile:
        for line in infile:
            a = line.split()
            #if int(a[-1])/(float(a[8])-float(a[7])) >= 0.9:
            if abs(int(a[1])-int(a[7])) <= 5 and abs(int(a[2])-int(a[8])) <= 5:
                validated_set[a[3]]=info[a[3]]['counts']

    subprocess.check_output('bedtools intersect -a tmp.bed -b human_brain_circRNA.hg38.bed -s -f 0.9 -wo > tmp',shell=True,stderr=subprocess.STDOUT)
    validated_set = {}
    with open('tmp','rU') as infile:
        for line in infile:
            a = line.split()
            #if int(a[-1])/(float(a[8])-float(a[7])) >= 0.9:
            if abs(int(a[1])-int(a[7])) <= 5 and abs(int(a[2])-int(a[8])) <= 5:
                validated_set[a[3]]=info[a[3]]['counts']

    subprocess.check_output('rm tmp.bed tmp',shell=True,stderr=subprocess.STDOUT)
    return len(validated_set)

def exonValidate(info,threshold):
    bedOverlap(info,threshold,'/mnt/shares2/annotations/hg38/refFlat_hg38.bed','-s -f 0.95 -wo')
    validated_set = {}
    with open('tmp','rU') as infile:
        for line in infile:
            a = line.split()
            starts = a[17].split(',')[:-1]
            ends = a[16].split(',')[:-1]
            matches = 0
            for i in range(len(starts)):
                starts[i] = int(starts[i])+int(a[7])
                if abs(starts[i]-int(a[1])) <= 5:
                    matches += 1
            for i in range(len(ends)):
                ends[i] = starts[i]+int(ends[i])
                if abs(ends[i]-int(a[2])) <= 5:
                    matches += 1
            if matches == 2:
                validated_set[a[3]]=info[a[3]]['counts']

    subprocess.check_output('rm tmp.bed tmp',shell=True,stderr=subprocess.STDOUT)
    return len(validated_set)

def getCircularNumbers(dirName,threshold):
    nums = []
    files = glob.glob(os.path.join(dirName,'*Chimeric.out.summary'))
    for filename in files:
        info=getInfo(filename)
        num1=0
        num2=0
        for key,value in info.iteritems():
            if value['counts'] >= threshold:
                if value['juncClass']=='circular':
                    num1 += 1
                else:
                    num2 += 1

        samplename = re.sub(r'Chimeric.out.summary','',os.path.basename(filename))
        nums.append({'sample': samplename,
                     'circular': num1,
                     'fusion': num2,
                     'circBase_validated': circBaseValidate(info,threshold),
                     'exon_validated': exonValidate(info,threshold),
                     'circBase_validated_pct': '{0:.4f}'.format(circBaseValidate(info,threshold)/float(num1)),
                     'exon_validated_pct': '{0:.4f}'.format(exonValidate(info,threshold)/float(num1))
                 })

        plotLenDist(os.path.join(dirName,samplename+'.lendist'),info,threshold)
        plotExonLenDist(os.path.join(dirName,samplename+'.exon'),info,threshold)
    return nums

def getOccurrences(threshold):
    occ = {}
    files = glob.glob('./*.summary')
    for filename in files:
        print filename
        info = getInfo(filename)
        counted = set([])
        for key,value in info.iteritems():
            a = key.split('|')
            if value['juncClass']=='circular' and value['counts'] >= threshold:
                found = False
                for key1,value1 in occ.iteritems():
                    b = key1.split('|')
                    if a[0]==b[0] and a[2]==b[2] and abs(int(a[1])-int(b[1])) <= 5 and abs(int(a[4])-int(b[4])) <= 5 and not key1 in counted:
                        occ[key1] = value1 + 1
                        found = True
                        counted.add(key1)
                        break

                if not found:
                    occ.update({key:1})

    return occ

def circRatio(info,threshold,destPrefix,estimateReadCounts):
    linear_junc = {}
    with open(destPrefix+'SJ.out.tab','rU') as infile:
        for line in infile:
            a = line.split()
            strand = '+'
            if a[3]=='2':
                strand = '-'
            key = '|'.join([a[0],a[1],a[2],strand])
            linear_junc[key] = int(a[6])

    ratio = {}
    linear_max_dict = {}
    bedOverlap(info,threshold,'/mnt/shares2/annotations/hg38/refFlat_hg38.bed','-s -f 0.95 -wo')
    with open('tmp','rU') as infile:
        for line in infile:
            a = line.split()
            starts = map(int,a[17].split(',')[:-1])
            ends = map(int,a[16].split(',')[:-1])
            start_match = -1
            end_match = -1
            for i in range(len(starts)):
                starts[i] = starts[i]+int(a[7])
                if abs(starts[i]-int(a[1])) <= 5:
                    start_match = i
            for i in range(len(ends)):
                ends[i] = starts[i]+ends[i]
                if abs(ends[i]-int(a[2])) <= 5: 
                    end_match = i

            if start_match == -1 or end_match == -1:
                continue

            linear_max = 0
            if start_match > 0:
                junc1 = '|'.join([a[0],str(ends[start_match-1]+1),str(starts[start_match]),a[5]])
                if junc1 in linear_junc and linear_junc[junc1] > linear_max:
                    linear_max = linear_junc[junc1]
            if end_match < len(ends)-1:
                junc2 = '|'.join([a[0],str(ends[end_match]+1),str(starts[end_match+1]),a[5]])
                if junc2 in linear_junc and linear_junc[junc2] > linear_max:
                    linear_max = linear_junc[junc2]

            ratio[a[3]] = (info[a[3]]['counts']+1)/float(linear_max+1)
            if a[3] in linear_max_dict:
                linear_max_dict[a[3]] += linear_max
            else:
                linear_max_dict[a[3]] = linear_max

    subprocess.check_output('rm tmp.bed tmp',shell=True,stderr=subprocess.STDOUT)

    #plot circRNA ratio distribution
    if estimateReadCounts==False:
        ratioList = []
        for value in ratio.values():
            ratioList.append(math.log(value,2))
        fig, ax = plt.subplots()
        plt.hist(ratioList,bins=20)
        plt.xlabel('Log 2 circular/linear ratio')
        plt.ylabel('Number of circRNAs')
        plt.title(destPrefix + '.circRNA.ratio')
        plt.savefig(destPrefix + '.circRNA.ratio.pdf')

        return ratio

    #estimate number of reads coming from circRNA
    out = open('overlap.bed','w')
    for key,value in linear_max_dict.iteritems():
        a = key.split('|')
        lower = a[1]
        upper = a[4]
        if int(a[1]) > int(a[4]):
            lower = a[4]
            upper = a[1]

        out.write('\t'.join([a[0],lower,upper,key,str(info[key]['counts']/float(info[key]['counts']+value)),a[5]])+'\n')
    out.close()

    subprocess.check_output('bedtools sort -i overlap.bed -faidx /mnt/shares2/annotations/hg38/STAR/GRCh38_Gencode24/chrNameLength.txt > overlap1.bed',shell=True,stderr=subprocess.STDOUT)
    subprocess.check_output('bedtools coverage -a overlap1.bed -b %s -sorted -nobuf -s -g /mnt/shares2/annotations/hg38/STAR/GRCh38_Gencode24/chrNameLength.txt > reads_counts.tab' % (destPrefix+'Aligned.sortedByCoord.out.bam'),shell=True,stderr=subprocess.STDOUT)
    circ_reads = 0.0
    with open('reads_counts.tab','rU') as infile:
        for line in infile:
            a = line.split()
            circ_reads += float(a[4])*float(a[6])
    print circ_reads
    '''
    subprocess.check_output('bedtools coverage -a refFlat_hg38_merged.bed -b %s -sorted -nobuf -s -g /mnt/shares2/annotations/hg38/STAR/GRCh38_Gencode24/chrNameLength.txt > reads_counts.tab' % sortedBamFile,shell=True,stderr=subprocess.STDOUT)
    total_reads = 0.0
    with open('reads_counts.tab','rU') as infile:
        for line in infile:
            a = line.split()
            total_reads += float(a[4])
    print total_reads
    subprocess.check_output('rm reads_counts.tab overlap*',shell=True,stderr=subprocess.STDOUT)
    '''
    total_reads = 0.0
    with open(destPrefix+'Aligned.sortedByCoord.out.bam.readDist.out') as infile:
        for line in infile:
            a = line.split()
            if re.search(r'_Exons',a[0]):
                total_reads += float(a[2])

    return circ_reads/total_reads


def main():
    '''
    occ = getOccurrences(10)
    sorted_occ = sorted(occ, key=occ.get, reverse=True)
    #sorted_occ = sorted(occ.items(), key=operator.itemgetter(1), reverse=True)

    with open('circRNA_sample_counts.tsv','w') as out:
        for key in sorted_occ:
            out.write(key+'\t'+str(occ[key])+'\n')

    plotRecurrence(occ)
    '''
    files = glob.glob(os.path.join(sys.argv[1],'*Chimeric.out.junction'))
    for filename0 in files:
        filename = re.sub(r'junction','summary',filename0)
        subprocess.check_output('cat %s | awk \'$1!="chrM" && $4!="chrM" && $1!~/ERCC-/ && $4!~/ERCC-/ && $7>0 && $8+$9<=5 {print $1,$2,$3,$4,$5,$6,$7,$8,$9}\' | sort | uniq -c | sort -k1,1rn > %s' % (filename0,filename),shell=True,stderr=subprocess.STDOUT)
        #info=getInfo(filename)
        #destPrefix = re.sub(r'Chimeric.out.junction','',filename0)
        #print filename,circRatio(info,1,destPrefix,True)

    runName = os.path.basename(sys.argv[1])
    if runName == '': 
        runName = os.path.basename(os.path.dirname(sys.argv[1]))

    summary = getCircularNumbers(sys.argv[1],10)
    with open(os.path.join(sys.argv[1],runName+'circRNA_summary.tsv'),'w') as out:
        writer = csv.DictWriter(out,fieldnames=summary[0].keys())
        writer.writeheader()
        for item in summary:
            writer.writerow(item)



if __name__=='__main__':
    main()
            
