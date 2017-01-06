#! /usr/bin/python

import re
import csv
import sys

def main():
    loci = []
    hap1 = {}
    hap2 = {}
    with open(sys.argv[1],'rU') as infile:
        lines = csv.reader(infile,delimiter=',')
        for line in lines:
            match = re.search(r'(\w+):(\d+)-(\d+)', line[0])
            if match:
                loci.append(':'.join([match.group(1),match.group(2)]))


    #with open(sys.stdin,'rU') as infile:
    for line in sys.stdin:
        if re.search(r'^#', line): continue

        fields = line.split('\t')
        coor = ':'.join([fields[0], fields[1]])
        if coor in loci:
            a = fields[9].split(':')
            match = re.search(r'(\d)\|(\d)',a[0])
            if match:
                if match.group(1)=='0':
                    hap1[coor] = fields[3]
                else:
                    hap1[coor] = fields[4]
                    
                if match.group(2)=='0':
                    hap2[coor] = fields[3]
                else:
                    hap2[coor] = fields[4]
                    
    #print hap1
    #print hap2

    with open(sys.argv[2],'rU') as infile:
        for line in infile:
            if not re.search(r'^chr',line): continue

            predictions = line.split('\t')
            loci_num = 0
            correct_num1 = 0
            correct_num2 = 0
            for pred in predictions[:-1]:
                a = pred.split(':')
                a[1] = int(a[1])+1
                coor = a[0]+':'+str(a[1])
                if coor in loci:
                    if a[2] == hap1[coor]:
                        correct_num1 += 1
                    if a[2] == hap2[coor]:
                        correct_num2 += 1
                    loci_num += 1

            print line[:-1]
            print '{:.4f}'.format(correct_num1/float(loci_num)),'{:.4f}'.format(correct_num2/float(loci_num))


if __name__ == "__main__":
    main()
