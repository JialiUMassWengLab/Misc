#! /usr/bin/bash

for i in [18]*.filtered.bam 
do 
    name=${i/.filtered.bam/} 
    bedtools intersect -a $i -b /mnt/shares/Users/jzhuang/2017Jul_MM/IG_TCR_genes.bed -wa > $name.IG_TCR.bam
    bam2fastq -o $name.IG_TCR#.fastq $name.IG_TCR.bam
    bam2fastq -o $name.unAligned#.fastq $name"unAligned.out.bam"
    cat $name.IG_TCR_1.fastq $name.unAligned_1.fastq | sed '/^@/!d;s//>/;N' >> $name.candidate.reads.fa
    cat $name.IG_TCR_2.fastq $name.unAligned_2.fastq | sed '/^@/!d;s//>/;N' >> $name.candidate.reads.fa    
done
