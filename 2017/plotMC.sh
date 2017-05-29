#! /usr/bin/bash

name=$1
gtfFile=$2
samtools view -H $name"Aligned.toTranscriptome.sorted.dedup.bam" | awk '{OFS="\t"; if ($0~/^@SQ/) {sub("SN:","",$2);sub("LN:","",$3);print $2,0,$3,"transcript"++i}}' > tx.bed
bedtools coverage -bed -a tx.bed -b $name"Aligned.toTranscriptome.sorted.dedup.bam" > tmp.bed
rm tx.bed

Rscript /mnt/pipeline/pipeline_scripts/plotMC.R $name $gtfFile
rm tmp.bed
