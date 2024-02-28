#! /usr/bin/bash

declare -a arr=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "X" "Y")
cd /mnt/data/

for i in "${arr[@]}"
do
    echo "chr$i"
    wget https://gnomad-public-us-east-1.s3.amazonaws.com/release/4.0/vcf/exomes/gnomad.exomes.v4.0.sites."chr$i".vcf.bgz
    wget https://gnomad-public-us-east-1.s3.amazonaws.com/release/4.0/vcf/exomes/gnomad.exomes.v4.0.sites."chr$i".vcf.bgz.tbi
    bcftools view -i 'INFO/AF_grpmax[0]>0.01' gnomad.exomes.v4.0.sites."chr$i".vcf.bgz | bcftools +split-vep -a vep -d -f '%CHROM\_%POS\_%REF\_%ALT\t%Gene\t%Consequence\t%IMPACT\t%LoF_flags\n' -s worst -i 'IMPACT="HIGH" || IMPACT="MODERATE" || LoF_flags="HC"' - > "chr$i".output.txt
    rm gnomad.exomes.v4.0.sites."chr$i".vcf.bgz*
done
