#! /usr/bin/bash

#for i in `seq 1 20`
for i in NIA_JG_1898_samples_GRM_WGS_b37_JointAnalysis01_2017-12-08_*_variants.vcf.gz
do
    anno=${i/variants.vcf.gz/variants.annotated.vcf.gz}
    name=${i/NIA_JG_1898_samples_GRM_WGS_b37_JointAnalysis01_2017-12-08_/}
    name=${name%.recalibrated_variants.vcf.gz}
    echo "chr$name"
    #plink --vcf NIA_JG_1898_samples_GRM_WGS_b37_JointAnalysis01_2017-12-08_$i.recalibrated_variants.vcf.gz --make-bed --double-id --keep-allele-order --keep-fam ROSMAP_wgs_patient.txt --vcf-min-gq 20 --maf 0.05 --out chr$i\_ROSMAP_maf0.05
    #plink --vcf $i --vcf-min-gq 20 --biallelic-only --make-bed --double-id --keep-allele-order --keep-fam ../Derived/ROSMAP_wgs_patient.txt --update-sex ../Output/ROSMAP_covariate4PLINK.txt 2 --maf 0.05 --out ../Derived/chr$name\_ROSMAP_maf0.05
    bcftools annotate -a $anno -c ID --collapse some -I +'%CHROM\_%POS\_%REF\_%FIRST_ALT' --threads 12 -Ov $i | plink --vcf /dev/stdin --vcf-min-gq 20 --biallelic-only --make-bed --double-id --keep-allele-order --keep-fam ../Derived/$1\_wgs_patient.txt --update-sex ../Output/$1\_covariate4PLINK.txt 2 --maf 0.001 --geno --out ../Derived/chr$name\_$1\_maf0.001
done

