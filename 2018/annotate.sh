#! /usr/bin/bash

for i in DESeq2*.tsv
do
    name=${i/.tsv/.sig}
    #join -a 1 -j 1 <(cut -f1,7,3 $i | sort -k1,1 | awk '{sub(/\.[0-9\.]+$/,"",$1);print $1,$2,$3}') <(zcat /usr/local/bioinfo/MS/tissueExpnDB/GTEx_Analysis_v6_RNA-seq_RNA-SeQCv1.1.8_gene_median_rpkm.gct.gz | cut -f1,2 | sort -k1,1 | awk '{sub(/\.[0-9\.]+$/,"",$1);print $1,$2}') | sort -k3,3g | awk '{OFS=","; if ($4=="") $4="NA";if (($3!="")&&($3<0.1)) print $1,$2,$3,$4}' > $name
    join -a 1 -j 1 <(cut -f1,7,3 $i | sort -k1,1 | awk '{sub(/\.[0-9\.]+$/,"",$1);print $1,$2,$3}') <(zcat /usr/local/bioinfo/MS/tissueExpnDB/GTEx_Analysis_v6_RNA-seq_RNA-SeQCv1.1.8_gene_median_rpkm.gct.gz | cut -f1,2 | sort -k1,1 | awk '{sub(/\.[0-9\.]+$/,"",$1);print $1,$2}') | sort -k3,3g | awk '{OFS=","; if ($4=="") $4="NA";if ($3!="") print $1,$2,$3,$4}' > $name

    name1=${i/DESeq2./}
    name1=${name1/.tsv/.upGenes.txt}
    name2=${i/DESeq2./}
    name2=${name2/.tsv/.downGenes.txt}
    awk -F "," '{if (($3<0.05)&&($2<0)) print $1}' $name > $name1
    awk -F "," '{if (($3<0.05)&&($2>0)) print $1}' $name > $name2
done
