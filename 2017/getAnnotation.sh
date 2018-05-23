#! /usr/bin/bash

awk '$3 == "transcript"' STAR/GRCm38_Gencode14_ERCC/gencode.vM14.primary_assembly.annotation.with_ERCC.gtf | cut -f1,4,5,7,9 | gawk '{OFS=","; match($0,/gene_id "([^"]+)"/,m); match($0,/gene_name "([^"]+)"/,n); match($0,/transcript_id "([^"]+)"/,l); match($0,/gene_type "([^"]+)"/,o); print l[1],m[1],n[1],o[1]}' > GeneInfo.mouse.csv
