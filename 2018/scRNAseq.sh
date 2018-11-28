cut -f1,30 -d"," runInfo 

~/bin/Salmon-0.8.2_linux_x86_64/bin/salmon quant -i /mnt/shares2/annotations/hg38/RSEMref/hg38.transcriptome.with_ERCC.transcripts.fa.salmon.idx/ -g /mnt/shares2/annotations/hg38/STAR/GRCh38_Gencode24_ERCC/gencode.v24.primary_assembly.annotation.with_ERCC.gtf -o Monocyte_rep4 -l IU -1 rep4/SRR1022982_1.fastq.gz rep4/SRR1022983_1.fastq.gz rep4/SRR1022984_1.fastq.gz -2 rep4/SRR1022982_2.fastq.gz rep4/SRR1022983_2.fastq.gz rep4/SRR1022984_2.fastq.gz

curl -OO sftp://bidas.kaist.ac.kr:4000/New_Collection/Choroid_plexus/N10/N10_{1,2}.fastq -u guest_180928_1 -k

wget ftp://ftp.sra.ebi.ac.uk/vol1/ERA159/ERA1598441/bam/P10A_Cerebella.bam



./Drop-seq_alignment2.sh -g mm10/STAR/ -r mm10/mm10.fasta -d /home/jzhuang@ms.local/bin/Drop-seq_tools-2.0.0 -o /mnt/nfs/sequencing_data/ClinChem_data/ -t /home/jzhuang@ms.local/ /mnt/nfs/sequencing_data/ClinChem_data/SRR6313172.unaligned.bam &
./DigitalExpression I=/mnt/nfs/sequencing_data/ClinChem_data/final.bam O=/mnt/nfs/sequencing_data/ClinChem_data/DGE.mat OUTPUT_LONG_FORMAT=/mnt/nfs/sequencing_data/ClinChem_data/DGE.long.mat NUM_CORE_BARCODES=5000 RARE_UMI_FILTER_THRESHOLD=0.01 TMP_DIR=/home/jzhuang@ms.local &
