#!/bin/bash -l
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --qos=normal
#SBATCH --time=0-5:00:00
#SBATCH --job-name=vep
#SBATCH --mail-type=end,fail
#SBATCH --mail-user=jiali.zhuang@takeda.com
#SBATCH -p queue1
 
#module load perl2/5.34.0
#module load vep/99.0
module load htslib
 
assembly=GRCh38
inputVCF=chr1..22_hard_filters_stringent_CR.95_44028samples.tidied.vep105.sitesonly.vcf.gz
dir=/mnt/efs/home/iet5740/Projects/GenesAndHealth

cd $dir
tabix -h $inputVCF chr21 > input.vcf
mkdir output_dir
chmod -R a+rwx output_dir

docker run --user root \
  -w $(readlink -f $(pwd)) \
  --mount type=bind,src=${dir},dst=${dir} \
  --mount type=bind,src=/mnt/efs/home/iet5740/.vep,dst=/root/.vep \
  --mount type=bind,src=/mnt/efs/home/iet5740/.vep/Plugins/loftee_data/${assembly},dst=/root/.loftee,readonly \
  yosuketanigawa/docker-ensembl-vep-loftee:release_105.0 \
  vep \
  --offline --cache \
  --allele_number --everything \
  --assembly ${assembly} --dir /root/.vep \
  -i input.vcf -o output_dir/GenesAndHealth.tidied.sitesonly.vep."chr21" \
  --plugin LoF,loftee_path:/opt/vep/src/loftee-master,human_ancestor_fa:/root/.loftee/human_ancestor.fa.gz,conservation_file:/root/.loftee/loftee.sql,gerp_bigwig:/root/.loftee/gerp_conservation_scores.homo_sapiens.GRCh38.bw

#docker run --user root -it -w $(readlink -f $(pwd)) --mount type=bind,src=/mnt/efs/home/iet5740/Projects/GenesAndHealth,dst=/mnt/efs/home/iet5740/Projects/GenesAndHealth --mount type=bind,src=/mnt/efs/home/iet5740/.vep,dst=/root/.vep --mount type=bind,src=/mnt/efs/home/iet5740/.vep/Plugins/loftee_data/GRCh38,dst=/root/.loftee,readonly yosuketanigawa/docker-ensembl-vep-loftee:release_101.0 vep --offline --cache --allele_number --everything --assembly GRCh38 -i input.vcf -o output_dir/GenesAndHealth.tidied.sitesonly.vep.chr21 --dir /root/.vep --plugin LoF,loftee_path:/opt/vep/src/loftee-master,human_ancestor_fa:/root/.loftee/human_ancestor.fa.gz,gerp_bigwig:/root/.loftee/gerp_conservation_scores.homo_sapiens.GRCh38.bw
