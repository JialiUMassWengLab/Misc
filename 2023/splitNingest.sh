#!/bin/bash -l
#SBATCH -J splitNingestVCF
#SBATCH --mail-type=end,fail
#SBATCH --mail-user=jiali.zhuang@takeda.com
#SBATCH -N 1
#SBATCH -p queue1
#SBATCH --mem=150G
#SBATCH --time=0-08:00:00
#SBATCH --qos=normal

echo "== Starting run at $(date)"
echo "== Job ID: ${SLURM_JOBID}"
echo "== Node list: ${SLURM_NODELIST}"
echo "== Submit dir. : ${SLURM_SUBMIT_DIR}"
echo ""

module load bcftools/1.9
module load htslib/1.9
dataset="FINAL_3041samples_2020_12_22.SCRAMBLED.100k_lines"

echo "=>Spliting VCF file..."
#bcftools +split $dataset.vcf.gz -Oz -o $dataset
for i in $dataset\/*.vcf.gz; do tabix -p vcf $i; done
echo "=>Completed"
echo

echo "=>Ingesting VCF into TileDB..."
. /mnt/efs/app/anaconda3/v2022.10/etc/profile.d/conda.sh
conda activate tiledb
python ingestVCF.py $dataset 500
echo "=>Completed"
