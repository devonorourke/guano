#!/bin/bash

#SBATCH -D /mnt/lustre/macmaneslab/devon/guano/NAU/OahuBird/clust/all
#SBATCH -p shared,macmanes
#SBATCH --job-name="ampDropAll"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --output=OahuBird_otudropd_all.log

module purge
module load anaconda/colsa
source activate amptk_env
PATH=/mnt/lustre/macmaneslab/devon/bin:$PATH

srun echo "      /\^._.^/\     starting process: `date`"

amptk drop \
--input /mnt/lustre/macmaneslab/devon/guano/NAU/OahuBird/clust/all/OahuBird.cluster.otus.fa \
--reads /mnt/lustre/macmaneslab/devon/guano/NAU/OahuBird/illumina/trimd_OahuBird.demux.fq \
--list OTU1730 OTU3130 OTU224 OTU1357 OTU949 OTU554 OTU264 OTU279 OTU2285 \
--out OahuBird_otudropd

echo "      /\^._.^/\     ending process: `date`"
