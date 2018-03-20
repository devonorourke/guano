#!/bin/bash

#SBATCH -D /mnt/lustre/macmaneslab/devon/guano/NAU/OahuBird/clust/all
#SBATCH -p shared,macmanes
#SBATCH --job-name="ampClustAll"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --output=OahuBird_clust_all.log

module purge
module load anaconda/colsa
source activate amptk_env
PATH=/mnt/lustre/macmaneslab/devon/bin:$PATH

srun echo "      /\^._.^/\     starting process: `date`"

amptk dada2 \
--fastq /mnt/lustre/macmaneslab/devon/guano/NAU/OahuBird/illumina/trimd_OahuBird.demux.fq \
--out OahuBird \
--length 180 \
--platform illumina \
--uchime_ref COI

echo "      /\^._.^/\     ending process: `date`"
