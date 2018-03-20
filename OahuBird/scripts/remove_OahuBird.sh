#!/bin/bash

#SBATCH -D /mnt/lustre/macmaneslab/devon/guano/NAU/OahuBird/illumina/
#SBATCH -p shared,macmanes
#SBATCH --job-name="rfiltAMPobird"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --output=OahuBird_ampRemove_all.log

module purge
module load anaconda/colsa
source activate amptk_env
PATH=/mnt/lustre/macmaneslab/devon/bin:$PATH

srun echo "      /\^._.^/\     starting process: `date`"

amptk remove \
--input OahuBird.demux.fq \
--threshold 1000 \
--out trimd_OahuBird.demux.fq

echo "      /\^._.^/\     ending process: `date`"
