#!/bin/bash

#SBATCH -D /mnt/lustre/macmaneslab/devon/guano/NAU/OahuBird/illumina/p10-2/readfiltd
#SBATCH -p shared,macmanes
#SBATCH --job-name="rfiltAMP102"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --output=OahuBird_ampRemove_p10-2.log

module purge
module load anaconda/colsa
source activate amptk_env
PATH=/mnt/lustre/macmaneslab/devon/bin:$PATH

srun echo "      /\^._.^/\     starting process: `date`"

amptk remove \
-i /mnt/lustre/macmaneslab/devon/guano/NAU/OahuBird/illumina/p10-2/trim_p102.demux.fq \
-t 50 \
-o readfiltd_p102.demux.fq

echo "      /\^._.^/\     ending process: `date`"
