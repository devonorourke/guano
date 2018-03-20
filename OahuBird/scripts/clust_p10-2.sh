#!/bin/bash

#SBATCH -D /mnt/lustre/macmaneslab/devon/guano/NAU/OahuBird/clust/p10-2
#SBATCH -p shared,macmanes
#SBATCH --job-name="ampClust2"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --output=OahuBird_clust_p10-2.log

module purge
module load anaconda/colsa
source activate amptk_env
PATH=/mnt/lustre/macmaneslab/devon/bin:$PATH

srun echo "      /\^._.^/\     starting process: `date`"

amptk dada2 \
--fastq /mnt/lustre/macmaneslab/devon/guano/NAU/OahuBird/illumina/p10-2/readfiltd/readfiltd_p102.demux.fq \
--out p10-2 \
--length 180 \
--platform illumina \
--uchime_ref COI

echo "      /\^._.^/\     ending process: `date`"
