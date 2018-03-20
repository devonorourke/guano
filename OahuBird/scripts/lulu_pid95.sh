#!/bin/bash

#SBATCH -D /mnt/lustre/macmaneslab/devon/guano/NAU/OahuBird/lulu/pid95
#SBATCH -p macmanes,shared
#SBATCH --job-name="ampLulu"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --output=lulu_pid95.log

module purge
module load linuxbrew/colsa
PATH=/mnt/lustre/macmaneslab/devon/bin:$PATH

srun echo "      /\^._.^/\     starting process: `date`"

amptk lulu \
-i /mnt/lustre/macmaneslab/devon/guano/NAU/OahuBird/filtd/all/OahuBird.final.txt \
-f /mnt/lustre/macmaneslab/devon/guano/NAU/OahuBird/filtd/all/OahuBird.filtered.otus.fa \
-o pid95 \
--min_match 95

echo "      /\^._.^/\     ending process: `date`"
