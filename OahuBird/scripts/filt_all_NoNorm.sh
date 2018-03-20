#!/bin/bash

#SBATCH -D /mnt/lustre/macmaneslab/devon/guano/NAU/OahuBird/filtd/all
#SBATCH -p macmanes,shared
#SBATCH --job-name="ampFiltallNoNorm"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --output=OahuBird_filter_all_NoNorm.log

module purge
module load linuxbrew/colsa
PATH=/mnt/lustre/macmaneslab/devon/bin:$PATH

srun echo "      /\^._.^/\     starting process: `date`"

amptk filter \
-i /mnt/lustre/macmaneslab/devon/guano/NAU/OahuBird/clust/all/OahuBird.cluster.otu_table.txt \
-f /mnt/lustre/macmaneslab/devon/guano/NAU/OahuBird/clust/all/OahuBird.cluster.otus.fa \
-b mockOahuBirdIM4 \
--delimiter csv \
--keep_mock \
--calculate all \
--subtract auto \
--mc /mnt/lustre/macmaneslab/devon/guano/mockFastas/CFMR_insect_mock4alt.fa \
--debug \
--threshold max \
-o all_NoNorm \
--normalize n

echo "      /\^._.^/\     ending process: `date`"
