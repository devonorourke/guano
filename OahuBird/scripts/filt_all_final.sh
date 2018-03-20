#!/bin/bash

#SBATCH -D /mnt/lustre/macmaneslab/devon/guano/NAU/OahuBird/filtd/all
#SBATCH -p macmanes,shared
#SBATCH --job-name="ampFiltallFinal"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --output=OahuBird_filter_allFinal.log

module purge
module load linuxbrew/colsa
PATH=/mnt/lustre/macmaneslab/devon/bin:$PATH

srun echo "      /\^._.^/\     starting process: `date`"

amptk filter \
-i /mnt/lustre/macmaneslab/devon/guano/NAU/OahuBird/clust/all/OahuBird_otudropd.cleaned.otu_table.txt \
-f /mnt/lustre/macmaneslab/devon/guano/NAU/OahuBird/clust/all/OahuBird_otudropd.cleaned.otus.fa \
-b mockOahuBirdIM4 \
--delimiter tsv \
--index_bleed 0.01 \
--subtract 50 \
--mc /mnt/lustre/macmaneslab/devon/guano/mockFastas/CFMR_insect_mock4alt.fa \
--debug \
--threshold max \
-o OahuBird \
--normalize n

echo "      /\^._.^/\     ending process: `date`"
