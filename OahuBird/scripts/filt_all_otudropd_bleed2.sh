#!/bin/bash

#SBATCH -D /mnt/lustre/macmaneslab/devon/guano/NAU/OahuBird/filtd/all
#SBATCH -p macmanes,shared
#SBATCH --job-name="ampFiltallOTUdropd_bleed2"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --output=OahuBird_filter_allOTUdropd_bleed2.log

module purge
module load linuxbrew/colsa
PATH=/mnt/lustre/macmaneslab/devon/bin:$PATH

srun echo "      /\^._.^/\     starting process: `date`"

amptk filter \
-i /mnt/lustre/macmaneslab/devon/guano/NAU/OahuBird/clust/all/OahuBird_otudropd.cleaned.otu_table.txt \
-f /mnt/lustre/macmaneslab/devon/guano/NAU/OahuBird/clust/all/OahuBird_otudropd.cleaned.otus.fa \
-b mockOahuBirdIM4 \
--delimiter csv \
--index_bleed 0.02 \
--subtract auto \
--keep_mock \
--mc /mnt/lustre/macmaneslab/devon/guano/mockFastas/CFMR_insect_mock4alt.fa \
--debug \
--threshold max \
-o allOTUdropd_bleed2 \
--normalize n

echo "      /\^._.^/\     ending process: `date`"
