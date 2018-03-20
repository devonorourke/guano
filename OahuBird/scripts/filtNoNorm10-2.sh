#!/bin/bash

#SBATCH -D /mnt/lustre/macmaneslab/devon/guano/NAU/OahuBird/filtd/p10-2
#SBATCH -p shared,macmanes
#SBATCH --job-name="amp101Filt"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --output=OahuBird_p10-2_filter_noNorm.log

module purge
module load linuxbrew/colsa
PATH=/mnt/lustre/macmaneslab/devon/bin:$PATH

srun echo "      /\^._.^/\     starting process: `date`"

amptk filter \
-i /mnt/lustre/macmaneslab/devon/guano/NAU/OahuBird/clust/p10-2/p10-2.cluster.otu_table.txt \
-f /mnt/lustre/macmaneslab/devon/guano/NAU/OahuBird/clust/p10-2/p10-2.cluster.otus.fa \
-b mockP10L2IM4 \
--delimiter csv \
--keep_mock \
--calculate all \
--subtract auto \
--mc /mnt/lustre/macmaneslab/devon/guano/mockFastas/CFMR_insect_mock4alt.fa \
--debug \
--threshold max \
-o 10-2_noNorm \
--normalize n

echo "      /\^._.^/\     ending process: `date`"
