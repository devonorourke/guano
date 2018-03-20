#!/bin/bash

#SBATCH -D /mnt/lustre/macmaneslab/devon/guano/NAU/OahuBird/filtd/p10-1
#SBATCH -p shared,macmanes
#SBATCH --job-name="amp101Filt"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --output=OahuBird_p10-1_filter_default.log

module purge
module load linuxbrew/colsa
PATH=/mnt/lustre/macmaneslab/devon/bin:$PATH

srun echo "      /\^._.^/\     starting process: `date`"

amptk filter \
-i /mnt/lustre/macmaneslab/devon/guano/NAU/OahuBird/clust/p10-1/p10-1.cluster.otu_table.txt \
-f /mnt/lustre/macmaneslab/devon/guano/NAU/OahuBird/clust/p10-1/p10-1.cluster.otus.fa \
-b mockIM4p10L1 \
--delimiter csv \
--keep_mock \
--calculate all \
--subtract auto \
--mc /mnt/lustre/macmaneslab/devon/guano/mockFastas/CFMR_insect_mock4alt.fa \
--debug \
--threshold max \
-o 10-1_default \
--normalize y

echo "      /\^._.^/\     ending process: `date`"
