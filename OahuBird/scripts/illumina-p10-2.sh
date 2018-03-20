#!/bin/bash

#SBATCH -D /mnt/lustre/macmaneslab/devon/guano/NAU/OahuBird/illumina/p10-2
#SBATCH -p shared,macmanes
#SBATCH --job-name="sAMPill2"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --output=OahuBird_ampIllumina-2.log

module purge
module load anaconda/colsa
source activate amptk_env
PATH=/mnt/lustre/macmaneslab/devon/bin:$PATH

srun echo "      /\^._.^/\     starting process: `date`"

amptk illumina \
-i /mnt/lustre/macmaneslab/devon/guano/NAU/OahuBird/fqraw/p10-2 \
-o trim_p10-2 \
--rescue_forward on \
--require_primer off \
--min_len 160 \
--full_length \
--read_length 300 \
-f GGTCAACAAATCATAAAGATATTGG \
-r GGWACTAATCAATTTCCAAATCC \
--cpus 24 \
--cleanup

echo "      /\^._.^/\     ending process: `date`"
