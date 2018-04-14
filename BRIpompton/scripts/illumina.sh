#!/bin/bash

#SBATCH -D /mnt/lustre/macmaneslab/devon/guano/Pompton/illumina
#SBATCH -p shared,macmanes
#SBATCH --job-name="PompIll"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --output=Pompton_ampIllumina.log

module purge
module load anaconda/colsa
source activate amptk_env
PATH=/mnt/lustre/macmaneslab/devon/bin:$PATH

srun echo "      /\^._.^/\     starting process: `date`"

amptk illumina \
-i /mnt/lustre/macmaneslab/devon/guano/Pompton/fqraw/catfiles \
-o trim_pomp \
--rescue_forward on \
--require_primer off \
--min_len 160 \
--full_length \
--read_length 250 \
-f GGTCAACAAATCATAAAGATATTGG \
-r GGWACTAATCAATTTCCAAATCC \
--cpus 24 \
--cleanup

echo "      /\^._.^/\     ending process: `date`"
