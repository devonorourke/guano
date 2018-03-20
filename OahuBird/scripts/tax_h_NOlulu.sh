#!/bin/bash

#SBATCH -D /mnt/lustre/macmaneslab/devon/guano/NAU/OahuBird/taxonomy/nolulu
#SBATCH -p shared,macmanes
#SBATCH --job-name="ampTAX.nolulu"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --output=OahuBird_NOlulu_tax.log

module purge
module load anaconda/colsa
source activate amptk_env
PATH=/mnt/lustre/macmaneslab/devon/bin:$PATH

srun echo "      /\^._.^/\     starting process: `date`"

amptk taxonomy \
-i /mnt/lustre/macmaneslab/devon/guano/NAU/OahuBird/filtd/all/OahuBird.final.txt \
--fasta /mnt/lustre/macmaneslab/devon/guano/NAU/OahuBird/filtd/all/OahuBird.filtered.otus.fa \
--out OahuBird_h \
--db COI \
--method hybrid \
--mapping_file /mnt/lustre/macmaneslab/devon/guano/NAU/OahuBird/taxonomy/FilteredMappingFile.txt

echo "      /\^._.^/\     ending process: `date`"
