#!/bin/bash

#SBATCH -D /mnt/lustre/macmaneslab/devon/guano/NAU/OahuBird/taxonomy/lulu
#SBATCH -p shared,macmanes
#SBATCH --job-name="ampTAX.lulu"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --output=OahuBird_lulu_tax.log

module purge
module load anaconda/colsa
source activate amptk_env
PATH=/mnt/lustre/macmaneslab/devon/bin:$PATH

srun echo "      /\^._.^/\     starting process: `date`"

amptk taxonomy \
-i /mnt/lustre/macmaneslab/devon/guano/NAU/OahuBird/lulu/pid95/pid95.lulu.otu_table.txt \
--fasta /mnt/lustre/macmaneslab/devon/guano/NAU/OahuBird/lulu/pid95/pid95.lulu.otus.fa \
--out OahuBird_lulu_h \
--db COI \
--method hybrid \
--mapping_file /mnt/lustre/macmaneslab/devon/guano/NAU/OahuBird/taxonomy/FilteredMappingFile.txt

echo "      /\^._.^/\     ending process: `date`"
