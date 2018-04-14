Starting within the `$HOME` directory of the project:
```
mkdir fqraw
cd fqraw
mkdir lane1 lane2
rsync -avz foster@cobb.unh.edu:/data/foster/devon/180302/*/*.gz .
mv *L001*.gz ./lane1
mv *L002*.gz ./lane2
cd lane1
for i in `ls`; do mv "$i" lane1-"$i"; done
cd ../lane2
for i in `ls`; do mv "$i" lane2-"$i"; done
mv *.gz ../
cd ../lane1
mv *.gz ../

mv lane1-mock-IM4p11_2_TGCGTCAA-GTCTAGTG_L001_R1_001.fastq.gz lane1-mock-IM4p11-2_TGCGTCAA-GTCTAGTG_L001_R1_001.fastq.gz
mv lane1-mock-IM4p11_2_TGCGTCAA-GTCTAGTG_L001_R2_001.fastq.gz lane1-mock-IM4p11-2_TGCGTCAA-GTCTAGTG_L001_R2_001.fastq.gz
mv lane2-mock-IM4p11_2_TGCGTCAA-GTCTAGTG_L002_R1_001.fastq.gz lane2-mock-IM4p11-2_TGCGTCAA-GTCTAGTG_L002_R1_001.fastq.gz
mv lane2-mock-IM4p11_2_TGCGTCAA-GTCTAGTG_L002_R2_001.fastq.gz lane2-mock-IM4p11-2_TGCGTCAA-GTCTAGTG_L002_R2_001.fastq.gz
```

# Creating the working environment
```
conda create -n amptk python=3.6 biopython natsort pandas numpy matplotlib seaborn python-edlib edlib biom-format psutil
source activate amptk

conda install -c bioconda vsearch
conda install r-base bioconductor-dada2
conda install r-base bioconductor-phyloseq
conda install r-tidyverse

R
install.packages('devtools')
library('devtools')
# did not run: install_github("tobiasgf/lulu")
q()

git clone https://github.com/nextgenusfs/amptk.git
cd $HOME/bin
ln -s /mnt/lustre/macmaneslab/devon/pkgs/amptk/bin/amptk .
```

# renaming Cobb raw.fq files
```
rename -n asc_ asc- *.gz
rename mock_IM3 mock-IM3 *.gz
find . -name "*.gz" | sort | sed -e 's|^\./||' > raw_filelist.txt
cut -d '_' -f 2- raw_filelist.txt > wantedc2.txt
## added wanted metadata names to wantedc0.txt
nano wantedc0.txt ## then pasted wanted metadata names from metadata.csv file
while read line; do for i in {1..2}; do echo "$line"; done; done < wantedc0.txt > wantedc1.txt
paste wantedc1.txt wantedc2.txt -d '_' > wantedlist.txt
paste raw_filelist.txt wantedlist.txt > final_filelist.txt
mv _CGAGAGTT-CGTTACTA_L002_R1_001.fastq.gz mockIM3_CGAGAGTT-CGTTACTA_L002_R1_001.fastq.gz
mv  p_CGAGAGTT-CGTTACTA_L002_R2_001.fastq.gz mockIM3_CGAGAGTT-CGTTACTA_L002_R2_001.fastq.gz
re name NTC_ NTC- *
```

# installing amptk
Installed v-1.1.3 `amptk` within new Conda environment:  

# running amptk
Illumina step:
```
#!/bin/bash

#SBATCH -D /mnt/lustre/macmaneslab/devon/guano/Pompton/illumina
#SBATCH -p shared,macmanes
#SBATCH --job-name="PompIll"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --output=Pompton_ampIllumina.log

module purge
module load anaconda/colsa
source activate amptk
PATH=/mnt/lustre/macmaneslab/devon/bin:$PATH

srun echo "      /\^._.^/\     starting process: `date`"

amptk illumina \
-i /mnt/lustre/macmaneslab/devon/guano/Pompton/fqraw \
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
```

This kept getting hung up on the demultiplexing step - the temporary directory that was generated with all the individual demultiplexed reads were present, so a single `cat` command sufficed to resolve the issue - note this resulted in an incomplete log file though.  
> perform this within the parent directory of where the temporary directory with all the individual files are:  
`cat trim_pomp/*.demux.fq | gzip > trim_pomp.demux.fq`
