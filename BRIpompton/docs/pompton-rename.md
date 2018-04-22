Reads from the sequencing center were demultiplexed by sample name and lane number, and placed within sample-named parent directories. Because the library was run among two independent lanes we generated technical replicates that were combined before processing in the `amptk` pipeline. Thus a single sample could have initially expected to have the following file structure:  
```
#parent directory named '.../Sample_pomp_4584/'
pomp-4584_TGAGTACG-TAGCGAGT_L001_R1_001.fastq.gz
pomp-4584_TGAGTACG-TAGCGAGT_L002_R1_001.fastq.gz
pomp-4584_TGAGTACG-TAGCGAGT_L001_R2_001.fastq.gz
pomp-4584_TGAGTACG-TAGCGAGT_L002_R2_001.fastq.gz
```

Files were transferred from the individual parent directories recursively and put into a new directory on Premise:
> The path to the parent directory for these raw .fastq files was: `/mnt/lustre/macmaneslab/devon/guano/Pompton/fqraw`

Files with matching base names were concatenated into (fake lane number) new file per forward/reverse read pair.   
```
cd fqraw
mkdir lane1 lane2
rsync -avz foster@cobb.unh.edu:/data/foster/devon/180302/*/*.gz .

for i in `ls | cut -d '_' -f 1,2 | sort -u`;
do
cat "$i"_L001_R1_001.fastq.gz "$i"_L002_R1_001.fastq.gz > "$i"_L003_R1_001.fastq.gz
cat "$i"_L001_R2_001.fastq.gz "$i"_L002_R2_001.fastq.gz > "$i"_L003_R2_001.fastq.gz
done
```

Using the example of files listed above, the output generated is a single read pair set:
```
pomp-4584_TGAGTACG-TAGCGAGT_L003_R1_001.fastq.gz
pomp-4584_TGAGTACG-TAGCGAGT_L003_R2_001.fastq.gz
```

The initial reads were removed, and these concatenated reads were then used as the raw data for the `amptk` workflow.



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
