# updating BOLD database
See the [bold R package](https://www.r-project.org/nosvn/pandoc/bold.html) site for install and usage details.

## install dependencies in R
In my system the developmental version could only be installed; I also then had to manually install an updated version of the `curl` R package _and had to load that package first before BOLD_:  

```
R
devtools::install_github("ropensci/bold")
install.packages("curl")
library('curl')
library('bold')
```

The goal is to create a tab-delimited file we can do a couple of things with:
1. export and use with Jon's python script to generate an fasta file formatted for UTAX;
2. use for our own purposes for data filtering, plotting, etc.

Jon's scripts require a number of fields to be retained; we want to keep all of these and a few more for our own purpose. See line 39 in [his Python script](https://github.com/nextgenusfs/amptk/blob/master/util/bold2utax.py) for the list of names he wants to retain.  

(1) First we're going to collect all the information in one chunk. This creates a master object (a data.frame) containing all specimen and sequence data for the BOLD record; we don't want all of that! (2) Next we're going to parse that data.frame object to grab just the data we want. (3) With that done we'll write the object as a text file to disk.  

```
#!/bin/bash

#SBATCH -D /mnt/lustre/macmaneslab/devon/guano/BOLDdb
#SBATCH -p shared,macmanes
#SBATCH --job-name="boldDBarth"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --output=boldRarth.log

module purge
module load linuxbrew/colsa
PATH=/mnt/lustre/macmaneslab/devon/bin:$PATH

srun echo "      /\^._.^/\     starting process: `date`"

library('bold',lib="/mnt/lustre/macmaneslab/devon/R/x86_64-pc-linux-gnu-library/3.4")
library('curl',lib="/mnt/lustre/macmaneslab/devon/R/x86_64-pc-linux-gnu-library/3.4")

tmpBOLD.df <- bold_seqspec(taxon='Osmia', sepfasta=FALSE)
allarth.df <- subset(tmpBOLD.df, select= c(processid,sampleid,bin_uri,genbank_accession,phylum_name,class_name,order_name,family_name,genus_name,species_name,lat,lon,country,province_state,marker_codes,identification_provided_by,nucleotides))
write.table(allarth.df, file ="arthBOLDdb.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

echo "      /\^._.^/\     ending process: `date`"
```
