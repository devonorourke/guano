# Introduction
All information created for this project is available at [this Github repo](https://github.com/devonorourke/guano/tree/master/BRIpompton). Please visit that page for more information regarding data tables, visualizations, and code used to complete this work.

## Molecular work  
Avian fecal samples were collected by Biodiversity Research Institute staff, led by Oksana Lane. Samples were placed directly in Qiagen/MoBio PowerFecal kit collection tubes to avoid potential cross contamination. Earlier work with avian samples suggested that this kit works well to extract DNA but potential PCR inhibitors remain in solution that prevent successful amplification of the arthropod COI target. Thus we followed the standard MoBio PowerFecal guide, eluting in 60 uL of Solution C6. The eluent was then subject to a 2x Ampure XP cleanup, and eluted in 20 uL of Nuclease Free Water (NFW).

COI amplicons were selectively amplified from the extracted and cleaned avian fecal DNA using custom primers adapted from Jusino 2017 (see preprint [here] (DOI: 10.7287/peerj.preprints.3184v1)); the key modification is that my design incorporates the entire Illumina adapter, barcode, pad, and linker, plus the COI primer sequence into a single oligo, rather than employing a 2-step process of COI PCR and subsequent adapter ligation. PCR products were pooled in equimolar fashion except when concentrations fell below 1.0 nM; those samples for which there was less than that mass of DNA were pooled with a maximum of 20 uL per sample. The subsequent pool was then filtered using the QiaQuick PCR cleanup spin column. The library was then submitted to Hubbard Center for Genome Studies at the University of New Hampshire.

## sequencing at UNH
The pooled library of COI amplicons were sequenced using a HiSeq 2500 platform following 250 bp PE sequencing using V3 chemistry set for a Rapid run on March 1, 2018. Raw numbers of reads and general run metrics are available to view from [this link](http://cobb.unh.edu/180302_orourke_P11_2_B_DemuxStats.html). Though it was not intended, the single library of pooled samples was spiked in on two separate lanes (we requested just one); thus there was about twice the expected number of sequences generated - reads of each sample among the two lanes were combined and renamed [described here](https://github.com/devonorourke/guano/blob/master/BRIpompton/docs/pompton-rename.md). These raw combined reads were then trimmed and quality filtered, clustered, and assigned taxonomy through the `amptk` pipeline described below.  

About 6.1 million reads were generated across 105 samples. The positive control (mock community) spike in was well balanced (containing about 30,000 reads), with DNA extraction negative controls showing among the lowest read counts (~ 2000) and the PCR negative demonstrating zero reads. The true samples comprised the majority of the sequence data, with about 98.7% of all data dedicated to bird fecal samples.  

## Creating the working environment for `amptk`
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

The current versions among the core programs used in the `amptk` pipleine are as follows:
- amptk v. 1.1.3-36d7eda
- usearch9 v9.2.64_i86linux32
- usearch10 v10.0.240_i86linux32
- vsearch  v2.7.0_linux_x86_64

Remaining python modules and R dependencies were installed via Conda (upgrade/updates with `pip` and/or `conda` performed 13-April-2018); install commands were:
```
pip install -U -I biopython natsort pandas numpy matplotlib seaborn edlib biom-format psutil
conda install r-base bioconductor-dada2
conda install r-base bioconductor-dada2
```  

The environment was activated as: `source activate amptk`. All bioinformatic processes described below occurred within this virtual environment.  

# amptk pipeline

[amptk](https://github.com/nextgenusfs/amptk) is a bioinformatic toolkit which performs all necessary tasks beginning with quality and adapter trimming of raw reads, clustering OTUs, denoising and chimera detection, through to assigning taxonomy to each identified cluster and generating (among other outputs) the list of per-sample taxa represented in the dataset. A full documentation of available parameters used for the program are [detailed here](http://amptk.readthedocs.io/en/latest/index.html).

As noted above, a similarly named virtual environment "`amptk`" was created when completing the installation process. Note that recent versions of `amptk` are now compatible with Python3, thus a virtual environment was created with that Python version to reflect the change.  

The following sections reflect the basic outlines of the scripts used to execute the `amptk` commands. However, as these jobs were submitted through a compute cluster they do not reflect the entirety of the script. See the [scripts](https://github.com/devonorourke/guano/tree/master/BRIpompton/scripts) directory for full details of each script employed.

Finally, though not implied in either the script nor in the code blocks below, child directories were created for each `amptk` process within the parent directory for this process (for example, the raw .fastq files are in `../project/fqraw`, the output of `amptk illumina` reside within `../project/illumina`, the output from the `amptk taxonomy` commands reside within `../project/taxonomy`, etc.).  

## adapter trimming and PE merging

The first step in the pipeline trims adapters (as a result of the insert length being less than the read length) and then uses USEARCH to merge paired end reads. Orphaned reads are discarded (this typically accounts for less then 2% of the overall number or reads in a sample).

```
amptk illumina \
-i /mnt/lustre/macmaneslab/devon/guano/NAU/Perlut/fqraw \
-o trim \
--rescue_forward on \
--require_primer off \
--min_len 160 \
--full_length \
--read_length 300 \
-f GGTCAACAAATCATAAAGATATTGG \
-r GGWACTAATCAATTTCCAAATCC \
--cpus 24 \
--cleanup
```

A single **.fastq.gz** file is output by concatenating all the individual paired reads with headers modified to specify the sample name. In addition, the **.amptk-demux.log** file documents the proportion of merged reads per sample. The library was fairly well balanced, with 76 samples having more than 5,000 reads each; all four of our DNA extraction negative controls had less than 5,000 reads but all had enough to warrant inclusion in the analysis to investigate low-level contamination. In fact, all samples were included in the subsequent filtering analyses to determine what minimal read threshold was required for inclusion in the clustering step.  

> note - to view the read counts from the `.amptk-demux.log` file, just run this little one-liner:
```
sed -n '/Found.*.barcoded samples/,$p' trim_pomp.amptk-demux.log | grep -v "^\[" | sed 's/^[[:blank:]]*//'
```

## dropping samples

There is an important tradeoff between the likelihood that a read is the result of index bleed versus a true representation of the amplicons in a sample; if one is to account and filter for index-bleed, then one is to likely reduce the number of reads in a sample. Because the mock community proportion of reads was well balanced in this run, the likelihood of index bleed is no higher among known community members in the positive control than any true sample (or negative control); we therefore run little risk in initially including all of our samples going into the preliminary clustering process. Nevertheless, we may revisit dropping samples after our first filtering estimations following clustering if it is determined that certain samples with low read numbers are likely containing just contaminant reads or those whose reads may be reduced due to our index-bleed calculations.  

Clustering was therefore performed initially on all reads without any samples dropped. One final note - it could be easy to just set an arbitrary threshold of 5000 reads, as this would then remove all NTC samples and about 20% of our true samples; we could then assert that all we'd need to filter thereafter was some proportion of index-bleed. However this approach doesn't leverage the fact that we can learn something about what OTUs may need to have their reads further reduced if they in fact are present in the negative controls. The power in keeping in those NTCs initially is that we can better assess which OTUs, and at what read depth, should be filtered.  

## clustering for OTUs

This is a two step process in which the **.fastq** file containing all reads is parsed first using the `DADA2` algorithm creating **Amplicon Sequence Variants (ASV)** (formerly termed **iSeqs**) candidate sequences. These unique sequences are then clustered to a specified similarity threshold (97%) using a traditional `UCLUST` approach. See Jon's documentation describing the differences [here](http://amptk.readthedocs.io/en/latest/clustering.html). In brief, **ASV** values are clustered at a 100% identity, whereas the resulting **OTUs** are clustered at 97% identity, meaning that the **ASV** sequences are more exclusive than the **OTUs**.  

In addition there is a chimera filtering step applied to the data; this requires the installation of the COI database provided by amptk:

```
amptk install -i COI
```

Then execute the clustering with the following code:

```
amptk dada2 \
--fastq /mnt/lustre/macmaneslab/devon/guano/Pompton/illumina/trim_pomp.demux.fq.gz \
--out rough \
--length 180 \
--platform illumina \
--uchime_ref COI
```

The output contains a pair of files which are applied in the next filtering strategy (for index bleed): the `.cluster.otu_table.txt` file which follows a traditional OTU matrix format, as well as the accompanying `.cluster.otus.fa` file which contains the OTU id in the header and the associated sequence. Each dataset is then filtered according to the following commands (see **Filtering** section below).  

We find that after quality and chimera filtering, about **86 % of reads** were used to identify unique sequence variants (ASVs) or clusters (OTUs). The following table summarizes information contained in the `Pompton_clust-rough.log` file generated in the `amptk clust` command:

|  | allSamples |
| --- | --- |
| # OTUs clustered | 2,152 |
| # reads mapped to OTUs | 2,595,284 |
| # ASVs clustered |  4,305 |
| # reads mapped to ASVs | 2,603,097 |


## filtering

Because a mock community was added to this project, the proportion of reads that are likely misassigned can be estimated on a per-OTU basis. In brief, we are certain of the OTUs likely to be present in mock community; any additional OTU is the result of index bleed (or potentially contamination, but this is unlikely given that this positive control is separately processed from the negative controls and true samples). By calculating the proportion of reads that are present in our mock sample which _shouldn't be there_ we can estimate what fraction of reads (on a per-OTU basis) should be subtracted from all true samples.

This process takes place by applying an initial filtering step that filters reads using the most strict criteria (taking the largest instance of an OTU bleed and applying that percentage to filter across true samples); intermediate files are kept to investigate how the index-bleed is distributed on a per-OTU basis. I have maintained a [separate document](https://github.com/devonorourke/guano/blob/master/Perlut/Perlut_filtering_notes.md) describing the detailed steps used to apply what I feel are the most appropriate filtering strategies for this dataset.   

In brief, this amounted to determining the proportion of index bleed _into the mock community_, the proportion of index bleed _from mock community into true samples_, and identifying any OTUs which were clustered yet not identified in the mock fasta file.  

The filtering results from the `dropd` and `trimd` datasets are as follows:  
- **For `trim`**, the initial dataset retained all **89** total samples, **3,543,867** reads (though > 1,245,00 were from the mock community alone), and **1912 OTUs**. Following filtering, we retained **2,198,389 reads** and **1,094 OTUs** (notably all mock reads have been removed).  
- **For `dropd`** the initial dataset was reduced to just **38** samples (including the mock community) prior to filtering, with **3,472,564 reads** retained among **1,792 OTUs**. Following filtering the dataset included **2,194,756 reads** and **1,254 OTUs**.     

A pair of output files after completing filtering steps are applied to then assign taxonomic information to our remaining reads. Following these steps, the `trim` dataset was used exclusively to assign taxonomy to the OTUs present.  

## taxonomy assignment
As described in the [amptk taxonomy](http://amptk.readthedocs.io/en/latest/taxonomy.html) section, the database used to assign taxonomy is derived from the Barcode of Life Database ([BOLD](http://v4.boldsystems.org/)). The sequences present in the database we're using are the result of two sequential clustering processes.
- BOLD's BIN data serve as the initial sequence material. These sequences themselves are initially derived from [a clustering process](http://v4.boldsystems.org/index.php/Public_BarcodeIndexNumber_Home).
- The BIN sequences are then clustered locally by amptk to 99% identity. These data are further processed to train the UTAX program which can be used in taxonomic assignment/prediction.  

> This database was updated as of 14-sept-2017, following the [release](https://github.com/nextgenusfs/amptk/releases/tag/1.0.0) of amptk v-1.0.0.  
> Complete database download is [available here](https://osf.io/4xd9r/files/)

Taxonomy was explored using both the _hybrid_ approach (default in amptk) as well as a _usearch only_ approach. In both instances, the `--method` reflected the given approach. See Jon's description of the steps used in his documentation at the link above.  

The following code was applied (in this example the method is the `usearch` approach):  

```
amptk taxonomy \
--i /mnt/lustre/macmaneslab/devon/guano/NAU/Perlut/filtd/finaltrim.final.txt \
--fasta /mnt/lustre/macmaneslab/devon/guano/NAU/Perlut/filtd/finaltrim.filtered.otus.fa \
--out Perlut_u \
--db COI \
--method usearch \
--mapping_file /mnt/lustre/macmaneslab/devon/guano/NAU/Perlut/illumina/dropd.mapping_file.txt
```

The output fasta sequence and OTU table with taxonomic information were uploaded to [the Github repo](https://github.com/devonorourke/guano/tree/master/Perlut).  

> Note that a value of 0 ("absence") could mean a variety of different things:  
> - it could be that the OTU is not truly in the sample of guano
> - it could be because an OTU was present but not amplified and sequenced  
> - it could be the OTU was sequenced but there wasn't enough reads to pass our filters (with `--index_bleed` and `--subtract` arguments in `amptk filter`))

 # Further analyses

 An R script was then used to manipulate the output `Perlut.otu_table.taxonomy.txt` file which includes both further data filtering, as well as the calculations for frequency tables and visualizations - [see here](https://github.com/devonorourke/guano/blob/master/Perlut/OTUanalysis.R).  

 > One such data filtering taking place here is the removal of reads associated with the mock community.

Please see the Perlut [Github repo](https://github.com/devonorourke/guano/blob/master/Perlut) for subseuqent data summaries and visualizations.  
