# Introduction
All information created for this project is available at [this Github repo](https://github.com/devonorourke/guano/tree/master/Mangan). Please visit that page for more information regarding data tables| visualizations| and code used to complete this work.

## Molecular work  
Guano was provided as both single sample and pooled samples. Samples were processed either in 96-well plates or as single tube extractions using Qiagen's PowerSoil kits. Negative controls were used in both plates and single tubes to evaluate potential for cross-contamination during extraction.
Arthropod COI sequences were amplified using custom primers which targets a 180 bp region of cytochrome oxidase subunit-1. The resulting PCR products were normalized using Applied Biosystem's SequalPrep kit. The resulting normalized PCR products were pooled in fixed volumes and submitted for sequencing.
> These primers contain unique forward and reverse barcode/indices as well as full length Illumina primer index and flow-cell binding sequence. The result of a single PCR run generates ready-to-sequence amplicons.

## sequencing at NAU
The pooled library of COI amplicons were sequenced on an Illumina MiSeq machine with 300 bp PE sequencing using V2 chemistry set for 600 cycles at Northern Arizona University's sequencing center on September 14, 2018. **14,490,123 raw reads** were generated

## Renaming
Renamed the default raw .fastq files to something easier to work with:

For the `CONTROL` files:
```
rename 's/-xx-xx-USA-2017-076-JF_//g' *.gz
rename 's/CONTROL-//g' *.gz
```

For the `NHCS` files:
```
rename 's/-xx-VE-USA-2017-076-JF_.*.L/_L/g' *.gz
rename 's/NHCS-//g' *.gz
```

This results in control files having names like:
```
ExtractionNTC9S97_L001_R1_001.fastq.gz
ExtractionNTC9S97_L001_R2_001.fastq.gz
```

and true samples having names like:
```
6212017EGA1_L001_R1_001.fastq.gz
6212017EGA1_L001_R2_001.fastq.gz
```

## Installation specifications


The current versions among the core programs used in the `amptk` pipleine are as follows:
- AMPtk v1.1.3-36d7eda
- usearch9 v9.2.64_i86linux32
- usearch10 v10.0.240_i86linux32
- vsearch  v2.7.0_linux_x86_64

# amptk pipeline

[amptk](https://github.com/nextgenusfs/amptk) is a bioinformatic toolkit which performs all necessary tasks beginning with quality and adapter trimming of raw reads| clustering OTUs| denoising and chimera detection| through to assigning taxonomy to each identified cluster and generating (among other outputs) the list of per-sample taxa represented in the dataset. A full documentation of available parameters used for the program are [detailed here](http://amptk.readthedocs.io/en/latest/index.html).

The following sections reflect the basic outlines of the scripts used to execute the `amptk` commands. Note that commands were submitted through a compute cluster and thus do not reflect the entirety of the script full script - we're showing just the arguments necessary to recapitulate the `amptk`-specific code.

Finally| though not implied in either the script nor in the code blocks below| subdirectories were created for each `amptk` process prior to executing the command. The directories were typically named for this process (for example| the raw .fastq files are in `../project/fqraw`| the output from the `amptk taxonomy` commands reside within `../project/taxonomy`| etc.).  

# adapter trimming and PE merging

The first step in the pipeline trims adapters (as a result of the insert length being less than the read length) and then uses USEARCH to merge paired end reads. Orphaned reads are discarded (this typically accounts for less then 2% of the overall number or reads in a sample).

> The following script was executed in `$PROJECTPATH/trimfq`

```
amptk illumina \
-i /mnt/lustre/macmaneslab/devon/guano/NAU/Mangan/rawfq \
--out trim \
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

A single **.fastq.gz** file is output by concatenating all the individual paired reads with headers modified to specify the sample name. In addition| the **.amptk-demux.log** file documents the proportion of merged reads per sample. All but 7 of the 308 samples had sufficient reads for analysis.

> note - to view the read counts from the `.amptk-demux.log` file| just run this little one-liner:
```
sed -n '/Found.*.barcoded samples/|$p' trim.amptk-demux.log | grep -v "^\[" | sed 's/^[[:blank:]]*//'
```

Notable concerns include substantial reads associated with negative control samples. We'll need to tack those samples and filter out possible contaminants from the entire dataset. For now| we'll remove the 7 samples with insufficient read depth:

# dropping samples

There is an important tradeoff with retaining all reads in all samples| versus retaining fewer reads in fewer samples: the likelihood of a false positive. In our dataset| false positives can occur in multiple ways: chimeric reads| mis-assigned barcodes| and PCR error. Retaining samples with very low read depth runs the risk of keeping reads which are near or below the index (barcode) switching error rate; we find it more appropriate to drop samples when there are fewer than a few thousand reads. We'll do that here - note there are a few ways in `amptk` to do this| but we'll just specify dropping any samples with less than 1000 reads (note that of the 7 samples we're dropping| the one with the _greatest_ number of reads is just 47):
> We are performing this task within the same path as the input: `$PROJECTPATH/trimfq`

```
amptk remove \
-i /mnt/lustre/macmaneslab/devon/guano/NAU/Mangan/trimfq/trim.demux.fq.gz \
-t 1000 \
-o dropd.demux.fq
```

We'll use the output `dropd.demux.fq` file for the next step in our pipeline: clustering OTUs from the available sequence data.

# clustering for OTUs

AMPTK's clustering pipeline (as we're using it) is a three step process in which:
1. The **dropd.demux.fq** file containing relevant reads is parsed first using `VSEARCH` to assess read quality
2. The subsequent quality-filtered reads are then clustered as OTUs. We're using the `UNOISE3` algorithm to cluster here| which "clusters" at 100% identity| meaning that we're keeping all unique sequence variants (with some built-in error modeling to merge variation likely due to sequencing). In addition to the clustering iteself| this process includes a chimera-filtering step driven by the COI database created with:
```
amptk install -i COI
```
3. The resulting _Zotus_ produced with UNOISE3 are then clustered at the conventional 97% identity with `UCLUST`.  

> See Jon's documentation describing the differences [here](http://amptk.readthedocs.io/en/latest/clustering.html) regarding other clustering options.

The command itself executed for clustering is:
```
~/bin/amptk unoise3 \
-i dropd.demux.fq \
-o unoise \
--uchime_ref COI
```

The clustering process results with two datasets consisting of the same kind of infromation - an OTU table and a fasta file. Each table/fasta pair relates to the clustering process mentioned above: one table/fasta contains 100% identity clusters| while the other contains 97% identity clusters. I've found that for the sake of performing simple diet analyses| clustering works fine. This is because we next append taxonomic information to these OTUs| and the resulting taxonomic identities of unique (at 100% identity) OTUs often collapse to the same taxonony| which would have been collapsed if we had used the 97% identity clustering.

We find that **11,544,724 reads** pass the quality filtering from the initial **13,882,162 reads** analyzed. These reads then are reduced to **1,517,626 unique sequences** which are then deduplicated within `UNOISE3`. These unique sequence variants are denoised to just **4,202** variants, of which **3,583** unique sequence variants pass chimera filtering (**520** (12.4%) chimeras, **3583** (85.3%) non-chimeras, and **99** (2.4%) borderline sequences in 4,202 total sequences. The **4,202** unique variants are then clustered into OTUs using a 97% identity threshold, resulting in **1,927** OTUs (Singletons: **1187**, 33.1% of seqs, 61.6% of clusters.

This last point is critical to notice, as it can greatly influence the interpretation in diversity metrics:
```
Singletons: **1187**, 33.1% of seqs, 61.6% of clusters.
```

This means that 33% of the data are sequences which occur only once in the dataset - they exist in only one sample. Likewise, over 60% of the OTUs (the clusters) are singletons. If we were to filter out singletons, we would be dropping the majority of our OTUs, and over 30% of our sequence data. The question is whether or not we want to analyze the rare variants or not.

We'll use  the 97% identity clustered output and apply that information to the next filtering step (for index bleed): the `.cluster.otu_table.txt` file which follows a traditional OTU matrix format| as well as the accompanying `.cluster.otus.fa` file which contains the OTU id in the header and the associated sequence.

## filtering

Because no positive control (mock community) was used with this project| the proportion of reads that are likely misassigned can not be estimated empirically _within our specific dataset_. Instead| we will assign the default filtering parameter that Jon has observed in his work in generating this pipeline (default for MiSeq Illumina runs is 0.5% - see [Jon's notes here](https://amptk.readthedocs.io/en/latest/filtering.html)). We sacrifice dropping a few reads and rare OTUs for the benefit of increased likelihood of assigning reads to their appropriate sample. One minor difference in this filtering approach compared to the default parameters of `amptk`: because this dataset was very well balanced we don't noramlize data (thus `--normalize n`).

```
amptk filter \
-i unoise.cluster.otu_table.txt \
-f unoise.cluster.otus.fa \
--index_bleed 0.005 \
--debug \
--out Mangan \
--normalize n
```

What we find is that we have reduced our initial number of samples because we've reduced the number of reads (on a per-OTU basis) by 0.5%; thus samples with very low read numbers - like negative controls - are often dropped. This is indeed the case as we lose three NTC's but no true samples by using this filtering approach. Interestingly there is no difference in the number of project-wide OTUs pre and post filtering (it remains at **1,927 OTUs**). What we do notice is the number of OTUs per sample is reduced quite a bit for most samples, and some negative controls retain certain OTUs. In all, we have reduced our dataset from **11,634,568** to **11,504,725 reads** . We're going to keep these in our analyses, assign taxonomy, then think about filtering these out at a later step in the pipeline.

For example, this table shows how true samples typically result in drops in numbers of unique OTUs before and after applying this 0.5% filtering step. In the case of some negative controls (_ExtractionNTC*_) we find that there are relatively few to start with; in three instances the entirety of available OTUs are reduced to zero which is why we lose those NTCs from further analyses. Notably, no true samples were dropped from analysis by applying this filter.

| sample | reads | preFilt OTUs | postFilt OTUs |
| --- | --- | --- | --- |
|9152017HBPoolB2|38344|93|56
|9152017HBPoolB3|38151|81|52
|9152017HBPoolB4|35090|102|58
|ExtractionNTC1S1|68391|67|47
|ExtractionNTC2S10|30104|16|11
|ExtractionNTC5S37|7|6|6
|ExtractionNTC9S97|15|10|10
|ExtractionNTC15S151|5|4|4
|blankS39|37761|59|41


## taxonomy assignment
As described in the [amptk taxonomy](http://amptk.readthedocs.io/en/latest/taxonomy.html) section| the database used to assign taxonomy is derived from the Barcode of Life Database ([BOLD](http://v4.boldsystems.org/)). The sequences present in the database we're using are the result of two sequential clustering processes.
- BOLD's BIN data serve as the initial sequence material. These sequences themselves are initially derived from [a clustering process](http://v4.boldsystems.org/index.php/Public_BarcodeIndexNumber_Home).
- The BIN sequences are then clustered locally by amptk to 99% identity. These data are further processed to train the UTAX program which can be used in taxonomic assignment/prediction.  

> This database was updated as of 14-sept-2017| following the [release](https://github.com/nextgenusfs/amptk/releases/tag/1.0.0) of amptk v-1.0.0.  
> Complete database download is [available here](https://osf.io/4xd9r/files/)

Taxonomy was assigned using the _hybrid_ approach (default in amptk). See [Jon's description](http://amptk.readthedocs.io/en/latest/taxonomy.html#amptk-taxonomy) of the method for complete details - in brief, this approach:
1. Performs a global alignment of all OTU sequences against the BOLD database using `USEARCH` to find the most taxonomic information up to a specified cutoff (we use default parameters which keep any hits with >= 70% identity)
2. Runs the `UTAX` and `SINTAX` classifiers to generate taxonomic strings with scores that pass a specified identity cutoff (here 80% for either)
3. The best hit, if found in USEARCH, is retained if above 97% identity; if less than that, then the best `UTAX` or `SINTAX` classifier is retained.
> When we filter our results later, we'll indicate which classifier was used

The code:
```
amptk taxonomy \
-i Mangan.final.txt \
--fasta Mangan.filtered.otus.fa \
--out Mangan \
--db COI \
--mapping_file Mangan.mappingFile.txt
```

> The `Mangan.mappingFile.txt` file was created by modifying the provided "_Collected Guano Samples_" Excel spreadsheet to fit the required format of the conversion. This involved ascribing each sample with the metadata information (each row a unique observation (sample), each column a variable like Site, Roost, Date, etc.)

The output fasta file, OTU table, and biom file with taxonomic information were uploaded to [the data section](https://github.com/devonorourke/guano/tree/master/Mangan/data) of this Github repo.  

 # Further analyses

 An [R script](https://github.com/devonorourke/guano/blob/master/Mangan/docs/Mangan.OTUanalysis.R) was then used to manipulate the [Mangan.otu_table.taxonomy.txt](https://github.com/devonorourke/guano/blob/master/Mangan/data/Mangan.otu_table.taxonomy.txt) file to both further filter potential contaminant reads as well as produce visualizations in the [data](https://github.com/devonorourke/guano/tree/master/Mangan/data) subdirectory of this repository.
 Note that a major difference between the `Mangan.otu_table.taxonomy.txt` input and the resulting `master.df` file which encompasses these data is the transformation of sequence reads in their relative counts (total number of filtered sequences) to presence/absence counts. "_Absence_" could mean a variety of different things:
 - it could be that the OTU is not truly in the sample of guano
 - it could be because an OTU was present but not amplified and sequenced  
 - it could be the OTU was sequenced but there wasn't enough reads to pass our filters (with the `--index_bleed` arguments in `amptk filter`)
