# Supplementary notes on filtering
Filtering is an often neglected portion of amplicon analyses, despite the well documented occurrence of amplicon artifacts which can lead to inflation of overall richness and diversity of OTUs perceived across a dataset. There is no one way to filter. What follows is a series of steps taken to find a set of empirically derived filters which can be applied to our data.  

# Filtering the `rough` dataset
We'll start with the default parameters established by `amptk` - normalizing data and using the maximum value for a single OTU to calculate index bleed. We'll carry though the entire filtering analysis for just the `rough` dataset first; once we've applied our necessary filters we may then drop samples, drop OTUs, and/or reduce numbers of sequences per sample. The consequence of these actions are then applied to our initial (`rough`) .demux.fq.gz dataset, and the remaining reads are reclustered. We finish by running a secondary filtering analysis to ensure that the intended actions we're going to define herein are appropriately applied to the expected final dataset that ultimately is assgined taxonomic identities.  

## Manipulating normalized data:
The following example shows an example command used for the `rough` dataset.  

> Note a new directory `filtd` was created prior to the execution of this code to retain output files. A subdirectory `rough` was also created and incorporates these initial analyses.

```
amptk filter \
-i /mnt/lustre/macmaneslab/devon/guano/Pompton/clust/rough/rough.cluster.otu_table.txt \
-f /mnt/lustre/macmaneslab/devon/guano/Pompton/clust/rough/rough.cluster.otus.fa \
-b mock-IM4p11-2 \
--delimiter csv \
--keep_mock \
--calculate in \
--subtract auto \
--mc /mnt/lustre/macmaneslab/devon/guano/mockFastas/CFMR_insect_mock4alt.fa \
--debug \
--threshold max \
-o mockIn \
--normalize y
```

The output generates the following summary statistics:  

```
OTU table contains 103 samples, 2,149 OTUs, and 2,595,284 reads counts
Index bleed, samples into mock: 6.376809%
Auto subtract filter set to 1747
mock-IM4p11-2 sample has 19 OTUS out of 25 expected; 0 mock variants; 0 mock chimeras; Error rate: 0.045%
Filtered OTU table contains 103 samples, 366 OTUs, and 2,209,160 read counts
```
The parameters above are identifying that:  
- A fairly high fraction of _normalized reads_ true samples are moving into the mock community (`Index bleed, samples into mock`); this is a high value which may be indicative of samples with low absolute read numbers being artificially inflated and influencing this calculation - this can be avoided by dropping the samples in the first place, but may also be investigated by looking at what happens to the `index-bleed` calculation when data isn't normalized (that's done next).  
- The `Auto subtract filter` identifies the number necessary to remove all unexpected reads in the mock sample. The value is extremely high, and in past experience generally indicates that a mock community sample is being misclassified. We'll look for that in a moment. We generally are looking for subtract values to be in the range of less than 50 or so.  
- The filtere OTU table resulting after this process drastically shrinks the number of OTUs, but doesn't really reduce the number of sequences nearly as much; this is often the case when applying the index-bleed filter, as it generally removes OTUs with very low sequence depth.  

The following commands are used to investigate what OTUs might be finding their way from true samples into the mock community, as well as finding how many reads (and from which OTUs) from the mock community are bleeding into the true samples - this uses the output of the filtered sample above.

```
sed -i 's/#OTU ID/OTUid/' mockIn.final.csv
sed -i 's/#OTU ID/OTUid/' mockIn.normalized.num.csv
awk -F ',' '{print NF; exit}' mockIn.final.csv
  ## there are 104 rows, thus we have 103 total samples - the mock community is the second line
cut mockIn.normalized.num.csv -d ',' -f 1,2,3 | sort -t ',' -k2,2nr | awk -F "," '$3 != "0.0" {print $0}'
```

Interpreting results from normalized data in which samples with low read abundances were included need to be taken with a grain of salt; this preliminary investigation is more concerned with what OTUs are found in the mock community, not necessarily their read abundances. We find that most of our mock community members are identified, and most reads are fairly high. Notably, other samples have higher proportions of mock reads in some instances - concerning, but understandable because there is a high likelihood these values are inflated from the normalization process. We'll compare these values again in a non-normalized analysis in a moment.  

Strikingly, one OTU was in much higher proportions than all non-mock members _within the positive control itself_: `OTU33`. To identify the sequences associated with each OTU of interest:  
```
grep "\\bOTU33\\b" mockIn.otus.counts.fa -A 1
```

Similarly, if you wanted to review the specific numbers of reads within each of those particular OTUs:
```
grep "\\bOTU33\\b" mockIn.normalized.num.csv
```

The output from the BLAST search is revealing.
- `OTU33` is a 100% identity match (with 100% coverage) for _Cellulosimicrobium funkei_, a gram positive bacterium; this sequence is possibly the result of our data being spiked in with another dataset (we only used about 15% of the lane so the remaining sequences are from other projects unaffiliated with this dataset) - if there were other samples included on the same lane that were present, especially if these were whole-genome sequencing projects, or say a 16-S microbiome project, it's possible this OTU was contaminated as part of that run. What's particularly strange is the distribution of this OTU: it's in everything. An additional BLAST-search limited to entries limited to just chordate references revealed no significant alignment matches; this OTU is most likely a good candidate to drop from the dataset entirely.  

Additionally high sequence depth samples within the mock community were familiar known contaminants. These need to be evaluated in the non-normalized dataset (next), but a few are worth noting for compairson between the normalized and non-noramlized analyses:  
- `OTU137` contains 100% coverage and 99% identity to a few species, but one is a known contaminant in other datasets: _Maccaffertium mediopunctatum_. It may be that this OTU is a true signature, but
- `OTU130` contains 100% coverage and 100% identity for _Chauilodes_, a mayfly know to be persistent in our datasets as a low-level contaminant.  


## Manipulating NON normalized data:

The same filtering approach was used as described above, with the `--normalize` flag switched from `y` to `n`. Output files were renamed to specify the fact we turned normalization off. I left **on** the `--subtract auto` option for the time being so we could make a direct comparison to the last filtering command:  

```
amptk filter \
-i /mnt/lustre/macmaneslab/devon/guano/Pompton/clust/rough/rough.cluster.otu_table.txt \
-f /mnt/lustre/macmaneslab/devon/guano/Pompton/clust/rough/rough.cluster.otus.fa \
-b mock-IM4p11-2 \
--delimiter csv \
--keep_mock \
--calculate in \
--subtract auto \
--mc /mnt/lustre/macmaneslab/devon/guano/mockFastas/CFMR_insect_mock4alt.fa \
--debug \
--threshold max \
-o mockIn \
--normalize n
```

By _not normalizing_, we see virtually no change in the overall index bleed rate, but the `subtract` filter value is reduced substantially (from 1747 top 601); this reflects the fact that the normalization step was inflating our reads contributing to the index-bleed filter. Additional details from the .log file are as follows:  

```
Index bleed, samples into mock: 6.375371%.
Auto subtract filter set to 601
mock-IM4p11-2 sample has 19 OTUS out of 25 expected; 0 mock variants; 0 mock chimeras; Error rate: 0.045%
Filtered OTU table contains 72 samples, 207 OTUs, and 2,118,635 read counts
```

Pretty interesting. We lose a bunch of samples - almost the exact amount that we would have dropped initially if we had just kept that minimum read threshold of about 5000 reads (that would have retained 76 samples total). Using non-normalized data significantly reduced the `subtract` value, but we need to further determine what OTUs are contributing to such a high value (this value is an extremely high value).  

```
sed -i 's/#OTU ID/OTUid/' mockIn.sorted.csv
cut mockIn.sorted.csv -d ',' -f 1,2,3,101,102,103,104 | sort -t ',' -k2,2nr | awk -F "," '$2 != "0" {print $0}'
```

What do we see? The top hitter is indeed the same OTU as before - `OTU33`, which an earlier BLAST search determined was that unexpected bacteria, _Cellulosimicrobium funkei_. We also notice that the (proportionally) high unexpected reads occur in the negative control samples, suggesting these OTUs are likely contaminants. To get a better sense of how these negative control-associated OTUs are distributed in our own dataset, we'll run `amptk filter` again, but add an additoinal argument `--negatives` as shown in the scrip below:
```
amptk filter \
-i /mnt/lustre/macmaneslab/devon/guano/Pompton/clust/rough/rough.cluster.otu_table.txt \
-f /mnt/lustre/macmaneslab/devon/guano/Pompton/clust/rough/rough.cluster.otus.fa \
-b mock-IM4p11-2 \
--delimiter csv \
--keep_mock \
--calculate in \
--subtract auto \
--mc /mnt/lustre/macmaneslab/devon/guano/mockFastas/CFMR_insect_mock4alt.fa \
--debug \
--threshold max \
-o mockIn \
--normalize n \
--negatives pomp-NTCa01 pomp-NTCb02 pomp-NTCc03 pomp-NTCd04
```

The .log file indicates that two of the four negative control samples did not contain enough reads post filtering to be used in the secondary filtering process (mapping negative control-associated OTUs against the entire dataset's OTU matrix); however the other two negative controls contained 188 OTUs which require further investigation. Rather than printing a list of each of these, let's look at the top reads associated with those two negative control samples: `pomp-NTCb02` and `pomp-NTCd04`:

```
sed -i 's/#OTU ID/OTUid/' mockIn.sorted.csv
cut mockIn.sorted.csv -d ',' -f 1,2,3,101,102,103,104 | sort -t ',' -k7,7nr
```


grep "\\bOTU435\\b" mockIn.sorted.csv

grep "\\bOTU424\\b" mockIn.otus.counts.fa -A 1


### off script R commands
```
setwd("~/Repos/guano/BRIpompton/data/amptk/")
## load in *.sorted.csv from `amptk filter` output:
df <- read.csv(file = "mockIn.sorted.csv")

## pull out only negative samples (my samples always have the prefix "NTC"):
ntc.df <- df[grepl("NTC", names(df))]
row.names(ntc.df) <- df$OTUid   ## remember what the OTUids are; use rowname so you can do math on data.frame

## find the highest value among NTC elements per row, then subtract that value from original df:

## first, find the maximum value in any NTC sample, per OTU:
ntcRowmax <- data.frame(apply(ntc.df, 1, max))
colnames(ntcRowmax) <- "counts"

## drop the OTUid column and keep as rowname
df1 <- df[,-1]  
rownames(df1) <- df$OTUid

## subtract the max value of reads in any NTC sample observed among for each OTU across all samples in the original matrix
filt.df <- sweep(df1,1,ntcRowmax$counts,"-")  ## for 'sweep' details see: https://bioinfomagician.wordpress.com/2014/08/12/my-favorite-commands-part3-sweep-function-in-r/
## and reduce values less than 0 to 0:
filt.df[filt.df<0] <- 0
```

Alternative filtering approach:
## drop the bacterial OTU (OTU33)
```
amptk drop \
-i /mnt/lustre/macmaneslab/devon/guano/Pompton/clust/rough/rough.cluster.otus.fa \
-r /mnt/lustre/macmaneslab/devon/guano/Pompton/illumina/trim_pomp.demux.fq \
-l OTU33 \
-o drop
```

Use the `drop.cleaned.otu_table.txt` and `drop.cleaned.otus.fa` for next stages.

Brief summary (so we've dropped just one OTU):
```
Loading 2152 OTUs
Dropping 1 OTUs
2,151 OTUs remaining
3,043,388 reads mapped to OTUs (100%)
```

## empty amptk
1. Run Jon's `amptk filter` program with empty parameters so that you can rename the OTUs that align to the mock community:

```
amptk filter \
-i /mnt/lustre/macmaneslab/devon/guano/Pompton/drop/drop.cleaned.otu_table.txt \
-f /mnt/lustre/macmaneslab/devon/guano/Pompton/drop/drop.cleaned.otus.fa \
-b mock.IM4p11.2 \
--delimiter csv \
--keep_mock \
--index_bleed 0 \
--subtract 0 \
--mc /mnt/lustre/macmaneslab/devon/guano/mockFastas/CFMR_insect_mock4alt.fa \
--debug \
-o nofilt \
--normalize n
```

A quick look at the .log file indicates that we have some filtering to do and confirms we haven't actually removed any reads:
```
mock-IM4p11-2 sample has 79 OTUS out of 25 expected; 4 mock variants; 41 mock chimeras; Error rate: 1.304%
Filtered OTU table contains 103 samples, 2,148 OTUs, and 2,551,644 read counts
```

Now run the R script to filter out reads according to max values observed per OTU among any negative control (NTC) sample.

## Run the NTCfilter R script
Use the `nofilt.sorted.csv` file and run the NTCfilter.R script to subtract reads by identifying the highest number of reads per negative control per OTU from all true samples.
```
rsync devon@premise.sr.unh.edu:/mnt/lustre/macmaneslab/devon/guano/Pompton/filt/rough/noFilt/noFilt.sorted.csv .
```

Run the R script:
```
setwd("~/Repos/guano/BRIpompton/data/amptk/")
df <- read.csv(file = "rough.cluster.otu_table.txt", sep = "\t")
ntc.df <- df[grepl("NTC", names(df))]
row.names(ntc.df) <- df$OTUid   ## remember what the OTUids are; use rowname so you can do math on data.frame
ntcRowmax <- data.frame(apply(ntc.df, 1, max))
colnames(ntcRowmax) <- "counts"
df1 <- df[,-1]
rownames(df1) <- df$OTUid
filt.df <- sweep(df1,1,ntcRowmax$counts,"-")  ## for 'sweep' details see: https://bioinfomagician.wordpress.com/2014/08/12/my-favorite-commands-part3-sweep-function-in-r/
filt.df[filt.df<0] <- 0
filt.df <- sweep(filt.df,1,filt.df$mock.IM4p11.2,"-")
filt.df[filt.df<0] <- 0   # set neg values back to zero again...
filt.df$OTUid <- df$X.OTU.ID  # adding back in the OTU column
filt.df <- filt.df[,c(104,1:103)] # positioning OTUid column back to first position
setwd("~/Repos/guano/BRIpompton/data/amptk/")
write.table(filt.df, file = "NTCreduced.otu_table.txt", quote = F, row.names = FALSE, sep = '\t')
```

Take the `NTCreduced.otu_table.txt` output from the R script and modify the header for the first field:
```
sed -i 's/OTUid/#OTU ID/' NTCreduced.otu_table.txt
```

Now plug back into the `amptk filter` program with the NTC removed table.

```
amptk filter \
-i /mnt/lustre/macmaneslab/devon/guano/Pompton/filt/rough/NTCfiltd/NTCreduced.otu_table.txt \
-f /mnt/lustre/macmaneslab/devon/guano/Pompton/drop/drop.cleaned.otus.fa \
-b mock.IM4p11.2 \
--delimiter csv \
--keep_mock \
--calculate in \
--subtract auto \
--mc /mnt/lustre/macmaneslab/devon/guano/mockFastas/CFMR_insect_mock4alt.fa \
--debug \
--threshold max \
-o NTCfiltd \
--normalize n
```

Review the `mockIn.sorted.csv` file and check out what OTUs remain in the mock:
```
sed -i 's/#OTU ID/OTUid/' NTCfiltd.sorted.csv
cut NTCfiltd.sorted.csv -d ',' -f 1,2,3 | sort -t ',' -k2,2nr | awk -F "," '$2 != "0" {print $0}'
```

Weirdly there is an OTU that's very closely aligning to a mock member, MOCKIM23_pident=98.3_OTU25 which seems to be triggering a really high subtract value (**1,741**). Let's look at each of these suspected members with high read numbers and see if they're just really mock community variants which we can cull in attempts at reducing the `subtract` rate:

For example:
```
grep "\\bOTU25\\b" NTCfiltd.otus.counts.fa -A 1
```

We find that OTU25 is a _Harmonia_ variant. Good to toss.
OTU179 is a likely chimera.
OTU156 is a likely contaminant (Caenis sp.)
OTU518 is a likely contaminant from (Maccaffertium)
OTU517 is a contaminant - a bat! (Lasiurus borealis)
OTU431 is a likely contaminant (Limoniidae)
... then it get's less clear from BLaST results alone ... investigate distributions of reads too
OTU532 is unclear - it's a click beetle and probably a contaminant
OTU790 is a cranefly - unclear if it's a contaminant
... then clear again ...
OTU137 was a previously suspected contaminant, a Maccaffertium/Stenonema bug
... then unclear ...
OTU819 is a tiger crane fly... possible contaminant but unclear
... then clear again...
OTU867 is a contaminant as we've seen it in other data sets... Dendroides canadensis is a click beetle

...

next step is to look at the distributions of these suspected contaminant reads...
See the R script `mockfilter.R`, which:
1. finds any OTUs with > 0 reads in the OTU table that aren't mock community members
2. finds the max value among any sample for that OTU
3. find the total number of samples that contain > 0 reads for a given OTU

It looks pretty clear that if an OTU is the result of index bleed from a true sample into a mock community, you get a  **max** value in the top-hitting true sample that is much higher than the mock sample. For example, a mock sample for `OTUx` might have 5 reads, but the highest number of reads in a true sample for `OTUx` might be something like _**40,000**_ reads. We can apply a `--subtract` filter that is low in this instance to appropriately filter out these low-level reads in our mock.  

It's not entirely clear what the signature for the OTUs that resulting from contamination look like (be it from PCR or sequencing chimeras, or from actual contamination during DNA extraction or PCR). We'll use a secondary filter in an R script later which will look to match OTUs by similar BOLDid's.  

## Summary of actions
It appears that among the mock OTUs not expected to be in the sample, those reads with between 50-100 reads are mostly suspected contaminants. We'll filter them out with another R script that compares BOLD identifiers between unrelated runs. Most likely then our index-bleed rate is rather low. We'll finsh up the `amptk filter` script by applying a generic (and somewhat conservative ) 2% index-bleed filter, and one additional parameter: any OTU only detected in a single sample will be dropped (so no singleton OTUs). This helps a lot in the diversity estimates and dissimilarity calculations downstream. The rationale is that if the OTU is only detected in one of about a hundred samples, it is likely not relevant. Note this may reduce our ability to detect the entire breadth of OTUs, but it's far more likely to just cut out chimeric sequences and other spurious OTUs.

```
amptk filter \
-i /mnt/lustre/macmaneslab/devon/guano/Pompton/filt/rough/NTCfiltd/NTCreduced.otu_table.txt \
-f /mnt/lustre/macmaneslab/devon/guano/Pompton/drop/drop.cleaned.otus.fa \
-b mock.IM4p11.2 \
--delimiter csv \
--index_bleed 0.02 \
--mc /mnt/lustre/macmaneslab/devon/guano/mockFastas/CFMR_insect_mock4alt.fa \
-o filtd \
--normalize n \
--min_samples_otu 2
```

We've dropped the four negative controls and eliminated our mock community sample from further analysis. All true samples remain for analysis. A summary of the output is as follows:  

```
Overwriting auto detect index-bleed, setting to 2.000000%
mock.IM4p11.2 sample has 60 OTUS out of 25 expected; 2 mock variants; 28 mock chimeras; Error rate: 1.107%
Dropped 1,354 OTUs found in fewer than 2 samples
Filtered OTU table contains 98 samples, 771 OTUs, and 1,353,398 read counts
Samples dropped: mock.IM4p11.2,pomp.NTCa01,pomp.NTCb02,pomp.NTCc03,pomp.NTCd04
```

That singleton OTU threshold is significant - we've dropped over a thousand OTUs because of their singleton status. It may be worthwhile investigating whether other clustering programs like LULU could potentially retain these singelton OTUs by identifying what other OTU they are derived from, but it's my feeling at the moment that there aren't enough examples clarifying exactly how well that program produces false positives.  

The resulting OTU table is then use to apply taxonomy next, as described in the [Pompton workflow](https://github.com/devonorourke/guano/blob/master/BRIpompton/docs/Pompton_workflow.md) document.  
