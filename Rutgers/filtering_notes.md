# Supplementary notes on filtering
Filtering is an often neglected portion of amplicon analyses, despite the well documented occurrence of amplicon artifacts which can lead to inflation of overall richness and diversity of OTUs perceived across a dataset. There is no one way to filter. What follows is a series of steps taken to find a set of empirically derived filters which can be applied to our data. Code applied is documented herein, while data tables referenced are linked through [this spreadsheet](https://docs.google.com/spreadsheets/d/1OQVuGjC5trpjTDsATPYikvr6J0B_bCLT1yoJc7i8NzQ/edit#gid=0). Note that there are multiple data tables, listed sequentially as S1, S2, etc., as tabs on the single spreadsheet document.

## Manipulating normalized data:

We'll start with the default parameters established by `amptk` - normalizing data and using the maximum value for a single OTU to calculate index bleed.  

```
amptk filter \
-i ../dropd.cluster.otu_table.txt \
-f ../dropd.cluster.otus.fa \
-b mockIM4p82redo \
--delimiter csv \
--keep_mock \
--calculate in \
--mc /mnt/lustre/macmaneslab/devon/guano/mockFastas/CFMR_insect_mock4alt.fa \
--debug \
--threshold max \
-o max_otu \
--normalize y
```

The output generates the following summary statistics:  

```
Index bleed, samples into mock: 4.798870%.
Auto subtract filter set to 161526
mockIM4p82redo sample has 8 OTUS out of 24 expected; 0 mock variants; 0 mock chimeras; Error rate: 0.000%
Filtering OTU table down to 13 OTUs and 3,669,013 read counts
```

This is a huge index bleed rate and suggests something about our parameters are wrong, or something about or sequence data itself is very problematic. To analyze the output files from this filtering process the following commands were executed:

```
sed -i 's/#OTU ID/OTUid/' filtMax.final.csv
sed -i 's/#OTU ID/OTUid/' filtMax.normalized.num.csv
awk -F ',' '{print NF; exit}' filtMax.final.csv
cut filtMax.normalized.num.csv -d ',' -f 1,2,69 | sort -t ',' -k2,2nr | awk -F "," '$2 != "0.0" {print $0}'
cut filtMax.final.csv -d ',' -f 1,2,69 | sort -t ',' -k2,2nr | awk -F  "," '$2 != "0" {print $0}'
```

 We see that OTU228 is the culprit. It contains 4748.0 normalized reads. The lowest mock OTU contains 714 reads, while all other non-mock OTUs are <= 13 reads (most are < 2). Thus the very large index bleed filter was applied to just one curious OTU. We'll want to check to determine how many other samples contain this same OTU:  

```
grep "\\bOTU228\\b" filtMax.normalized.num.csv
```
The output suggests that several samples contain this OTU in high amounts, making it a likely candidate for either index-bleed or contamination. However, because the mock community was not amplified with these samples, contamination via extraction or PCR is impossible. Thus index bleed is the only culprit. To determine the likely taxa contributing to this contamination:

```
grep "\\bOTU228\\b" filtMax.filtered.otus.fa -A 1
```

Using NCBI BLAST, this OTU aligns with 98% identity (and full query length) to _Harmonia axyridis_, a member of the mock community. This demonstrates that the OTU228 is potentially just a misassigned read which wasn't included in the mock community fasta file which all OTUs are aligned against when calculating index bleed. Index bleed _into_ the mock community is very small among all other OTUs. However, there are two complicating factors at work:  

  1. This same OTU is present in large proportions in several other samples; because these data were normalized, it may be that there is significant index bleed from mock into true samples, but it is also possible that this proportion is greatly inflated because of the normalization process used to estimate these values. We'll explore the same dataset _without_ normalizing reads next.  
  2. This same OTU was a model organism used in the lab in which PCR (but not extraction) was performed. Thus the somewhat bimodal distribution of reads for this OTU (several high, several very low) in true samples (and negative controls) could likely be explained by potential contamination at the PCR step.

## Manipulating NON normalized data:

The same filtering approach was used as described above, with the `--normalize` flag switched from `y` to `n`.  This generated the following summary statistics:  

```
Index bleed, samples into mock: 4.795240%.
Auto subtract filter set to 4748
mockIM4p82redo sample has 8 OTUS out of 24 expected; 0 mock variants; 0 mock chimeras; Error rate: 0.000%
Filtering OTU table down to 133 OTUs and 6,391,405 read counts
```

Interestingly, using non-normalized data doesn't significantly reduce the estimation of index-bleed. We next evaluate whether a single OTU is causing this high index bleed as before:  

```
sed -i 's/#OTU ID/OTUid/' filtMaxNoNorm.final.csv
sed -i 's/#OTU ID/OTUid/' filtMaxNoNorm.sorted.csv
awk -F ',' '{print NF; exit}' filtMaxNoNorm.final.csv
cut filtMaxNoNorm.normalized.csv -d ',' -f 1,2,69 | sort -t ',' -k2,2nr | awk -F "," '$2 != "0.0" {print $0}'
cut filtMaxNoNorm.final.csv -d ',' -f 1,2,69 | sort -t ',' -k2,2nr | awk -F  "," '$2 != "0" {print $0}'
```

Indeed the same OTU is identified as with normalized data (which is expected). Because data isn't normalized, the lowest expected mock OTU has a high read depth (24291 reads), while next highest non-mock contains just 433 reads. This confirms that the proportion of index bleed outside of the single OTU228 is minimal; if 433 represents the floor, and the mock community contained 3401801 total reads, then index bleed is 1/100th of the original calculation. The question remains as to whether or not OTU228 is a result of contamination or index bleed. To determine how the proportion of normalized vs. non-nomalized reads compare, the same code was executed on this output file:  

```
grep "\\bOTU228\\b" filtMaxNoNorm.sorted.csv
```

Information from this table was combined with the original data table and can [be viewed on Table S1 here](https://docs.google.com/spreadsheets/d/1OQVuGjC5trpjTDsATPYikvr6J0B_bCLT1yoJc7i8NzQ/edit#gid=0). Unfortunately it doesn't immediately provide any definitive answer as to whether OTU228 is present in a sample due to index bleed or contamination. There are instances in which relatively low read depth contains high proportion of its overall reads to this one sample (which could be explained either by contamination or index bleed), while other samples with low read depth don't contain any reads for this OTU. When looking at non-normalized data, what is clear is the majority of samples don't contain significant amounts of this OTU - 44 of 68 samples have less than 5 reads each.

## Manipulating data with modified mock fasta

Filtering parameters were same as described in the initial script. The difference is I added OTU228's sequence into the fasta file; this will effectively filter out that Harmonia read and recalculate index bleed _into_ mock.

Output is identical in terms of reads attributed to all non-mock OTUs (nothing changes in terms of read numbers or OTUs; all that's changed is the previous OTU288 is now identified as a mock sample). However, the summary statistics portray a different index bleed estimation:  

```
Index bleed, samples into mock: 0.047002%
mockIM4p82redo sample has 35 OTUS out of 25 expected; 0 mock variants; 0 mock chimeras; Error rate: 0.039%  
Filtering OTU table down to 932 OTUs and 9,330,163 read counts
```

This demonstrates that if you assume this OTU is part of the mock, then the index bleed _into_ mock is reduced to 1/100th of the original estimation - **this suggests that whether this OTU is a source of contamination or index bleed, we are better off removing this OTU by either assigning it to the mock community or ignoring it when caculating index bleed**.

I think our best approach is to remove any reads attributed to Harmonia once we've assigned taxonomy to our dataset, but it's worth mentioning in our analysis that it is possible these could be potential dietary sources in some instance. I looked into whether samples with significant numbers of Harmonia reads were associated with some sort of similar barcode; this doesn't appear to be the case. Of the 24 samples with >= 100 Harmonia reads (per sample), there were 10 forward and 10 reverse barcodes used. It is likely that similar barcodes are not the culprit for misidentification; however, index bleed is still as possibility to be explored such that the mock community Harmonia reads were bleeding into the true samples. We'll calculate that next.

For what it's worth, the non-normalized data were also compared with this modified fasta file, with the `--normalize` flag switched from `y` to `n`. As with the normalized data, the output looks identical to the non-normalized data in terms of raw reads assigned to OTUs that should/shouldn't be in the mock.

```
Index bleed, samples into mock: 0.050620%.
mockIM4p82redo sample has 67 OTUS out of 25 expected; 0 mock variants; 0 mock chimeras; Error rate: 0.041%
Filtering OTU table down to 953 OTUs and 9,324,140 read counts
```

In summary, whether or not we normalize data will not affect our calculation for index bleed **assuming we account for the OTU228** issue described above (there are multiple ways to deal with it). Index bleed with this single OTU accounted for reduces our estimation very close to the `amptk` default of 0.005.

## Calculating bleed from mock into true samples

Because the mock community had such a high proportion of reads relative to the overall library, the rate of index bleed from the mock _into_ the samples should be addressed. `amptk` caclulates this by default and we'll pass the same argument as above but modify just one item (turning altering the `--calculate` function to `all` from `in`):

Code for normalized estimation was as follows:

```
amptk filter \
-i ../dropd.cluster.otu_table.txt \
-f ../dropd.cluster.otus.fa \
-b mockIM4p82redo \
--delimiter csv \
--keep_mock \
--calculate all \
--mc /mnt/lustre/macmaneslab/devon/guano/mockFastas/CFMR_insect_mock4alt.fa \
--debug \
--threshold max \
-o mockIn \
--normalize y
```

The normalized values show a massive 95% index bleed of mock into true samples.

```
Index bleed, mock into samples: 95.355336%.  Index bleed, samples into mock: 0.047002%.
mockIM4p82redo sample has 6 OTUS out of 25 expected; 1 mock variants; 4 mock chimeras; Error rate: 16.429%
Filtering OTU table down to 932 OTUs and 3,374,623 read counts
```

This calculation may be largely due to the normalization of reads and having a disproportionally small number of actual reads in some true samples relative to the mock sample. We'll run a non-normalized estimate to compare. As before for non-normalized estimation `--normalize` was altered from `y` to `n`. We observe a huge difference:

```
Index bleed, mock into samples: 8.307496%.  Index bleed, samples into mock: 0.050620%.
mockIM4p82redo sample has 45 OTUS out of 25 expected; 5 mock variants; 14 mock chimeras; Error rate: 0.040%
Filtering OTU table down to 953 OTUs and 8,839,573 read counts
```

In this case, using non-normalized data demonstrates a dramatically different estimation of index bleed from the mock into the true samples: just ~8.4%. And that likely represents the extreme end. To determine which OTUs are contributing to the highest proportion of index bleed (we'd assume the big contributors would be the mock sequences given the higher throughput of mock sequences than any true reads) the analysis was carried out using a bit of R code. The output file `mockIn_noNorm.sorted.csv` was used for this analysis.  

```
## read in data
reads.df <- read_csv("~/Desktop/guano/Rutgers/mockIn_noNorm.sorted.csv")

## create matrix from data.frame, creating row.names from first column of data.frame:
reads.mat <- as.matrix(reads.df[,-1])
namelist <- reads.df$otuID
rownames(reads.mat) <- namelist

## create output
maxRead.df <- data.frame(row.names = colnames(reads.mat),
                MaxVal = apply(reads.mat, 2, max),
                WhichMax = apply(reads.mat, 2, which.max))

## this creates a "WhichMax" value which is the row number, not row name. Working a bit more to fix that:
maxRead.df$SampleID <- rownames(maxRead.df)
counter = (1:953)
swap.df <- data.frame(counter, namelist)
colnames(swap.df) <- c("WhichMax", "OTUid")

Final.df <- merge(maxRead.df, swap.df)
Final.df$WhichMax <- NULL
```

The goal of this analysis was to search through each sample and identify the OTU which contained the greatest number of reads. The output lists the OTU identity, the number of reads, and the sample name. The output can be viewed on [tableS2](https://docs.google.com/spreadsheets/d/1OQVuGjC5trpjTDsATPYikvr6J0B_bCLT1yoJc7i8NzQ/edit#gid=0).  

The raw output from the R script only generates three fields. Two additional fields were created and are explained below. Note that the **OTUid** field indicates extensive numbers of suspected mock chimeras - this is a bug in the software which prints this value out when performing `--calculate all` in the `amptk filter` script _and you keep your `--keep mock` flag_. For further details about the naming convention Jon applies, see lines 185-189 in [his python script](https://github.com/nextgenusfs/amptk/blob/master/bin/amptk-filter.py).  

The data is structured as follows:  
- **SampleID** represents the uniquely sequenced guano sample  
- **TotalNonNormReads** represent the total number of reads passing our initial `amptk illumina` and `amptk cluster` filters, summed on a per-sample basis  
- **MaxVal** is identified by finding the OTU which contains the most number of reads on a per sample basis  
- **percReads** represents the proportion of reads represented by the **MaxVal** value for a given sample. In other words, this identifies how much one OTU is contributing to the overall proportion of reads in a sample  
- **OTUid** represents an unclassified (but unique) sequence as determined by the `amptk cluster` scripts  

For example, the first observation occurs for the mock community sample, which identifies a Mock OTU as the OTU with the most reads associated to that sample (which is 223034 reads) - that's exactly what we'd expect (the mock community should contain only mock sequences). The next observation is that of sample `rut16-1382`; the OTU for this sample which has the greatest number of reads is OTU2, which contains 21124 reads total.  

Highlighted **OTUid** fields are those cases which are of concern, and suggest that in some, but not all cases, we have mock community sequences producing significant index bleed into the overall proportion of our samples. 19 of 68 of these samples have a mock community member contributing to the highest proportion of their overall reads - this is a likely a direct consequence of the high proportion of mock community sequences relative to sample sequences. The only way to account for these sequences are to remove them - this eliminates the false positive situation, but importantly, does nothing to address the fact that these may be samples in which that represents a real dietary component. Notably, these sequences are **not removed** from our filtering process, but remain in place even after taxonomy is assigned. However, when examining dietary trends among our samples, these mock samples will be filtered out.  

What is just as revealing is that 2/3 of the samples don't contain a mock OTU as the dominant sequence type. This is encouraging because it suggests that if index bleed is happening, it's happening randomly and is not universally extensive (that is, even if every sample had some amount of mock OTU sequence in it, it's not in the majority). What is impossible to calculate is the extent with which a true sample is bleeding into another true sample.  

## Summary observations
The proportion of index bleed _into_ the mock community is trivial - less than 0.05%. This is expected given the large proportion of reads in the mock community relative to any other single sample in the library. However, the proportion of mock community reads _into a non-mock sample_ varies; in about 1/3 of the overall instances the mock sequence represents the highest proportion of reads, yet this generally represents less than 20% of the overall number of total reads in a given sample. This would appear to set the upper ceiling of index bleed at about 20% (an enormous value which would drastically reduce our dataset).  

The tradeoff is as follows:  
- **Keep majority of data, lose specific OTUS**. We can greatly reduce the observed index bleed if we eliminate mock sequences entirely from our analysis. However, the only rate with which we can empirically derive a value for index bleed if doing so would be from the `--calculate in` function which estimated a low index bleed rate of about 0.05%. The result in doing so would eliminate our potential to include the mock OTUs in any true sample, even if it were true that bats are eating one of those ~20 species. The benefit is that the vast majority of OTUs are retained; I would likely apply a slightly more conservative index bleed of about 0.1% (doubling the empirically derived observation).
- **Keep all OTUs, lose majority of data**. This would drastically cut down the total number of reads remaining in our dataset, but would retain instances in which mock OTUs are identified in non-mock samples. This does not seem like the most prudent option because you are retaining potential false positives (which you can't empirically differentiate from true signal short of resequencing or targeting via qPCR), and you lose the majority of the overall read numbers in your data set.  

## (almost) Final filtering protocol  
I am going to apply an index-bleed filter of 1% across all reads, drop the mock community sample from our analysis, and carry forward with applying taxonomic information to the OTUs in the dataset. In addition, I will use the modified mock fasta file which added OTU228 as an additional Harmonia sequence.  

This is performed with the following code:  

```
amptk filter \
-i ../dropd.cluster.otu_table.txt \
-f ../dropd.cluster.otus.fa \
-b mockIM4p82redo \
--delimiter csv \
--index_bleed 0.01 \
--mc /mnt/lustre/macmaneslab/devon/guano/mockFastas/CFMR_insect_mock4alt.fa \
--threshold max \
-o noNorm_bleed_1.0 \
--normalize n
```

We observe the following summary statistics:  

```
Index bleed, mock into samples: 8.307496%.  Index bleed, samples into mock: 0.050620%.
mockIM4p82redo sample has 56 OTUS out of 25 expected; 10 mock variants; 20 mock chimeras; Error rate: 0.041%
Filtering OTU table down to 941 OTUs and 5,677,849 read counts
```

See **tableS3** for this output - a table of non-normalized number of reads per OTU per sample. Of note, this table suggests that among the OTUs detected in the mock sample the top 3 OTUs (in terms of number of reads) are **216, 29, and 11**. All other OTUs identified in the mock sample contain less than 10 reads. This is suggests that we may set a threshold for inclusion in our binary "present" or "absent" data table for future analyses in which somewhere between 10-217 reads are required to be considered "present", rather than just 1 or 2. This is ultimately the last major filtering consideration - what number of reads should any OTU, per sample, have _at a minimum_?  

This _read minimum threshold_ can be applied in the same `amptk filter` script. When we append the above final filtering script to include `--min_reads_otu` and a value (either 216 or 29, for example) we see the following summary outputs:  

** for a minimum number of reads > 29 **

```
mockIM4p82redo sample has 35 OTUS out of 25 expected; 4 mock variants; 7 mock chimeras; Error rate: 0.038%
Filtering OTU table down to 703 OTUs and 5,673,534 read counts
```

** for a minimum number of reads > 216 **

```
mockIM4p82redo sample has 25 OTUS out of 25 expected; 1 mock variants; 0 mock chimeras; Error rate: 0.036%
Filtering OTU table down to 442 OTUs and 5,626,434 read counts
```

We find that there are no unintended OTUs present in the mock community and have thus resolved any traces of index bleed _into the mock sample_ with the more stringend **217** minimum read threshold. As expected, the additional filtering parameter of requiring a read minimum of 217 reads per sample don't significantly reduce the total number of reads from the library (it's a net loss of less than 1%). What is significant is the number of OTUs which are discarded - we lose more than half in the case of a minimum read depth set to **217**, but lose only about a quarter of reads when set to **30**. To be most conservative I'm inclined to enact a stricter threshold for read inclusion and will set the minimum number of reads to that higher value to ensure we have removed all trace of index bleed into our mock community.  


The final filtering command is as follows:  

```
amptk filter \
-i ../dropd.cluster.otu_table.txt \
-f ../dropd.cluster.otus.fa \
-b mockIM4p82redo \
--mc /mnt/lustre/macmaneslab/devon/guano/mockFastas/CFMR_insect_mock4alt.fa \
--index_bleed 0.01 \
--threshold max \
--subtract 217 \
-o fullFilt \
--delimiter csv \
--normalize n
```

One final note: the filtering threshold is applied to the dataset which contained all samples. Because some samples contained less than the minimum read number on a per-OTU basis
