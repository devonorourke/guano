# Supplementary notes on filtering
Filtering is an often neglected portion of amplicon analyses, despite the well documented occurrence of amplicon artifacts which can lead to inflation of overall richness and diversity of OTUs perceived across a dataset. There is no one way to filter. What follows is a series of steps taken to find a set of empirically derived filters which can be applied to our data.  

# Filtering the initial clustered datasets
We'll start with the default parameters established by `amptk` - normalizing data and using the maximum value for a single OTU to calculate index bleed. Because we have two datasets which were independently clustered, we'll need to apply the same code to both libraries of data. We'll carry though the entire filtering analysis for just the `p10-1` dataset first, then switch to the `p10-2` data in the second section.  

## Manipulating normalized data for `p10-1`:
The following example shows an example command used for the `p10-1` dataset.  

> Note a new directory `filtd` was created prior to the execution of this code to retain output files  

```
amptk filter \
-i path/to/{filename}.cluster.otu_table.txt \
-f path/to/{filename}.cluster.otus.fa \
-b une-mockIM4 \
--delimiter csv \
--keep_mock \
--calculate all \
--subtract auto \
--mc /mnt/lustre/macmaneslab/devon/guano/mockFastas/CFMR_insect_mock4alt.fa \
--debug \
--threshold max \
-o mockIn \
--normalize y
```

The `.log` file includes summary information regarding the resulting filtered reads:  
```
Index bleed, mock into samples: 70.050345%.  
Index bleed, samples into mock: 0.048002%.
Auto subtract filter set to 17
mockIM4p10L1 sample has 23 OTUS out of 25 expected; 0 mock variants; 0 mock chimeras; Error rate: 0.084%
Filtering OTU table down to 2,768 OTUs and 4,216,429 read counts
```

The parameters above are identifying that:  
- A large number of reads expected to be present in the mock community are mapping to true samples (`Index bleed, mock into samples`); this isn't surprising given that nearly 1/3 of the lane consisted of the mock community. This was due to an improperly balanced library. As such, reads that are present in the mock community should be filtered out from the real data; we should take note that nearly all our expected OTUs are in the mock; it's often the case that the missing pair of OTUs (we have 23 of 25 identified) are actually clustered together from our `amptk cluster` script because they are so similar in sequence identity)  
- A fairly small fraction of true samples are moving into the mock community (`Index bleed, samples into mock`); in addition, because we've normalized our dataset in this initial filtering analysis it's very likely that there has been limited index bleed of any true sample into our OTUs (that's a good thing).  
- The `Auto subtract filter` identifies the number necessary to remove all unexpected reads in the mock sample. The value isn't particularly high, however we see that we drop out a few _expected_ OTUs from our mock community. We'll want to next look into what the OTUs are that are causing this "subtract" filter to be activated in the first place (ie. the OTUs present in the mock sample that shouldn't be there)  
- The final OTU table and number of reads are sharply reduced after filtering the `p10-1` data; we lose about a half of our overall reads, and cut the number of OTUs by about 80%. This is because of the absurdly high mock bleed in rate, and the subtract filter working together. We're going to modify those values next.  

The following commands are used to investigate what OTUs might be finding their way from true samples into the mock community, as well as finding how many reads (and from which OTUs) from the mock community are bleeding into the true samples - this uses the output of the filtered sample above.

```
sed -i 's/#OTU ID/OTUid/' 10-1_default.sorted.csv
awk -F ',' '{print NF; exit}' 10-1_default.sorted.csv   # returns a value of 386 (the number of samples in our OTU table (385), plus the first column listing our OTUs)
cut 10-1_default.sorted.csv -d ',' -f 1,2,386 | sort -t ',' -k3,3nr | awk -F "," '$3 != "0" {print $0}'
```

It's clear from the output of the third line of the code above that our mock community looks good - here's a portion of that output:  

```
MockIM33_pident=100.0_OTU49,0,36402
MockIM49_pident=100.0_OTU56,0,28017
OTU1256_suspect_mock_variant,0,270
OTU2336_suspect_mock_variant,0,99
OTU1159_suspect_mock_chimera,0,61
OTU707_suspect_mock_chimera,0,37
OTU36;_suspect_mock_variant,1,34
OTU163_suspect_mock_chimera,0,33
OTU524_suspect_mock_chimera,0,30
OTU205_suspect_mock_variant,0,29
```
 - there are at least 28,000 reads in every mock OTU, and the next highest number of reads of an unexpected read has just 270 reads To deduce a taxonomic identity for those top potential contaminant OTUs we could also perform a [BLAST search](https://blast.ncbi.nlm.nih.gov/Blast.cgi) through NCBI; here's how to get the sequence to BLAST:  

```
grep "OTU1256" 10-1_default.otus.counts.fa -A 1
grep "OTU2336" 10-1_default.otus.counts.fa -A 1
grep "OTU1159" 10-1_default.otus.counts.fa -A 1
grep "OTU707" 10-1_default.otus.counts.fa -A 1
grep "OTU36;" 10-1_default.otus.counts.fa -A 1
grep "OTU163;" 10-1_default.otus.counts.fa -A 1
grep "OTU524" 10-1_default.otus.counts.fa -A 1
grep "OTU205;" 10-1_default.otus.counts.fa -A 1
```

Similarly, if you wanted to review the specific numbers of reads within each of those particular OTUs:
```
grep "OTU1256" 10-1_default.sorted.csv
grep "OTU2336" 10-1_default.sorted.csv
# ... repeat for all other OTUs of interest listed above
```

The output from each of the three BLAST searchers coupled with the distribution of reads per sample are revealing:  
- `OTU1256` is a 97% identity match (with 100% coverage) for _Harmonia axyridis_, a mock community member itself. However it only exists in a fraction of samples - just 7 of 385 samples (see `grep "OTU2336" 10-1_default.sorted.csv | tr ',' '\n' | uniq -c` for a little one-liner to do such a count). We'll likely chalk this up to a mutation that was propagated via PCR either during clustering (Illumina step) or during conventional PCR (with our arthropod primers). Either way, the OTU will be eliminated entirely as it's part of our expected mock community, albeit a little mutated from the expected sequence.  
- `OTU2336` is a 97% match (with 99% coverage) for some type of moth (there are a few BLAST hits with different species listed). As with the first OTU, it's only present in a tiny fraction of samples - just 2 - one in our mock, and one in a true sample (and there are just 10 reads in the true sample). This OTU will be dropped from the dataset.  
- `OTU1159` has a 100% match with 100% coverage to a certain beetle species. This has been identified in other sequencing projects and is likely a minor contaminant either in the mock community directly or in the PCR mix (as those are the only two reagents consistently shared among all libraries). Because this sample occurs _only_ in the mock community sample, it's likely either a minor contaminant in the mock community directly, or this sequence is similar enough to another true mock community member that a mutation during PCR is being identified as this other beetle. Either way, the OTU is going to be dropped from further analysis.  
- `OTU707` is a 100% match with 100% coverage to a species of moth commonly seen in other sequencing projects. Unlike the earlier contaminant OTUs, this one appears in 31 different samples, yet the highest value in any of these samples if just 58 total reads. This is thus quite likely a ver low-level contaminant possibly carried over in the PCR mix and will be eliminated as an OTU from analysis.  
- `OTU36` has several possible BLAST matches at 98% identity and 100% coverage, yet unlike the other OTUs in question here, this OTU is clearly real. Nevertheless, the proportion of reads associated with this OTU among samples varies wildly - for example, the top 10 reads among samples with this OTU (none of which are the mock) are:  47065,16208,1673,199,115,114,107,87,74,72. In other words, just a handful of samples likely contain high proportions of this OTU, but this is clearly a real OTU and not a contaminant. This may ultimately be the value in which our `--subtract` filter is set.
- `OTU136` is a common house fly; it's a contaminant seen in many sequencing runs. Unlike the previous OTU (`OTU36`), there are no samples with highly abundant read numbers, though these read depths are not inconsequential: 642,516,430,348,280,232,223,214,180,151 reads are identified in our top 10 samples (though the mock has just 33). The fact this sample has some amount of reads detected in 131 different samples suggests it is a low level-contaminant that can be removed (it it was naturally part of the diet there would be at least a few samples with elevated levels).  
- `OTU524` is detected in 49 different samples, but like `OTU707` there are no samples with even moderate read depths - the maximum value in any one sample is 58 total reads. As with other OTUs described here the alignment result has been identified in other sequencing projects as a low-level contaminant. This OTU will be dropped.  
- `OTU205` follows a similar pattern to other OTUs described here with it being frequently observed among samples but always in very low read numbers. This OTU will be dropped from further analysis.

The remaining concern is to what to do about the high index bleed of the mock community members into our true samples. As mentioned previously, this is to be expected because of the unbalanced library that was sequenced: there were far more reads dedicated to the single mock community sample than to any true sample. Thus there is a far greater likelihood of these OTUs present in the mock community to have found there way into the true samples. However, I wondered if the incidence of bleed into the mock community was actually inflated by the fact that we were normalizing our reads in the first place. Thus the next step was to run the same command filtering our reads as above, but passing setting the `--normalize` argument to `n` ("no" normalizing); this is explained in the next section.  

## Manipulating NON normalized data:

The same filtering approach was used as described above, with the `--normalize` flag switched from `y` to `n`. Output files were renamed to specify the fact we turned normalization off. I left **on** the `--subtract auto` option for the time being so we could make a direct comparison to the last filtering command:  

```
amptk filter \
-i path/to/{filename}.cluster.otu_table.txt \
-f path/to/{filename}.cluster.otus.fa \
-b une-mockIM4 \
--delimiter csv \
--keep_mock \
--calculate all \
--subtract auto \
--mc /mnt/lustre/macmaneslab/devon/guano/mockFastas/CFMR_insect_mock4alt.fa \
--debug \
--threshold max \
-o mockIn \
--normalize n
```

By _not normalizing_, we see a monumental difference in the overall index bleed rate, and this is entirely attributed to the calculation of mock reads bleeding into the true samples:  

```
Index bleed, mock into samples: 1.842435%.  
Index bleed, samples into mock: 0.051200%.
Auto subtract filter set to 270
mockIM4p10L1 sample has 24 OTUS out of 25 expected; 0 mock variants; 0 mock chimeras; Error rate: 0.080%
Filtering OTU table down to 417 OTUs and 7,480,186 read counts
```

Pretty interesting. Using non-normalized data significantly reduced the estimation of index-bleed overall, yet the subtract filter now jumped up. We know that if the `--subtract` filter is set to 270, then it's misrepresenting those mutant _Harmonia_ reads that's setting that value - we'll fix that in the next script. The main observation here is that by _not normalizing_ our data, we have a very different picture of how often mock reads are finding their way into our community: it's tiny - just 1.8 percent, which is right in line with what Jon and the default `amptk` program generally suggests (2% is the default). This suggests we shouldn't be normalizing this data, given how skewed our distribution of reads in the mock community is.  

Another takeaway: we have a lot fewer OTUs, but have retained many more reads. This is expected: a drop in OTUs should occur because we've observed an _increase_ in the **subtract** value (551 versus the previous 45); likewise, a _decrease_ in the index bleed should retain more overall reads.  The next question is to determine whether or not that **subtract** increase is warranted:  

```
sed -i 's/#OTU ID/OTUid/' 10-1_noNorm.sorted.csv
cut 10-1_noNorm.sorted.csv -d ',' -f 1,2,386 | sort -t ',' -k3,3nr | awk -F "," '$3 != "0" {print $0}'
```

What do we see? The same contaminant OTUs as before:  

```
# Top OTUs from our mock samples:
MockIM32_pident=100.0_OTU10,0,98650
MockIM20_pident=100.0_OTU13,0,91076
MockIM42_pident=100.0_OTU11,0,86046
...
# Bottom OTUs from our mock samples:
MockIM47_pident=99.4_OTU47,0,40336
MockIM33_pident=100.0_OTU49,0,36402
MockIM49_pident=100.0_OTU56,0,28017
...
# Top OTUs from our non-mock samples
OTU1256_suspect_mock_variant,0,270    ## Harmonia
OTU2336_suspect_mock_variant,0,99
OTU1159_suspect_mock_chimera,0,61
OTU707_suspect_mock_chimera,0,37
OTU36_suspect_mock_variant,1,34
OTU163_suspect_mock_chimera,0,33
OTU524_suspect_mock_chimera,0,30
OTU205_suspect_mock_variant,0,29
...
# Remaining in list have less than 20 reads
```

## Final steps for `p10-1` data
My takeaways regarding the `p10-1` dataset and best filtering practice: there is no need to apply as large a `--subtract` filter in this dataset as initially calculated, because these top hits are from expected _Harmonia_ sequences. We _should apply an index bleed of **2.1%**_ as a relatively conservative estimate of index-bleed. Finally, there are specific OTUs which we should drop out directly, and that could be done by executing an `amptk drop` command below. **Note we DO NOT actually execute this command** because we have a pair of libraries which must be clustered collectively before any OTUs are dropped. These OTUs (and associated taxonomic identities) are to be compared with the other `p10-2` library; following those comparisons we'll determine which OTUs to drop collectively, what the `--subtract` value should be set to, and what the `--index_bleed` rate will be globally.  

Nevertheless, here's how you'd drop a series of OTUs you don't want:  

```
amptk drop \
--input path/to/{filename}.cluster.otus.fa \
--reads path/to/{filename}.demux.fq \
--list OTU3 OTU14 OTU15 OTU92 OTU653 \
--out dropdOTUs
```

This would generate a pair of new files to use in the final filtering steps: **dropdOTUs.cleaned.otus.fa** and **dropdOTUs.cleaned.otu_table.txt**. We'll generate those only _after_ we analyse `p10-2`, which is the next thing we'll do.  

## Filtering `p10-2` data
A few things hold true between both datasets:
1. We'd predict that at least some of the same contaminants to be present in our `p10-1` data as in our `p10-2` data, because we assumed some of these contaminants are present in the mock community (which was the same sample, added to each library separately); in addition, these low-level contaminants may possibly have been derived from the PCR mix itself, applied to both libraries - we'll probably find similar OTUs.  
2. We may find additional contaminants not present in the `p10-1` set; thus we'll need to start with a default filtering strategy to identify all possible contaminants.  
3. We will **not** normalize data - this would only exacerbate the problems explained above (principally that it is artificially inflating the `--index_bleed` calculation. Because of this, we can jump directly into the second part of the filtering, and apply `--normalize n` to our filter script with the `--subtract auto` argument passed.  

Let's begin with evaluating the dataset when we filter it without normalizing:  

```
amptk filter \
-i /path/to/p10-2.cluster.otu_table.txt \
-f path/to/p10-2.cluster.otus.fa \
-b mockP10L2IM4 \
--delimiter csv \
--keep_mock \
--calculate all \
--subtract auto \
--mc /mnt/lustre/macmaneslab/devon/guano/mockFastas/CFMR_insect_mock4alt.fa \
--debug \
--threshold max \
-o 10-2_noNorm \
--normalize n
```

Great. Similar to what we observed with `p10-1`, when we aren't normalizing the data the proportion of mock reads bleeding into the true samples, or the true samples bleeding into the mock community are lower than the expected 2%. The `--subtract` value is much higher than we want, however, so we'll investigate which OTUs are causing this number to spike as with the previous library - often these are the same contaminant issues, though note that the OTU numbers are not necessarily the same (as the OTU id for each cluster is independently assigned between each library).  

```
Index bleed, mock into samples: 0.452564%.  
Index bleed, samples into mock: 0.049501%.
Auto subtract filter set to 223
mockP10L2IM4 sample has 23 OTUS out of 25 expected; 0 mock variants; 0 mock chimeras; Error rate: 0.043%
Filtering OTU table down to 348 OTUs and 10,055,324 read counts
```

Let's figure out which OTUs are causing this `--subtract` value to be so high (we're suspecting it's probably that _Harmonia_ sample just like in `p10-1`):  

```
sed -i 's/#OTU ID/OTUid/' 10-2_noNorm.sorted.csv
awk -F ',' '{print NF; exit}' 10-2_noNorm.sorted.csv    ## prints out 244 (there are 243 samples)
cut 10-2_noNorm.sorted.csv -d ',' -f 1,2,244 | sort -t ',' -k3,3nr | awk -F "," '$3 != "0" {print $0}'
```

Output suggests a similar result - a handful of contaminant OTUs with a modest number of reads (less than 300 in any one sample), and a smattering of OTUs with less than 50 reads in the mock community which aren't expected.

```
# Bottom expected mock OTUs
MockIM7_pident=100.0_OTU40,0,81955
MockIM49_pident=100.0_OTU44,0,67726
...
# Top unexpected OTUs in mock
OTU469_suspect_mock_variant,0,223
OTU456_suspect_mock_chimera,0,89
OTU113_suspect_mock_variant,0,88
OTU346_suspect_mock_variant,0,80
OTU110_suspect_mock_chimera,0,66
OTU700_suspect_mock_variant,0,62
OTU1_suspect_mock_chimera,0,49
OTU93_suspect_mock_chimera,0,49
...
additional OTUs with < 34 reads
```

We next identify the identify these OTUs listed above. Note that because clustering processes were independent, the same sequence string may be assigned to a totally different OTU name:  

```
grep "\\bOTU469\\b" 10-2_noNorm.otus.counts.fa -A 1
grep "\\bOTU456\\b" 10-2_noNorm.otus.counts.fa -A 1
grep "\\bOTU113\\b" 10-2_noNorm.otus.counts.fa -A 1
grep "\\bOTU346\\b" 10-2_noNorm.otus.counts.fa -A 1
grep "\\bOTU110\\b" 10-2_noNorm.otus.counts.fa -A 1
grep "\\bOTU700\\b" 10-2_noNorm.otus.counts.fa -A 1
grep "\\bOTU1\\b" 10-2_noNorm.otus.counts.fa -A 1
grep "\\bOTU93\\b" 10-2_noNorm.otus.counts.fa -A 1
```

We also add information about how these contaminants are distributed in terms of numbers of reads per sample for each OTU above (only the two of eight command shown below):  

```
grep "OTU469_" 10-2_noNorm.sorted.csv
grep "OTU346_" 10-2_noNorm.sorted.csv
```

We find the following observations:
- `OTU469` was not detected in our top contaminant OUTs in the `p10-1` dataset, however it is only present in two samples: the mock community with 223 reads, our top contaminant OTU in that sample, yet just 2 reads in a single true sample. This is highly likely to be a mutant read from the _Agrotis ipsilon_ mock DNA sequence. It will be discarded from further analyses.
- `OTU456` is the same as `OTU135` from the `p10-1` dataset (the house fly, _Fannia_). It will be discarded.
- `OTU113` was not detected in our top contaminant OUTs in the `p10-1` dataset. It has a peculiar distribution, with the top sample having 3405 reads (from sample `OahuBird-646`), then the next most amount of reads is just 88 (in the mock sample), then the next highest read count is just 6. This strikes me as an instance of a real OTU that is finding its way via index bleed into a a few samples; further, a BLAST search reveals that it's a march grass moth which certainly could be part of the diet of a bird. We'll keep track of this specific sequence when we cluster our combined libraries - perhaps it's present in more samples than expected. It may have to be filtered but we may also have a way to subtract it's value out from the mock dataset. Stay tuned.  
- `OTU346` is a moth with several potentially high % identity/coverage options. As with `OTU113` above, there are just a few samples with any appreciable number of reads: the top value has 307 reads in a true sample, followed by another true sample with 95 reads, then our mock community sample with 80 reads, then the next highest has just 3 reads. This might be a PCR contaminant, but might also be a PCR error, or perhaps indicative on index-bleed. We'll flag this OTU and look for it in our post-clustered table when we combine libraries next.  
- `OTU110` doesn't have any specific BLAST match with high coverage, yet is broadly identified with moderate amounts of reads in many samples: 12 samples have > 100 reads, and 23 samples have > 50. If it's a contaminant, it's pervasive and must be removed; if it's a chimera, it may have been missed by the software. If it's a true sample, we'd expect to see it in both datasets, but we only seem to detect it in `p10-2`. We'll keep an eye on this one and see how many times it occurs in `p10-1` and `p10-2` when we cluster our libraries collectively.
- `OTU700` is a perfect match for a moth species found in continental North America. We see it in one and only one sample: our mock community It would be dropped from final analysis.  
- `OTU1` is found in many samples in high abundance. This is a classic example of index-bleed. It's a fruit fly found in Hawaii and thus likely not a contaminant from earlier sequencing runs. It would be retained, but the value identified in our mock may drive what we ultimately set our `--subtract` filter to. However, that may be a bit high of a value to set such a filter, given that this OTU would represent a situation perfect for a _high_ amount of index bleed because true sample shave such high read depths.  
- `OTU93` is found in the mock community and in many other samples with moderate read depth; I think this is another example of PCR contamination carrying over into multiple samples.  

## `p10-1` and `p10-2` filtering summary
Given that all we have detected in both `p10-1` and `p10-2` are minor levels of contamination in the mock, but potentially moderate levels of contamination in our true samples (likely from PCR mix) we'll take the following approach to analyze the dataset as one singular pile:  

- the `index_bleed` value isn't particularly high when we don't normalize our datasets; we'll set this value to optimize the removal of suspect OTUs from our mock without dropping too many reads from our true samples (explained in the final filtering section below)  
- the `--subtract` filter will need to be used if we find that there are suspect OTUs remaining after the `--index_bleed` filter has been applied; it's unclear what this value will be, but it may result in quite a few OTUs being dropped  
- there are many low-level contaminants which will need to be individually BLASTed and confirmed whether they are likely indicative of true samples or contaminants on a case-by-case basis. We'll explain each of those decisions in the single-dataset clustering output below.

# Filtering considerations for combined dataset, `OahuBird`

Reads from the two libraries, `p10-1` and `p10-2` were combined as described in the [workflow.md document](https://github.com/devonorourke/guano/blob/master/OahuBird/docs/workflow.md#combining-datasets-and-renaming-mock-community-headers). In addition, samples with less than 1000 reads were removed from the dataset. The sequences remaining from the samples that survived filtering were then clustered. The filtering approach discussed here applies to the resulting cluster of the filtered, joint libraries of `p10-1` and `p10-2` termed `OahuBird`.  

We begin by applying the default parameters within `amptk filter`, including normalizing our reads.

```
amptk filter \
-i /path/to/{filename}.otu_table.txt \
-f /path/to/{filename}.otus.fa \
-b mockOahuBirdIM4 \
--delimiter csv \
--keep_mock \
--calculate all \
--subtract auto \
--mc /mnt/lustre/macmaneslab/devon/guano/mockFastas/CFMR_insect_mock4alt.fa \
--debug \
--threshold max \
-o all_default \
--normalize y
```

And we'll look into the output files to check a few things:
- Do we see a similar index_bleed as observed in the two independent libraries as before?  
- Do we see a similar level of suspected contaminant OTUs in terms of read depth?  

```
sed -i 's/#OTU ID/OTUid/' all_default.sorted.csv
awk -F ',' '{print NF; exit}' all_default.sorted.csv    ## prints out 397 (there are 396 samples including the one mock)
cut all_default.sorted.csv -d ',' -f 1,2,397 | sort -t ',' -k3,3nr | awk -F "," '$3 != "0" {print $0}'
```

From the `.log` file it's apparent that our index-bleed calculations appear just like the other libraries in which normalization was performed prior to read filtering:  
```
Index bleed, mock into samples: 71.344836%.  
Index bleed, samples into mock: 0.054010%.
Auto subtract filter set to 15
mockOahuBirdIM4 sample has 23 OTUS out of 25 expected; 0 mock variants; 0 mock chimeras; Error rate: 0.092%
Filtering OTU table down to 3,673 OTUs and 12,613,743 read counts
```
Note that we lose about a third of our reads with these default parameters, due to the extremely high index-bleed calculation. The output looks as expected in terms of the distribution of reads of expected mock community OTUs relative to unexpected OTUs in the mock sample:

```
## Top 3 expected mock OTUs
MockIM5_pident=100.0_OTU18,0,229709
MockIM23_pident=100.0_OTU19,0,228075
MockIM29_pident=100.0_OTU23,0,226421
...
## Bottom 3 expected mock OTUs
MockIM7_pident=100.0_OTU42,0,126121
MockIM40_pident=100.0_OTU41,0,123435
MockIM49_pident=100.0_OTU47,0,95743
...
## Top UNexpected mock OTUs
OTU1730_suspect_mock_variant,0,650
OTU3130_suspect_mock_variant,0,255
OTU224_suspect_mock_chimera,0,122
OTU1357_suspect_mock_chimera,0,95
OTU949_suspect_mock_chimera,0,88
OTU554_suspect_mock_variant,0,80
OTU264_suspect_mock_variant,0,78
```

Now let's see what the results look like if we don't normalize our data prior to filtering:
> Not shown: we executed the same filtering commands as described for the normalized data above, but switched the `--normalize` flag from `y` to `n`    

```
Index bleed, mock into samples: 0.925361%.  
Index bleed, samples into mock: 0.067208%.
Auto subtract filter set to 650
mockOahuBirdIM4 sample has 24 OTUS out of 25 expected; 0 mock variants; 0 mock chimeras; Error rate: 0.087%
Filtering OTU table down to 425 OTUs and 16,892,039 read counts
```

Like with the independent libraries before, we see that by not normalizing our data we increase the number of reads retained, but decrease the number of OTUs present in our data. It's a drastic difference: the normalized data has more than eight times the number of OTUs, and it's highly likely these are not representative of distinct taxonomic units (it's more likely they're the result of PCR error, contamination, sequencing error, etc.). Let's check to see if we see the same _UNexpected_ OTUs in the non-normalized data:  

```
sed -i 's/#OTU ID/OTUid/' all_NoNorm.sorted.csv
cut all_NoNorm.sorted.csv -d ',' -f 1,2,397 | sort -t ',' -k3,3nr | awk -F "," '$3 != "0" {print $0}'
```

We find that the number of reads per OTU may have changed, but the order and their relative proportions are the same - this indicates that as we look into filtering these OTUs we can likely retain more reads by dropping a few particular contaminant OTUs, setting a modest index-bleed filter, and reducing the `--subract` value to something between what is expected between the normalized and non-normalized data. Here's a similar summary of selected expected and unexpected OTUs in the mock community:  

```
## Top 3 expected mock OTUs - same as in Normalized data
MockIM5_pident=100.0_OTU18,0,229709
MockIM23_pident=100.0_OTU19,0,228075
MockIM29_pident=100.0_OTU23,0,226421
## Bottom 3 expected mock OTUs - same as Normalized
MockIM7_pident=100.0_OTU42,0,126121
MockIM40_pident=100.0_OTU41,0,123435
MockIM49_pident=100.0_OTU47,0,95743
## Top UNexpected mock OTUs (more exist below what is shown here)
OTU1730_suspect_mock_variant,0,650
OTU3130_suspect_mock_variant,0,255
OTU224_suspect_mock_chimera,0,122
OTU1357_suspect_mock_chimera,0,95
OTU949_suspect_mock_chimera,0,88
OTU554_suspect_mock_variant,0,80
OTU264_suspect_mock_variant,0,78
```

We'll now look at those top unexpected OTUs to figure out whether or not they're the same taxa we observed in our independent libraries as before by BLASTing these OTUs in question and looking at their read distributions:  

```
## BLAST the resulting sequence (just the first few OTUs shown; many more analyzed in detail below)
grep "\\bOTU1730\\b" all_NoNorm.otus.counts.fa -A 1
grep "\\bOTU3130\\b" all_NoNorm.otus.counts.fa -A 1
grep "\\bOTU224\\b" all_NoNorm.otus.counts.fa -A 1

## Examine the read counts per sample (just one example shown)
grep "OTU700" all_NoNorm.sorted.csv
```

We find the following:
- `OTU1730` is _Harmonia axyridis_. We'll drop this OTU.
- `OTU3130` matches `OTU2336` from `p10-1`; it occurred just twice in that library (one of which is the mock) and that stays true here - this is clearly a contaminant and will be dropped from further analysis.  
- `OTU224` is our favorite house fly; `OTU456` from `p10-2` and `OTU135` from `p10-1`. It occurs in _lots_ of samples - 110 including the mock - yet just 18 of these samples have more than 100 reads total (the highest single read count being 654). This is very likely a contaminant and will be dropped from analysis.
- `OTU1357` matches `OTU1159` from `p10-1`, and like with `OTU3130` above, it remains in just two samples (one of which is the mock). It is a contaminant which will be dropped from analysis.
- `OTU949` matches `OTU707` from `p10-1`, it was detected in 28 samples, yet all at low frequency. It is another OTU which will be dropped.
- `OTU554` matches `OTU700` from `p10-2`. It was detected 31 times, yet no sample had more than 200 reads, and just three of these had more than 100 reads. It is a likely contaminant OTU which will be discarded.  
- `OTU264` matches `OTU205` from `p10-1`. It is found in low abundance in several samples and will be discarded as well.  
- `OTU279` matches `OTU113` from `p10-2`. It's a strange example because there are just two samples with more than 10 reads per sample: our mock (with 88) and a single true sample, **OahuBird-646**, with 3045 reads. That's a lot of reads to be a contaminant. Moreover, the sample itself doesn't seem particularly strange - it contains OTUs which are present in many other samples. In this particular case the tradeoff is dropping an OTU only in one sample (which contain many other OTUs) and lowering out `--subtract` value quite a bit, or keeping that OTU and leaving our `subtract` value at **88**. I'm inclined to just drop the OTU, as lowering the `subtract` value will retain more reads and more OTUs overall.
- `OTU69`matches `OTU346` from `p10-2`. This OTU is identified in seven samples with more than 100 reads, two of which have more than 10,000 reads. It's likely a real signal and shouldn't be dropped yet. Because of this, the `--subtract` value can likely be set to the number of reads detected in our mock for this OTU.  
- `OTU167` matches `OTU110` from `p10-2`. It's detected in over 200 samples, yet no single sample has more than 1500 reads (however, 31 samples have more than 100 reads, and 67 samples have more than 50 reads). If this is a contaminant, it's a highly abundant contaminant (we don't usually see these); if it's a chimera or true signal, it's unclear why this sequence isn't producing more reads in at least a few samples. What's curious is that the BLAST search generates a low % identity for any match (top hit just 89% over 100 % coverage); this makes me suspicious that it's in fact a chimeric sequence. However it could also be that it's a true signal which doesn't have good representation in our database. We'll retain this read for now.  
- `OTU1` is our top hit - _D. suzukii_ - and is clearly a real signal. It's a known Hawaiian fruit fly. If we had dropped `OTU69`, this likely would be where we'd set our `subtract` value. This OTU is retained.
- `OTU2285` is a match for `OTU700` in `p10-2`. It is identified in low read abundance in many samples , but the highest read number per sample is 55 (and that's in the mock!); this OTU will be removed.  

It's clear that several OTUs should be dropped from further analysis. Let's do this first and then determine how this influences our subsequent filtering (using default parameters and any modified ones thereafter). To drop the OTUs discussed above:  
```
amptk drop \
--input /mnt/lustre/macmaneslab/devon/guano/NAU/OahuBird/clust/all/OahuBird.cluster.otus.fa \
--reads /mnt/lustre/macmaneslab/devon/guano/NAU/OahuBird/illumina/trimd_OahuBird.demux.fq \
--list OTU1730 OTU3130 OTU224 OTU1357 OTU949 OTU554 OTU264 OTU279 OTU2285 \
--out OahuBird_otudropd
```

After compleing this process we've removed 9 OTUs. The main questions remaining to resolve concern what value we should set the `--index_bleed` to as well as the `--subtract` filter.  One option would be to apply an `index_bleed` value calculated from either the _bleed-in_ observed (~ 0.5%) when we don't normalize, or what we have observed for the _mock-in_ values (~1%); though note we would expect in an index bleed of ~2% in typical MiSeq run. Let's check what happens when we apply a 2% `index_bleed` on non-normalized data, and compare it with an `index_bleed` value calculated at 1% (close to the default estimate); we'll set the `--subtract` value to auto for the moment, then investigate which OTUs remain in our mock at modest values:  

>code not shown for the 1% estimate; change the `--index_bleed` value from **0.02** to **0.01**  

```
amptk filter \
-i /mnt/lustre/macmaneslab/devon/guano/NAU/OahuBird/clust/all/OahuBird_otudropd.cleaned.otu_table.txt \
-f /mnt/lustre/macmaneslab/devon/guano/NAU/OahuBird/clust/all/OahuBird_otudropd.cleaned.otus.fa \
-b mockOahuBirdIM4 \
--delimiter csv \
--index_bleed 0.02 \
--subtract auto \
--keep_mock \
--mc /mnt/lustre/macmaneslab/devon/guano/mockFastas/CFMR_insect_mock4alt.fa \
--debug \
--threshold max \
-o all_NoNorm_bleed2 \
--normalize n
```

We see there isn't a big difference between the two `index_bleed` values. The same number of OTUs are detected in the mock, the same number of OTUs remain, and there are just a few more sequences remaining in the lower `index_bleed` calculation (note that this value depends on whether we're counting the mock reads; if we remove them (by not including the `--keep_mock` argument above) we get a different output of OTUs and reads remaining; this is because that value includes (or doesn't) the mock reads depending upon if we pass that argument. For our purposes we should be comparing what the effect on read number and OTU number is when the mock is removed, as that is what our library will reflect).  
> Note that our final read counts and OTUs **include mock community reads and OTUs**; these will not be included in the final table  

```
## For `index_bleed` = 0.01
0.925364%.  Index bleed, samples into mock: 0.037660%.
Overwriting auto detect index-bleed, setting to 1.000000%
Auto subtract filter set to 89
mockOahuBirdIM4 sample has 24 OTUS out of 25 expected; 0 mock variants; 0 mock chimeras; Error rate: 0.087%
Filtering OTU table down to 1,148 OTUs and 18,144,460 read counts

## For `index_bleed` = 0.02
Index bleed, mock into samples: 0.925364%.  Index bleed, samples into mock: 0.037660%.
Overwriting auto detect index-bleed, setting to 2.000000%
Auto subtract filter set to 89
mockOahuBirdIM4 sample has 24 OTUS out of 25 expected; 0 mock variants; 0 mock chimeras; Error rate: 0.087%
Filtering OTU table down to 1,148 OTUs and 18,048,470 read counts
```

However, I'm concerned we're dropping a lot of OTUs which contain reads which aren't particularly high in abundance (per sample). What is driving the `subract` filter to be set to **89**?  

```
sed -i 's/#OTU ID/OTUid/' allOTUdropd_default.sorted.csv
cut allOTUdropd_default.sorted.csv -d ',' -f 1,2,397 | sort -t ',' -k3,3nr | awk -F "," '$3 != "0" {print $0}'
```

It's strange - we see that there are a different number of sequences mapping to these reads between our `.sorted.csv` files before and after dropping OTUs. In the case of the post-OTU dropped data, it's **OTU69** driving the `subtract` value - with a value of **89**. What's strange is that there were only **76** reads identified in the pre-OTU dropped data for this OTU. There must be a bit of stochasticity in the read mapping process itself which is driving this variation; while this doesn't seem all that big a difference, the `subtract` function really can dramatically influence the resulting dataset that is retained and/or discarded. For example, if we filter according to an `--index_bleed 0.01` and specify `--subtract 50` (rather than the auto-detected **89** value as above), we retain a lot more information:  
> code not shown. replace values specified in above statement using similar filtering code as usual  
```
Index bleed, mock into samples: 0.925364%.  Index bleed, samples into mock: 0.037660%.
Overwriting auto detect index-bleed, setting to 1.000000%
Subtracting 50 from OTU table
mockOahuBirdIM4 sample has 25 OTUS out of 25 expected; 0 mock variants; 1 mock chimeras; Error rate: 0.087%
Filtering OTU table down to 1,553 OTUs and 18,410,121 read counts
```

What we find is that we retain over **400 more OTUs**, yet those gains are made up in less than 2% of the total sequences available. If we look at the dataset this `--subtract` filter doesn't quite remove all the OTUs from our mock that shouldn't be in there (because we collapse two similar mock members in the clustering process, we really only have 24 expected members):  

```
sed -i 's/#OTU ID/OTUid/' allOTUdropd_sub50.normalized.csv
cut allOTUdropd_sub50.normalized.csv -d ',' -f 1,2,397 | sort -t ',' -k3,3nr | awk -F "," '$3 != "0" {print $0}'
```

Here we find that **OTU167** retains a total of 16 reads in our mock, while all other contaminant OTUs are filtered out. I think this is a reasonable tradeoff: we have a very, very minor amount of contaminant read in our positive control, yet we've retained many more OTUs than if we had applied a `subtract` value by our default value. To create the final filtered table, we'll re-run the `amptk filter` command as above with our OTU-dropped table/fasta files, as well as specify the subtract filter to a value of `50`.

```
amptk filter \
-i /mnt/lustre/macmaneslab/devon/guano/NAU/OahuBird/clust/all/OahuBird_otudropd.cleaned.otu_table.txt \
-f /mnt/lustre/macmaneslab/devon/guano/NAU/OahuBird/clust/all/OahuBird_otudropd.cleaned.otus.fa \
-b mockOahuBirdIM4 \
--delimiter csv \
--index_bleed 0.01 \
--subtract 50 \
--mc /mnt/lustre/macmaneslab/devon/guano/mockFastas/CFMR_insect_mock4alt.fa \
--debug \
--threshold max \
-o OahuBird \
--normalize n
```

In total we retain **1,529 OTUs** and **13,950,987 reads** (we've lost a lot of reads because of the high read count in our mock community. We'll next use the `OahuBird.filtered.otus.fa` and `OahuBird.final.csv` files in the taxonomy assignment portion of the workflow, as described in the [workflow document](https://github.com/devonorourke/guano/blob/master/OahuBird/docs/workflow.md).  
