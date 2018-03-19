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
## Top UNexpected mock OTUs
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
## BLAST the resulting sequence (just the first three OTUs shown)
grep "\\bOTU1730\\b" all_NoNorm.otus.counts.fa -A 1
grep "\\bOTU3130\\b" all_NoNorm.otus.counts.fa -A 1
grep "\\bOTU224\\b" all_NoNorm.otus.counts.fa -A 1
grep "\\bOTU1357\\b" all_NoNorm.otus.counts.fa -A 1
grep "\\bOTU949\\b" all_NoNorm.otus.counts.fa -A 1
grep "\\bOTU554\\b" all_NoNorm.otus.counts.fa -A 1
grep "\\bOTU264\\b" all_NoNorm.otus.counts.fa -A 1

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

**STARTHERE.. comment about which OTUs to drop, what `subtract` filter to set, what index bleed to set...**
**likely start by dropping certain OTUs, then applying a 2% index bleed filter, then looking at what OTUs remain, then applying _that_ `subtract` value. that should give us confidence in retaining that OTUs that are left**

Which produces a similar output as with out sanity check above, except the filtering produces fewer reads and OTUs (because we've now discarded anything associated with our mock community sample):  

```
Filtering OTU table down to 1,277 OTUs and 3,412,100 read counts
```

We carry forward the `finaldropd.filtered.otus.fa` and `finaldropd.final.txt` files into taxonomy assignment. The next portion of the filtering analysis focuses on the other dataset - the `trim` dataset which did not drop any samples. The same principles apply in how we investigate what OTUs to drop, what value to assign for the `--subtract` argument, and what the `--index_bleed` values should be. The code applied is documented below, but the results and explanations used are shortened, as they follow the similar ideas already described for the `dropd` dataset.  

**!!!!!!!!!!!!!!!!!**






```

- A bunch of OTUs are only present in the mock community sample: `OTU526`, `OTU770`,
- `OTU184` is present just in two samples (the same as in the `dropd` dataset)
- `OTU166` is present in 7 samples including the mock community; this makes me think it's probably not a contaminant at this point and will likely not be dropped from the final analysis. If that's the case, we'll have to increase the `--subtract` option a bit higher to eliminate it from the mock sample.  
- `OTU1155` was present in four samples, the highest being 23 total reads to the mock community. This will be dropped from analysis because all other true samples with any read number are less than 23, and we're setting the subtract filter higher than that value.  

There was one outlier OTU which didn't match up with our earlier `dropd` dataset: `OTU1181` in the `trimd` dataset. This was a beetle that a BLAST search only aligned 95% identity (across 100% coverage), thus the species-level assignment is likely incorrect. Nevertheless, it's clearly a beetle of some sort. However, this was detected in only one true sample in addition to the mock community, thus it would be dropped as a singleton OTU in our final processing. We'll eliminate it along with all the other OTUs as follows:  

```
amptk drop \
--input /mnt/lustre/macmaneslab/devon/guano/NAU/Perlut/clust/trim.cluster.otus.fa \
--reads /mnt/lustre/macmaneslab/devon/guano/NAU/Perlut/illumina/trim.demux.fq.gz \
--list OTU527 OTU1181 OTU770 OTU184 OTU1155 \
--out trim_dropdOTUs
```

We then apply a default filter once more to determine how dropping these OTUs influences the index-bleed calculation (we expect the subtract value to equal 24 because we've retained `OTU166`):  

```
amptk filter \
-i /mnt/lustre/macmaneslab/devon/guano/NAU/Perlut/dropd/trim_dropdOTUs.cleaned.otu_table.txt \
-f /mnt/lustre/macmaneslab/devon/guano/NAU/Perlut/dropd/trim_dropdOTUs.cleaned.otus.fa \
-b une-mockIM4 \
--delimiter csv \
--keep_mock \
--calculate all \
--subtract auto \
--mc /mnt/lustre/macmaneslab/devon/guano/mockFastas/CFMR_insect_mock4alt.fa \
--debug \
--threshold max \
-o trim_OTUdropd \
--normalize n
```

We find that the _subtract_ value was reduced to **24** as expected, and the _index bleed_ value remained at 3.8%.

```
Index bleed, mock into samples: 3.790306%.  
Index bleed, samples into mock: 0.012968%.
Auto subtract filter set to 24
une-mockIM4 sample has 24 OTUS out of 25 expected; 0 mock variants; 0 mock chimeras; Error rate: 0.081%
Filtering OTU table down to 1,117 OTUs and 3,416,600 read counts
```

The index bleed value is expected to be higher than the `dropd` data set because of the default way it's calculated: it sums up the entirety of mock-associated reads in true samples, and divides that sum by the number of reads in the actual mock sample (on a per OTU basis). That percentage is then applied to a `--threshold` value, which can be calculated a variety of ways (the sum of all OTUs, the single max OTU, the top5 % OTU values, etc.), but the way the percent is determined is always the same. Because we have about 2x the number of samples in the `trimd` data set, we have a greater changes that additional mock read will be part of the sum of the non-mock sample counts, thus there's an intrinsic tradeoff here: more samples to investigate means a higher index bleed.  

> Side note - if you want to see what these index-bleed calculations should look like, try running this little R script:
```
# written: 8-March-2018
# author: devon o'rourke
# Part 1 - Motivation:
## Jon palmer described the index-bleed being calculated as follows:

#   
#    So how it is actually being calculated is as follows:
#      1) OTU table is sliced for the mock community sample
#      2) Total number of reads in mock community is calculated
#      3) the total number of reads from non-mock community OTUs are then summed
#      4) index-bleed into the mock community is then calculated by dividing total non-mock reads by the total reads in the mock
#
## We want to know what OTUs are causing the index-bleed to be elevated:

# Part 2 - script
#runonce:
install.packages('data.table')
## load packages
library(data.table)
## read in data
reads.df <- fread('https://raw.githubusercontent.com/devonorourke/guano/master/Perlut/trim_OTUdropd.sorted.csv')
## create matrix from data.frame, creating row.names from first column of data.frame:
reads.mat <- as.matrix(reads.df[,-1])
## let's first drop the last column from the matrix of values to calculate the top non-mock read number present in every OTU (row):
nonmock.mat <- reads.mat[,1:87]
## now we'll calculate the sum of all those non-mock reads per OTU
non_mock_sum <- apply(nonmock.mat, 1, sum)
## pull out a vector of just the mock community values:
mock_vals <- reads.mat[,88]
## create a data.frame with those two vectors, then make a new column that follows Jon's logic above:
df <- data.frame(non_mock_sum, mock_vals)
options(scipen = 999)   ## this converts the forthcoming `df$perc` values to be printed as intergers rather than in scientific notation
df$sumperc <- (df$non_mock_sum / df$mock_vals * 100)
df <- subset(df, mock_vals != 0)
```

I tend to trust this higher value because of the fact that we have so many more mock reads than any other single true sample. We're still retaining a lot of OTUs (and my guess is many of these will be dropped once we remove singleton OTUs in the next R script), but the benefit of increasing this index-bleed value is that we can now retain more individual guano samples. Even if that means we only detect a few OTUs per sample, we have many more samples.  

The final filter applied is similar to the `dropd` set, except we're going to increase our index bleed calculation to 3.8%, and set the subtract filter to 24.  

```
amptk filter \
-i /mnt/lustre/macmaneslab/devon/guano/NAU/Perlut/dropd/trim_dropdOTUs.cleaned.otu_table.txt \
-f /mnt/lustre/macmaneslab/devon/guano/NAU/Perlut/dropd/trim_dropdOTUs.cleaned.otus.fa \
--mc /mnt/lustre/macmaneslab/devon/guano/mockFastas/CFMR_insect_mock4alt.fa \
-b une-mockIM4 \
--delimiter tsv \
--index_bleed 0.038 \
--subtract 24 \
-o finaltrim \
--normalize n
```

This produces the following output after mock reads and OTUs are removed:  

```
Index bleed, mock into samples00: 3.790306%.  
Index bleed, samples into mock: 0.012968%.
Subtracting 24 from OTU table
une-mockIM4 sample has 24 OTUS out of 25 expected; 0 mock variants; 0 mock chimeras; Error rate: 0.081%
Filtering OTU table down to 1,094 OTUs and 2,198,389 read counts
```

We'll carry the `trim` final output file pair (`finaltrim.filtered.otus.fa` and `finaltrim.final.txt`) into the final taxonomy analysis.
