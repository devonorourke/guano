# Supplementary notes on filtering
Filtering is an often neglected portion of amplicon analyses, despite the well documented occurrence of amplicon artifacts which can lead to inflation of overall richness and diversity of OTUs perceived across a dataset. There is no one way to filter. What follows is a series of steps taken to find a set of empirically derived filters which can be applied to our data.  

# Filtering the `dropd` dataset 
We'll start with the default parameters established by `amptk` - normalizing data and using the maximum value for a single OTU to calculate index bleed. Because we have two datasets which were independently clustered, we'll need to apply the same code to both `dropd` and `trim` data. We'll carry though the entire filtering analysis for just the `dropd` dataset first, then switch to the `trimd` data thereafter. 

## Manipulating normalized data:
The following example shows an example command used for the `dropd` dataset.  

> Note a new directory `filtd` was created prior to the execution of this code to retain output files  

```
amptk filter \
-i /mnt/lustre/macmaneslab/devon/guano/NAU/Perlut/clust/dropd.cluster.otu_table.txt \
-f /mnt/lustre/macmaneslab/devon/guano/NAU/Perlut/clust/dropd.cluster.otus.fa \
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

The output generates the following summary statistics:  

```
Index bleed, mock into samples: 71.011751%.  
Index bleed, samples into mock: 0.103005%
Auto subtract filter set to 45
une-mockIM4 sample has 19 OTUS out of 25 expected; 0 mock variants; 0 mock chimeras; Error rate: 0.125%
Filtering OTU table down to 877 OTUs and 2,887,946 read counts
```
The parameters above are identifying that:  
- A large number of reads expected to be present in the mock community are mapping to true samples (`Index bleed, mock into samples`); this isn't surprising given that nearly 1/3 of the lane consisted of the mock community. This was due to an improperly balanced library. As such, reads that are present in the mock community should be filtered out from the real data; we should take note of what these ~20 expected OTUs are in the mock becuase there is a chance that some of these OTUs are likely to be in our true samples (this is always on a case-by-case basis with every differnt project)  
- A fairly small fraction of true samples are moving into the mock community (`Index bleed, samples into mock`); in addition, because we've normalized our dataset in this initial filtering analysis it's very likely that there has been limited index bleed of any true sample into our OTUs (that's a good thing).  
- The `Auto subtract filter` identifies the number necessary to remove all unexpected reads in the mock sample. The value isn't particularly high, however we see that we drop out a few _expected_ OTUs from our mock community. We'll want to next look into what the OTUs are that are causing this "subtract" filter to be activated in the first place (ie. the OTUs present in the mock sample that shouldn't be there)  
- The final OTU table and number of reads are sharply reduced after filtering the `dropd` data; we lose about a third of our overall reads, and cut the number of OTUs in half. This is because of the absurdly high mock bleed in rate, and the subtract filter working together.   

The following commands are used to investigate what OTUs might be finding their way from true samples into the mock community, as well as finding how many reads (and from which OTUs) from the mock community are bleeding into the true samples - this uses the output of the filtered sample above. 

```
sed -i 's/#OTU ID/OTUid/' mockIn.final.csv
sed -i 's/#OTU ID/OTUid/' mockIn.normalized.num.csv
awk -F ',' '{print NF; exit}' mockIn.final.csv
cut mockIn.normalized.num.csv -d ',' -f 1,2,38 | sort -t ',' -k3,3nr | awk -F "," '$3 != "0.0" {print $0}'
```

It's clear from the output of the fourth line of the code above that our mock community looks good - there are at least 2000 reads in every mock OTU, and the next highest number of reads of an unexpected read has just 45 (normalized) reads In fact, there are only three potential contaminant OTUs which are identified in our mock which we are suspicious of. We could also perform a blast search to identify the taxonomy that could be assigned to each OTU:

```
grep "\\bOTU1004\\b" mockIn.otus.counts.fa -A 1
grep "\\bOTU924\\b" mockIn.otus.counts.fa -A 1
grep "\\bOTU190\\b" mockIn.otus.counts.fa -A 1
```

Similarly, if you wanted to review the specific numbers of reads within each of those particular OTUs:
```
grep "\\bOTU1004_suspect_mock_variant\\b" mockIn.normalized.num.csv
grep "\\bOTU924_suspect_mock_variant\\b" mockIn.normalized.num.csv
grep "\\bOTU190_suspect_mock_variant\\b" mockIn.normalized.num.csv
```

The output from each of the three BLAST searchers are revealing. `OTU1004` is a 98% identity match (with 100% coverage) for _Harmonia axyridis_, a mock community member itself. `OTU924` is a 98% match (with 99% coverage) for _Harmonia axyridis_ also. These two initial results are likely pointing to the fact that we have a mutation in our mock community member that is driving the false representation of contamination; these OTUs can be ignored and thus the `--subtract` value we should be cautious of is that of the remaining OTU.  

The last BLAST search for `OTU190` was a 100% identity match with 100% coverage for two _Gretchena_ species (_amatana_ and _deludana_); a third match at 99% identity was associated with _Gretchena watchungana_. This last search indicates the potential for contamination into our mock is very low - the total number of reads that this one OTU was bleeding into our mock community? **2** reads. We can therefore proceed without adding the `--subtract` filter at all.  

The remaining concern is to what to do about the high index bleed of the mock community members into our true samples. As mentioned previously, this is to be expected because of the unbalanced library that was sequenced: there were far more reads dedicated to the single mock community sample than to any true sample. Thus there is a far greater likelihood of these OTUs present in the mock community to have found there way into the true samples. However, I wondered if the incidence of bleed into the mock community was actually inflated by the fact that we were normalizing our reads in the first place. Thus the next step was to run the same command filtering our reads as above, but passing setting the `--normalize` argument to `n` ("no" normalizing); this is explained in the next section.  

## Manipulating NON normalized data:

The same filtering approach was used as described above, with the `--normalize` flag switched from `y` to `n`. Output files were renamed to specify the fact we turned normalization off. I left **on** the `--subtract auto` option for the time being so we could make a direct comparison to the last filtering command:  

```
amptk filter \
-i /mnt/lustre/macmaneslab/devon/guano/NAU/Perlut/clust/dropd.cluster.otu_table.txt \
-f /mnt/lustre/macmaneslab/devon/guano/NAU/Perlut/clust/dropd.cluster.otus.fa \
-b une-mockIM4 \
--delimiter csv \
--keep_mock \
--calculate all \
--subtract auto \
--mc /mnt/lustre/macmaneslab/devon/guano/mockFastas/CFMR_insect_mock4alt.fa \
--debug \
--threshold max \
-o mockIn.noNorm \
--normalize n
```


By _not normalizing_, we see a monumental difference in the overall index bleed rate, and this is entirely attributed to the calculation of mock reads bleeding into the true samples:  

```
Index bleed, mock into samples: 2.067058%.  
Index bleed, samples into mock: 0.106595%.
Auto subtract filter set to 551
une-mockIM4 sample has 24 OTUS out of 25 expected; 0 mock variants; 0 mock chimeras; Error rate: 0.101%
Filtering OTU table down to 184 OTUs and 3,262,931 read counts
```

Pretty interesting. Using non-normalized data significantly reduced the estimation of index-bleed overall, yet the subtract filter now jumped up. We'll evaulate what's behind the subtract filter being so high once more (though I suspect it's going to be those _Harmonia_ reads causing the headache) in a minute. The first point to stress is that by _not normalizing_ our data, we have a very different picture of how often mock reads are finding their way into our community: it's tiny - just 2.1, which is right in line with what Jon and the default `amptk` program generally suggests (2% is the default). This suggests we shouldn't be normalizing this data, given how skewed our distribution of reads in the mock community is.  

Another takeaway: we have a lot fewer OTUs, but have retained many more reads. This is expected: a drop in OTUs should occur because we've observed an _increase_ in the **subtract** value (551 versus the previous 45); likewise, a _decrease_ in the index bleed should retain more overall reads.  The next question is to determine whether or not that **subtract** increase is warranted:  q 

```
sed -i 's/#OTU ID/OTUid/' mockIn.noNorm.sorted.csv
cut mockIn.noNorm.sorted.csv -d ',' -f 1,2,38 | sort -t ',' -k3,3nr | awk -F "," '$3 != "0" {print $0}'
```

What do we see? The top hitter is indeed the same OTU as before - `OTU1004`, which an earlier BLAST search determined was _Harmonia axyridis_. Likewise, if you look at the output from the above `cut` command, you'll notice a that there is a huge difference in the expected number of reads from our mock members, and about a 100-fold drop in read depth for all other 'contaminant' OTUs present in our mock sample - I've annotated the output below.

```
# Top OTUs from our mock samples: 
MockIM32_pident=100.0_OTU7,103,70571
MockIM20_pident=100.0_OTU8,194,66401
MockIM27_pident=100.0_OTU9,337,63916
...
# Bottom OTUs from our mock samples:
MockIM47_pident=99.4_OTU34,55,32224
MockIM40_pident=100.0_OTU33,83,32166
MockIM49_pident=100.0_OTU38,135,24965
...
# Top OTUs from our non-mock samples
OTU1004_suspect_mock_variant,1,551    ## Harmonia (from earlier BLAST search)
OTU487_suspect_mock_variant,0,302
OTU924_suspect_mock_variant,1,243     ## Harmonia (from earlier BLAST search)
OTU701_suspect_mock_chimera,0,49
OTU175_suspect_mock_variant,0,26
OTU190_suspect_mock_chimera,105,24    ## Gretchena (from earlier BLAST search)
OTU1079_suspect_mock_chimera,0,23
...
# Remaining in list have less than 15 reads
```

So we see that the top hit for a 'contaminant' OTU is the same Harmonia as before, but what about `OTU487` and others not identified? Performing the same BLAST search with the earlier commands to identify the sequence to search for:  

```
grep "\\bOTU487\\b" mockIn.noNorm.otus.counts.fa -A 1
grep "\\bOTU701\\b" mockIn.noNorm.otus.counts.fa -A 1
grep "\\bOTU175\\b" mockIn.noNorm.otus.counts.fa -A 1
grep "\\bOTU1079\\b" mockIn.noNorm.otus.counts.fa -A 1
```

Also as with before, we can determine the read depth per sample for each OTU in question:  

```
sed -i 's/#OTU ID/OTUid/' mockIn.noNorm.sorted.csv
grep "OTU487" mockIn.noNorm.sorted.csv
grep "OTU701" mockIn.noNorm.sorted.csv
grep "OTU175" mockIn.noNorm.sorted.csv
grep "OTU1079" mockIn.noNorm.sorted.csv
```

Observations from each OTU:
- `OTU487` didn't have a clear strong winner in the BLAST search. The highest match was 96% identity across 100% coverage (pretty good!), but that species is found in Australia (so unless these sparrows travel to South Wales??...). What's more likely is we're within the right family (and possibly genus). However, the number of samples for which this OTU is present is just 3 (out of 88!), and the reads per sample are distributed as 302, 97, and 1. In fact, the sample with 302 reads is our mock community! This strongly suggest to me that this is a chimeric read and can be discarded from analyses entirely.  
- `OTU701` was more convincing, as it has 99% identity across 100% coverage in the BLAST search. However, there were a total of just 49 reads in the entire library, and all of these were present in the mock community sample (not present in any true sample); this again may reflect some sort of chimera forming, though I'm doubtful of that given that we have a perfect alignment. Perhaps this species was present in very low abundance in the mock sample as a contaminant? Either way, because the mock sample isn't included in our analyses downstream, we can elimiate this OTU entirely.  
- `OTU175` was unequivocally matched as _Psilocorsis quercicella_ (oak leaftier moth). It's an OTU we see in our bat guano samples and I wouldn't be surprised to find in the bird samples either. However, the fact that this OTU is present in very low abundance (just 26 reads in our OTU), and identified in only one other sample (with 606 reads) makes me suspicious that such a read is of value to retain. We likely would discard it from our downstream analyses because we often discard singleton OTUs anyway (as it can really throw a wrench in diversity estiamtes). My inclination is to drop this OTU from further analysis entirely
- `OTU1079` is a classic problem of trying to determine if something present in our lab could also be present out in nature: the sample is a clear match for _Fannia canicularis_ (common housefly). This same OTU has been detected in many other projects, including bat and bird guano from the northeast, western American regions as well as Central American countries. My guess is that this contaminant was present in the master mix used because it's the common component shared across all projects (water was changed, primers were changed, multiple extraction buffers were used). It's a particularly low-level contaminant, with only three samples having any identifiable reads (with 23, 19, and 16 reads represented in the mock community and two true samples, respectively). Given it's low read abundance and low read depth, I would discard this OTU entirely from analysis. 

## Final filtering steps for `dropd` data
My takeaways regarding the `dropd` dataset and best filtering practice: there is no need to apply as large a `--subtract` filter in this dataset as initially calculated, because these top hits are from expected _Harmonia_ sequences. We _should apply an index bleed of **2.1%**_ as a relatively conservative estimate of index-bleed. Finally, there are specific OTUs which we should drop out directly, and that can be done by executing the following `amptk drop` command:  

```
amptk drop \
--input /mnt/lustre/macmaneslab/devon/guano/NAU/Perlut/clust/dropd.cluster.otus.fa \
--reads /mnt/lustre/macmaneslab/devon/guano/NAU/Perlut/dropd/dropd.demux.fq \
--list OTU1004 OTU487 OTU924 OTU701 OTU175 OTU190 OTU1079 \
--out dropdOTUs
```

This generate a pair of new files to use in the final filtering steps: **dropdOTUs.cleaned.otus.fa** and **dropdOTUs.cleaned.otu_table.txt**. We'll use these files for our final filtering. First, perform a sanity check to look through the filtered files and ensure that the appropriate OTUs are present in the dataset:  

```
amptk filter \
-i /mnt/lustre/macmaneslab/devon/guano/NAU/Perlut/dropd/dropdOTUs.cleaned.otu_table.txt \
-f /mnt/lustre/macmaneslab/devon/guano/NAU/Perlut/dropd/dropdOTUs.cleaned.otus.fa \
-b une-mockIM4 \
--delimiter csv \
--keep_mock \
--calculate all \
--subtract auto \
--mc /mnt/lustre/macmaneslab/devon/guano/mockFastas/CFMR_insect_mock4alt.fa \
--debug \
--threshold max \
-o dropdOTUs \
--normalize n
```

Great. So what we see here is that our subtract filter has been reduced to just 15 (for the remaining OTUs which crept into our mock sample), our index bleed was set to 2.1%, and we have no mock variants or chimeras in our positive control (mock community).  

```
Index bleed, mock into samples: 2.067033%.  
Index bleed, samples into mock: 0.010760%.
Auto subtract filter set to 15
une-mockIM4 sample has 24 OTUS out of 25 expected; 0 mock variants; 0 mock chimeras; Error rate: 0.102%
Filtering OTU table down to 1,277 OTUs and 3,412,100 read counts
```

Let's finish up this process by applying the final filter, this time removing all the temporary files and discarding the reads associated with the mock community:  

```
amptk filter \
-i /mnt/lustre/macmaneslab/devon/guano/NAU/Perlut/dropd/dropdOTUs.cleaned.otu_table.txt \
-f /mnt/lustre/macmaneslab/devon/guano/NAU/Perlut/dropd/dropdOTUs.cleaned.otus.fa \
--mc /mnt/lustre/macmaneslab/devon/guano/mockFastas/CFMR_insect_mock4alt.fa \
-b une-mockIM4 \
--delimiter tsv \
--index_bleed 0.021 \
--subtract 15 \
-o finaldropd \
--normalize n
```

Which produces a similar output as with out sanity check above, except the filtering produces fewer reads and OTUs (because we've now discarded anything associated with our mock community sample):  

```
Filtering OTU table down to 1,277 OTUs and 3,412,100 read counts
```

We carry forward the `finaldropd.filtered.otus.fa` and `finaldropd.final.binary.txt` files into taxonomy assignment. The next portion of the filtering analysis focuses on the other dataset - the `trim` dataset which did not drop any samples. The same principles apply in how we investigate what OTUs to drop, what value to assign for the `--subtract` argument, and what the `--index_bleed` values should be. The code applied is documented below, but the results and explanations used are shortened, as they follow the similar ideas already described for the `dropd` dataset.  

# Filtering the `trim` dataset
A few things hold true between both datasets: 
1. We expect_at least_ the same contaminants to be present in our `trim` data as in our `dropd` data, because the `trim` set is just a subset of the `dropd`.  
2. We may find additional contaminants not present in the `trim` set; thus we'll need to start with a default filtering strategy to identify all possible contaminants.  
3. We will **not** normalize data - this would only exacerbate the problems explained above (principally that it is artificially inflating the `--index_bleed` calculation. Because of this, we can jump directly into the second part of the filtering, and apply `--normalize n` to our filter script with the `--subtract auto` argument passed.  

## Default filtering 

```
amptk filter \
-i /mnt/lustre/macmaneslab/devon/guano/NAU/Perlut/clust/trim.cluster.otu_table.txt \
-f /mnt/lustre/macmaneslab/devon/guano/NAU/Perlut/clust/trim.cluster.otus.fa \
-b une-mockIM4 \
--delimiter csv \
--keep_mock \
--calculate all \
--subtract auto \
--mc /mnt/lustre/macmaneslab/devon/guano/mockFastas/CFMR_insect_mock4alt.fa \
--debug \
--threshold max \
-o trim_default \
--normalize n
```

Produces this output:
```
Index bleed, mock into samples: 3.791080%.  
Index bleed, samples into mock: 0.049073%.
Auto subtract filter set to 302
une-mockIM4 sample has 24 OTUS out of 25 expected; 0 mock variants; 0 mock chimeras; Error rate: 0.081%
Filtering OTU table down to 267 OTUs and 3,300,018 read counts
```

Apparently with a marginally higher index bleed we see a big drop in the number of OTUs retained relative to our final `dropd` filtering argument set `--index_bleed` to 2.1%. We'll now identify who those contaminants might be:  

