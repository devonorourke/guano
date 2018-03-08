# Supplementary notes on filtering
Filtering is an often neglected portion of amplicon analyses, despite the well documented occurrence of amplicon artifacts which can lead to inflation of overall richness and diversity of OTUs perceived across a dataset. There is no one way to filter. What follows is a series of steps taken to find a set of empirically derived filters which can be applied to our data. Code applied is documented herein, while data tables referenced are linked through [this spreadsheet](https://docs.google.com/spreadsheets/d/19LlFd7W81gD3AIJEnKlICdZQh-2EwNxpmJkzRuA9kMA/edit#gid=0). Note that there are multiple data tables, listed sequentially as S1, S2, etc., as tabs on the single spreadsheet document.

## Manipulating normalized data:

We'll start with the default parameters established by `amptk` - normalizing data and using the maximum value for a single OTU to calculate index bleed. Because we have two datasets which were independently clustered, we'll need to apply the same code to both `dropd` and `trim` data. The following example shows an example command used for the `dropd` dataset.  

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

### My takeaways regarding the `dropd` dataset and best filtering practice
There is no need to apply a highly conservative `--subtract` filter in this dataset, however we _should apply an index bleed of **2.1%**_ as a relatively conservative estimate of index-bleed. This will certainly reduce the number of OTUs we've identified in our raw OTU table, but it ultimately is just cutting out all the OTUs with very low read depth that are just as likely to be DNA present in a true sample at low concentrations as they are to be the product of tag-switching. We will _not_ normalized the dataset for the `dropd` samples, but we must make sure to remove any OTUs which are identified as matches for our mock community samples (this is done in the downstream R script). Finally, there are specific OTUs which we should drop out directly, and that can be done in the final filtering script.  

The final filtering command is as follows:  

```
amptk filter \
-i ../dropd.cluster.otu_table.txt \
-f ../dropd.cluster.otus.fa \
-b mockIM4p82redo \
--mc /mnt/lustre/macmaneslab/devon/guano/mockFastas/CFMR_insect_mock4alt.fa \
--index_bleed 0.01 \
--threshold max \
--calculate in \
--subtract 217 \
-o DropdFilt \
--delimiter tsv \
--normalize n
```

See **tableS4** for output of the final number of reads associated per sample per OTU.  

One final note: the filtering threshold applied to this dataset could also reasonably be applied to the other dataset which contained _all samples_; note, however, that because some samples contained less _total reads_ than the minimum read number on a per-OTU basis, some of these would be dropped anyway again. Because we based our filtering from index bleed into the mock community using samples with the greatest number of reads it's highly likely that the filtering parameters would have changed by including additional samples which contained even fewer reads. This was applied as a comparison, with the `-i` and `-f` flags reflecting the changed input files from `dropd` to `trim`. There were a few more OTUs preserved, a few more reads preserved, and more samples.  

See **tableS5** for output when using the `trimd` dataset.
