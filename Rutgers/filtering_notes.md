I. Dropd data
OTU table contains 953 OTUs and 9,335,697 read counts

  A. Not normalized, MAX filter
Index bleed, samples into mock: 4.798870%.
Auto subtract filter set to 161526
mockIM4p82redo sample has 8 OTUS out of 24 expected; 0 mock variants; 0 mock chimeras; Error rate: 0.000%
Filtering OTU table down to 13 OTUs and 3,669,013 read counts

  B. Normalized, MAX filter
Index bleed, samples into mock: 4.795240%.
Auto subtract filter set to 4748
mockIM4p82redo sample has 8 OTUS out of 24 expected; 0 mock variants; 0 mock chimeras; Error rate: 0.000%
Filtering OTU table down to 133 OTUs and 6,391,405 read counts

  C. Normalized, sum filter, auto subtract off...
Index bleed, samples into mock: 4.795240%.
mockIM4p82redo sample has 18 OTUS out of 24 expected; 0 mock variants; 0 mock chimeras; Error rate: 0.047%
Filtering OTU table down to 932 OTUs and 7,045,341 read counts

  D. Normalized, top5 filter, autosubtract off...
Index bleed, samples into mock: 4.795240%.
mockIM4p82redo sample has 31 OTUS out of 24 expected; 0 mock variants; 0 mock chimeras; Error rate: 0.091%
Filtering OTU table down to 932 OTUs and 9,105,181 read counts

  B. Normalized, max filter, added OTU228 to mock.fasta
Index bleed, samples into mock: 0.047002%
mockIM4p82redo sample has 35 OTUS out of 25 expected; 0 mock variants; 0 mock chimeras; Error rate: 0.039%  
Filtering OTU table down to 932 OTUs and 9,330,163 read counts

## Manipulating normalized data:
```
sed -i 's/#OTU ID/OTUid/' filtMax.final.csv
sed -i 's/#OTU ID/OTUid/' filtMax.normalized.num.csv
awk -F ',' '{print NF; exit}' filtMax.final.csv
cut filtMax.normalized.num.csv -d ',' -f 1,2,69 | sort -t ',' -k2,2nr | awk -F "," '$2 != "0.0" {print $0}'
cut filtMax.final.csv -d ',' -f 1,2,69 | sort -t ',' -k2,2nr | awk -F  "," '$2 != "0" {print $0}'
```

  - We see that OTU228 is the culprit. It contains 4748.0 normalized reads. The lowest mock OTU contains 714 reads, while all other non-mock OTUs are <= 13 reads (most are < 2).

```
grep "\\bOTU228\\b" filtMax.normalized.num.csv
```
  - Table suggests that several samples contain this OTU in high amounts, making it a likely candidate for either index-bleed or contamination. However, because the mock community was not amplified with these samples, contamination via extraction or PCR is impossible. Thus index bleed is the only culprit.

To determine the likely taxa contributing to this contamination:

```
grep "\\bOTU228\\b" filtMax.filtered.otus.fa -A 1
```

  - This OTU blasts with 98% identity to Harmonia axyridis, a member of the mock community. This demonstrates that the OTU228 is a member of the mock community. Index bleed _into_ the mock community is very small among all other OTUs. However, there are two complicating factors at work:  

  1. This same OTU is present in large proportions in several other samples; because these data were normalized, it may be that there is significant index bleed from mock into true samples, but it is also possible that this proportion is greatly inflated because of the normalization process used to estimate these values. We'll explore the same dataset _without_ normalizing reads next.  
  2. This same OTU was a model organism used in the lab in which PCR (but not extraction) was performed. Thus the somewhat bimodal distribution of reads for this OTU (several high, several very low) in true samples (and negative controls) could likely be explained by potential contamination at the PCR step.

## Manipulating NON normalized data:

```
sed -i 's/#OTU ID/OTUid/' filtMaxNoNorm.final.csv
sed -i 's/#OTU ID/OTUid/' filtMaxNoNorm.sorted.csv
awk -F ',' '{print NF; exit}' filtMaxNoNorm.final.csv
cut filtMaxNoNorm.normalized.csv -d ',' -f 1,2,69 | sort -t ',' -k2,2nr | awk -F "," '$2 != "0.0" {print $0}'
cut filtMaxNoNorm.final.csv -d ',' -f 1,2,69 | sort -t ',' -k2,2nr | awk -F  "," '$2 != "0" {print $0}'
```

- Same OTU as with normalized data (which is expected). Because data isn't normalized, lowest mock read is much higher at 24291, while next highest non=mock is 433 reads. This confirms that the proportion of index bleed is minimal; if 433 represents the floor, and the mock community contained 3401801 total reads, then index bleed is less than 0.01 % (rather than the estimated ~5%)

The question remains as to whether or not OTU228 is a result of contamination or index bleed.

```
grep "\\bOTU228\\b" filtMaxNoNorm.sorted.csv
```

- The output of this table doesn't provide any definitive answer as to whether OTU228 is present in a sample due to index bleed or contamination. There are instances in which relatively low read depth contains high proportion of its overall reads to this one sample (which could be explained either by contamination or index bleed), while other samples with low read depth don't contain any reads for this OTU. When looking at non-normalized data, what is clear is the majority of samples don't contain significant amounts of this OTU - 44 of 68 samples have less than 5 reads each.

## Manipulating data with modified mock fasta

Using the same parameters as above, but I added OTU228 into the fasta file; this will effectively filter out that Harmonia read and recalculate index bleed _into_ mock at first. We'll need to explore index bleed from mock into other samples separately.  

The command:  

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

And the output filtering parameters:

```
sed -i 's/#OTU ID/OTUid/' max_otu.final.csv
sed -i 's/#OTU ID/OTUid/' max_otu.normalized.num.csv
awk -F ',' '{print NF; exit}' max_otu.final.csv
cut max_otu.normalized.num.csv -d ',' -f 1,2,69 | sort -t ',' -k2,2nr | awk -F "," '$2 != "0.0" {print $0}'
cut max_otu.final.csv -d ',' -f 1,2,69 | sort -t ',' -k2,2nr | awk -F  "," '$2 != "0" {print $0}'
```

Output is identical in terms of reads attributed to all non-mock OTUs (nothing changes in terms of read numbers or OTUs; all that's changed is the previous OTU288 is now identified as a mock sample). However, this also demonstrates that if you assume this OTU is part of the mock, then the index bleed _into_ mock is reduced to 1/100th of the original estimation - **this suggests that whether this OTU is a source of contamination or index bleed, we are better off removing this OTU by either assigning it to the mock community or ignoring it when caculating index bleed**.

I think our best approach is to remove any reads attributed to Harmonia once we've assigned taxonomy to our dataset, but it's worth mentioning in our analysis that it is possible these could be potential dietary sources in some instance. I looked into whether samples with significant numbers of Harmonia reads were associated with some sort of similar barcode; this doesn't appear to be the case. Of the 24 samples with >= 100 Harmonia reads (per sample), there were 10 forward and 10 reverse barcodes used. It is likely that similar barcodes are not the culprit for misidentification; however, index bleed is still as possibility to be explored such that the mock community Harmonia reads were bleeding into the true samples. We'll calculate that next.

For what it's worth, the non-normalized data were also compared with this modified fasta file, with the `--normalize` flag switched from `y` to `n`. As with the normalized data, the output looks identical to the non-normalized data in terms of raw reads assigned to OTUs that should/shouldn't be in the mock. In summary, whether or not we normalize data will not affect our calculation for index bleed **assuming we account for the OTU228** issue described above (there are multiple ways to deal with it). Index bleed with this single OTU accounted for reduces our estimation very close to the `amptk` default of 0.005.

## Calculating bleed from mock into true samples

Because the mock community had such a high proportion of reads relative to the overall library, the rate of index bleed from the mock _into_ the samples should be addressed. `amptk` caclulates this by default and we'll pass the same argument as above but modify just one item (turning off the `--calculate in` function):

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

---




sed -i 's/#OTU ID/OTUid/' top5.final.csv
sed -i 's/#OTU ID/OTUid/' top5.normalized.num.csv
awk -F ',' '{print NF; exit}' top5.final.txt
cut top5.normalized.num.txt -d "," -f 1,2,69 | sort -k2,2n | awk '$2 != "0.0" {print $0}'
cut top5.final.txt -f 1,2,69 | sort -k2,2n | awk '$2 != "0" {print $0}'


grep "\\bOTU228\\b" top5.normalized.num.csv



II. Trimd data
  A. No normalized
