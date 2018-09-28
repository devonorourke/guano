A brief description of the following datasets:

- **Mangan.mappingFile.txt** - Contains sample ID's and associated metadata (ex. collection dates, site names, roost IDs). Used to process `.biom` file.
- **Mangan.otus.taxonomy.fa** - The fasta file of sequences for OTUs used in this analysis which were index-bleed and chimera-filtered and then taxonomy was assigned. Includes some OTUs which were filtered by the R script used to identify potential contaminant reads and bat OTUs.
- **Mangan.final.txt** - An OTU table. Specifically, a matrix with columns as OTUs (first column) and Sample ID's (columns 2 - (N-1)), and rows as counts of sequence reads.
- **Mangan.biom** - A `.biom` (JSON-encoded) file containing the same taxonomy and sequence abundance data that can be visualized online with programs like [phinch.org](phinch.org)
- **Mangan.chordata.csv** - A comman-seprated text file containing the frequencies of unique Chordata removed from the dataset (including bat, mule deer, fish, and reptiles). Data includes BOLDid's and OTUid's for potentially more updated alignment queries.
