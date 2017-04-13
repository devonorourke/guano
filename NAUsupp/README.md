# Overview of NAU contents

Supplementary data tables are available as .csv files, with brief summaries below. In addition, code executed using the AMPtk pipeline (used to take raw Illumina reads and ultimately assign taxonomy to unique OTUs) is available as a Jupyter Notebook - see **nau_amptk_pipeline.ipynb**. A single R script was used to generate data tables and figures shared - see **otu_analysis.R**.  

Please contact me should you have any questions.

devon.orourke@gmail.com

# Per sample explanations:

**BOLDonly.df.csv**
This is the master data table containing all metadata and taxonomic information ascribed to each sample. It is formatted with each row representing a unique observation per sample, per OTU. That is, a single guano sample may have 5 different OTU sequences present, each with their own taxonomic identities. These are represented as five independent rows, sharing identical metadata information, and just varying by the corresponding taxonomy classification and OTU info. Most additional data tables are derived from this master table. Importantly, this is not the entirety of the dataset. Several other OTUs not present here exist and are not identified here because there was no match for that OTU sequence with the reference (BOLD) database. It is possible to search for additional taxonomic information among those OTUs using the file **nauall.otus.taxonomy.fa**.  

**BarPlotArtFamilies.PerBatSpecies.PerSite.pdf**
A faceted bar plot showing the number of counts of OTUs categorized within particular taxonomic Families identified per site, per bat species. OTU detection is binary, such that the presence of the OTU was determined by meeting some minimal threshold of reads (in the case of this dataset, that minimum was 48 total).

**BarPlotArtOrders.PerBatSpecies.PerSite.pdf**
A faceted bar plot showing the number of counts of OTUs categorized within particular taxonomic Orders identified per site, per bat species. As above, each count indicates the presence of that OTU in a particular sample and is not an indication of relative numbers of reads within a sample. Note that the color scheme used in this depiction is identical to a separate figure **treemap.OrdersPerBatSpecies.SitesAggregated.pdf**.

**BatSpeciesTable.BOLDonly.csv**
A frequency table indicating the number of guano samples containing a particular bat species. **PTPA** is more than twice as prevalent as any other bat species in this analysis.

**EUMAdf.csv**
A subset of the larger **BOLDonly.df.csv** file containing just the information related to the single guano sample from EUMA.

**EUMAonly.OTUreads.csv**
A table demonstrating the total number of sequence reads detected on a per-OTU basis for the one EUMA sample. This is particularly useful when comparing the total number of OTUs present in the **EUMAdf.csv** table above with the relative numbers of sequence reads associated with each OTU. Notably there are just three OTUs accounting for 99% of all reads among the nine total OTUs detected.

**OTUcounts.perBatAndSite.csv**
The total number of OTUs detected for each bat species at each site. Data was used to generate the **OTUsPerBatSpeciesAndSite.pdf** file.

**OTUhistogram.pdf**
Summarizes the distribution of OTUs detected in a single guano sample across the entire study. Most samples contain between 5-15 OTUs.

**OTUsPerBatSpeciesAndSite.pdf**
A graphical representation of the frequency with which an OTU is detected per bat species at a given site. With nearly 1/6 of all samples associating with one bat species, we see that *PTPA* accounts for a large proportion of all OTUs detected at both *Escamecca* and *Escamaquita* sites.
