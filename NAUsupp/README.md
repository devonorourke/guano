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

**PerSampleOTUcounts.bySiteandBatSpecies.csv**
The distribution of OTUs detected on a per-sample basis among all samples in the study.

**SiteTable.BOLDonly.csv**
A frequency table indicating the number of guano samples collected at a particular site. *ESCAMAQUITA* and *NICARAGUA* account for half of all samples analyzed.  

**Taxa.perBatSpecies.byFamily.csv** and **Taxa.perBatSpecies.byFamily.csv**
Data tables contining proportions of arthropod taxa identified for each bat species at the Family and Order level respectively. Note that if a bat species is present at multiple sites, those data are aggregated together in these tables.

**Taxa.perSite.byFamily.csv** and **Taxa.perSite.byOrder.csv**
Data tables contining proportions of arthropod taxa identified for each Site at the Family and Order level respectively. Note that if an arthropod taxa is preyed upon by multiple bat species at a given site those data are aggregated together in these tables.  

**batOTUs.df.csv**
Subset of the **BOLDonly.df.csv** file specific to any OTU classified as a bat. See explanation in R script for additional comments and cautions about interpreting these data (lines 251 - 253).

**masterdf.csv**
The complete data table containing all OTUs identified by all three classification programs: alignment of OTUs with USEARCH (to the BOLD database), and naive classification with UTAX and SINTAX. All OTUs and information assigned with UTAX and SINTAX were eliminated from any data tables and figures presented within these files. Further refinement of OTUs using another database such as Genbank's 'nr' database may provide additional taxonomic information not identified using these programs.  

**metadf.csv**
The metadata file containing Site, Date, BatSpecies, Sex, and BatTagID values for each guano sample analyzed.  

**nau_amptk_pipeline.ipynb**
A notebook explaining the methods and tools employed to process Illumina read data, cluster unique OTUs, and assigne taxonomic information to each OTU. The file contains explanations and motivations behind the entire workflow and includes links when possible to the programs used for each step.  

**nauall.otus.taxonomy.fa**
The fasta file of all uniquely clustered OTUs in the dataset. A useful file to perform secondary alignments against non-BOLD databases for both confirmation of existing taxonomic classifications as well as corrections or improvements when lacking in our current BOLD-dependent calls.  

**otu_analysis.R**
The R script used to generate most of the data tables and figures presented herein.  

**percFamily.bySiteandBatSpecies.csv** and **percOrderbySiteandBatSpecies.csv**
Similar to the **Taxa.perSite.by{Order|Family}.csv** files, but further subset the data into *both bat species and site, rather than aggregate for all sites or all bat species. Useful if you want to determine if site-specific or bat species-specific trends exist in your dataset.  

**treemap.OrdersPerBatSpecies.SitesAggregated.pdf**
A figure that depicts the area of a square in relative proportions to others based on total number of counts of that observation. In this instance, the larger boxes are separated by bat species. Boxes nested within each bat species visualize the proportion of a particular arthropod Order identified (on a per-count basis) for that bat, while accounting for the total number of counts for any particular Order between bats. Larger boxes, more frequently detected. Boxes are colored by taxonomic Order as indicated in the legend. The same coloring scheme is used in **BarPlot.ArtOrders.PerBatSpecies.PerSite.pdf**.
