# Overview of contents

Supplementary data tables are available as .csv files, with brief summaries below. In addition, code executed using the AMPtk pipeline (used to take raw Illumina reads and ultimately assign taxonomy to unique OTUs) is available as a Jupyter Notebook - see **sean_amptk_pipeline.ipynb**. A single R script was used to generate data tables and figures shared - see **otu_analysis.R**.  

Please contact me should you have any questions.

devon.orourke@gmail.com

# Per sample explanations:

**BOLDonly.df.csv**
This is the master data table containing all metadata and taxonomic information ascribed to each sample. It is formatted with each row representing a unique observation per sample, per OTU. That is, a single guano sample may have 5 different OTU sequences present, each with their own taxonomic identities. These are represented as five independent rows, sharing identical metadata information, and just varying by the corresponding taxonomy classification and OTU info. Most additional data tables are derived from this master table. Importantly, this is not the entirety of the dataset. Several other OTUs not present here exist and are not identified here because there was no match for that OTU sequence with the reference (BOLD) database. It is possible to search for additional taxonomic information among those OTUs using the file **sean.otus.taxonomy.fa**.  

**BOLDonlyMinRead110.df.csv**
Same information as with **BOLDonly.df.csv** (above), however the data has been subset to include only observations with at least 110 reads. Note these are the post-filtered, non-normalized values following the 'AMPtk filter' step as indicated in the **sean_amptk_pipeline.ipynb** file.

**OTUidTaxaClassifications.csv**
Table of the output of the **taxifying_blastout.R** script. Contains potential alignment matches which may improve upon the taxonomic resolution in the existing **BOLDonlyMinRead110.df.csv** file. Note that multiple hits may exist for a given OTU.

**OTUlist.110minReads.txt**
A single filed text file listing all OTUs present in the **min110Reads.fa** file. Used for the subsequent BLAST search used in the **taxifying_blastout.R** script.

**cleanBlastout.txt**
The partially filtered BLAST output. See Part 5 within **sean_amptk_pipeline.ipynb** for details. Information was used in the **taxifying_blastout.R** script to generate taxonomic assignments for each OTU present in the **OTUlist.110minReads.txt** file.

**masterdf.csv**
The complete data table containing all OTUs identified by all three classification programs: alignment of OTUs with USEARCH (to the BOLD database), and naive classification with UTAX and SINTAX. All OTUs and information assigned with UTAX and SINTAX were eliminated from any data tables and figures presented within these files. Further refinement of OTUs using another database such as Genbank's 'nr' database may provide additional taxonomic information not identified using these programs.  

**metadf.csv**
The metadata file containing Site, Date, BatSpecies, Sex, and BatTagID values for each guano sample analyzed.

**min110Reads.fa**
Same information as in **sean.otus.taxonomy.fa** but containing only OTUs which were part of the **BOLDonlyMinRead110.df.csv** subset. In other words, it's a reduced fasta file with only OTUs with at least 110 reads represented in every sample containing that OTU.

**otu_analysis.R**
The R script used to generate most of the data tables and figures presented herein.  

**taxifying_blastout.R**
R script which takes a fasta file of unknown taxonomic classification and assigns identity using Genbanks 'nr' database.


**sean_amptk_pipeline.ipynb**
A notebook explaining the methods and tools employed to process Illumina read data, cluster unique OTUs, and assigne taxonomic information to each OTU. The file contains explanations and motivations behind the entire workflow and includes links when possible to the programs used for each step.  

**sean.final.txt**
Tab-delimited file consisting of the total number of filtered reads mapped to each OTU called from the DADA2 pipeline on a per-sample basis. Very useful table when looking at the binary 'presence/absence' taxonomic information and determining whether or not an OTU (and subsequent classification) is something to further pursue, as well as deterimine the range of sequences per OTU per sample (ie. is there an equal number of sequences per OTU within a sample, or do just a few OTUs dominate all reads).

**sean.otu_table.taxonomy.txt**
Output file from Part 4 in the AMPtk pipeline (see **sean_amptk_pipeline.ipynb** ); contains a matrix of binary presence/absence of a given OTU for each sample. The file serves as the input for generating the **BOLDonly.df.csv** file (and derivatives) - see **otu_analysis.R** for details.

**sean.otus.taxonomy.fa**
The fasta file of all uniquely clustered OTUs in the dataset. A useful file to perform secondary alignments against non-BOLD databases for both confirmation of existing taxonomic classifications as well as corrections or improvements when lacking in our current BOLD-dependent calls.  


