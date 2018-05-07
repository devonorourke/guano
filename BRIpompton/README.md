## Background

All files created for joint project with the Foster and Biodiversity Research Institute (BRI). See [BRIpompton_workflow.md](https://github.com/devonorourke/guano/blob/master/BRIpompton/docs/Pompton_workflow.md) for details describing molecular and bioinformatic processes applied to sequence data. Contact [devon.orourke@gmail.com](mailto:devon.orourke@gmail.com) with questions.  

---

## File description
- **README.md** - this message  

Within the `docs` directory:  
- **OTUanalysis.R** - R script used to generate all non amptk-derived outputs  
- **Pompton_filtering_notes.md** - details elucidating filtering strategies used within `amptk` pipeline  
- **Pompton_workflow.md** - description and example code provided for steps used in `amptk` pipeline  
- **pompton-rename.md** - supplementary info describing how raw .fq were renamed from NAU files  

Within the `data` directory, there are two main subdirectories for each of the major bioinformatic process: `amptk` and `Routput`.  
- **Within the `amptk` subdirectory:  
  - **noFilt.final.csv** is the unfiltered output from the `amptk` taxonomy process; a matrix type file with rows as OTUs and sample names as columns; elements within the matrix represent reads per sample/OTU
  - **NTCreduced.otu_table.txt** is a similar table as above, but removes potential contaminants detected in the negative control samples  
  - **Pompton_h.otu_table.taxonomy.txt** - raw output from `amptk taxonomy`  
  - **Perlut_h.otus.taxonomy.fa** - updated fasta file from `amptk taxonomy` converting unlabeled OTUs to include as complete taxonomic information as possible  
  - **rawCounts.txt** - the total number of reads detected per Sample (inclusive to any OTU per sample)
  - **rough.cluster.otu_table.txt** - the raw output from the `amptk clust` command; no taxonomy assigned

- **Within the `Routput` subdirectory:
  - **BRIPompton_FullFilteredOTUtable_noMock.csv** - the mock-removed OTU table which serves as the foundation for most metadata anlyses  
  - **BRIPompton_FullFilteredOTUtable_withMock.csv** - as above, but with mock reads/OTUs included  
  - **BRIPompton_rawOTUtable.csv** - essentially just the `Pompton_h.otu_table.taxonomy.txt` file but converted from a matrix to a single-observation type data table - the foundational object from which the `OTUanalysis.R` filtering is applied to 
  - **master.csv** the final object which incorporates all filtering steps for OTUs per sample and includes metadata; data is organized such that each row is an observation of some OTU for a single sample. Columns represent distinct variables. 
  - **OTU_per_Site_and_Species.csv** - data table analyzing the number of OTUs detected for each bird species at each unique Site 
  - **OTU_per_Site.csv** - data table analyzing the number of OTUs detected for each Site (all bird species merged)  
  - **OTU_per_Territory.csv** - data table analyzing the number of OTUs detected for each Territory (all bird species merged)
  - **OTU_counts.csv** - data table with number of detections per OTU
  - **species_frq_table.csv** - data table with number of detections per OTU which included species name information; not inclusive to entirey data set (only about 25% of all detected OTUs had information which included species name)... infer with caution 
  - **suspectedContaminants.csv** - a data table of OTUs which were removed from dataset because they appeared in other analyses from projects with which similar OTUs should not be shared due to geographic or host species distinctiveness 
  
- Within the `scripts` directory are all shell scripts used to execute `amptk functions`. Only one example (`illumina.sh`) currently provided, but I can certainly include any other ones should these be of interest. The commands are also generally demonstrated in the `Pompton_workflow.md` and `Pompton_filtering_notes.md` files, but do not include the additional features specific to our SLURM job submission commands required for the computer cluster with which these processes were executed. 

