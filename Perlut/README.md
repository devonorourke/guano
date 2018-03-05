## Background

All files created for joint project with Foster and Maslo labs. See [Perlut_workflow.md](https://github.com/devonorourke/guano/tree/master/Perlut) for details describing molecular and bioinformatic processes applied to sequence data. Contact [devon.orourke@gmail.com](mailto:devon.orourke@gmail.com) with questions.  

---

## File description

- **OTU_per_Site.csv** - data table derived from _master.csv_ which calculates numbers of OTUs detected per site  
- **OTU_per_SiteWeek.csv** - data table derived from _master.csv_   which calculates numbers of OTUs detected per site per week  
- **OTUanalysis.R** - R script used to generate all non   amptk-derived outputs  
- **README.md** - this message  
- **filtering_notes.md** - workflow describing details of filtering strategies for `amptk filter` process  
- **genusPestMatch.csv** - data table derived by matching _master.csv_ `genus_name` observations with those from [pestlist.csv](https://github.com/devonorourke/guano/blob/master/pestlist.csv) file (a list of pests from Perdue, Canada Natural Resources, and the USDA)  
- **master.csv** - the fully filtered data set of all OTUs detected with all associated metadata  
- **metadata.txt** - data table with sample name, date of collection, and location information   
- **rut16_h.otu_table.taxonomy.txt** - raw output from `amptk taxonomy` script in binary format  
- **rut16_h.otus.taxonomy.fa** - updated fasta file from `amptk taxonomy` converting unlabeled OTUs to include as complete taxonomic information as possible  
- **rutgers_workflow.md** - description and example code provided for steps used in `amptk` pipeline
- **speciesPestMatch.csv** - data table derived by matching _master.csv_ `species_name` observations with those from _pestlist.csv_  
- **species_frq_table.csv** - data table describing frequency of detection of unique species, grouping all sites and weeks  
