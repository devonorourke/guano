## Background

All files created for joint project with Foster and Maslo labs. See [BRIpompton_workflow.md](https://github.com/devonorourke/guano/blob/master/BRIpompton/BRIpompton_workflow.md) for details describing molecular and bioinformatic processes applied to sequence data. Contact [devon.orourke@gmail.com](mailto:devon.orourke@gmail.com) with questions.  

---

## File description
- **README.md** - this message  

Within the `docs` directory:  
- **Pompton_filtering_notes.md** - details elucidating filtering strategies used within `amptk` pipeline  
- **Pompton_workflow.md** - description and example code provided for steps used in `amptk` pipeline  
- **OTUanalysis.R** - R script used to generate all non amptk-derived outputs  
- **renaming_scheme.md** - supplementary info describing how raw .fq were renamed from NAU files  

Within the `data` directory:  
- **Pompton_h.otu_table.taxonomy.txt** - raw output from `amptk taxonomy`  
- **Perlut_h.otus.taxonomy.fa** - updated fasta file from `amptk taxonomy` converting unlabeled OTUs to include as complete taxonomic information as possible  

Within the `scripts` directory:  
