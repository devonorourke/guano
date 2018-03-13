## purpose: amptk OTU analysis
## written: 28-jan-2018
## author: devon orourke

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # 
               ######     Part 1 - data wrangling     ######     
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #

## load libraries:
library(data.table)
library(reshape2)
library(tidyr)

## read in data:
w1.df <- fread('https://raw.githubusercontent.com/devonorourke/guano/master/Rutgers/rut16_h.otu_table.taxonomy.txt')

## reformat matrix into data.frame, split strings into appropriate columns
w1.df[w1.df == 0] <- NA       # we loaded a binary matrix; all zero's indicate absence of data, all 1's indicate presence of an amplicon
tmp.df <- melt(w1.df)           # converting from wide format to long format (useful for plots)
rm(w1.df)
colnames(tmp.df) <- c("OTUid", "Taxonomy", "SampleID", "Presence")
tmp.df <- tmp.df[complete.cases(tmp.df),]
tmp.df$Presence <- NULL
tmp.df <- separate(data = tmp.df, 
                    col = Taxonomy, 
                    into = c("TaxMethod", "AlignScore", "BOLDid", "kingdom_name", "phylum_name", "class_name", 
                             "order_name", "family_name", "genus_name", "species_name"), 
                    sep = "\\;|,|\\|")

## reformat taxonomy columns to discard given prefix
tmp.df$kingdom_name <- sub("k:", "", tmp.df$kingdom_name)
tmp.df$phylum_name <- sub("p:", "", tmp.df$phylum_name)
tmp.df$class_name <- sub("c:", "", tmp.df$class_name)
tmp.df$order_name <- sub("o:", "", tmp.df$order_name)
tmp.df$family_name <- sub("f:", "", tmp.df$family_name)
tmp.df$genus_name <- sub("g:", "", tmp.df$genus_name)
tmp.df$species_name <- sub("s:", "", tmp.df$species_name)

## remove Chiroptera reads; remove reads assigned to Mock Community sequences
tmp.df <- subset(tmp.df, order_name != "Chiroptera")
tmp.df <- tmp.df[grep("^MockIM", tmp.df$OTUid, invert = TRUE),]
tmp.df <- tmp.df[grep("^Mock-Harmonia", tmp.df$OTUid, invert = TRUE),]

## merge with metadata information
meta.df <- fread('https://raw.githubusercontent.com/devonorourke/guano/master/Rutgers/metadata.txt')
meta.df$SampleID <- sub("\\.", "-", meta.df$SampleID)
master.df <- merge(tmp.df, meta.df)
rm(tmp.df, meta.df)

## adding hyperlink to link BOLD BIN value to website
master.df$onclick <- paste("http://v4.boldsystems.org/index.php/Public_BarcodeCluster?clusteruri=",
           as.character(master.df$BOLDid), sep = "")
master.df$onclick <- gsub('.{2}$', '', master.df$onclick)     # had to remove last 2 characters for link to work
master.df <- master.df[grepl("None", BOLDid), onclick := "no_link_available"];    # when BOLDid not available, removed broken link

setwd("~/Desktop/guano/Rutgers/")
write.csv(master.df, "master.csv", row.names = F, quote = F)

## Notrun: write.table(master.df, "PHINCHmaster.txt", row.names = F, quote = F, sep = "\t")


# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # 
                 ######     Part 2 - data analyses     ######     
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #

## load libraries and set working directory to print out data tables:
library(plyr)
setwd("~/Desktop/guano/Rutgers/")

## following our filtering, how many samples remain with at least 1 OTU? 2 OTUs? 10 OTUs?
OTUperSample = count(master.df, vars = c("SampleID"))   # There are 92 remaining samples; 5 of these are negative controls
sum(OTUperSample$freq > 1)    # There are 74 samples with at least 2 OTUs (includes all 5 negative controls)
sum(OTUperSample$freq > 4)    # There are 41 samples with at least 5 OTUs (includes 2 negative controls)
sum(OTUperSample$freq > 9)    # There are 26 samples with at least 10 OTUs (includes 1 negative control)


## how many observations of OTUs contain complete information (ie. include 'species_name')... a.k.a. species frequency table
speciesOnly.df <- na.omit(master.df)
freq_species <- as.data.frame(table(speciesOnly.df$species_name))     # frequency table of species detected
colnames(freq_species) <- c("species_name", "counts")
write.csv(freq_species, "species_frq_table.csv", row.names = F, quote = F)   
sum(freq_species$counts > 1)  # note 172 species identified, but almost all rare (just 55 OTUs detected more than once)


## how many OTUs are called per site?
OTUperSite = count(master.df, vars = c("Location"))
OTUperSite$Location <- sub("^$", "-control", OTUperSite$Location)
write.csv(OTUperSite, "OTU_per_Site.csv", row.names = F, quote = F)   


## how many OTUs are called per site per week?
OTUperSiteWeek = count(master.df, vars = c("Location", "WeekOfYear"))
OTUperSiteWeek$Location <- sub("^$", "-control", OTUperSiteWeek$Location)
write.csv(OTUperSiteWeek, "OTU_per_SiteWeek.csv", row.names = F, quote = F)   


## what pests are detected?
pestlist.df <- fread('https://raw.githubusercontent.com/devonorourke/guano/master/pestlist.csv')

## create a list of unique names from 'pestlist.df': unique species as well as unique genera
pestSpeciesNames.df <- as.data.frame(unique(pestlist.df$LinneanName))     ## creates non-redundant list
colnames(pestSpeciesNames.df) <- c("genus_name")
pestSpeciesNames.df$genus_name <- as.character(pestSpeciesNames.df$genus_name)
pestSpeciesNames.df$species_name <- pestSpeciesNames.df$genus_name
pestSpeciesNames.df$genus_name <- gsub(" .*","", pestSpeciesNames.df$genus_name)  ## creates genus list, but not unique
pestGenusNames.df <- data.frame(unique(pestSpeciesNames.df$genus))      ## creates a unique data.frame of unique genera
colnames(pestGenusNames.df) <- "genus_name"
pestGenusNames.df$genus_name <- as.character(pestGenusNames.df$genus_name)
pestGenusNames.df$status <- "pest"
pestSpeciesNames.df$genus_name <- NULL                                         ## creates a unique data.frame of unique species
pestSpeciesNames.df$status <- "pest"
rm(pestlist.df)

## what 'species_name'values in 'master.df' match our 'pestSpeciesName.df' list?
library(tidyverse)
speciesPestMatch.df <- master.df %>% inner_join(pestSpeciesNames.df)

## what 'genus_name'values in 'master.df' match our 'pestSpeciesName.df' list?
    ## note these are not exact matches; 
    ## these merely represent insects in 'master.df' which share a genus with pests listed in the same genera
genusPestMatch.df <- master.df %>% inner_join(pestGenusNames.df)

## how often was each species detected?
genusPestMatch.table <- data.frame(table(genusPestMatch.df$species_name))
length(unique(genusPestMatch.df$species_name))    ## 26 unique matches (not all unique matches named to species level)

setwd("~/Desktop/guano/Rutgers/")
write.csv(speciesPestMatch.df, "speciesPestMatch.csv", row.names = F, quote = F)
write.csv(genusPestMatch.df, "genusPestMatch.csv", row.names = F, quote = F)

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # 
              ######     Part 3 - data visualization     ######     
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #

