## purpose: OTU analysis for M. Mangan and J. Foster - M. sodalis diets
## written: 28-sept-2018
## author: devon orourke

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
######     Part 1 - data wrangling     ######
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #

## install packages (run once) - commented out by default
# install.packages('data.table')
# install.packages('reshape2')
# install.packages('tidyr')
# install.packages('ggplot2')
# install.packages('dplyr')
# install.packages('plyr')
#   source("https://bioconductor.org/biocLite.R")
#   biocLite("biomformat")

## load libraries:
library(reshape2)
library(tidyr)
library(ggplot2)
library(dplyr)
library(biomformat)

## read in data:
setwd("~/Desktop")

otutable <- read.csv('https://raw.githubusercontent.com/devonorourke/guano/master/Mangan/data/Mangan.otu_table.taxonomy.txt', sep="\t")
colnames(otutable)[1] <- ""

## reformat matrix into data.frame, split strings into appropriate columns
otutable[otutable == 0] <- NA       # we loaded a binary matrix; all zero's indicate absence of data, all 1's indicate presence of an amplicon
tmp.df <- melt(otutable)           # converting from wide format to long format (useful for plots)
rm(otutable)
colnames(tmp.df) <- c("OTUid", "Taxonomy", "SampleID", "CountReads")
tmp.df <- tmp.df[complete.cases(tmp.df),]
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

## reformat `SampleID` column to remove the "X" positioned because it hates numbers as sample names
tmp.df$SampleID <- sub("^X", "", tmp.df$SampleID)

## remove the BOLDid `.`
tmp.df$BOLDid <- sub("\\..*", "", tmp.df$BOLDid)

## remove the `blankS39` sample, which represents a known contaminated single-tube sample
tmp.df <- tmp.df %>% 
  filter(!grepl('blankS39', SampleID))

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
######     Part 2 - data filtering     ######
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #

## 1a. Remove all chordates:
arth.df <- subset(tmp.df, phylum_name == "Arthropoda")
## 1b. Generate a table of the frequencies with which unique chordate OTU/BOLD ID's occur:
chordata.df <- subset(tmp.df, phylum_name == "Chordata") %>%
  group_by(OTUid, BOLDid, order_name, family_name, genus_name, species_name) %>%
  count()
write.csv(chordata.df, file="Mangan.chordata.csv", row.names = FALSE, quote = FALSE)

## how many unique OTUs remain?
length(unique(arth.df$OTUid)) ## lots: 1,906
## how many unique BOLDid's?
length(unique(arth.df$BOLDid)) ## 1,499; this is why you assign taxonomy...

## 2. Among any remaining negative control samples, what OTUs are present and how are the distributed among our true samples?
## subset the data for all 'SampleID' values that match "ExtractionNTC" string
toMatch <- c('ExtractionNTC*')    # make a list of the regex terms to match
matches <- unique (grep(paste(toMatch,collapse="|"),
                        arth.df$SampleID, value=TRUE))   # make a vector of every one of those matches fitting the regex above
NTCs.df <- arth.df[arth.df$SampleID %in% matches,]         # capture all matches of regex in 'arth.df' object

## how diverse is the NTC dataset among those 9 samples?
length(unique(NTCs.df$OTUid))     # there are 86 unique OTUs (out of a 1,906 )

## How many total reads are there among all these NTCs?
sum(NTCs.df$CountReads) # there are 313,924 total reads among all NTCs idenitified (out of 11,388,887 remaining)
                        # this amounts to about 3% of all the remaining data

## Identify OTUs in NTCs also found in true samples
NTC_otu.list <- NTCs.df$OTUid
## Create a data frame that counts the frequency of occurances of each OTU among contaminant samples, and calculate read depth for each OTU across NTC samples
contamcounts.df <- NTCs.df %>% 
  group_by(OTUid, BOLDid) %>%
  summarise(NTCreads = sum(CountReads), NTChits = n())
## Create a data frame that sums read depth across all matches AND counts frequency of occurance of each OTU among overlaps betweeen negative control and true samples
contam.df <- arth.df[arth.df$OTUid %in% NTC_otu.list,] %>%
  group_by(OTUid, BOLDid) %>%
  summarise(ALLreads = sum(CountReads), ALLhits = n())
## add in the NTC info from `contamcounts.df`
contam.df <- merge(contamcounts.df, contam.df, by = c("OTUid", "BOLDid"))
## get differnces in which read depth sums and frequencies are subtracted from entire dataset (to define the sums/depths associated exclusively with true samples)
contam.df$TRUEhits <- contam.df$ALLhits - contam.df$NTChits
contam.df$TRUEreads <- contam.df$ALLreads - contam.df$NTCreads
## calculate percent of reads and freqencies that NTC comprise for a given OTU
contam.df$NTChitFrac <- contam.df$NTChits / contam.df$ALLhits
contam.df$NTCreadFrac <- contam.df$NTCreads / contam.df$ALLreads

## save to disk
write.csv(contam.df, file = "Mangan.contam.csv", row.names = FALSE, quote = FALSE)

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
#   We're going to remove OTUs that have >= 20% for either 
#   1. `NTChitFRac` or :
#   2. `NTCreadFRac`
#   Essentially, we're dropping OTUs that are more associated with the contaminants than true samples
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #

## create a drop.list object to pull subset out suspected contaminant OTUs
hitsOTUdrops.df <- subset(contam.df, NTChitFrac > 0.2)
drop.list.tmp1 <- as.character(hitsOTUdrops.df$OTUid)
readsOTUdrops.df <- subset(contam.df, NTCreadFrac > 0.2)
drop.list.tmp2 <- as.character(readsOTUdrops.df$OTUid)
drop.list <- (unique(c(drop.list.tmp1, drop.list.tmp2)))

## create filtered dataset removing these OTUs from the `arth.df` object:
filt.df <- arth.df[!arth.df$OTUid %in% drop.list,]
length(unique(filt.df$OTUid))   ## 1,883 OTUs remain
length(unique(filt.df$BOLDid))  ## 1,487 BOLDid's remain

## save to disk
write.csv(filt.df, file = "Mangan.contamFilt.csv", row.names = FALSE, quote = FALSE)

## cleanup
rm(arth.df, chordata.df, contam.df, contamcounts.df, hitsOTUdrops.df, NTCs.df, OTUdropslist.hits, readsOTUdrops.df, tmp.df)
rm(drop.list, drop.list.tmp1, drop.list.tmp2, matches, NTC_otu.list, toMatch)


# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
######     Part 3 - incorporating metadata     ######
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #

## merge with metadata information
#not run: filt.df <- read.csv('https://raw.githubusercontent.com/devonorourke/guano/master/Mangan/data/Mangan.contamFilt.csv', header = TRUE)
meta.df <- read.csv('https://raw.githubusercontent.com/devonorourke/guano/master/Mangan/data/Mangan.metadata.csv')
master.df <- merge(filt.df, meta.df)
#rm(filt.df, meta.df)

# write file to disk:
write.csv(master.df, "Mangan.master.csv", row.names = F, quote = F)


# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
######     Part 4 - data analyses     ######
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #

