## purpose: amptk OTU analysis
## written: 28-jan-2018
## edited: 13-mar-2018
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

## load libraries:
library(data.table)
library(reshape2)
library(tidyr)
library(ggplot2)
library(dplyr)

## read in data:
h_otutable.df <- fread('https://raw.githubusercontent.com/devonorourke/guano/master/Perlut/data/amptk/Perlut_h.otu_table.taxonomy.txt')

## reformat matrix into data.frame, split strings into appropriate columns
h_otutable.df[h_otutable.df == 0] <- NA       # we loaded a binary matrix; all zero's indicate absence of data, all 1's indicate presence of an amplicon
tmp.df <- melt(h_otutable.df)           # converting from wide format to long format (useful for plots)
rm(h_otutable.df)
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

## save this raw (unfiltered) data frame
setwd("~/Repos/guano/Perlut/data/Routput/")
write.csv(tmp.df, "Perlut_rawOTUtable.csv", quote = F)

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
            ######     Part 2 - data filtering     ######
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #

## A bunch of considerations here:
# 1. What is the lowest number of reads we should retain (should we keep all low reads)?
# 2. What OTUs are creeping into our Extraction NTCs?
# 3. What OTUs are clearly a result of contamination (ex. Chiropteran DNA)

## To address issue #1, let's plot the distribution of counts of reads per OTU looks like...
## note we're restricting our views to read counts <= 50,000 (anything bigger is going to be vastly rare)
## After plotting this first bit of code it's clear there are usually less than 10,000 reads per OTU...
histReads <- ggplot(data = tmp.df,
                    aes(x = CountReads)) +
  geom_histogram(binwidth = 100) +
  scale_x_continuous(limits = c(0, 50000))
histReads

## if we alter the x-axis to view just those between 0 - 1000 reads and change the 'binwidth'...
## to just 25 reads, and it's quite noticable that most of our data are coming from reads with less than 50 total sequences...
histReads_1000 <- ggplot(data = tmp.df,
                    aes(x = CountReads)) +
  geom_histogram(binwidth = 25) +
  scale_x_continuous(limits = c(0, 1000))
histReads_1000
  # it's clear we shouldn't cut out anything more than 50, but not clear if we should include less than 50

## What OTUs are occurring from our extraction NTCs?
## first, subset the data for all 'SampleID' values that match "extBlank" string
toMatch <- c("^extBlank", "nosample", "^solutBlank")    # make a list of the regex terms to match
matches <- unique (grep(paste(toMatch,collapse="|"),
                        tmp.df$SampleID, value=TRUE))   # make a vector of every one of those matches fitting the regex above
NTCs.df <- tmp.df[tmp.df$SampleID %in% matches]         # capture all matches of regex in 'tmp.df' object
  # there are 44 total OTUs (out of a possible 1734), or ~ 2.3% of our total OTUs among all samples...

## Let's a series of questions:
# 1. How many total reads are there among all these NTCs?
sum(NTCs.df$CountReads)
  # there are just 7,093 total reads among all NTCs idenitified to an OTU (recall we generated over 2 million reads in this dataset)
# 2. How many OTUs are associated for each unique NTC?
# 3. How many reads are associated for each NTC?
tmp_counts <- NTCs.df %>%
  group_by(SampleID) %>%
  count()
tmp_sums <- NTCs.df %>%
  group_by(SampleID) %>%
  summarise(TotalCounts = sum(CountReads))
NTC_summary <- merge(tmp_counts, tmp_sums)
rm(tmp_counts, tmp_sums)
  # These results suggest there isn't a single NTC sample driving the contamination
#4. By looking through the  'NTC.df' data frame we don't really see any single OTU driving contamination either...
# but we might be interested in finding out how many reads are associated per OTU, rather than per sample
NTC_OTU_sums <- NTCs.df %>%
  group_by(OTUid) %>%
  summarise(TotalCounts = sum(CountReads))
  # these results suggest there isn't really a single OTU associated with the majority of the contamination...
#5. The last thing we might want to find out is how many different samples (true and NTCs) share OTUs present in the NTCs:
NTC_otu.list <- NTCs.df$OTUid     # generate the list of all unique OTUs present in the NTCs
contam.df <- tmp.df[tmp.df$OTUid %in% NTC_otu.list,]    # find all matches from the above list in our original 'tmp.df' data frame object
  # we find there are 121 matches; we know 39 of these are from the NTCs themselves, so there must be 72 instances in which potential contaminant OTUs are present in true samples
#5...
## how many total reads are there?
sum(contam.df$CountReads)
  # 48,393 total reads...
## how many differnt samples have these suspected contaminant OTUs, and how many reads are there per sample?
tmp_counts <- contam.df %>%
  group_by(SampleID) %>%
  count()
tmp_sums <- contam.df %>%
  group_by(SampleID) %>%
  summarise(TotalCounts = sum(CountReads))
contam_summary <- merge(tmp_counts, tmp_sums)
rm(tmp_counts, tmp_sums)
  # some observations from that 'contam_summary' object:
  # 1. We see a lot of samples with some potential contaminants (34 true samples have at last one suspected contaminant)
  # 2. Nearly all of these true samples with some contaminant have relatively few unique OTUs (typically 1-3), with one exception (une-31)
  # 3. Of the contaminant-containing true samples, the total number of reads aren't negligible; given our amptk filtering with index-bleed and subtracting, these OTUs are probably present in those samples
  # 4. Several of these contaminant samples are identified in earlier bat guano; making me suspect these are carried over from the PCR mix... thinking it might be useful to filter against a list of the previous guano projects
  # 5. Sample 'une-31' has vastly more suspected OTU contaminants; unclear what happened there, but given the extent of contamination in that one sample it's safer to remove from analysis

## Let's figure out how to address problem #4 - are any of these things matches from earlier work?
## first, get the data from the old project
nauProject.df <- fread('https://raw.githubusercontent.com/devonorourke/guano/master/NAUsupp/nauProject_BINlist.txt', header = F)
colnames(nauProject.df) <- c("BOLDid", "Taxonomy")

## now clean up the data frame
nauProject.df <- separate(data = nauProject.df,
                   col = Taxonomy,
                   into = c("kingdom_name", "phylum_name", "class_name",
                            "order_name", "family_name", "genus_name", "species_name"),
                   sep = ",")
nauProject.df$kingdom_name <- sub("k:", "", nauProject.df$kingdom_name)
nauProject.df$phylum_name <- sub("p:", "", nauProject.df$phylum_name)
nauProject.df$class_name <- sub("c:", "", nauProject.df$class_name)
nauProject.df$order_name <- sub("o:", "", nauProject.df$order_name)
nauProject.df$family_name <- sub("f:", "", nauProject.df$family_name)
nauProject.df$genus_name <- sub("g:", "", nauProject.df$genus_name)
nauProject.df$species_name <- sub("s:", "", nauProject.df$species_name)

## append an alternate BOLDid column to the 'tmp.df' objecte to make matches with the old taxonomy naming convention in `amptk`
tmp.df$BOLDalt <- tmp.df$BOLDid
tmp.df <- separate(tmp.df, col = BOLDalt, into = c("BOLDalt", "delete"), sep = "\\.")
tmp.df$delete <- NULL

## get a list of those BOLDids... note we're probably missing some given that the old data set didn't follow the exact same taxonomy script (it was done over a year ago; things have changed within 'amptk')
nauProject.list <- nauProject.df$BOLDid

## Using that 'nauProject.list' object, query that list against the 'tmp.df$BOLDalt' vector
nauProjectMatches.df <- tmp.df[tmp.df$BOLDalt %in% nauProject.list,]
  # how many reads here?
  sum(nauProjectMatches.df$CountReads)    # 628,109 ... a big number, but > 500,000 are from just a one OTU (Dicranomyia (Idiopyga) halterella)
  # how many unique OTUs are identified?
  length(unique(nauProjectMatches.df$OTUid))    # 41 unique matches identified
## Finally, use that 'nauProjectMatches.df' object and determine how many reads and how frequently are these matches identified in our real data?
  tmp_counts <- nauProjectMatches.df %>%
    group_by(BOLDid) %>%
    count()
  tmp_sums <- nauProjectMatches.df %>%
    group_by(BOLDid) %>%
    summarise(TotalCounts = sum(CountReads))
  nauMatch_summary <- merge(tmp_counts, tmp_sums)
  rm(tmp_counts, tmp_sums)

  ## the 'nauMatch_summary' object makes it pretty clear: just a couple of OTUs appear to be repeatedly identified in our output, while the ...
  ## vast majority of contaminants are only every identified once or twice...
  ## however, it's clear that we can't conclusively determine whether these OTUs are true or artifacts from PCR mix contamination ...
  ## thus we'll remove all of these OTUs from further analysis; but we'll note a few of these highly occurring taxa and print out this data.frame

  setwd("~/Repos/guano/Perlut/data/Routput/")
  write.csv(nauMatch_summary, "suspectedContaminants.csv", quote = F)
  write.csv(tmp.df, "Perlut_rawOTUtable.csv", quote = F)

## Now let's finally remove those OTUs...
  # make a list of OTUs intended to be removed
  OTUdrop.list <- nauMatch_summary$BOLDid
  # then write a little function and remove them from 'tmp.df' object
  '%!in%' <- function(x,y)!('%in%'(x,y))
  tmpfilt.df <- tmp.df[tmp.df$BOLDid %!in% OTUdrop.list,]


## Now let's look back at this filtered table and see how many extraction blanks remain:
  ## What OTUs are occurring from our extraction NTCs?
  ## first, subset the data for all 'SampleID' values that match "extBlank" string
  NTCs.df <- tmpfilt.df[tmpfilt.df$SampleID %in% matches]
  # there are 43 total OTUs (out of a possible 1654)

  ## Let's a the same series of questions as before:
  # 1. How many total reads are there among all these NTCs?
  sum(NTCs.df$CountReads)
  # there are just 7,054 total reads among all NTCs idenitified to an OTU; clearly our NTC OTUs aren't the same as our previous run OTUs we filtered out
  # 2. How many OTUs are associated for each unique NTC?
  # 3. How many reads are associated for each NTC?
  tmp_counts <- NTCs.df %>%
    group_by(SampleID) %>%
    count()
  tmp_sums <- NTCs.df %>%
    group_by(SampleID) %>%
    summarise(TotalCounts = sum(CountReads))
  NTC_summary <- merge(tmp_counts, tmp_sums)
  rm(tmp_counts, tmp_sums)
  # These results suggest there isn't a single NTC sample driving the contamination
  #4. By looking through the  'NTC.df' data frame we don't really see any single OTU driving contamination either...
  # but we might be interested in finding out how many reads are associated per OTU, rather than per sample
  NTC_OTU_sums <- NTCs.df %>%
    group_by(OTUid) %>%
    summarise(TotalCounts = sum(CountReads))
  # these results suggest there isn't really a single OTU associated with the majority of the contamination...
  #5. The last thing we might want to find out is how many different samples (true and NTCs) share OTUs present in the NTCs:
  NTC_otu.list <- NTCs.df$OTUid     # generate the list of all unique OTUs present in the NTCs
  contam.df <- tmpfilt.df[tmpfilt.df$OTUid %in% NTC_otu.list,]    # find all matches from the above list in our original 'tmp.df' data frame object
  # we find there are 135 matches; let's remove these OTUs from out dataset as well
  #5...
  ## how many total reads are there?
  sum(contam.df$CountReads)
  # 48,177 total reads...
  ## how many differnt samples have these suspected contaminant OTUs, and how many reads are there per sample?
  tmp_counts <- contam.df %>%
    group_by(SampleID) %>%
    count()
  tmp_sums <- contam.df %>%
    group_by(SampleID) %>%
    summarise(TotalCounts = sum(CountReads))
  contam_summary <- merge(tmp_counts, tmp_sums)
  rm(tmp_counts, tmp_sums)
  # As observed before, most true samples have only 1 or 2 suspected contaminants; we'll remove all of these in our final filtered list

## remove these contaminants from our data frame:
tmpfilt2.df <- tmpfilt.df[tmpfilt.df$OTUid %!in% NTC_otu.list,]

## What's left to filter?
## remove the remaining OTUs which clearly are contaminants:
# 1. Chiroptera reads and reads assigned to Mock Community sequences
# 2. Other chordates from fish, reptiles, and amphibians (keeping only Arthropod information)
filtdOTUs.df <- subset(tmpfilt2.df, phylum_name == "Arthropoda")

setwd("~/Repos/guano/Perlut/data/Routput/")
write.csv(filtdOTUs.df, "filteredOTUtable.csv", quote = F)
write.csv(contam.df, "contaminantTable.csv", quote = F)

# cleanup...
rm(tmp.df, tmpfilt.df, tmpfilt2.df,
   contam_summary, contam.df,
   nauMatch_summary, nauProject.df,
   NTC_OTU_sums, NTC_summary, NTCs.df,
   matches, nauProject.list, NTC_otu.list, OTUdrop.list, toMatch)

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
          ######     Part 3 - incorporating metadata     ######
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #


## merge with metadata information
# meta.df <- fread('../Perlut_metadata.csv', header = T)
meta.df <- fread('https://raw.githubusercontent.com/devonorourke/guano/master/Perlut/data/Perlut_metadata.csv', header = T)
master.df <- merge(filtdOTUs.df, meta.df)
rm(filtdOTUs.df, meta.df)

## adding hyperlink to link BOLD BIN value to website
master.df$onclick <- paste("http://v4.boldsystems.org/index.php/Public_BarcodeCluster?clusteruri=",
           as.character(master.df$BOLDalt), sep = "")
master.df <- master.df[grepl("None", BOLDid), onclick := "no_link_available"];    # when BOLDid not available, removed broken link

# write file to disk:
setwd("~/Repos/guano/Perlut/data/Routput/")
write.csv(master.df, "master.csv", row.names = F, quote = F)

## Notrun: write.table(master.df, "PHINCHmaster.txt", row.names = F, quote = F, sep = "\t")

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
                 ######     Part 4 - data analyses     ######
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #

## load libraries and set working directory to print out data tables:
library(plyr)
setwd("~/Repos/guano/Perlut/data/Routput/")

## following our filtering, how many samples remain with at least 1 OTU? 2 OTUs? 10 OTUs?
OTUperSample = count(master.df, vars = c("SampleID"))   # There are 66 remaining true samples (all NTCs have been filtered out by dropping related OTUs from dataset)
sum(OTUperSample$freq > 1)    # There are 63 samples with at least 2 OTUs
sum(OTUperSample$freq > 4)    # There are 51 samples with at least 4 OTUs
sum(OTUperSample$freq > 9)    # There are 32 samples with at least 10 OTUs

## how many observations of OTUs contain complete information (ie. include 'species_name')... a.k.a. species frequency table
speciesOnly.df <- na.omit(master.df)
freq_species <- as.data.frame(table(speciesOnly.df$species_name))     # frequency table of species detected
colnames(freq_species) <- c("species_name", "counts")
write.csv(freq_species, "species_frq_table.csv", row.names = F, quote = F)
sum(freq_species$counts > 1)  # note 408 species identified, but most are not abundant (no OTU ID'd more than 11 samples)

## how many OTUs are called per site?
OTUperSite = count(master.df, vars = c("NestBox"))
colnames(OTUperSite) <- c("NestBox", "NumberOfDetections")
write.csv(OTUperSite, "OTU_per_Site.csv", row.names = F, quote = F)


## how many OTUs are called per site per week?
OTUperSiteWeek = count(master.df, vars = c("NestBox", "Date"))
colnames(OTUperSiteWeek) <- c("NestBox", "Date", "NumberOfDetections")
write.csv(OTUperSiteWeek, "OTU_per_SiteWeek.csv", row.names = F, quote = F)


# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
              ######     Part 5 - data visualization     ######
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #

# ... To be added soon
