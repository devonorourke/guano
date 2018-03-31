## purpose: amptk OTU analysis for A. Scinto and J. Foster OahuBird diets
## written: 28-jan-2018
## edited: 21-mar-2018
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
setwd("~/Repos/guano/OahuBird/")               
               
h_otutable.df <- fread('https://raw.githubusercontent.com/devonorourke/guano/master/OahuBird/data/amptk/OahuBird_lulu_h.otu_table.taxonomy.txt')
colnames(h_otutable.df)[1] <- ""

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

## append any of the 'SampleID' values which were actually NTCs (but not labeled as such):  
## import target and replacement data.frame:
replace.df <- fread('https://raw.githubusercontent.com/devonorourke/guano/master/OahuBird/data/Rrenamelist.df')

## then use dplyr to make a new column for the replacement matches
tmp2.df <- dplyr::left_join(tmp.df,replace.df, by = "SampleID")
## then substitute any `N/A` value with the original value for those values that didn't have a NTC match:
tmp.df <- within(tmp2.df, X <- ifelse(is.na(replacements), SampleID, replacements))
tmp.df$SampleID <- NULL
tmp.df$replacements <- NULL
colnames(tmp.df)[13] <- "SampleID"
rm(tmp2.df, replace.df)

## save this raw (unfiltered) data frame
setwd("~/Repos/guano/OahuBird/data/Routput/")
write.csv(tmp.df, "OahuBird_rawOTUtable.csv", quote = FALSE, row.names = FALSE)

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
            ######     Part 2 - data filtering     ######
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #

## One way to address potential sources of contamination is to identify the OTUs present in samples across multiple libraries...
## ...with the idea being that if you see an OTU present in a diverse set of libraries, it's likely this is a contaminant
## We're going to grab data from three separate projects: 
# 1. the `nau.df` project was from a bat guano collected in Central America
# 2. the `rut.df` project was from bat guano collected in New Jersey
# 3. the `perl.df` project was from bird guano collected in Maine
## We don't expect any of the overlapping OTUs to be present in Hawaiian samples 

rut.df <- fread('https://raw.githubusercontent.com/devonorourke/guano/master/Rutgers/master.csv', header = TRUE)
perl.df <- fread('https://raw.githubusercontent.com/devonorourke/guano/master/Perlut/data/Routput/master.csv', header = TRUE)
nau.df <- fread('https://raw.githubusercontent.com/devonorourke/guano/master/NAUsupp/masterdf.csv', header = TRUE)
colnames(nau.df)[3] <- c("BOLDalt")

## Before we can match we need to append the `BOLDid` vectors to remove the "." delimiter in the `tmp.df` and `rut.df` objects: 
tmp.df$BOLDalt <- tmp.df$BOLDid
tmp.df <- separate(tmp.df, col = BOLDalt, into = c("BOLDalt", "delete"), sep = "\\.")
tmp.df$delete <- NULL

rut.df$BOLDalt <- rut.df$BOLDid
rut.df <- separate(rut.df, col = BOLDalt, into = c("BOLDalt", "delete"), sep = "\\.")
rut.df$delete <- NULL


## Get a list of all non-redundant `BOLDalt` elements among all three projects:
# because 'none' is not a unique identifier we're removing it from this list
# we're also removing 'CFMR:IM4', 'SINTAX', and 'UTAX' as these were earlier relics of an amptk naming scheme we're not using
string <- unique(c(rut.df$BOLDalt, nau.df$BOLDalt, perl.df$BOLDalt))
find.list <- list("None", "UTAX", "SINTAX", "CFMR:IM4")
find.string <- paste(unlist(find.list), collapse = "|")
tmpx.list <- gsub(find.string, replacement = "", x = string)
allmatch.list <- tmpx.list[tmpx.list != ""]
rm(string, find.list, find.string, tmpx.list)


## Using that 'allmatch.list' object, query that list against the 'tmp.df$BOLDalt' vector to find common BOLD(alt) id's
ProjectMatches.df <- tmp.df[tmp.df$BOLDalt %in% allmatch.list,]
# how many reads here?
sum(ProjectMatches.df$CountReads)    # 5,295,342 ... a huge number; represents about 2/5 of our overall data (13,950,987 reads)
# how many unique OTUs are identified?
length(unique(ProjectMatches.df$OTUid))    # 501 unique matches identified... (our data originally had 1416 OTUs)
## Finally, use that 'nauProjectMatches.df' object and determine how many reads and how frequently are these matches identified in our real data?
tmp_counts <- ProjectMatches.df %>%
  group_by(BOLDalt) %>%
  count()
tmp_sums <- ProjectMatches.df %>%
  group_by(BOLDalt) %>%
  summarise(TotalCounts = sum(CountReads))
Match_summary <- merge(tmp_counts, tmp_sums)
rm(tmp_counts, tmp_sums)

x <- tmp.df[,c(1,14)]         ## grab just the OTUid and BOLDalt vectors from the `tmp.df` object
y <- x[!duplicated(x[1:2]),]  ## obtain only unique combinations of OTUid and BOLDids (note there can be the same BOLDid for multiple OTUid's)
z <- merge(Match_summary, y)  ## paste these unique OTUid's for every BOLDid in the `Match_summary` object to identify trends in OTU read abundance and frequency of detection

## now group the datasets from our `Match_summary` list to include just the "BOLDalt" and "SampleID" columns, then make a frequency table from that:
tiny.nau = nau.df[,c(1,3)]
colnames(tiny.nau) <- c("SampleID", "BOLDalt")
tiny.nau$library <- "nau"

tiny.rut = rut.df[,c(1,17)]
colnames(tiny.rut) <- c("SampleID", "BOLDalt")
tiny.rut$library <- "rut"

tiny.perl = perl.df[,c(1,14)]
colnames(tiny.perl) <- c("SampleID", "BOLDalt")
tiny.perl$library <- "perl"

tmpMatch.df <- rbind(tiny.nau, tiny.perl, tiny.rut)  ## this is used for our frequency counts
## drop the 'SINTAX', 'UTAX', 'None', and 'CFMR:IM4' labels as they won't produce any meaningful matches:
tmp1.df <- subset(tmpMatch.df, BOLDalt != "SINTAX")
tmp2.df <- subset(tmp1.df, BOLDalt != "UTAX")
tmp3.df <- subset(tmp2.df, BOLDalt != "CFMR:IM4")
match.df <- subset(tmp3.df, BOLDalt != "None")
rm(tmp1.df, tmp2.df, tmp3.df)


## Let's build a "naughty.list" that would include any BOLDalt value that:
# Test-A: is detected in more than two samples in any library at least once (discard this)
# Test-B: passes Test-A, but is detected in at least 2 libraries (discard this)

## this generates a list counting the number of times a 'BOLDalt' item is identified, not the number of unique "library" values present...
library(dplyr)
x <- match.df %>%
  group_by(library, BOLDalt) %>%
  tally() 

## Test-A: create list of values for all non-singleton occurrences of a BOLDalt entry per library:
y <- subset(x, n <= 2)
remove1 <- subset(x, n > 2)

## Test-B: filter out duplicate entries by "BOLDalt" as these are matches present in multiple libraries
remove2 <- y[duplicated(y$BOLDalt),]
z <- y[!duplicated(y$BOLDalt),]

## combine 'remove1' and 'remove2' entries into a vector to exclude from final list:
tmpdrop.list <- unique(c(remove1$BOLDalt, remove2$BOLDalt))
## There was a single exception to the above filter - for samples containing 'BOLD:AAO4373'...
## ...which contained vastly more reads in the Hawaiian dataset than other samples and was likely thus a contaminant into those datasets from Hawaii
tmpx.list <- sub("BOLD:AAO4373", "", tmpdrop.list)
drop.list <- tmpx.list[tmpx.list != ""]

## Now let's finally remove those OTUs...
## First write a little function and remove them from 'tmp.df' object
'%!in%' <- function(x,y)!('%in%'(x,y))
tmpfilt.df <- tmp.df[tmp.df$BOLDalt %!in% drop.list,]
# note we drop from 10,720 observations to 5,564 observations

# ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ 

## Further filtering: 
## 1. Remove all non-arthropod chordates:
tmpfilt2.df <- subset(tmpfilt.df, phylum_name == "Arthropoda")
  ## how many unique OTUs remain?
  length(unique(tmpfilt2.df$OTUid)) ## lots: 1,176
  ## how many unique BOLDid's?
  length(unique(tmpfilt2.df$BOLDalt)) ## 994; interesting degree of redundancy here...

## 2. Among any remaining NTC samples, what OTUs are present and how are the distributed among our true samples?
## subset the data for all 'SampleID' values that match "extBlank" string
toMatch <- c("^NTC")    # make a list of the regex terms to match
matches <- unique (grep(paste(toMatch,collapse="|"),
                        tmpfilt2.df$SampleID, value=TRUE))   # make a vector of every one of those matches fitting the regex above
NTCs.df <- tmpfilt2.df[tmpfilt2.df$SampleID %in% matches,]         # capture all matches of regex in 'tmp.df' object
  # there are 5 total OTUs (out of a 1162 )

## How many total reads are there among all these NTCs?
sum(NTCs.df$CountReads) # there are just 7,485 total reads among all NTCs idenitified (out of 11,782,862 remaining)

## Create an object listing only the OTUs present in out NTCs found in true samples too:
NTC_otu.list <- NTCs.df$OTUid     # generate the list of all unique OTUs present in the NTCs
contam.df <- tmpfilt2.df[tmpfilt2.df$OTUid %in% NTC_otu.list,]    # find all matches from the above list in our original 'tmp.df' data frame object
# we find there are 95 matches; we know 5 of these are from the NTCs themselves
# each of these OTUs seem like real contaminants from the Hawaiian samples themselves - they appear in multiple samples with high read counts
# we'll retain each of these reads, but note that they were identified in a single contaminant sample after filtering

## Save these data:
setwd("~/Repos/guano/OahuBird/data/Routput/")
write.csv(Match_summary, "suspectedContaminants.csv", quote = FALSE)
write.csv(tmp.df, "OahuBird_rawOTUtable.csv", quote = FALSE)
write.csv(tmpfilt2.df, "FilteredOTUs.csv", quote = FALSE, row.names = FALSE)  

## cleanup:
rm(x,y,z,
   nau.df, rut.df, perl.df,
   tiny.nau, tiny.perl, tiny.rut,
   tmp.df, tmpfilt.df, tmpMatch.df, NTCs.df,
   Match_summary, match.df, contam.df, ProjectMatches.df,
   remove1, remove2)

rm(allmatch.list, drop.list, matches, NTC_otu.list, tmpdrop.list, tmpx.list, toMatch)

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
          ######     Part 3 - incorporating metadata     ######
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #

## merge with metadata information
#not run: tmpfilt2.df <- fread('~/Repos/guano/OahuBird/data/Routput/FilteredOTUs.csv', header = TRUE)
meta.df <- fread('https://raw.githubusercontent.com/devonorourke/guano/master/OahuBird/data/OahuBird_metadata.csv', header = TRUE)
meta.df <- meta.df[,c(4:6)]
colnames(meta.df) <- c("SampleID", "DateSampled", "BirdSpecies")
master.df <- merge(tmpfilt2.df, meta.df)
#rm(tmpfilt2.df, meta.df)

# write file to disk:
setwd("~/Repos/guano/OahuBird/data/Routput/")
write.csv(master.df, "master.csv", row.names = F, quote = F)

## Notrun: write.table(master.df, "PHINCHmaster.txt", row.names = F, quote = F, sep = "\t")

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
                 ######     Part 4 - data analyses     ######
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #

detach("package:dplyr", unload=TRUE)
detach("package:data.table", unload=TRUE)
detach("package:reshape2", unload=TRUE)
detach("package:tidyr", unload=TRUE)

##not run: master.df <- fread('https://raw.githubusercontent.com/devonorourke/guano/master/OahuBird/data/Routput/master.csv', header = T)

library(plyr)
setwd("~/Repos/guano/OahuBird/data/Routput/")

## following our filtering, how many samples remain with at least 1 OTU? 2 OTUs? 10 OTUs?
OTUperSample = count(master.df, vars = c("SampleID"))   # There are 99 remaining true samples (all bit one NTC has been filtered out)
sum(OTUperSample$freq > 1)    # There are 90 samples with at least 2 OTUs
sum(OTUperSample$freq > 4)    # There are 72 samples with at least 4 OTUs
sum(OTUperSample$freq > 9)    # There are 51 samples with at least 10 OTUs

## how many observations of OTUs contain complete information (ie. include 'species_name')... a.k.a. species frequency table
speciesOnly.df <- na.omit(master.df)
freq_species <- as.data.frame(table(speciesOnly.df$species_name))     # frequency table of species detected
colnames(freq_species) <- c("species_name", "counts")
write.csv(freq_species, "species_frq_table.csv", row.names = F, quote = F)
sum(freq_species$counts > 1)  # note 142 species identified, but most are not abundant (no OTU ID'd more than 11 samples)

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
