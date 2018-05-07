## purpose: amptk OTU analysis for BRI Pompton project
## written: 28-jan-2018
## edited: 22-apr-2018
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
library(reshape2)
library(tidyr)
library(ggplot2)
library(dplyr)

## read in data:
setwd("~/Repos/guano/BRIpompton/")

h_otutable.df <- read.csv('https://raw.githubusercontent.com/devonorourke/guano/master/BRIpompton/data/amptk/Pompton_h.otu_table.taxonomy.txt', sep = '\t')
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

## save this raw (unfiltered) data frame
setwd("~/Repos/guano/BRIpompton/data/Routput/")
write.csv(tmp.df, "Pompton_rawOTUtable.csv", quote = FALSE, row.names = FALSE)

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
            ######     Part 2 - data filtering     ######
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #

## One way to address potential sources of contamination is to identify the OTUs present in samples across multiple libraries...
## ...with the idea being that if you see an OTU present in a diverse (in terms of geographic sampling origins) set of libraries, it's likely this is a contaminant
## We're going to grab data from two separate projects:
# 1. the `nau.df` project was from a bat guano collected in Central America
# 2. the `oahu.df` project was from bird guano collected in Hawaii
## We don't expect any of the overlapping OTUs to be present in the samples from New Jersey

oahu.df <- read.csv('https://raw.githubusercontent.com/devonorourke/guano/master/OahuBird/data/Routput/master.csv', header = TRUE)
nau.df <- read.csv('https://raw.githubusercontent.com/devonorourke/guano/master/NAUsupp/masterdf.csv', header = TRUE)
colnames(nau.df)[3] <- c("BOLDalt")

## Before we can match we need to append the `BOLDid` vectors to remove the "." delimiter in the `tmp.df` and `rut.df` objects:
tmp.df$BOLDalt <- tmp.df$BOLDid
tmp.df <- separate(tmp.df, col = BOLDalt, into = c("BOLDalt", "delete"), sep = "\\.")
tmp.df$delete <- NULL

oahu.df$BOLDalt <- oahu.df$BOLDid
oahu.df <- separate(oahu.df, col = BOLDalt, into = c("BOLDalt", "delete"), sep = "\\.")
oahu.df$delete <- NULL


## Get a list of all non-redundant `BOLDalt` elements among the two projects:
# because 'none' is not a unique identifier we're removing it from this list
# we're also removing 'CFMR:IM4', 'SINTAX', and 'UTAX' as these were earlier relics of an amptk naming scheme we're not using
string <- unique(c(oahu.df$BOLDalt, nau.df$BOLDalt))
find.list <- list("None", "UTAX", "SINTAX", "CFMR:IM4")
find.string <- paste(unlist(find.list), collapse = "|")
tmpx.list <- gsub(find.string, replacement = "", x = string)
allmatch.list <- tmpx.list[tmpx.list != ""]
rm(string, find.list, find.string, tmpx.list)


## Using that 'allmatch.list' object, query that list against the 'tmp.df$BOLDalt' vector to find common BOLD(alt) id's
ProjectMatches.df <- tmp.df[tmp.df$BOLDalt %in% allmatch.list,]

# how many reads here?
sum(ProjectMatches.df$CountReads)    # 85,184 ... not a massive number, but something to take notice of; in addition, we have 231 observations that match!
# how many unique OTUs are identified?
length(unique(ProjectMatches.df$OTUid))    # 59 unique matches identified... (our tmp.df dataset has 3725 OTUs)
## Finally, use 'ProjectMatches.df' to determine how many reads and how frequently are these matches identified in our tmp.df object
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


## Observation: The majority of these are suspected mock chimeras.
## To keep as conservative an estimate as possible we'll remove each of these from the dataset.

## Let's build a "naughty.list" that would include any BOLDalt value that is detected in the ProjectMatches.df object:

drop.list <- unique(ProjectMatches.df$BOLDalt)

## Now let's finally remove those OTUs...
## First write a little function and remove them from 'tmp.df' object
'%!in%' <- function(x,y)!('%in%'(x,y))
tmpfilt.df <- tmp.df[tmp.df$BOLDalt %!in% drop.list,]
# note we drop from 3,280 observations to 3,049 observations

# ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~

## Further filtering:
## 1. Remove all non-arthropod chordates:
tmpfilt2.df <- subset(tmpfilt.df, phylum_name == "Arthropoda")
  ## how many unique OTUs remain?
  length(unique(tmpfilt2.df$OTUid)) ## lots: 677
  ## how many unique BOLDid's?
  length(unique(tmpfilt2.df$BOLDalt)) ## 543; interesting degree of redundancy here...

## 2. Remove the big contaminators we've seen from previous datasets:
## This is a subjective decision, but from many previous sequencing runs I've noticed that the following inverts are often detected at low levels...
## In this dataset, none of these BOLDID's had a single sample with high read depths but were present in low levels in many samples, indicative of conaminants
othercontams <- c("BOLD:AAH3593", "BOLD:AAO2062", "BOLD:AAC4559", "BOLD:AAC3145", "BOLD:AAE8479", "BOLD:AAG2645", "BOLD:AAA6614")
## These were: 'Chauliodes pectinicornis', 'Chauliodes rastricornis', 'Maccaffertium mediopunctatum', 'Maccaffertium terminatum', 'Maccaffertium mexicanum', 'Hexagenia limbata', 'Psilocorsis reflexella'

## now let's remove these other contaminants:
tmpfilt3.df <- tmpfilt2.df[tmpfilt2.df$BOLDalt %!in% othercontams,]
  ## I'm not a huge fan of this subjective approach; so let's explore one final one...

## 3. Final option - and sort of most "nuclear" option. Look back to the Oahu project contaminant list...
## This list was derived from similar datasets as the filters here: NAU and Oahu, but also included Rutgers...
## The danger is that Rutgers datasets overlap geographically with this Pompton one, and it's possible that the birds and bats are eating similar inverts
## However, after applying this filter, notice that we only reduce about 25% of our overall obserations; this would indicate we have done everything possible to remove contaminants, but haven't lost most of our data
## what about searching through another bigger list from the Oahu project?
bigcontamlist <- read.csv(file = 'https://raw.githubusercontent.com/devonorourke/guano/master/BRIpompton/data/oahubigcontamlist.csv', header = FALSE, stringsAsFactors = FALSE)
biglist <- bigcontamlist$V1
tmpfilt4.df <- tmpfilt2.df[tmpfilt2.df$BOLDalt %!in% biglist,]


## 4. Last item: we did not remve the "Mock" OTUs present in the dataset yet...
## This is problematic because there may be instances in which the OTUs present in that mock community really may be consumed by the birds
## We'll split our dataset into a final filtered library in which the mock remained, versus one in which the mock is removed:

## drop the "_suspect_mock_chimera" tag in the 'OTUid' field:
tmpfilt4.df$OTUid <- gsub("_.*$","",tmpfilt4.df$OTUid)
  # We'll keep the `tmpfilt4.df` object as the data.frame which retains those Mock reads

## drop any rows which match a mock OTU:
tomatch <- c("MockIM")
tmpfilt5.df <- subset(tmpfilt4.df, !grepl(paste(tomatch, collapse= "|"), OTUid))


## Save these data:
setwd("~/Repos/guano/BRIpompton/data/Routput/")
write.csv(Match_summary, "suspectedContaminants.csv", quote = FALSE)
write.csv(tmp.df, "BRIPompton_rawOTUtable.csv", quote = FALSE)
write.csv(tmpfilt4.df, "BRIPompton_FullFilteredOTUtable_withMock.csv", quote = FALSE, row.names = FALSE)
write.csv(tmpfilt5.df, "BRIPompton_FullFilteredOTUtable_noMock.csv", quote = FALSE, row.names = FALSE)

## cleanup:
rm(x,y,z,
   nau.df, oahu.df,
   tmp.df, tmpfilt.df, tmpfilt2.df, tmpfilt3.df, tmpfilt4.df, tmpMatch.df, NTCs.df, bigcontamlist,
   Match_summary, match.df, contam.df, ProjectMatches.df, allmatch.list,
   remove1, remove2)

rm(allmatch.list, drop.list, matches, NTC_otu.list, tmpdrop.list, tmpx.list, tomatch, biglist, othercontams)

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
          ######     Part 3 - incorporating metadata     ######
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #

## merge with metadata information
#not run: tmpfilt5.df <- read.csv('https://raw.githubusercontent.com/devonorourke/guano/master/BRIpompton/data/Routput/BRIPompton_FullFilteredOTUtable_noMock.csv', header = TRUE)
meta.df <- read.csv('https://raw.githubusercontent.com/devonorourke/guano/master/BRIpompton/data/BRIPomptonMeta.tsv',
                    sep="\t", header = TRUE)
meta.df$Alias <- gsub("-","\\.",meta.df$Alias)
colnames(meta.df) <- c("SampleID", "BandNum", "SampleDate", "SiteName", "Territory", "RecapStatus", "BirdSpecies", "BirdAge", "BirdSex")

master.df <- merge(tmpfilt5.df, meta.df)
rm(tmpfilt5.df, meta.df)

# write file to disk:
setwd("~/Repos/guano/BRIpompton/data/Routput/")
write.csv(master.df, "master.csv", row.names = F, quote = F)

## Notrun: write.table(master.df, "PHINCHmaster.txt", row.names = F, quote = F, sep = "\t")

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
                 ######     Part 4 - data analyses     ######
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #

detach("package:dplyr", unload=TRUE)
detach("package:reshape2", unload=TRUE)
detach("package:tidyr", unload=TRUE)

##not run:
master.df <- read.csv('https://raw.githubusercontent.com/devonorourke/guano/master/BRIpompton/data/Routput/master.csv', header = TRUE)

library(plyr)
setwd("~/Repos/guano/BRIpompton/data/Routput/")

## following our filtering, how many samples remain with at least 1 OTU? 2 OTUs? 10 OTUs?
OTUperSample = count(master.df, vars = c("SampleID"))   # There are 95 true samples remaining from our initial 99 submitted
sum(OTUperSample$freq > 1)    # There are 89 samples with at least 2 OTUs
sum(OTUperSample$freq > 4)    # There are 81 samples with at least 4 OTUs
sum(OTUperSample$freq > 9)    # There are 62 samples with at least 10 OTUs

## how many observations of OTUs contain complete information (ie. include 'species_name')... a.k.a. species frequency table
speciesOnly.df <- na.omit(master.df)    ## 575 observations remain (about 1/3 of our dataset is named to the species level...)
freq_species <- as.data.frame(table(speciesOnly.df$species_name))     # frequency table of species detected... so we see there are 173 unique species named
colnames(freq_species) <- c("species_name", "counts")
sum(freq_species$counts > 1)  # note 156 species identified, but most are not abundant (only 15 species ID'd in >= 10 samples)
write.csv(freq_species, "species_frq_table.csv", row.names = F, quote = F)


## what is the frequency of each OTU?
OTUcounts = count(master.df, vars = c("OTUid"))
colnames(OTUcounts) <- c("OTUid", "NumberOfDetections")
write.csv(OTUcounts, "OTUcounts.csv", row.names = F, quote = F)
  ## most abundant OTU (OTU361) present in 35 samples - can't read too much into this because it's never detected in high amounts
  ## probably the most frequent high read winner is OTU20 - a gypsy moth - as this has many samples with high numbers of reads detected

## how many OTUs are called per site?
OTUperSite = count(master.df, vars = c("SiteName"))
colnames(OTUperSite) <- c("SiteName", "NumberOfDetections")
write.csv(OTUperSite, "OTU_per_Site.csv", row.names = F, quote = F)
  ## decent balance of detections per site

## how many OTUs are called per Territory?
OTUperTerritory = count(master.df, vars = c("Territory"))
colnames(OTUperSite) <- c("Territory", "NumberOfDetections")
write.csv(OTUperSite, "OTU_per_Territory.csv", row.names = F, quote = F)
  ## as with SiteName, decent number of observations per Territory to make comparisons across all elements of this variable

## what about a few multivariate questions...
## what OTUs do we observe per Site per Species?
OTUperSiteperSpecies = count(master.df, vars = c("SiteName", "BirdSpecies"))
colnames(OTUperSiteperSpecies) <- c("SiteName", "BirdSpecies", "NumberOfDetections")
write.csv(OTUperSiteperSpecies, "OTU_per_Site_and_Species.csv", row.names = F, quote = F)

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
              ######     Part 5a - taxa sampled viz     ######
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #

setwd("~/Repos/guano/OahuBird/data/Routput/")
master.df <- read.csv('https://raw.githubusercontent.com/devonorourke/guano/master/BRIpompton/data/Routput/master.csv', header = TRUE)
plot.df <- master.df
rm(master.df)

## First plot will be to generate a couple of stacked bar plots with frequency of detections by Taxonomic Order per:
## 1. Bird Species
## 2. Site
## 3. Bird Species, partitioned (faceted) by Site

## It's worth generating a table to examine the frequency of detection for each bird per Taxonomic Order, grouped with Taxonomic Class:
## We wont' use this in the plotting directly, but it's good to refer to when figuring out how to lay out the figure
library(plyr)
detach("package:dplyr", unload=TRUE)
freqOrders <- count(plot.df, vars = c("class_name", "order_name"))

## to figure out a coloring scheme for this legend, let's figure out how many taxonomic orders there are in our dataset
## note we could choose to exclude any value missing information at this Taxonomic level (it appears as "NA")...
#  but I think it's better to show the viewer how much of our information is missing

## generate a color palette where we group each taxonomic Class a single color, and use a variety of hues to discriminate among the orders within those classes:
length(unique(plot.df$class_name))  # we have 8 unique colors to make
length(unique(plot.df$order_name))  # we have 33 unique orders, including the one 'NA' value associated to the "insecta"

## What's the order ggplot will print things? It's automatically alphabetical...
## don't forget that it will also include "NA" at the end (not within the alphabetical order)
sort(unique(plot.df$order_name))

## so we can make a color palette that is grouped by color, then hue'd by taxonomic Class:
## helpful site: https://www.w3schools.com/colors/colors_picker.asp

## what are all the unique class/order taxonomic level pairings?
bcolours <- unique(plot.df[c("class_name", "order_name")])
bcolours <- bcolours[with(bcolours, order(class_name)), ]
  ## brutal. so we have a lot of different numbers of hues here...
  ## we're going to assigne the following color spectrum into 'http://tools.medialab.sciences-po.fr/iwanthue/' to get what we need:

## in that `bcolours` object you'll see there are:
##    (6)  Arachnida - red
##    (1)  Branchiopoda - slate blue
##    (3)  Collembola - green
##    (2)  Diplopoda - orange/yellow
##    (1)  Hexanauplia - light gray
##    (15) Insecta - blue/violet
##    (3)  Malacostraca - brown
##    (2)  Maxillopoda - dark grey

mypal <- c("#e6ccb3", "#ff4d4d", "#732e74", "#C0C0C0", 
           "#4fc2e4", "#DCDCDC", "#cc9966", "#000000", 
           "#c743df", "#79d279", "#5d94d1", "#708090", 
           "#6739ce", "#bea0e4", "#86592d", "#ffa31a", 
           "#462287", "#e18cdb", "#990000", "#332460", 
           "#d968db", "#51589c", "#206020", "#ffd633", 
           "#a232a6", "#4d0000", "#ff9999", "#ffe6e6", 
           "#d9f2d9", "#627ee0", "#a664ad", "#7462da", "#777777")


## if you want to make a similar plot that doesn't plot NA values, use this line instead
## b <- ggplot(data = (plot.df[!is.na(plot.df$order_name),]),

b <- ggplot(data = plot.df,
            aes(BirdSpecies, ..count.., fill = order_name)) +
  geom_bar() +
  scale_fill_manual(values = mypal) +
  labs(title = "Relative taxonomic Orders detected by avian guano - BRI Pompton project",
       subtitle = "Taxonomic order defined using approaches described in 'amptk' bioinformatic pipeline\n using Barcode of Life Database references",
       x = "Bird Species",
       y = "Number of detections") +
  guides(fill = guide_legend(title = "Taxonomic Order", ncol = 2)) +
  theme(legend.position = "right", axis.text.y=element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))
b
ggsave(filename = "TaxOrder_byBirdSpecies.png", 
       b,
       device = "png", 
       path = "~/Repos/guano/BRIpompton/data/Routput/plots/")


s <- ggplot(data = plot.df,
            aes(SiteName, ..count.., fill = order_name)) +
  geom_bar() +
  scale_fill_manual(values = mypal) +
  labs(title = "Relative taxonomic Orders detected by avian guano - BRI Pompton project",
       subtitle = "Taxonomic order defined using approaches described in 'amptk' bioinformatic pipeline\n using Barcode of Life Database references",
       x = "Site Name",
       y = "Number of detections") +
  guides(fill = guide_legend(title = "Taxonomic Order", ncol = 2)) +
  theme(legend.position = "right", axis.text.y=element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))
s

ggsave(filename = "TaxOrder_bySite.png", 
       s,
       device = "png", 
       path = "~/Repos/guano/BRIpompton/data/Routput/plots/")

bs <- ggplot(data = plot.df,
                aes(BirdSpecies, ..count.., fill = order_name)) +
  geom_bar() +
  scale_fill_manual(values = mypal) +
  labs(title = "Relative taxonomic Orders detected by avian guano - BRI Pompton project",
       subtitle = "Taxonomic order defined using approaches described in 'amptk' bioinformatic pipeline\n using Barcode of Life Database references",
       x = "Bird Species",
       y = "Number of detections") +
  guides(fill = guide_legend(title = "Taxonomic Order", ncol = 2)) +
  theme(legend.position = "right", axis.text.y=element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~ SiteName)
bs

ggsave(filename = "TaxOrder_bySite_andBirdSpecies.png", 
       bs,
       device = "png", 
       path = "~/Repos/guano/BRIpompton/data/Routput/plots/")
