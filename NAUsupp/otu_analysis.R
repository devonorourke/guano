#amptk_2017_OTUtable_analysis.R
#written 10.12.16
#modified 4.12.17
#Devon O'Rourke
#amptk analysis for NAU samples from various regions

########################################################################################
############################    SECTION 0: metadata files  #############################
########################################################################################

#################### process started: creating metadata file  ####################
date()

#import collection data
setwd("~/Documents/Lab.Foster/guano/")
library(readr)
coldat16 <- read_csv("~/Documents/Lab.Foster/guano/collectiondata2016.csv",
                               col_types = cols(PlateNumber = col_character()))

#rename and combine first two fields to create a "tagID" field
#rename sampleID values so they are all 4-digit values (for the traditional 'xxxx' numbering scheme)...
coldat16$SampleID <- sprintf("%04s", coldat16$SampleID)  
#this leaves the non-traditional variants unchanged (ie. "01B02" is unchanged because it's more than four digits)

#create new 'tagID' field by concantenating project name and sampleID values
coldat16$tagID <- paste(coldat16$StudyID, ".", coldat16$SampleID, sep = "")

#drop the StudyID and SampleID fields
coldat16$StudyID <- NULL
coldat16$SampleID <- NULL

#reorder the columns to put this new 'tagID' field near the front
#how many total fields are there?
length(names(coldat16))
#in this example there are 14, so you want to move the last field ('tagID') to the first position
#reorder data.frame columns:
coldat16 <- coldat16[,c(14,1:13)]

#relabel the first column to "SampleID"
colnames(coldat16)[1] <- "SampleID"

#look at just the nau data
naudat = coldat16[grepl("nau", coldat16$SampleID),]

#rename that one EUMA sample so it matches how it's named in the sequencing dataset:
naudat$SampleID <- gsub('_', '', naudat$SampleID)

#import list of SampleID names completed through AMPtk pipeline for NAU samples data
setwd("~/Documents/Lab.Foster/guano/guano_nau/all/")
NAUnames <- read_delim("~/Documents/Lab.Foster/guano/guano_nau/all/nauNames.txt", 
                       "\t", escape_double = FALSE, col_names = FALSE, 
                       trim_ws = TRUE)
colnames(NAUnames)[1] <- "SampleID"
#fix the same EUMA name once more
NAUnames$SampleID <- gsub('_', '', NAUnames$SampleID)

#merge so that only samples listed in the 'NAUnames' data.frame are pulled from the 'naudat' data.frame
onlyseqd <- merge(NAUnames, naudat, by = "SampleID", all.x = T)
  #noticed that three samples had incomplete entries; this was because they were from DNA extracted from NAU instead of me, so there was no record in the CollectionData.csv form. I'm going to export this 'onlyseqd' file and update it with the three samples metadata.
setwd("~/Documents/Lab.Foster/guano/guano_nau/all/")
write.csv(onlyseqd, file = "updateme.csv", row.names = F, quote = F)
  #entered relevant metadata to this file and will now update that file as a different name 'temp_meta.csv' data.frame to replace old version... this is the metadata file you're going to use with a little modification:
library(readr)
meta_df <- read_csv("~/Documents/Lab.Foster/guano/guano_nau/all/temp_meta.csv")

#check for unique records by SampleID:
#the unique ones:
nrow(meta_df[!duplicated(meta_df$SampleID),])
#or to store as data.frame:
#meta_df_u <- (meta_df[!duplicated(meta_df$SampleID),])
#the dupilcates:
nrow(meta_df[duplicated(meta_df$SampleID),])
#or to store as data.frame:
#meta_df_d <- meta_df[duplicated(meta_df$SampleID),]
#if all clear, just remove those two files, otherwise you have some work to do to remove what you need to:
#rm(meta_df_d, meta_df_u)

#remove fields that are unnecessary for this project:
meta_df$WeekOfYear <- NULL
meta_df$Status <- NULL
meta_df$PlateNumber <- NULL
meta_df$WellNumber <- NULL
meta_df$Age <- NULL
meta_df$ReproStat <- NULL
meta_df$Roost <- NULL
meta_df$Notes <- NULL

#substitute 'NA' cells for certain terms in specific columns:
meta_df[is.na(meta_df)] <- "noRecord"

#rename fields as needed:
colnames(meta_df)[c(2,3)] <- c("Site", "Date")

#make a few frequency tables of various meta_dfdata to ensure you have the fields/values labeled correctly:
library(plyr)
Site.table <- count(meta_df, Site)
  #good. nothing out of the ordinary.
Date.table <- count(meta_df, Date)
  #dates are as exptected; single value with no record to perhaps revisit
BatSpecies.table <- count(meta_df, BatSpecies)
  #they have a lot of bat species. might want to check on that single 'EPBR?' value...
Sex.Table <- count(meta_df, Sex)
  #may want to find a way to update records for the 5 missing instances

#the above tables suggest that our meta_dfdata doesn't appear to need further modifications or updates.

#final cleanup
rm(coldat16, naudat, NAUnames, onlyseqd)
rm(BatSpecies.table, Date.table, Sex.Table, Site.table)

#write meta_dfdata to disk for aca.samples
setwd("~/Documents/Lab.Foster/guano/guano_nau/all/R_OTUcall_output/")
write.csv(meta_df, file = "metadf.csv", row.names = F, quote = F)

#################### process ended: file "meta_df" created  ####################

########################################################################################
############################    SECTION 1: data wrangline  #############################
########################################################################################

#################### process started: OTU table import and renaming  ####################
date()

#import packages
library(tidyr)
library(reshape2)
library(data.table)

#set working directory
setwd("~/Documents/Lab.Foster/guano/guano_nau/all/")

#import taxonomy-classified OTU Table:
amptk_otu_tableA <- read.delim("nauall.otu_table.taxonomy.txt", stringsAsFactors = F)
#rename that first column
colnames(amptk_otu_tableA)[1] <- "OTUid"
#rename the column that had the inappropriate name with the 'euma' sample:
colnames(amptk_otu_tableA)[60] <- "nau.euma1bandersonNM"

#transform data
amptk_otu_tableA[amptk_otu_tableA == 0 ] <- NA
amptk_otu_tableA$OTUnTax <- paste(amptk_otu_tableA$OTUid, amptk_otu_tableA$Taxonomy, sep = ";")
amptk_otu_tableA$OTUid <- NULL
amptk_otu_tableA$Taxonomy <- NULL
amptk_dfA <- melt(amptk_otu_tableA, id = "OTUnTax")
colnames(amptk_dfA) <- c("OTUnTax", "tagID", "binaryPresence")
amptk_dfA <- subset(amptk_dfA, binaryPresence >= 1)
amptk_dfA$binaryPresence <- NULL

#separate 'Taxonomy' field into respective taxa levels:
amptk_dfA <- separate(data = amptk_dfA, col = OTUnTax, into = c("OTUid", "OTU_identifier", "kingdom_name", "phylum_name", "class_name", "order_name", "family_name", "genus_name", "species_name"), sep = "\\;|,")


#remove the prefixes leading into the taxa names:
amptk_dfA$kingdom_name <- sub("k:", "", amptk_dfA$kingdom_name)
amptk_dfA$phylum_name <- sub("p:", "", amptk_dfA$phylum_name)
amptk_dfA$class_name <- sub("c:", "", amptk_dfA$class_name)
amptk_dfA$order_name <- sub("o:", "", amptk_dfA$order_name)
amptk_dfA$family_name <- sub("f:", "", amptk_dfA$family_name)
amptk_dfA$genus_name <- sub("g:", "", amptk_dfA$genus_name)
amptk_dfA$species_name <- sub("s:", "", amptk_dfA$species_name)


#add in a designated field to show which data.frame belongs to which library (oro##-1 or oro##-2)
#amptk_dfA$library <- "library01"

#make list of data frames to merge
l = list(amptk_dfA)

#use 'rbindlist' function in 'data.table' package to join two data frames
amptk_dfA <- rbindlist(l, use.names = TRUE, fill = TRUE)

#remove redundant calls
#this will cut out a lot of instances in which you have poorly described taxonomies for a sample using UTAX as the 'OTU_identifier'
#as a precaution, what you can say is that you have confidence in eliminating reads in which we can't associate it to a unique 'phlya',...
#so what's the point in keeping it anyway?
#it will keep a single instance in which some value is described at just the phylum level, yet this may actually be cutting out things like our bat reads
#so we'll keep these reads that are cut in the "dt_dupilc" data frame; you can then use the OTUId from each of these and go back and search for those in the 
#original fasta file, and see if any of them are things like bat DNA sequences

#using 'data.table'
setkey(amptk_dfA, OTU_identifier, tagID, species_name, genus_name, family_name, order_name, class_name)
amptk_dt_unique <- amptk_dfA[!duplicated(amptk_dfA)]   
#we're pulling out repeated calls just one situation with this dataset:
##1. you've named something where there are multiple unique OTUs assigned to an identical taxonomic description; there isn't a point ...
## in double counting these because we're collapsing a taxonomic description to a binary presence/absence (without this filter you...
## end up adding duplicate counts for a single species per sample)
amptk_dt_duplic <- amptk_dfA[duplicated(amptk_dfA)]    
#these represent redundant calls explained above. in some datasets (not this one) these can also represent situations where there may be MULTIPLE DISTINCT...
## taxa you could assign IF the taxonomic description is only classified through to a "phylum" level (and "class", "order", .... "species" are all "NA")...
## you could look into identifying if these are truly different, but if they are so poorly classified to begin with it's likely you won't have enough...
## percent identity in an alignment score to definitively improve the calls. the only way it'll improve is if your existing databases just missed the read...
## You can look at the unique calls in which there is no 'phylum_name' known to the database using the following (comment muted) commands:
##nophyla_df <- amptk_dt_unique[is.na(amptk_dt_unique$phylum_name),]
##setkey(nophyla_df, OTUid)
##nophyla_df_unique <- nophyla_df[!duplicated(nophyla_df)]   
##these represent all the unique per OTUId between both libraries; note that the same call may occur between the two libraries but...
##contain different OTUId's because of how they were independently assigned from the amptk pipeline

#CLEANUP
rm(amptk_dfA, amptk_otu_tableA, l, amptk_dt_duplic)

#################### process complete: OTU table import and renaming  ####################

#################################################################################################
#####################   SECTION 2:  Merge Metadata and OTUtable     #############################
#################################################################################################

#################### process started: OTU table import and metadata merging  ####################
date()

#first create a list of the OTU's used in this study - just make a new dataframe from the 'tagID' column
#then search through the 'collection_data' file in Google Drive and pull out the samples that match
#that new file, called "{something}...metadata.txt" is what we're importing to start here...

#import metadata
setwd("~/Documents/Lab.Foster/guano/guano_nau/all/R_OTUcall_output/")
meta_df = read.csv(file = "metadf.csv", header = T, stringsAsFactors = F)

#rename that first column
colnames(meta_df)[1] <- "tagID"

#replace the hyphen with a period in the "tagID" column (may not be necessary - check both data.frames for consistency in delimiter)
#metadata$tagID <- sub("-", "\\.", metadata$tagID)

#merge OTU calls with metadata by 'tagID' (really just the sample ID)
masterdf = merge(amptk_dt_unique, meta_df, by = "tagID")

#write to disk:
setwd("~/Documents/Lab.Foster/guano/guano_nau/all/R_OTUcall_output/")
write.csv(masterdf, file = "masterdf.csv", row.names = F, quote = F)

#to keep things conservative for this first go of things, we're going to limit our observations excusively to BOLD-called OTUs (we may improve upon their description though); we are not going to include the SINTAX and UTAX terms here, though these may play a role in subsequent analyses. These are going to be removed at this point.

#this gets rid of the SINTAX and UTAX assigned samples
BOLDonly.df = masterdf[!grepl("SINTAX", masterdf$OTU_identifier),]
  #drops to 363 observations (you lose about half of all calls)
BOLDonly.df = BOLDonly.df[!grepl("UTAX", BOLDonly.df$OTU_identifier),]
  #drops to 351; not nearly as many UTAX assignments as SINTAX

#and this gets rid of the samples that are likely the result of the mock community bleeding into the true samples (another index-bleed happenstance)
BOLDonly.df = BOLDonly.df[!grepl("CFMR", BOLDonly.df$OTU_identifier),]
  #just drops one more observation, so this wasn't particularly problematic.

#rename a few columns so we're explicit about what species is being described
colnames(BOLDonly.df)[c(1, 12, 13)] <- c("SampleID", "Date", "BatTagID")

#it's worth renaming empty (N/A) values for taxonomic classes with "unknown" for future plotting purposes (we'll want to show in our data how many OTUs weren't explicitly described relative to those that were)
BOLDonly.df[is.na(BOLDonly.df)] <- "unknown"
  #careful with this term, as it will convert any 'N/A' in the dataset, regardless of what field it belongs to, to "unknown"

#we're also going to make a separate file for the bat-specific reads; a quick analysis suggests these are working alright when we compare the known bats to the predicted names, but they're by no means perfect. for example:
#1. there are only 26 total calls out of 59 potenteial samples with bat calls
#2. of those that worked many are defined as Mollosus despite the bat being identified by a researcher as something else
#3. in one case there are two bat species associated with a single sample
batOTUs.df = BOLDonly.df[grepl("Chiroptera", BOLDonly.df$order_name),]

#now remove the Bat OTUs from the 'BOLDonly.df' data.frame
BOLDonly.df = BOLDonly.df[!grepl("Chiroptera", BOLDonly.df$order_name),]

#write object 'masterdf' to disk
setwd("~/Documents/Lab.Foster/guano/guano_nau/all/R_OTUcall_output/")
write.csv(BOLDonly.df, file = "BOLDonly.df.csv", row.names = F, quote = F)

#write the 'batOTUs.df' file as well:
write.csv(batOTUs.df, file = "batOTUs.df.csv", row.names = F, quote = F)

#CLEANUP
rm(meta_df, amptk_dt_unique, masterdf, batOTUs.df)  

#################### process ended: OTU table import and metadata merging  ####################
############### BOLDonly.df file contains all information for subsequent work  ################

##################################################################################################
###################     SECTION 3: Calculate Frequencies for plots           #####################
##################################################################################################

################## process started: data tables created for plots and analysis  ##################

rm(BOLDonly.df)

#load file if necessary:
setwd("~/Documents/Lab.Foster/guano/guano_nau/all/R_OTUcall_output/")
library(readr)
BOLDonly.df <- read_csv("~/Documents/Lab.Foster/guano/guano_nau/all/R_OTUcall_output/BOLDonly.df.csv", 
                        col_types = cols(BatTagID = col_character(), 
                                         Sex = col_character()))

#need to generate a suite of data tables for each unique plotting scenario: 
## 1. you want to know the relative frequency (in percentage terms) of each 'family_name'(or some other taxonomic unit) when you aggregate all counts for:
##  a. All bats at each unique WOY or Site, or
##  b. All WOY or Sites for each unique BatSpecies
## 2. What'll likely get used more for analyses: the relative frequency of each 'family_name' identified when aggregating counts for:
### a. each unique BatSpecies per Site (aggregating all dates, even across years)
#In addition, we're going to generate a few other frequency tables that may prove useful:
## 3a. a frequency table showing how many total OTUs were identified (regardless of taxonomic class) for each BatSpecies at each Site
## 3b. a frequency table showing how many unique OTUs were identified (by taxonomic 'family') for each sample, for each BatSpecies at each Site


######## START Part 1 ########
#first, make a table with the number of counts per 'family_name' for all our metadata (Site, WOY, BatSpecies, Sex, Age, and ReproStat) that's unique:
library(plyr)
Familycounts_table <- ddply(BOLDonly.df, .(BOLDonly.df$Site, BOLDonly.df$Date, BOLDonly.df$BatSpecies, BOLDonly.df$Sex, BOLDonly.df$family_name), nrow)
names(Familycounts_table) <- c("Site", "Date", "BatSpecies", "Sex", "ArthropodFamily", "RelativeCounts")

#next, repeat for 'order_name' across each BatSpecies and WOY that's unique
Ordercounts_table <- ddply(BOLDonly.df, .(BOLDonly.df$Site, BOLDonly.df$Date, BOLDonly.df$BatSpecies, BOLDonly.df$Sex, BOLDonly.df$order_name), nrow)
names(Ordercounts_table) <- c("Site", "Date", "BatSpecies", "Sex", "ArthropodOrder", "RelativeCounts")

#write to disk:
setwd(dir = "~/Documents/Lab.Foster/guano/guano_nau/all/R_OTUcall_output/")
write.csv(Familycounts_table, file="Familycounts_table.csv", quote = F, row.names = F)
write.csv(Ordercounts_table, file="Ordercounts_table.csv", quote = F, row.names = F)

#next, sum up the frequency values given a unique Site (aggregating across all dates samples were collected):
##1. For each unique Site, grouping all BatSpecies (ie. how many times did we see a family/order of arthropod get eaten at a given Site, regardless of BatSpecies?)
library(dplyr)
FamilyCount.bySite <- Familycounts_table %>% 
  group_by(Site, ArthropodFamily) %>% 
  summarise(Freq = sum(RelativeCounts))
names(FamilyCount.bySite) <- c("Site", "ArthropodFamily", "RelativeCounts")

TotalFamilyCount.perSite <- Familycounts_table %>% 
  group_by(Site) %>% 
  summarise(Freq = sum(RelativeCounts))
names(TotalFamilyCount.perSite) <- c("Site", "TotalCounts")

#and repeat for the "order" table
OrderCount.bySite <- Ordercounts_table %>% 
  group_by(Site, ArthropodOrder) %>% 
  summarise(Freq = sum(RelativeCounts))
names(OrderCount.bySite) <- c("Site", "ArthropodOrder", "RelativeCounts")

TotalOrderCount.perSite <- Ordercounts_table %>% 
  group_by(Site) %>% 
  summarise(Freq = sum(RelativeCounts))
names(TotalOrderCount.perSite) <- c("Site", "TotalCounts")


#finally, merge the two tables (repeat for each taxonimic group) by "Site"
##first, for family
Taxa.perSite.byFamily <- merge(FamilyCount.bySite, TotalFamilyCount.perSite, by = "Site")
##then, for order
Taxa.perSite.byOrder <- merge(OrderCount.bySite, TotalOrderCount.perSite, by = "Site")

#calculate the relative frequency per Arthropod by taxa level per Site
##first, for family
Taxa.perSite.byFamily$percent_taxa <- (Taxa.perSite.byFamily$RelativeCounts / Taxa.perSite.byFamily$TotalCounts * 100)
##then, for order
Taxa.perSite.byOrder$percent_taxa <- (Taxa.perSite.byOrder$RelativeCounts / Taxa.perSite.byOrder$TotalCounts * 100)

###these final two tables (by taxoomic Family and Order) allow you to investigate the relative proportions of each unique taxa; because some Families or Orders are going to be far more represented DUE TO SAMPLING AT ONE Site MORE THAN ANOTHER (not because of diet preferences), this table reflects the RELATIVE proportions of a given taxonomic Family or Order at each Site so as to reduce the sampling skew.

#write these tables to disk:
getwd()
setwd("/Users/devonorourke/Documents/Lab.Foster/guano/guano_nau/all/R_OTUcall_output/")
write.csv(Taxa.perSite.byFamily, "Taxa.perSite.byFamily.csv", quote = F, row.names = F)  
write.csv(Taxa.perSite.byOrder, "Taxa.perSite.byOrder.csv", quote = F, row.names = F)  
#cleanup
rm(FamilyCount.bySite, TotalFamilyCount.perSite, OrderCount.bySite, TotalOrderCount.perSite)

##1b. For each unique BatSpecies, grouping all Sites (ie. how many times did we see an taxonomic group of arthropods get eaten by a particular BatSpecies, regardless of Site?). The work follows the same format as above but replaces the "Site" with "BatSpecies" factor where appropriate.

#first by Family
FamilyCount.byBatSpecies <- Familycounts_table %>% 
  group_by(BatSpecies, ArthropodFamily) %>% 
  summarise(Freq = sum(RelativeCounts))
names(FamilyCount.byBatSpecies) <- c("BatSpecies", "ArthropodFamily", "RelativeCounts")

TotalFamilyCount.perBatSpecies <- Familycounts_table %>% 
  group_by(BatSpecies) %>% 
  summarise(Freq = sum(RelativeCounts))
names(TotalFamilyCount.perBatSpecies) <- c("BatSpecies", "TotalCounts")

#then by Order
OrderCount.byBatSpecies <- Ordercounts_table %>% 
  group_by(BatSpecies, ArthropodOrder) %>% 
  summarise(Freq = sum(RelativeCounts))
names(OrderCount.byBatSpecies) <- c("BatSpecies", "ArthropodOrder", "RelativeCounts")

TotalOrderCount.perBatSpecies <- Ordercounts_table %>% 
  group_by(BatSpecies) %>% 
  summarise(Freq = sum(RelativeCounts))
names(TotalOrderCount.perBatSpecies) <- c("BatSpecies", "TotalCounts")


#finally, merge the two tables (repeat for each taxonimic group) by "BatSpecies"
##first, for family
Taxa.perBatSpecies.byFamily <- merge(FamilyCount.byBatSpecies, TotalFamilyCount.perBatSpecies, by = "BatSpecies")
##then, for order
Taxa.perBatSpecies.byOrder <- merge(OrderCount.byBatSpecies, TotalOrderCount.perBatSpecies, by = "BatSpecies")

#calculate the relative frequency per Arthropod by taxa level per Site
##first, for family
Taxa.perBatSpecies.byFamily$percent_taxa <- (Taxa.perBatSpecies.byFamily$RelativeCounts / Taxa.perBatSpecies.byFamily$TotalCounts * 100)
##then, for order
Taxa.perBatSpecies.byOrder$percent_taxa <- (Taxa.perBatSpecies.byOrder$RelativeCounts / Taxa.perBatSpecies.byOrder$TotalCounts * 100)

###these final two tables (by taxoomic Family and Order) allow you to investigate the relative proportions of each unique taxa; because some Families or Orders are going to be far more represented DUE TO SAMPLING OF MORE OF ONE BatSpecies THAN ANOTHER (not because of diet preferences), this table reflects the RELATIVE proportions of a given taxonomic Family or Order for each BatSpecies so as to reduces the sampling skew.

#write to disk:
setwd("/Users/devonorourke/Documents/Lab.Foster/guano/guano_nau/all/R_OTUcall_output/")
write.csv(Taxa.perBatSpecies.byFamily, "Taxa.perBatSpecies.byFamily.csv", quote = F, row.names = F)
write.csv(Taxa.perBatSpecies.byOrder, "Taxa.perBatSpecies.byOrder.csv", quote = F, row.names = F)
#cleanup:
rm(FamilyCount.byBatSpecies, TotalFamilyCount.perBatSpecies, OrderCount.byBatSpecies, TotalOrderCount.perBatSpecies)
######## END Part 1 ########


######## START Part 2 ########
#For each unique BatSpecies at each unique Site (for both taxonomic levels of arthropods)

### Part 2 ###
#2. For Sote and BatSpecies: the relative counts for each unique pair of BatSpecies and Sote classes:
#first for Family:
FamilyCount.bySiteandBatSpecies <- Familycounts_table %>% 
  group_by(BatSpecies, Site, ArthropodFamily) %>% 
  summarise(Freq = sum(RelativeCounts))
names(FamilyCount.bySiteandBatSpecies) <- c("BatSpecies", "Site", "ArthropodFamily", "RelativeCounts")
#then for Order:
OrderCount.bySiteandBatSpecies <- Ordercounts_table %>% 
  group_by(BatSpecies, Site, ArthropodOrder) %>% 
  summarise(Freq = sum(RelativeCounts))
names(OrderCount.bySiteandBatSpecies) <- c("BatSpecies", "Site", "ArthropodOrder", "RelativeCounts")

##Next, tally the total counts of BatSpecies split by Site. It's the same regardless of Family/Order taxa level:  
SiteandBatSpecies.TotalCounts <- Familycounts_table %>% 
  group_by(BatSpecies, Site) %>% 
  summarise(Freq = sum(RelativeCounts))
names(SiteandBatSpecies.TotalCounts) <- c("BatSpecies", "Site", "TotalCounts")

#merge the original '{Family/Order}counts_table' object with the "SiteandBatSpecies.TotalCounts" object with R package 'dplyr'
library(dplyr)
#for Family taxa:
percFamily.bySiteandBatSpecies <- FamilyCount.bySiteandBatSpecies %>% inner_join(SiteandBatSpecies.TotalCounts)
#and calculate relative frequency per insect Family per Site AND BatSpecies
percFamily.bySiteandBatSpecies$percent_taxa <- (percFamily.bySiteandBatSpecies$RelativeCounts / percFamily.bySiteandBatSpecies$TotalCounts * 100)

#for Order taxa:
percOrder.bySiteandBatSpecies <- OrderCount.bySiteandBatSpecies %>% inner_join(SiteandBatSpecies.TotalCounts)
percOrder.bySiteandBatSpecies$percent_taxa <- (percOrder.bySiteandBatSpecies$RelativeCounts / percOrder.bySiteandBatSpecies$TotalCounts * 100)

#These tables represent the relative proportions for a Family/Order that was detected for a given BatSpecies at a specific Site.

#write tables to disk:
setwd("~/Documents/Lab.Foster/guano/guano_nau/all/R_OTUcall_output/")
write.csv(percFamily.bySiteandBatSpecies, "percFamily.bySiteandBatSpecies.csv", quote = F, row.names = F)
write.csv(percOrder.bySiteandBatSpecies, "percOrder.bySiteandBatSpecies.csv", quote = F, row.names = F)  
#these tables partition Site and BatSpecies factors into unique classes, thus separating out the total and relative counts into smaller proportions while being able to attribute both Site and BatSpecies differences for a given Family or Order abundance.

#cleanup:
rm(FamilyCount.bySiteandBatSpecies, OrderCount.bySiteandBatSpecies, SiteandBatSpecies.TotalCounts)
##########    END Part 2B     ##########

######## START Part 3 ########
##3a. this calculates the number of times an OTU (at any taxonomic level) was detected for each BatSpecies at each Site.
#note you may want to alter this so that you only include those fully identified through to some specific level (say at least "Order" or "Family") so that you don't include many calls where the OTU could be described only to a Phylum level (that is, you only know that it's an arthropod and can't even distinguish if it's an insect or spider according to this initial search with our BOLD database)

#3a - for Site
library(plyr)
OTUcounts.perBatAndSite <- ddply(BOLDonly.df, .(BOLDonly.df$BatSpecies, BOLDonly.df$Site), nrow)
names(OTUcounts.perBatAndSite) <- c("BatSpecies", "Site", "Counts")

##3b. this calculates the number of OTUs called per sample of guano while also retaining the BatSpecies and Site information:
library(plyr)
#for Site
PerSampleOTUcounts.bySiteandBatSpecies <- ddply(BOLDonly.df, .(BOLDonly.df$SampleID, BOLDonly.df$BatSpecies, BOLDonly.df$Site), nrow)
names(PerSampleOTUcounts.bySiteandBatSpecies) <- c("SampleID", "BatSpecies", "Site", "Counts")
#what's the mean?
mean(PerSampleOTUcounts.bySiteandBatSpecies$Counts)
#about 5.9
#what's the standard deviation?
sd(PerSampleOTUcounts.bySiteandBatSpecies$Counts)
#about 4.3
#what's the median?
median(PerSampleOTUcounts.bySiteandBatSpecies$Counts)
#it's 5  

#write all to disk:
setwd("~/Documents/Lab.Foster/guano/guano_nau/all/R_OTUcall_output/")
write.csv(OTUcounts.perBatAndSite, "OTUcounts.perBatAndSite.csv", quote = F, row.names = F)
write.csv(PerSampleOTUcounts.bySiteandBatSpecies, "PerSampleOTUcounts.bySiteandBatSpecies.csv", quote = F, row.names = F)

########################################################################################
############################    SECTION 4: data analysis  ##############################
########################################################################################

#########################################################################################
########### PART A - Broad analyses - OTUs per WOY, BatSpecies, or Sample ##############

#very basic: how frequently are taxa are called per Site for a given {Family/Order} among all birds?
head(Taxa.perSite.byFamily)
head(Taxa.perSite.byOrder)
#if BatSpecies is unknown:
#OTUcounts.perBatAndWOY$BatSpecies <- as.character(OTUcounts.perBatAndWOY$BatSpecies)    
#unnecessary if object imported with 'stringsAsFactors = F'
#OTUcounts.perBatAndWOY[is.na(OTUcounts.perBatAndWOY)] <- "undetermined"

#how frequently are taxa called to a given BatSpecies among all Sites?
head(Taxa.perBatSpecies.byFamily)
head(Taxa.perBatSpecies.byOrder)

#how frequently is a given family identified at a unique Site and unique BatSpecies?
#see SECTION 2
head(percFamily.bySiteandBatSpecies)
head(percOrder.bySiteandBatSpecies)

#how many total OTUs were identified for each BatSpecies at each Site?
#see object: "OTUcounts_perSample" from SECTION 3, #3a.
head(OTUcounts.perBatAndSite)

#what kind of range existed for OTU detection on a per-sample basis across every BatSpecies at each Location? 
#see object: "OTUcounts_perSample" from SECTION 3, #3b.
head(PerSampleOTUcounts.bySiteandBatSpecies)
###############################   END PART A    #########################################
#########################################################################################


#########################################################################################
######### PART B - Factor-specific analyses: Sex, Age, and Reproductive Status ############

setwd("~/Documents/Lab.Foster/guano/guano_nau/all/R_OTUcall_output/")

#make a new table of numbers of two factors: Site and BatSpecies
BOLDonly.dfSamples <- subset(BOLDonly.df, !duplicated(SampleID))
#keep only relevant fields
BOLDonly.dfSamples <- BOLDonly.dfSamples[,c("Site", "BatSpecies")]

#how many of each Site are there?
SiteTable.BOLDonly = as.data.frame(table(BOLDonly.dfSamples$Site))
colnames(SiteTable.BOLDonly) <- c("Site", "counts")
write.csv(SiteTable.BOLDonly, "SiteTable.BOLDonly.csv", quote = F, row.names = F)
  #note the single NewMexico sample is the unique EUMA sample

#how many of each BatSpecies are there?
BatSpeciesTable.BOLDonly = as.data.frame(table(BOLDonly.dfSamples$BatSpecies))
colnames(BatSpeciesTable.BOLDonly) <- c("BatSpecies", "counts")
write.csv(BatSpeciesTable.BOLDonly, "BatSpeciesTable.BOLDonly.csv", quote = F, row.names = F)
#a lot of singleton BatSpecies samples...
###############################   END PART B    #########################################
#########################################################################################

##################################################################################################
########################      SECTION 5: Making the plots     ####################################
##################################################################################################

#########################################################################################
##########################           PART A - ggplots        ############################

#load ggplot then play around with formatting
library(ggplot2)
require(graphics)
require(RColorBrewer)

##############################  Historgram plots  ##############################  
##########  plot-h1: range of OTUs detected for samples ##########  

#see here for useful tips when constructing histograms: 
# http://t-redactyl.io/blog/2016/02/creating-plots-in-r-using-ggplot2-part-7-histograms.html

#load in the data and observe default distribution
dat.h1 = read.csv(file = "~/Documents/Lab.Foster/guano/guano_nau/all/R_OTUcall_output/PerSampleOTUcounts.bySiteandBatSpecies.csv")
h1 <- ggplot(dat.h1, aes(x=Counts))
h1 + geom_histogram()

#set colours ahead of time
barfill <- "#4271AE"
barlines <- "#1F3552"

#make your plot
h1 + geom_histogram(binwidth = 5, colour = barlines, fill = barfill) +
  labs(title="Distribution of OTU frequency detections for all NAU samples") +
  labs(x="Number of OTUs detected", y="Number of Samples") +
  scale_x_continuous(breaks = seq(0, 40, 5),
                     limits = c(0, 40))

#as indicated above, most samples contain between ~5-15 OTUs.



################################    Bar plots   ################################  
###for most plots you're going to want to tailor the color scheme... there's a bunch of places to read up on choices:
#see here for examples: http://www.r-bloggers.com/choosing-colour-palettes-part-ii-educated-choices/)
#for example, if you wanted to create a palette of 20 colors (because you had 20 factors to differentiate) from a palette containing a total of 12 unique colors, you would enter the following command:
# Pal20cols = colorRampPalette(brewer.pal(12, "Paired"))(20)
#another palette option... go here: http://tools.medialab.sciences-po.fr/iwanthue/
#use the "HEX json" output to copy and paste the relevant values to create the palette you want
#for example, if you had 13 unique WOYs to colour, you could use that WOY then copy that text into a new palette like this:
#WOYpal <- c("#ffacf4", "#5cff6f","#5d2ed9","#3a5500","#78cf78","#01fbf2","#ff5645","#00b9e4","#922700","#e9d287","#b6b9b2","#a60064","#00313f")


##################################################################################################
#####   Plot-b1: How many times is an OTU detected per BatSpecies   ######
#load in the data
dat.b1 = read.csv(file = "~/Documents/Lab.Foster/guano/guano_nau/all/R_OTUcall_output/OTUcounts.perBatAndSite.csv")

#how many bat species?
length(unique(dat.b1$BatSpecies))
  #there are 24 BatSpecies. 

#chose your palette with the right number of colors. Because there are 24 unique Bat species, it may be easier to group by hue within a genus. I'm not going to do that here, but it could be done pretty easily: 
b1pal = c("#547a9f",
          "#78da4c",
          "#653ec6",
          "#cfd04a",
          "#c94fd3",
          "#79db96",
          "#b9408f",
          "#5c9139",
          "#542e7a",
          "#cf8c39",
          "#7478d3",
          "#d64c34",
          "#70d7d0",
          "#cc486b",
          "#528a71",
          "#d494cb",
          "#333e28",
          "#a4c1dc",
          "#793528",
          "#ccd09c",
          "#3c243e",
          "#d49c8b",
          "#7d7137",
          "#8b6372")

#set colours ahead of time for this basic example (no special fill required by some other factor)
barfill <- "#4271AE"
barlines <- "#1F3552"

#plot data 
b1 <- ggplot(dat.b1, aes(x = Site, y = Counts, fill = BatSpecies, label = BatSpecies)) +
  geom_bar(stat = "identity", colour=barlines) +
  scale_fill_manual(values = b1pal, name="Bat Species") + 
  ylab("Number of times an OTU is detected") +
  xlab("Site") + 
  theme(axis.text.y = element_text(size = 8)) +
  theme(axis.text.x = element_text(size = 10)) +
  labs(title = "Number of OTUs detected per Bat Species and Site using DADA2/AMPtk pipeline") +
  theme(plot.title = element_text(face = "bold")) +
  #theme(legend.key.size = unit(4, "mm")) +
  theme(legend.position = "none") +
  geom_text(size = 3, position = position_stack(vjust = 0.5), angle = 45) +
  coord_flip()
b1

## so clearly PTPA accounts for a big chunk of all OTUs detected along with EPFUR, EPBR?, RHBI, and PTPE. Lots of detections at low levels for many other species, which is more an indication of the few numbers of samples for those species, not anything to do with diet detection bias per Bat species.
## you can also see that Finca_Guadelupe has very few OTUs and would probably be worth tossing if you wanted to do a Site by Site comparison...
##################################################################################################

#####   Plots-b2a and b2b: For a given BatSpecies, what proportions of arthropod Orders are they eating at a given Site?   ######
#load in the data
dat.b2 = read.csv("~/Documents/Lab.Foster/guano/guano_nau/all/R_OTUcall_output/percOrder.bySiteandBatSpecies.csv")
  #this breaks up BatSpecies within each site

#the problem with these data is that we have an unequal distribution of OTUs per Site; i'm going to break up the New Mexico sample to be it's own plot, and keep the Central American plots as their own group:

#subset just the EUMA
dat.b2euma <- dat.b2[grepl("EUMA", dat.b2$BatSpecies),]
#subset everything else not 'EUMA'
dat.b2cams <- dat.b2[!grepl("EUMA", dat.b2$BatSpecies),]

#to create the appropriate palette you need to determine how many factors (Arthropod orders) there are for each plot:

#for the larger subsets (everything but EUMA)
length(unique(dat.b2cams$ArthropodOrder))
#there are 11
#for the small subset (EUMA only)
length(unique(dat.b2euma$ArthropodOrder))
#there are 2

#Of note: there are 11 total Orders, but we should keep colors consistent between both graphs. I like using darker hues for the majority of colors which show low count totals for the larger plot, then use bright and distinct colors for the more frequent ones (in this case there are four Orders that are well represented, and 7 that aren't, the frequent ones being Coleoptera, Diptera, Hemiptera, and Lepidoptera). However, the way the bigger plot is going to look you're not going to see a lot of all those distinct orders next toe ach other each time, so it's best to have as much contrast as possible among all orders. So we're just going for 11 different colors.

#these are the 11 colors for the all Orders, These are present in alphabetical order:
palb2cams <- c("#657800",
               "#e86bff",
               "#84ff99",
               "#013aab",
               "#df8700",
               "#79c7ff",
               "#ff4557",
               "#003a05",
               "#e3d1ff",
               "#760032",
               "#9a9a9b")

#for b2euma we'll stay consistent and just use the two Orders detected: Diptera and Lepidoptera
palb2euma <- c("#84ff99", "#e3d1ff")

#set colour for lines around bars as usual:  
barlines <- "#1F3552"

#order data.frame by Arthropod Order name for both sets
dat.b2cams <- dat.b2cams[with(dat.b2cams, order(ArthropodOrder)), ] 
dat.b2euma <- dat.b2euma[with(dat.b2euma, order(ArthropodOrder)), ]

#facet-names - used to relabel the facets to show samples sizes contributing to the resulting illustration
#not used. Though you could check out sample size as follows: 
  # b2cams.facetnames <- c(CAPE = "CAPE, n = 2", EPBR? = "EPBR?, n = 1", etc...)
  # see "BatSpeciesTable.BOLDonly.csv' file for full counts for all BatSpecies

b2cams <- ggplot(dat.b2cams, aes(x = ArthropodOrder, y = RelativeCounts, fill = ArthropodOrder)) + 
  geom_bar(stat = "identity", colour=barlines) +
  scale_fill_manual(values = palb2cams, name="Arthropod Order") + 
  guides(col = guide_legend(nrow = 1)) + 
  theme(legend.position = "none") +
  theme(axis.text.y = element_text(size = 8)) +
  theme(axis.text.x = element_text(size = 10)) +
  labs(title = "Total Counts of Arthropod Order detected among Central American BatSpecies at all Central American Sites using BOLD database") +
  labs(y="Total Counts of each Arthropod Order Detected") +
  theme(legend.key.size = unit(4, "mm")) +
  #facet_grid(~BatSpecies, labeller=labeller(BatSpecies=b2a.facetnames)) +
  #facet_grid(~BatSpecies) +
  ## facet_grid(~BatSpecies, scales = "free", space = "free") +
  facet_wrap(~BatSpecies + Site, ncol = 8, scales = "free_x") +
  scale_y_continuous(breaks = c(0,5,10,20,30,40,50)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
b2cams


#or for b2euma
b2euma <- ggplot(dat.b2euma, aes(x = ArthropodOrder, y = RelativeCounts, fill = ArthropodOrder)) + 
  geom_bar(stat = "identity", colour=barlines) +
  scale_fill_manual(values = palb2euma, name="Arthropod Order") + 
  guides(col = guide_legend(nrow = 1)) + 
  theme(legend.position = "none") +
  theme(axis.text.y = element_text(size = 8)) +
  theme(axis.text.x = element_text(size = 10)) +
  labs(title = "Total Counts of Arthropod Order detected among single EUMA sample using BOLD database") +
  labs(y="Total Counts of each Arthropod Order Detected") +
  theme(legend.key.size = unit(4, "mm")) +
  facet_wrap(~BatSpecies + Site) +
  #scale_y_continuous(breaks = c(0,5,10,20,30,40,50)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
b2euma
  #This is a stupid graph for EUMA as there is so little diversity at the Order level...
  #Even plotting at the Family level is a bit silly because you have just a total of 9 observations. Data Table is probably more appropriate here...

##################################### FAMILY PLOTS   ########################################
#####   Plot-b3: How many times is an Arthropod Family detected per BatSpecies per Site?   ######
#load data:
dat.b3 = read.csv("~/Documents/Lab.Foster/guano/guano_nau/all/R_OTUcall_output/percFamily.bySiteandBatSpecies.csv")
#how many Family factors for ArthropodFamily class?

#subset Central American species (no EUMA)
dat.b3cams <- dat.b3[!grepl("EUMA", dat.b3$BatSpecies),]

#to create the appropriate palette you need to determine how many factors (Arthropod Families) there are for each plot:

#for the larger subsets (everything but EUMA)
length(unique(dat.b3cams$ArthropodFamily))
#there are 30

#Of note: there are 11 total Orders, but there are 30 Families within those Orders. Probably the best way to plot this would be to take hues of each Order color and nest them for the individual Families. I'm not going to do that here, but it would follow the same principle.

#generate your pal3bcams list here:
  #palb3cams <- c(your list of colors)

#set colour for lines around bars as usual:  
barlines <- "#1F3552"

#order data.frame by Arthropod Order name for both sets
dat.b3cams <- dat.b3cams[with(dat.b3cams, order(ArthropodFamily)), ] 

#b3cams <- ggplot(dat.b3cams, aes(x = ArthropodFamily, y = RelativeCounts, fill = ArthropodFamily)) + 
b3cams <- ggplot(dat.b3cams, aes(x = ArthropodFamily, y = RelativeCounts)) +
  geom_bar(stat = "identity", colour=barlines) +
  #scale_fill_manual(values = palb3cams, name="Arthropod Family") + 
  guides(col = guide_legend(nrow = 1)) + 
  theme(legend.position = "none") +
  theme(axis.text.y = element_text(size = 8)) +
  theme(axis.text.x = element_text(size = 10)) +
  labs(title = "Total Counts of Arthropod Family detected among Central American BatSpecies at all Central American Sites using BOLD database") +
  labs(y="Total Counts of each Arthropod Family Detected") +
  theme(legend.key.size = unit(4, "mm")) +
  facet_wrap(~BatSpecies + Site, ncol = 8, scales = "free_x") +
  #scale_y_continuous(breaks = c(0,5,10,20,30,40,50)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
b3cams


##################################################################################################
##############################      Part B: Treemaps     #########################################
#load package
library(treemap)
require(graphics)

#load file if necessary:
setwd("~/Documents/Lab.Foster/guano/guano_nau/all/R_OTUcall_output/")
library(readr)
BOLDonly.df <- read_csv("~/Documents/Lab.Foster/guano/guano_nau/all/R_OTUcall_output/BOLDonly.df.csv", 
                        col_types = cols(BatTagID = col_character(), 
                                         Sex = col_character()))

treemap_df <- BOLDonly.df
treemap_df <- as.data.frame(treemap_df)
treemap_df$counts <- 1

#subset everything else not 'EUMA'
treemap_df <- treemap_df[!grepl("EUMA", treemap_df$BatSpecies),]

#recall how many arthropod Orders are in this dataset:
length(unique(treemap_df$order_name))
#11 total

#we're going to use the "palb2cams" palette described in the ggplot setup above. This is the palette:
palb2cams <- c("#657800", "#e86bff", "#84ff99", "#013aab", "#df8700",
               "#79c7ff", "#ff4557", "#003a05", "#e3d1ff", "#760032", "#9a9a9b")

#order the data.frame by Arthropod Order ('order_name'), then by Arthropod Family
treemap_df <- treemap_df[with(treemap_df, order(order_name, family_name)), ] 

#option2 - look at ORDER eaten among BatSpecies
treemap(treemap_df, 
        index=c("BatSpecies", "order_name"),
        vSize="counts",
        vColor="order_name",
        type="categorical",
        title = "Total Counts of Arthropod Order by BatSpecies for all Central American samples",
        title.legend = "Taxonomic Order",
        position.legend = "right",
        align.labels=list(c("center", "top"), c("left", "bottom")), 
        fontsize.labels=c(14,6),
        lowerbound.cex.labels=0.2, 
        overlap.labels=0.4,
        palette=palb2cams
)

###############################   END PART B    #########################################
#########################################################################################

##############################################################################
########################## The EUMA bit  ##########################
##############################################################################

EUMAdf <- BOLDonly.df[grepl("EUMA", BOLDonly.df$BatSpecies),]
write.csv(EUMAdf, file = "~/Documents/Lab.Foster/guano/guano_nau/all/R_OTUcall_output/EUMAdf.csv")
