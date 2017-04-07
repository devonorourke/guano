#amptk_2017_OTUtable_analysis.R
#written 10.12.16
#modified 3.29.17
#Devon O'Rourke
#amptk analysis for Acadia samples from BRI

########################################################################################
############################    SECTION 0: metadata files  #############################
########################################################################################

#################### process started: creating metadata file  ####################
date()

#import collection data
setwd("~/Documents/Lab.Foster/guano/")
coldat16 <- read.csv(file = "collectiondata2016.csv", stringsAsFactors = F)
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

#loook at just the bri data
bridat = coldat16[grepl("bri16", coldat16$SampleID),]

#import list of names from p1 run associated with BRI data
setwd("~/Documents/Lab.Foster/guano/guano_bri/catp1p2/")
p1names <- read.delim(file = "p1names.txt", col.names = F)
p2names <- read.delim(file = "p2names.txt", col.names = F)
brinames <- rbind(p1names, p2names)
colnames(brinames) <- "SampleID"
brinames$SampleID <- gsub('-', '.', brinames$SampleID)

#merge so that only samples listed in the 'brinames' data.frame are pulled from the 'bridat' data.frame
onlyseqd <- merge(brinames, bridat, by = "SampleID")

#pull out just the Acadia samples
aca.df <- subset(onlyseqd, onlyseqd$LocationName == "ACA")

#which samples are from Acadia in the entire collection (not just those sequenced thus far)?
aca.coll <- subset(coldat16, coldat16$LocationName == "ACA")

#this data.frame will serve as your metadata file:
meta_df <- aca.coll
rm(aca.coll)

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

#make new field called "SiteWeek" which may help identify the nature of the outliers in question
meta_df$SiteWeek <- paste(meta_df$LocationName, meta_df$WeekOfYear)

#remove fields that are unnecessary for this project:
meta_df$PlateNumber <- NULL
meta_df$WellNumber <- NULL
meta_df$Status <- NULL
meta_df$Roost <- NULL
meta_df$Notes <- NULL

#finally, substitute blank cells for certain terms in specific columns:
meta_df$BatSpecies[meta_df$BatSpecies==""]<-"unknown"
meta_df$BatID[meta_df$BatID==""]<-"no-record"
meta_df$Sex[meta_df$Sex==""]<-"no-record"
meta_df$Age[meta_df$Age==""]<-"no-record"
meta_df$ReproStat[meta_df$ReproStat==""]<-"no-record"

#make a few frequency tables of various meta_dfdata to ensure you have the fields/values labeled correctly:
library(plyr)
Location.table <- count(meta_df, 'LocationName')
  #good. only one location, and all are for Acadia
WeekOfYear.table <- count(meta_df, 'WeekOfYear')
  #we have eight weeks worth of data, though most was collected in Weeks 24, 35, and 36
BatSpecies.table <- count(meta_df, 'BatSpecies')
  #almost all samples are from MYLE and MYLU; only a single PESU, two EPFUs, and three MYSE's
Sex.Table <- count(meta_df, 'Sex')
  # 49 males : 27 females
Age.Table <- count(meta_df, 'Age')
  # 60 adults : 15 juvenilles; only 4 records without info
ReproStat.Table <- count(meta_df, 'ReproStat')
  # majority of samples are 'NR' == non-reproductive (likely from so many males)

#the above tables suggest that our meta_dfdata doesn't appear to need further modifications or updates.

#final cleanup
rm(aca.df, brinames, onlyseqd, p1names, p2names, coldat16, bridat)
rm(Age.Table, BatSpecies.table, Location.table, ReproStat.Table, Sex.Table, WeekOfYear.table)

#write meta_dfdata to disk for aca.samples
setwd("~/Documents/Lab.Foster/guano/guano_bri/catp1p2/R_OTUcall_output/")
write.table(meta_df, file = "meta_dfACAonly.txt", col.names = T, row.names = F, sep = "\t", quote = F)

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
setwd("~/Documents/Lab.Foster/guano/guano_bri/catp1p2/")

#import taxonomy-classified OTU Table:
amptk_otu_tableA <- read.delim("bri.otu_table.taxonomy.txt", stringsAsFactors = F)

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
setwd("~/Documents/Lab.Foster/guano/guano_bri/catp1p2/R_OTUcall_output/")
meta_df = read.delim(file = "meta_dfACAonly.txt", header = T, stringsAsFactors = F)

#rename that first column
colnames(meta_df)[1] <- "tagID"

#replace the hyphen with a period in the "tagID" column (may not be necessary - check both data.frames for consistency in delimiter)
#metadata$tagID <- sub("-", "\\.", metadata$tagID)

#merge OTU calls with metadata by 'tagID' (really just the sample ID)
masterdf = merge(amptk_dt_unique, meta_df, by = "tagID")

#write to disk:
setwd("~/Documents/Lab.Foster/guano/guano_bri/catp1p2/R_OTUcall_output/")
write.table(masterdf, file = "masterdf.txt", col.names = T, row.names = F, sep = "\t", quote = F)

#########################################################
#we don't expect the same number of observations in the 'amptk_dt_unique' table here because the ACA samples represent just a fraction of the overall number of samples that were processed in the OTU-calling AMPtk pipeline (there are also samples from VT, NY, MA, etc.). The following scripts are useful only if you expect the same number, and can be ignored in this specific instance.

#side note if you're getting weird results after merging (not same number of observations pre and post merge from above)
  #uniqA <- unique(amptk_dt_unique$tagID)
  #uniqB <- unique(meta_df$tagID)

#what's in both?
  #both <- intersect(uniqA, uniqB)
  #length(both)
    #how many?

#whats in A not in B?
  #AnotB <- setdiff(uniqA, uniqB)
    #how many characters? if zero, then all of A in B.

#whats in B not in A?
  #BnotA <- setdiff(uniqB, uniqA)
    #how many characters? if more than zero, that's your issue.
#########################################################

#drop out any reads you don't want to include in subsequent analyses:

#to keep things conservative for this first go of things, we're going to limit our observations excusively to BOLD-called OTUs (we may improve upon their description though); we are not going to include the SINTAX and UTAX terms here, though these may play a role in subsequent analyses. These are going to be removed at this point.

#this gets rid of the parasitic mites
BOLDonly.df = masterdf[!grepl("SINTAX", masterdf$OTU_identifier),]
BOLDonly.df = BOLDonly.df[!grepl("UTAX", BOLDonly.df$OTU_identifier),]

#we're also going to make a separate file for the bat-specific reads; a quick analysis suggests these aren't worth using:
  #1. there are only 15 total calls out of 79 potenteial samples
  #2. of the few that worked and they are clearly wrong in some instances (Mollosus?)
  #3. of these, some times there are two bat species associated with a single sample
batOTUs.df = BOLDonly.df[grepl("Chiroptera", BOLDonly.df$order_name),]

#now remove the Bat OTUs from the 'BOLDonly.df' data.frame
BOLDonly.df = BOLDonly.df[!grepl("Chiroptera", BOLDonly.df$order_name),]

#rename a few columns so we're explicit about what species is being described
colnames(BOLDonly.df)[11:14] <- c("Site", "Date", "WOY", "BatTagID")

#it's worth renaming empty (N/A) values for taxonomic classes with "unknown" for future plotting purposes (we'll want to show in our data how many OTUs weren't explicitly described relative to those that were)
BOLDonly.df[is.na(BOLDonly.df)] <- "unknown"
  #careful with this term, as it will convert any 'N/A' in the dataset, regardless of what field it belongs to, to "unknown"

#write object 'masterdf' to disk
setwd("~/Documents/Lab.Foster/guano/guano_bri/catp1p2/R_OTUcall_output/")
write.table(BOLDonly.df, file = "BOLDonly.df.txt", col.names = T, row.names = F, sep = "\t", quote = F)

#CLEANUP
rm(meta_df, amptk_dt_unique, masterdf, batOTUs.df)  

#################### process ended: OTU table import and metadata merging  ####################
############### BOLDonly.df file contains all information for subsequent work  ################

##################################################################################################
###################     SECTION 3: Calculate Frequencies for plots           #####################
##################################################################################################

################## process started: data tables created for plots and analysis  ##################

#load file if necessary:
setwd("~/Documents/Lab.Foster/guano/guano_bri/catp1p2/R_OTUcall_output/")
BOLDonly.df <- read.delim("BOLDonly.df.txt")

#need to generate a suite of data tables for each unique plotting scenario: 
## 1. you want to know the relative frequency (in percentage terms) of each 'family_name'(or some other taxonomic unit) when you aggregate all counts for:
##  a. All bats at each unique WOY or Site, or
##  b. All WOY or Sites for each unique BatSpecies
## 2. What'll likely get used more for analyses: the relative frequency of each 'family_name' identified when aggregating counts for:
### a. each unique BatSpecies at each unique WOY or Site

#In addition, we're going to generate a few other frequency tables that may prove useful:
## 3a. a frequency table showing how many total OTUs were identified (regardless of taxonomic class) for each BatSpecies at each WOY or Site
## 3b. a frequency table showing how many unique OTUs were identified (by taxonomic 'family') for each sample, for each BatSpecies at each WOY or Site


######## START Part 1 ########
#first, make a table with the number of counts per 'family_name' for all our metadata (Site, WOY, BatSpecies, Sex, Age, and ReproStat) that's unique:
library(plyr)
Familycounts_table <- ddply(BOLDonly.df, .(BOLDonly.df$Site, BOLDonly.df$WOY, BOLDonly.df$BatSpecies, BOLDonly.df$Sex, BOLDonly.df$Age, BOLDonly.df$ReproStat, BOLDonly.df$family_name), nrow)
names(Familycounts_table) <- c("Site", "WOY", "BatSpecies", "Sex", "Age", "ReproStat", "ArthropodFamily", "RelativeCounts")

#next, repeat for 'order_name' across each BatSpecies and WOY that's unique
Ordercounts_table <- ddply(BOLDonly.df, .(BOLDonly.df$Site, BOLDonly.df$WOY, BOLDonly.df$BatSpecies, BOLDonly.df$Sex, BOLDonly.df$Age, BOLDonly.df$ReproStat, BOLDonly.df$order_name), nrow)
names(Ordercounts_table) <- c("Site", "WOY", "BatSpecies", "Sex", "Age", "ReproStat", "ArthropodOrder", "RelativeCounts")

##Note the old way of doing things was to make a count table for each factor (just for "Site", or "WOY", or "Sex", or "Age", or "ReproStat"):
  #Ordercounts_table <- ddply(BOLDonly.df, .(BOLDonly.df$BatSpecies, BOLDonly.df$WOY, BOLDonly.df$order_name), nrow)
  #names(Ordercounts_table) <- c("BatSpecies", "WOY", "ArthropodOrder", "RelativeCounts")

#write to disk:
setwd(dir = "~/Documents/Lab.Foster/guano/guano_bri/catp1p2/R_OTUcall_output/")
write.table(Familycounts_table, file="Familycounts_table.txt", quote = F, sep = "\t", row.names = F, col.names = T)
write.table(Ordercounts_table, file="Ordercounts_table.txt", quote = F, sep = "\t", row.names = F, col.names = T)

#next, sum up the frequency values given one of the two situations described above:
##1a. For each unique WOY, grouping all BatSpecies (ie. how many times did we see a family/order of arthropod get eaten at a given WOY, regardless of BatSpecies?)
library(dplyr)
FamilyCount.byWOY <- Familycounts_table %>% 
  group_by(WOY, ArthropodFamily) %>% 
  summarise(Freq = sum(RelativeCounts))
names(FamilyCount.byWOY) <- c("WOY", "ArthropodFamily", "RelativeCounts")

TotalFamilyCount.perWOY <- Familycounts_table %>% 
  group_by(WOY) %>% 
  summarise(Freq = sum(RelativeCounts))
names(TotalFamilyCount.perWOY) <- c("WOY", "TotalCounts")

#and repeat for the "order" table
OrderCount.byWOY <- Ordercounts_table %>% 
  group_by(WOY, ArthropodOrder) %>% 
  summarise(Freq = sum(RelativeCounts))
names(OrderCount.byWOY) <- c("WOY", "ArthropodOrder", "RelativeCounts")

TotalOrderCount.perWOY <- Ordercounts_table %>% 
  group_by(WOY) %>% 
  summarise(Freq = sum(RelativeCounts))
names(TotalOrderCount.perWOY) <- c("WOY", "TotalCounts")


#finally, merge the two tables (repeat for each taxonimic group) by "WOY"
  ##first, for family
  Taxa.perWOY.byFamily <- merge(FamilyCount.byWOY, TotalFamilyCount.perWOY, by = "WOY")
  ##then, for order
  Taxa.perWOY.byOrder <- merge(OrderCount.byWOY, TotalOrderCount.perWOY, by = "WOY")

#calculate the relative frequency per Arthropod by taxa level per WOY
  ##first, for family
  Taxa.perWOY.byFamily$percent_taxa <- (Taxa.perWOY.byFamily$RelativeCounts / Taxa.perWOY.byFamily$TotalCounts * 100)
  ##then, for order
  Taxa.perWOY.byOrder$percent_taxa <- (Taxa.perWOY.byOrder$RelativeCounts / Taxa.perWOY.byOrder$TotalCounts * 100)
  
###these final two tables (by taxoomic Family and Order) allow you to investigate the relative proportions of each unique taxa; because some Families or Orders are going to be far more represented DUE TO SAMPLING AT ONE WOY MORE THAN ANOTHER (not because of diet preferences), this table reflects the RELATIVE proportions of a given taxonomic Family or Order at each WOY so as to reduce the sampling skew.

#write these tables to disk:
getwd()
setwd("/Users/devonorourke/Documents/Lab.Foster/guano/guano_bri/catp1p2/R_OTUcall_output/")
write.table(Taxa.perWOY.byFamily, "Taxa.perWOY.byFamily.txt", quote = F, sep = "\t", row.names = F, col.names = T)  
write.table(Taxa.perWOY.byOrder, "Taxa.perWOY.byOrder.txt", quote = F, sep = "\t", row.names = F, col.names = T)  
#cleanup
rm(FamilyCount.byWOY, TotalFamilyCount.perWOY, OrderCount.byWOY, TotalOrderCount.perWOY)


#repeat for "Site" instead of 'WOY' variable
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


#finally, merge the two tables (repeat for each taxonimic group) by "WOY"
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
setwd("/Users/devonorourke/Documents/Lab.Foster/guano/guano_bri/catp1p2/R_OTUcall_output/")
write.table(Taxa.perSite.byFamily, "Taxa.perSite.byFamily.txt", quote = F, sep = "\t", row.names = F, col.names = T)  
write.table(Taxa.perSite.byOrder, "Taxa.perSite.byOrder.txt", quote = F, sep = "\t", row.names = F, col.names = T)  
#cleanup
rm(FamilyCount.bySite, TotalFamilyCount.perSite, OrderCount.bySite, TotalOrderCount.perSite)



##1b. For each unique BatSpecies, grouping all WOYs (ie. how many times did we see an taxonomic group of arthropods get eaten by a particular BatSpecies, regardless of WOY?). The work follows the same format as above but replaces the "Site" with "BatSpecies" factor where appropriate.

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
  
###these final two tables (by taxoomic Family and Order) allow you to investigate the relative proportions of each unique taxa; because some Families or Orders are going to be far more represented DUE TO SAMPLING OF MORE OF ONE BIRD SPECIES THAN ANOTHER (not because of diet preferences), this table reflects the RELATIVE proportions of a given taxonomic Family or Order for each BatSpecies so as to reduces the sampling skew.

#write to disk:
setwd("/Users/devonorourke/Documents/Lab.Foster/guano/guano_bri/catp1p2/R_OTUcall_output/")
write.table(Taxa.perBatSpecies.byFamily, "Taxa.perBatSpecies.byFamily.txt", quote = F, sep = "\t", row.names = F, col.names = T)
write.table(Taxa.perBatSpecies.byOrder, "Taxa.perBatSpecies.byOrder.txt", quote = F, sep = "\t", row.names = F, col.names = T)
#cleanup:
rm(FamilyCount.byBatSpecies, TotalFamilyCount.perBatSpecies, OrderCount.byBatSpecies, TotalOrderCount.perBatSpecies)
######## END Part 1 ########


######## START Part 2 ########
#For each unique BatSpecies at each unique WOY or Site (for both taxonomic levels of arthropods)

### Part 2a ###
#2a. For WOY and BatSpecies: the relative counts for each unique pair of BatSpecies and WOY classes:
#first for Family:
FamilyCount.byWOYandBatSpecies <- Familycounts_table %>% 
  group_by(BatSpecies, WOY, ArthropodFamily) %>% 
  summarise(Freq = sum(RelativeCounts))
names(FamilyCount.byWOYandBatSpecies) <- c("BatSpecies", "WOY", "ArthropodFamily", "RelativeCounts")
#then for Order:
OrderCount.byWOYandBatSpecies <- Ordercounts_table %>% 
  group_by(BatSpecies, WOY, ArthropodOrder) %>% 
  summarise(Freq = sum(RelativeCounts))
names(OrderCount.byWOYandBatSpecies) <- c("BatSpecies", "WOY", "ArthropodOrder", "RelativeCounts")

##Next, tally the total counts of BatSpecies split by WOY. It's the same regardless of Family/Order taxa level:  
WOYandBatSpecies.TotalCounts <- Familycounts_table %>% 
  group_by(BatSpecies, WOY) %>% 
  summarise(Freq = sum(RelativeCounts))
names(WOYandBatSpecies.TotalCounts) <- c("BatSpecies", "WOY", "TotalCounts")

#merge the original '{Family/Order}counts_table' object with the "WOYandBatSpecies.TotalCounts" object with R package 'dplyr'
library(dplyr)
#for Family taxa:
  percFamily.byWOYandBatSpecies <- FamilyCount.byWOYandBatSpecies %>% inner_join(WOYandBatSpecies.TotalCounts)
  #and calculate relative frequency per insect Family per WOY AND BatSpecies
  percFamily.byWOYandBatSpecies$percent_taxa <- (percFamily.byWOYandBatSpecies$RelativeCounts / percFamily.byWOYandBatSpecies$TotalCounts * 100)

#for Order taxa:
  percOrder.byWOYandBatSpecies <- OrderCount.byWOYandBatSpecies %>% inner_join(WOYandBatSpecies.TotalCounts)
  percOrder.byWOYandBatSpecies$percent_taxa <- (percOrder.byWOYandBatSpecies$RelativeCounts / percOrder.byWOYandBatSpecies$TotalCounts * 100)

#These tables represent the relative proportions for a Family/Order that was detected for a given BatSpecies at a specific WOY.

#write tables to disk:
setwd("~/Documents/Lab.Foster/guano/guano_bri/catp1p2/R_OTUcall_output/")
write.table(percFamily.byWOYandBatSpecies, "percFamily.byWOYandBatSpecies.txt", quote = F, sep = "\t", row.names = F, col.names = T)
write.table(percOrder.byWOYandBatSpecies, "percOrder.byWOYandBatSpecies.txt", quote = F, sep = "\t", row.names = F, col.names = T)  
#these tables partition WOY and BatSpecies factors into unique classes, thus separating out the total and relative counts into smaller proportions while being able to attribute both WOY and BatSpecies differences for a given Family or Order abundance.

#cleanup:
rm(FamilyCount.byWOYandBatSpecies, OrderCount.byWOYandBatSpecies, WOYandBatSpecies.TotalCounts)
##########    END Part 2A     ##########


### Part 2a ###
#2a. For WOY and BatSpecies: the relative counts for each unique pair of BatSpecies and WOY classes:
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
setwd("~/Documents/Lab.Foster/guano/guano_bri/catp1p2/R_OTUcall_output/")
write.table(percFamily.bySiteandBatSpecies, "percFamily.bySiteandBatSpecies.txt", quote = F, sep = "\t", row.names = F, col.names = T)
write.table(percOrder.bySiteandBatSpecies, "percOrder.bySiteandBatSpecies.txt", quote = F, sep = "\t", row.names = F, col.names = T)  
#these tables partition Site and BatSpecies factors into unique classes, thus separating out the total and relative counts into smaller proportions while being able to attribute both Site and BatSpecies differences for a given Family or Order abundance.

#cleanup:
rm(FamilyCount.bySiteandBatSpecies, OrderCount.bySiteandBatSpecies, SiteandBatSpecies.TotalCounts)
##########    END Part 2B     ##########

######## START Part 3 ########
##3a. this calculates the number of times an OTU (at any taxonomic level) was detected for each BatSpecies at each Site (1) or WOY (2).
  #note you may want to alter this so that you only include those fully identified through to some specific level (say at least "Order" or "Family") so that you don't include many calls where the OTU could be described only to a Phylum level (that is, you only know that it's an arthropod and can't even distinguish if it's an insect or spider according to this initial search with our BOLD database)

#3a1 - for Site
library(plyr)
OTUcounts.perBatAndSite <- ddply(BOLDonly.df, .(BOLDonly.df$BatSpecies, BOLDonly.df$Site), nrow)
names(OTUcounts.perBatAndSite) <- c("BatSpecies", "Site", "Counts")

#3a2 - for WOY
library(plyr)
OTUcounts.perBatAndWOY <- ddply(BOLDonly.df, .(BOLDonly.df$BatSpecies, BOLDonly.df$WOY), nrow)
names(OTUcounts.perBatAndWOY) <- c("BatSpecies", "WOY", "Counts")

##3b. this calculates the number of OTUs called per sample of guano while also retaining the BatSpecies and Site/WOY information:
library(plyr)
#for WOY
PerSampleOTUcounts.byWOYandBatSpecies <- ddply(BOLDonly.df, .(BOLDonly.df$tagID, BOLDonly.df$BatSpecies, BOLDonly.df$WOY), nrow)
names(PerSampleOTUcounts.byWOYandBatSpecies) <- c("tagID", "BatSpecies", "WOY", "Counts")
  #what's the mean?
  mean(PerSampleOTUcounts.byWOYandBatSpecies$Counts)
    #about 13
  #what's the standard deviation?
  sd(PerSampleOTUcounts.byWOYandBatSpecies$Counts)
    #about 8
  #what's the median?
  median(PerSampleOTUcounts.byWOYandBatSpecies$Counts)
    #it's "10"

#for Site
PerSampleOTUcounts.bySiteandBatSpecies <- ddply(BOLDonly.df, .(BOLDonly.df$tagID, BOLDonly.df$BatSpecies, BOLDonly.df$Site), nrow)
names(PerSampleOTUcounts.bySiteandBatSpecies) <- c("tagID", "BatSpecies", "Site", "Counts")
  #what's the mean?
    mean(PerSampleOTUcounts.bySiteandBatSpecies$Counts)
  #about 13
  #what's the standard deviation?
    sd(PerSampleOTUcounts.bySiteandBatSpecies$Counts)
  #about 8
  #what's the median?
    median(PerSampleOTUcounts.bySiteandBatSpecies$Counts)
  #it's "10"  
      
#write all to disk:
setwd("~/Documents/Lab.Foster/guano/guano_bri/catp1p2/R_OTUcall_output/")
write.table(OTUcounts.perBatAndWOY, "OTUcounts.perBatAndWOY.txt", quote = F, sep = "\t", row.names = F, col.names = T)
write.table(OTUcounts.perBatAndSite, "OTUcounts.perBatAndSite.txt", quote = F, sep = "\t", row.names = F, col.names = T)
write.table(PerSampleOTUcounts.bySiteandBatSpecies, "PerSampleOTUcounts.bySiteandBatSpecies.txt", quote = F, sep = "\t", row.names = F, col.names = T)
write.table(PerSampleOTUcounts.byWOYandBatSpecies, "PerSampleOTUcounts.byWOYandBatSpecies.txt", quote = F, sep = "\t", row.names = F, col.names = T)


########################################################################################
############################    SECTION 4: data analysis  ##############################
########################################################################################


########################################################################################
################# PART A - if you want to get very species specific ####################
### will generate:
### 1. a data frame containing only those rows with taxa classified completely to species name ('BOLDspecies.df')
### 2a. a list of those species
### 2b. write that list to disk
### 3a. a table listing the frequency of each species (data set wide, not WOY or species or sex, etc. specific)
### 3b. write that table to disk

#subset out only cases in which species is known
BOLDspecies.df <- BOLDonly.df[!grepl("unknown", BOLDonly.df$species_name),]

#generate a df of just the 'species_name' and 'sampleID' fields
BOLDarthspeciesList <- subset(BOLDspecies.df, select = c("tagID", "species_name"))

#remove all duplicates to generate a list of the species
BOLDarthspeciesList <- BOLDarthspeciesList[!duplicated(BOLDarthspeciesList$species_name),]
BOLDarthspeciesList$tagID <- NULL
#save as text file
setwd("~/Documents/Lab.Foster/guano/guano_bri/catp1p2/R_OTUcall_output/")
write.table(BOLDarthspeciesList, "ACAonlyBOLDarthspeciesList.txt", row.names = F, quote = F, col.names = F)

#create a frequency table to determine the occurence of each species
BOLDarthspecies.counts <- as.data.frame(table(BOLDspecies.df$species_name))
colnames(BOLDarthspecies.counts) <- c("Arthropod_name", "counts")
#add in the underscore if you want to: species.counts <- as.data.frame(sapply(species.counts,gsub,pattern=" ",replacement="_"))

#save as text file:
setwd("~/Documents/Lab.Foster/guano/guano_bri/catp1p2/R_OTUcall_output/")
write.table(BOLDarthspecies.counts, "BOLDarthspecies.counts.txt", row.names = F, quote = F, col.names = F, sep = "\t")

#CLEANUP
rm(BOLDarthspeciesList)
rm(BOLDarthspecies.counts)
###############################   END PART A    #########################################
#########################################################################################


#########################################################################################
########### PART B - Broad analyses - OTUs per WOY, BatSpecies, or Sample ##############

#very basic: how frequently are taxa are called per WOY for a given {Family/Order} among all birds?
head(Taxa.perWOY.byFamily)
head(Taxa.perWOY.byOrder)
  #if BatSpecies is unknown:
  #OTUcounts.perBatAndWOY$BatSpecies <- as.character(OTUcounts.perBatAndWOY$BatSpecies)    
  #unnecessary if object imported with 'stringsAsFactors = F'
  #OTUcounts.perBatAndWOY[is.na(OTUcounts.perBatAndWOY)] <- "undetermined"

#how frequently are taxa called to a given BatSpecies among all WOYs?
#see SECTION 3, #1b
head(Taxa.perBatSpecies.byFamily)
head(Taxa.perBatSpecies.byOrder)

#how frequently is a given family identified at a unique WOY and unique bird species?
#see SECTION 2
head(percFamily.byWOYandBatSpecies)
head(percOrder.byWOYandBatSpecies)

#how many total OTUs were identified for each BatSpecies at each WOY?
#see object: "OTUcounts_perSample" from SECTION 3, #3a.
head(OTUcounts.perBatAndWOY)

#what kind of range existed for OTU detection on a per-sample basis across every BatSpecies at each Location? 
#see object: "OTUcounts_perSample" from SECTION 3, #3b.
head(PerSampleOTUcounts.bySiteandBatSpecies)


###############################   END PART B    #########################################
#########################################################################################


#########################################################################################
######### PART C - Factor-specific analyses: Sex, Age, and Reproductive Status ############

setwd("~/Documents/Lab.Foster/guano/guano_bri/catp1p2/R_OTUcall_output/")

#first thing to do is generate a new df of just unique tagIDs to determine counts for things like Sex, Age, and Reproductive Status
BOLDonly.dfSamples <- subset(BOLDonly.df, !duplicated(tagID))
#keep only relevant fields
BOLDonly.dfSamples <- BOLDonly.dfSamples[,c("tagID", "WOY", "WOY", "BatSpecies", "Sex", "Age", "ReproStat")]

#how many of each Sex are there?
SexTable.BOLDonly = as.data.frame(table(BOLDonly.dfSamples$Sex))
colnames(SexTable.BOLDonly) <- c("sex", "counts")
write.table(SexTable.BOLDonly, "SexTable.BOLDonly.txt", quote = F, sep = "\t", row.names = F, col.names = T)
  # 32 Males to 17 Females; likely enought for a meaningful comparison

#how many of each Age are there?
AgeTable.BOLDonly = as.data.frame(table(BOLDonly.dfSamples$Age))
colnames(AgeTable.BOLDonly) <- c("age", "counts")
write.table(AgeTable.BOLDonly, "AgeTable.BOLDonly.txt", quote = F, sep = "\t", row.names = F, col.names = T)
  #don't discard the "no-record" single sample, but if any Age-analysis is to be conducted make sure to drop that sample (bri.0765)

#how many of each Species are there?
BatSpeciesTable.BOLDonly = as.data.frame(table(BOLDonly.dfSamples$BatSpecies))
colnames(BatSpeciesTable.BOLDonly) <- c("BatSpecies", "counts")
write.table(BatSpeciesTable.BOLDonly, "BatSpeciesTable.BOLDonly.txt", quote = F, sep = "\t", row.names = F, col.names = T)
  #only single record for PESU, EPFU, and MYSE...
  #mostly MYLE (n=26), and MYLU (n=20)

#how many discrete sampling weeks ("WOY" == week of year) were recorded?
WOYTable.BOLDonly = as.data.frame(table(BOLDonly.dfSamples$WOY))
colnames(WOYTable.BOLDonly) <- c("WOY", "counts")
write.table(WOYTable.BOLDonly, "WOYTable.BOLDonly.txt", quote = F, sep = "\t", row.names = F, col.names = T)
#few records for earlier weeks (just n=9 samples for weeks 28, 29, and 30); most records (n=40) for weeks 35,36,37

#how many of each Reproductive Status characters are there?
ReproStatTable.BOLDonly = as.data.frame(table(BOLDonly.dfSamples$ReproStat))
colnames(ReproStatTable.BOLDonly) <- c("ReproStat", "counts")
write.table(ReproStatTable.BOLDonly, "ReproStatTable.BOLDonly.txt", quote = F, sep = "\t", row.names = F, col.names = T)
  #almost all are 'NR'... not likely worth while to investigate in a statistical way, but certainly coudl qualitatively compare; however there may be species-effects that confound this observation, and the sample sizes are likely too low if they are distributed between more than one bat species


###############################   END PART C    #########################################
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
dat.h1 = read.delim(file = "~/Documents/Lab.Foster/guano/guano_bri/catp1p2/R_OTUcall_output/PerSampleOTUcounts.bySiteandBatSpecies.txt", header = TRUE)
h1 <- ggplot(dat.h1, aes(x=Counts))
h1 + geom_histogram()

#set colours ahead of time
barfill <- "#4271AE"
barlines <- "#1F3552"

#make your plot
h1 + geom_histogram(binwidth = 5, colour = barlines, fill = barfill) +
  labs(title="Distribution of OTU frequency detections for AcadiaNP samples") +
  labs(x="Number of OTUs detected", y="Frequency") +
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
dat.b1 = read.delim(file = "~/Documents/Lab.Foster/guano/guano_bri/catp1p2/R_OTUcall_output/OTUcounts.perBatAndWOY.txt", header = TRUE)
#convert "WOY" variable to character class
dat.b1$WOY <- as.character(dat.b1$WOY)

#with only one Site, no need to create a specific fill setting. however if there are multiple Sites, you'd run the following code to determine the number of colours you'd need to set the fill setting (or, to create the appropriate palette):
  # length(unique(OTUcounts.perBatAndWOY$Site))

#chose your palette with the right number of colors. Because there are 6 unique WOYs here but only 1 site, we can use a stacked bar graph to show how collections vaired over time for a given bat species; the color options are using just two basic colors (red and blue) but altering hues to represent the weeks because there are really two splits in collection weeks for these samples: 
b1pal = c("#70cf54",
          "#879b60",
          "#cbd356",
          "#306072",
          "#92c5d9",
          "#50a9d1")

#set colours ahead of time for this basic example (no special fill required by some other factor)
barfill <- "#4271AE"
barlines <- "#1F3552"

#plot data 
b1 <- ggplot(dat.b1, aes(x = BatSpecies, y = Counts, fill = WOY)) +
  geom_bar(stat = "identity", colour=barlines) +
  scale_fill_manual(values = b1pal, name="Week Of Year") + 
  ylab("Number of times an OTU is detected") +
  xlab("Bat species") + 
  theme(axis.text.y = element_text(size = 8)) +
  theme(axis.text.x = element_text(size = 10)) +
  labs(title = "Number of OTUs detected per Bat Species at AcadiaNP using DADA2/AMPtk pipeline") +
  theme(legend.key.size = unit(4, "mm"))
b1

## so clearly, most of the OTUs we're detecting are in the two bat species we're identifying.
## in addition, you can see that the bulk of data was collecetd in weeks 35-37, though quite a few OTUs were detected for MYLE earlier in the season (especially Week 29)
##################################################################################################


#####   Plots-b2a and b2b: For a given bat species, what proportions of arthropod orders are they eating?   ######
#load in the data
dat.b2 = read.delim(file = "~/Documents/Lab.Foster/guano/guano_bri/catp1p2/R_OTUcall_output/percOrder.bySiteandBatSpecies.txt", header = TRUE)

#the problem with these data is that we have an unequal distribution of BatSpecies samples:
  #just one of MYSE, EPFU, and PESU each
  #about 20 each of MYLE and MYLU
#we're going to split up these data into two separate plots so that you can view the snapshot of the singleton samples, while also then making more meaningful comparisons of the two larger BatSpecies (MYLE and MYLU)

#subset just the MYLE and MYLU
dat.b2a <- dat.b2[!grepl("EPFU|MYSE|PESU", dat.b2$BatSpecies),]
#subset just the MYSE, PESU, and EPFU
dat.b2b <- dat.b2[grepl("EPFU|MYSE|PESU", dat.b2$BatSpecies),]
                   
#to create the appropriate palette you need to determine how many factors (Arthropod orders) there are for each plot:

#for the large subset (MYLE and MYLU)
length(unique(dat.b2a$ArthropodOrder))
  #there are 17!  

#for the small subset (EPFU, MYSE, and PESU)
length(unique(dat.b2b$ArthropodOrder))
  #there are 7
  
#A few precautions
#1. There are 18 total Orders, so there is a unique order in 'b2b' that isn't in 'b2a'. Be mindful when assigning this color.
#2. We need to keep colors consistent between both graphs, but we also want to have color patterns that make sense (we should easily demonstrate the larger counts for each graph).
#My current solution is to darken the hue for the majority of colors which show low count totals for either large or small figure, then use bright and distinct colors for the more frequent ones (think beetle, moth, and dipteran).
  #these are the initial 18 developed:
  bigpal <- c("#a68741","#00aaf5","#97d271","#ff584e","#046357","#29101e","#403e21","#e4d6ff","#ba6300","#596fac","#6e0016","#598294","#7d3226","#462570","#996c5e","#491f3c","#9b5e8e","#32414e")
    #notably, just three colors ("#00aaf5","#97d271","#ff584e") are particularly bright relative to all others; these should be assigned to the predominant Orders detected in our datasets (Diptera, Coleoptera, Lepidoptera).
  

  
  #we're going to construct our palettes for each graph:
palb2a <- c("#fffe9f", 
            "#627a84", 
            "#97d271", 
            "#ff584e", 
            "#627a84", 
            "#627a84", 
            "#d7a257",
            "#627a84", 
            "#8aade4", "#627a84", "#627a84", "#627a84", "#627a84","#627a84","#627a84","#ea9de9","#627a84")

palb2b <- c("#fffe9f", 
            "#627a84", 
            "#97d271", 
            "#ff584e", 
            "#8aade4", "#32414e", "#627a84")
  #note that color ("#32414e") must be assigned to the unique Order 'dat.b2a' (=="Megaloptera")
  #note that colors () must be assigned to Orders Coleoptera, Diptera, and Lepidoptera, respectively.

#set colour for lines around bars as usual:  
barlines <- "#1F3552"

#facet-names - used to relabel the facets to show samples sizes contributing to the resulting illustration
b2a.facetnames <- c(MYLE = "MYLE, n = 26", MYLU = "MYLU, n = 20")

b2a <- ggplot(dat.b2a, aes(x = ArthropodOrder, y = RelativeCounts, fill = ArthropodOrder)) + 
  geom_bar(stat = "identity", colour=barlines) +
  scale_fill_manual(values = palb2a, name="Arthropod Order") + 
  guides(col = guide_legend(nrow = 1)) + 
  theme(legend.position = "none") +
  theme(axis.text.y = element_text(size = 8)) +
  theme(axis.text.x = element_text(size = 10)) +
  labs(title = "Total Counts of Arthropod Order detected from MYLE and MYLU at AcadiaNP using BOLD database") +
  labs(y="Total Counts of each Arthropod Order Detected") +
  theme(legend.key.size = unit(4, "mm")) +
  facet_grid(~BatSpecies, labeller=labeller(BatSpecies=b2a.facetnames)) +
  scale_y_continuous(breaks = c(0,10,20,50,100,150)) +
  theme(axis.title.x = element_blank()) +
  #theme(axis.text.x=element_blank() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
b2a

#so for MYLE and MYLU, it's clear that Dipterans and Lepidopterans dominate their feeding strategy


### what about the singleton EPFU, MYSE, and PESU samples?

#facet-names - used to relabel the facets to show samples sizes contributing to the resulting illustration:
b2b.facetnames = c(EPFU = "EPFU, n=1", MYSE = "MYSE, n=1", PESU = "PESU, n=1")

b2b <- ggplot(dat.b2b, aes(x = ArthropodOrder, y = percent_taxa, fill = ArthropodOrder)) + 
  geom_bar(stat = "identity", colour=barlines) +
  scale_fill_manual(values = palb2b, name="Arthropod Order") + 
  guides(col = guide_legend(nrow = 1)) + 
  theme(legend.position = "none") +
  theme(axis.text.y = element_text(size = 8)) +
  theme(axis.text.x = element_text(size = 10)) +
  labs(title = "Percent Arthropod Order detected from single representative sample of \n EPFU, MYSE, and PESU individuals at AcadiaNP using BOLD database") +
  labs(y="Percent of each Arthropod Order Detected") +
  theme(legend.key.size = unit(4, "mm")) +
  facet_grid(~BatSpecies, labeller=labeller(BatSpecies=b2b.facetnames)) +
  scale_y_continuous(breaks = c(0,10,20,50,100,150)) +
  theme(axis.title.x = element_blank()) +
  #theme(axis.text.x=element_blank() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
b2b

###################################        END Plot 1        #####################################
##################################################################################################
#load data
setwd("~/Documents/Lab.Foster/guano/guano_bri/catp1p2/R_OTUcall_output/")
Ordercounts_table <- read.delim("Ordercounts_table.txt", header = TRUE)

#subset only MYLE and MYLU samples for this plot:
dat.b3a <- Ordercounts_table[!grepl("EPFU|MYSE|PESU", Ordercounts_table$BatSpecies),]

#facet-names - used to relabel the facets to show samples sizes contributing to the resulting illustration
topfacetnames <- c(MYLE = "MYLE, n = 26", MYLU = "MYLU, n = 20")

#set colour for lines around bars as usual:  
barlines <- "#1F3552"

#Check there are the same number of arthropod orders as in b2a plots earlier:
length(unique(dat.b3a$ArthropodOrder))
  #17 still. So use the same "palb2a" to be consistent when filling in color for each column.

palb2a <- c("#fffe9f", "#627a84", "#97d271", "#ff584e", "#627a84", "#627a84", "#d7a257","#627a84", "#8aade4", "#627a84", "#627a84", "#627a84", "#627a84","#627a84","#627a84","#ea9de9","#627a84")


b3a <- ggplot(dat.b3a, aes(x = ArthropodOrder, y = RelativeCounts, fill = ArthropodOrder)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = palb2a, name="Arthropod Order") + 
  guides(col = guide_legend(nrow = 1)) + 
  theme(legend.position = "none") +
  theme(axis.text.y = element_text(size = 8)) +
  theme(axis.text.x = element_text(size = 10)) +
  labs(title = "Total Counts per Arthropod Order detected from MYLE and MYLU \n for each Week at AcadiaNP using BOLD database") +
  labs(y="Counts of each Arthropod Order Detected") +
  theme(legend.key.size = unit(4, "mm")) +
  facet_grid(WOY~BatSpecies, labeller=labeller(BatSpecies=topfacetnames)) +
  scale_y_continuous(breaks = c(0,20,40,60,80, 100)) +
  theme(axis.title.x = element_blank()) +
  #theme(axis.text.x=element_blank() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
b3a


#####################################################################################################################

##################################################################################################
##############################      Part B: Treemaps     #########################################

library(treemap)

treemap_df <- BOLDonly.df
treemap_df <- as.data.frame(treemap_df)
treemap_df$counts <- 1

#recall how many arthropod Orders are in this dataset:
length(unique(treemap_df$order_name))
  #18 total

#we're going to use the original "bigpal" palette described in the ggplot setup above. This is the palette:
bigpal <- c("#a68741","#00aaf5","#97d271","#ff584e","#046357","#29101e","#403e21","#e4d6ff","#ba6300","#596fac","#6e0016","#598294","#7d3226","#462570","#996c5e","#491f3c","#9b5e8e","#32414e")

#option2a - look at ORDER eaten among BatSpecies
treemap(treemap_df, 
        index=c("BatSpecies", "order_name"),
        vSize="counts",
        vColor="order_name",
        type="categorical",
        title = "Total Counts of Arthropod Order by BatSpecies for all dates",
        title.legend = "Taxonomic Order",
        position.legend = "right",
        align.labels=list(c("center", "top"), c("left", "bottom")), 
        fontsize.labels=c(14,6),
        lowerbound.cex.labels=0.2, 
        overlap.labels=0.4,
        palette=bigpal
)

#further split into individual weeks (need to fix the number as it's kind of messy)
treemap(treemap_df, 
        index=c("BatSpecies", "order_name", "WOY"),
        vSize="counts",
        vColor="order_name",
        type="categorical",
        title = "Total Counts of Arthropod Order by BatSpecies \n for each Week of Year sampled",
        title.legend = "Taxonomic Order",
        position.legend = "right",
        align.labels=list(c("center", "top"), c("left", "bottom")), 
        fontsize.labels=c(12,1,9),
        fontcolor.labels = c("#272627", "#272627", "#272627"),
        bg.labels = "#00000000",
        lowerbound.cex.labels=0.2, 
        overlap.labels=0.2,
        palette=bigpal,
        fontsize.title = 20
)

  #note that the "bg.labels" sets the background color to completely transparent!

#option2b - look at FAMILY eaten among bats grouped by location
#first, reorder:
treemap_df <- treemap_df[with(treemap_df, order(family_name, order_name)), ] 

treemap(treemap_df, 
        index=c("WOY", "family_name"),
        vSize="counts",
        vColor="order_name",
        type="categorical",
        title = "Numbers of Family counts per WOY colored by taxonomic Order",
        title.legend = "Taxonomic Order",
        position.legend = "right",
        align.labels=list(c("center", "top"), c("left", "bottom")), 
        fontsize.labels=c(14,6),
        lowerbound.cex.labels=0.2, 
        overlap.labels=0.4,
        palette=bigpal,
        height = 200
)
###############################   END PART B    #########################################
#########################################################################################
