#amptk_2017_OTUtable_analysis.R
#written 10.12.16
#modified 4.10.17
#Devon O'Rourke
#amptk analysis for mud samples from Montana

########################################################################################
############################    SECTION 0: metadata files  #############################
########################################################################################

#################### process started: creating metadata file  ####################
date()

#import collection data
setwd("~/Documents/Lab.Foster/guano/mudMT/")
meta_df <- read.csv(file = "sampleinfo.csv", stringsAsFactors = F)

#drop the WellNumber field
meta_df$WellNumber <- NULL

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

#write updated meta_df data to disk
setwd("~/Documents/Lab.Foster/guano/mudMT/p2_data/")
write.table(meta_df, file = "meta_df.txt", col.names = T, row.names = F, sep = "\t", quote = F)

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
setwd("~/Documents/Lab.Foster/guano/mudMT/p2_data/")

#import taxonomy-classified OTU Table:
amptk_otu_tableA <- read.delim("mud.otu_table.taxonomy.txt", stringsAsFactors = F)
#rename that first column
colnames(amptk_otu_tableA)[1] <- "OTUid"

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
setwd("~/Documents/Lab.Foster/guano/mudMT/p2_data/")
meta_df = read.delim(file = "meta_df.txt", header = T, stringsAsFactors = F)

#rename that first column
colnames(meta_df)[1] <- "tagID"

#replace the hyphen with a period in the "tagID" column (may not be necessary - check both data.frames for consistency in delimiter)
#metadata$tagID <- sub("-", "\\.", metadata$tagID)

#merge OTU calls with metadata by 'tagID' (really just the sample ID)
masterdf = merge(amptk_dt_unique, meta_df, by = "tagID")

#write to disk:
setwd("~/Documents/Lab.Foster/guano/mudMT/p2_data/")
write.table(masterdf, file = "masterdf.txt", col.names = T, row.names = F, sep = "\t", quote = F)

#########################################################
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

#this gets rid of the SINTAX and UTAX assigned samples
BOLDonly.df = masterdf[!grepl("SINTAX", masterdf$OTU_identifier),]
BOLDonly.df = BOLDonly.df[!grepl("UTAX", BOLDonly.df$OTU_identifier),]
#and this gets rid of the samples that are likely the result of the mock community bleeding into the true samples (another index-bleed happenstance)
BOLDonly.df = BOLDonly.df[!grepl("CFMR", BOLDonly.df$OTU_identifier),]
#this gets rid of any unwanted other taxa, such as bat reads, which shouldn't be there...
BOLDonly.df = BOLDonly.df[!grepl("Mammalia", BOLDonly.df$class_name),]

#rename the first column
colnames(BOLDonly.df)[1] <- "SampleID"

#it's worth renaming empty (N/A) values for taxonomic classes with "unknown" for future plotting purposes (we'll want to show in our data how many OTUs weren't explicitly described relative to those that were)
BOLDonly.df[is.na(BOLDonly.df)] <- "unknown"
  #careful with this term, as it will convert any 'N/A' in the dataset, regardless of what field it belongs to, to "unknown"

#write object 'masterdf' to disk
setwd("~/Documents/Lab.Foster/guano/mudMT/p2_data/")
write.table(BOLDonly.df, file = "BOLDonly.df.txt", col.names = T, row.names = F, sep = "\t", quote = F)

#CLEANUP
rm(meta_df, amptk_dt_unique, masterdf, batOTUs.df)  

#################### process ended: OTU table import and metadata merging  ####################
############### BOLDonly.df file contains all information for subsequent work  ################