#ufits_2016_OTUtable_analysis.R
#written 10.12.16
#modified 1.7.17
#Devon O'Rourke
#ufits analysis from 2016 data but with modified BLAST output

########################################################################################

#import packages
library(tidyr)
library(reshape2)
library(data.table)

#set working directory
setwd("~/Documents/Lab.Foster/guano/guano_bri/bri-chunk1/")

#import sheet

ufits_otu_tableA <- read.delim("modOTU_table.taxonomy.txt", stringsAsFactors = F)

########################################################################################
############################    SECTION 1: data wrangline  #############################
########################################################################################


#transform data
ufits_otu_tableA[ufits_otu_tableA == 0 ] <- NA

ufits_otu_tableA$OTUnTax <- paste(ufits_otu_tableA$OTUid, ufits_otu_tableA$Taxonomy, sep = ";")

ufits_otu_tableA$OTUid <- NULL

ufits_otu_tableA$Taxonomy <- NULL

ufits_dfA <- melt(ufits_otu_tableA, id = "OTUnTax")

colnames(ufits_dfA) <- c("OTUnTax", "tagID", "binaryPresence")

ufits_dfA <- subset(ufits_dfA, binaryPresence >= 1)

ufits_dfA$binaryPresence <- NULL


#separate 'Taxonomy' field into respective taxa levels:
ufits_dfA <- separate(data = ufits_dfA, col = OTUnTax, into = c("OTUid", "OTU_identifier", "kingdom_name", "phylum_name", "class_name", "order_name", "family_name", "genus_name", "species_name"), sep = "\\;|,")


#remove the prefixes leading into the taxa names:
ufits_dfA$kingdom_name <- sub("k:", "", ufits_dfA$kingdom_name)
ufits_dfA$phylum_name <- sub("p:", "", ufits_dfA$phylum_name)
ufits_dfA$class_name <- sub("c:", "", ufits_dfA$class_name)
ufits_dfA$order_name <- sub("o:", "", ufits_dfA$order_name)
ufits_dfA$family_name <- sub("f:", "", ufits_dfA$family_name)
ufits_dfA$genus_name <- sub("g:", "", ufits_dfA$genus_name)
ufits_dfA$species_name <- sub("s:", "", ufits_dfA$species_name)


#add in a designated field to show which data.frame belongs to which library (oro##-1 or oro##-2)
#ufits_dfA$library <- "library01"

#make list of data frames to merge
l = list(ufits_dfA)

#use 'rbindlist' function in 'data.table' package to join two data frames
ufits_dfA <- rbindlist(l, use.names = TRUE, fill = TRUE)

#remove redundant calls
  #this will cut out a lot of instances in which you have poorly described taxonomies for a sample using UTAX as the 'OTU_identifier'
  #as a precaution, what you can say is that you have confidence in eliminating reads in which we can't associate it to a unique 'phlya',...
    #so what's the point in keeping it anyway?
  #it will keep a single instance in which some value is described at just the phylum level, yet this may actually be cutting out things like our bat reads
    #so we'll keep these reads that are cut in the "dt_dupilc" data frame; you can then use the OTUId from each of these and go back and search for those in the 
      #original fasta file, and see if any of them are things like bat DNA sequences

#using 'data.table'
setkey(ufits_dfA, OTU_identifier, tagID, species_name, genus_name, family_name, order_name, class_name)
ufits_dt_unique <- ufits_dfA[!duplicated(ufits_dfA)]   
  #we're pulling out repeated calls just one situation with this dataset:
  ##1. you've named something where there are multiple unique OTUs assigned to an identical taxonomic description; there isn't a point ...
  ## in double counting these because we're collapsing a taxonomic description to a binary presence/absence (without this filter you...
  ## end up adding duplicate counts for a single species per sample)
ufits_dt_duplic <- ufits_dfA[duplicated(ufits_dfA)]    
  #these represent redundant calls explained above. in some datasets (not this one) these can also represent situations where there may be MULTIPLE DISTINCT...
  ## taxa you could assign IF the taxonomic description is only classified through to a "phylum" level (and "class", "order", .... "species" are all "NA")...
  ## you could look into identifying if these are truly different, but if they are so poorly classified to begin with it's likely you won't have enough...
  ## percent identity in an alignment score to definitively improve the calls. the only way it'll improve is if your existing databases just missed the read...
  ## You can look at the unique calls in which there is no 'phylum_name' known to the database using the following (comment muted) commands:
    ##nophyla_df <- ufits_dt_unique[is.na(ufits_dt_unique$phylum_name),]
    ##setkey(nophyla_df, OTUid)
    ##nophyla_df_unique <- nophyla_df[!duplicated(nophyla_df)]   
  ##these represent all the unique per OTUId between both libraries; note that the same call may occur between the two libraries but...
  ##contain different OTUId's because of how they were independently assigned from the UFITS pipeline

#CLEANUP
rm(ufits_dfA, ufits_otu_tableA, l)


#################################################################################################
###########################   SECTION 2:  Add in metadata     ###################################
#################################################################################################

#first create a list of the OTU's used in this study - just make a new dataframe from the 'tagID' column
#then search through the 'collection_data' file in Google Drive and pull out the samples that match
#that new file, called "{something}...metadata.txt" is what we're importing to start here...

#import metadata
metadata = read.csv(file = "~/Documents/Lab.Foster/guano/guano_metadata_2016.csv",
                    header = T, stringsAsFactors = F)

#merge OTU calls with metadata by 'tagID' (really just the sample ID)
metaPlusOTUdf = merge(ufits_dt_unique, metadata, by = "tagID")

#drop out any reads you don't want to include in subsequent analyses:
  #this gets rid of the parasitic mites
  metaPlusOTUdf = metaPlusOTUdf[!grepl("Trombidiformes", metaPlusOTUdf$order_name),]
  #this gets rid of the bat reads
  metaPlusOTUdf = metaPlusOTUdf[!grepl("Chiroptera", metaPlusOTUdf$order_name),]

  #note that we're keeping the 'Unknown' values with the "NA" for the 'order_name' field ... we should substitute those NA values though...
  metaPlusOTUdf$order_name <- as.character(metaPlusOTUdf$order_name)
  metaPlusOTUdf[is.na(metaPlusOTUdf)] <- "unknown"
  
  #if necessary, you're going to want to get rid of the instances in which a bat species isn't listed (was 'NA', now is 'unknown'):
    #metaPlusOTUdf = metaPlusOTUdf[!grepl("unknown", metaPlusOTUdf$BatSpecies),]

#write object 'metaPlusOTUdf' to disk
    write.table(metaPlusOTUdf, file = "metaPlusOTUdf.txt", col.names = T, row.names = F, sep = "\t", quote = F)

#CLEANUP
rm(metadata, ufits_dt_unique)  

##################################################################################################
###################     SECTION 3: Calculate Frequencies for plots           #####################
##################################################################################################

#need to generate a suite of data tables for each unique plotting scenario: 
## 1. you want to know the relative frequency (in percentage terms) of each order eaten when you aggregate all counts for:
##  a. All bats at each unique site, or
##  b. All sites for each unique bat species
## 2. What'll likely get used more for analyses: the relative frequency of each order eaten when aggregating counts for:
### a. each unique species of bat at each unique site

#In addition, we're going to generate a few other frequency tables that may prove useful:
## 3a. a frequency table showing how many total OTUs were identified (regardless of taxonomic class) for each BatSpecies at each LocationName
## 3b. a frequency table showing how many unique OTUs were identified (by taxonomic 'order') for each sample, for each BatSpeacies at each LocationName

#first, make a table with the number of counts per order name across each bat species and site that's unique
library(plyr)
counts_table <- ddply(metaPlusOTUdf, .(metaPlusOTUdf$BatSpecies, metaPlusOTUdf$LocationName, metaPlusOTUdf$order_name), nrow)
names(counts_table) <- c("BatSpecies", "LocationName", "order_name", "Freq")

#next, sum up the frequency values given one of the two situations described above:
##1a. For each unique site, grouping all bats (ie. how many times did we see an order of insect get eaten at a given site, regardless of bat?)
library(dplyr)
all_bats_orderCounts <- counts_table %>% 
  group_by(LocationName, order_name) %>% 
  summarise(Freq = sum(Freq))
names(all_bats_orderCounts) <- c("LocationName", "order_name", "relCounts")

allOTUs_per_LocationName <- counts_table %>% 
  group_by(LocationName) %>% 
  summarise(Freq = sum(Freq))
names(allOTUs_per_LocationName) <- c("LocationName", "TotalCounts")

#then merge the two tables by "LocationName"
relative_OTUs_perLocation <- merge(all_bats_orderCounts, allOTUs_per_LocationName, by = "LocationName")
#calculate the relative frequency per insect order per Location
relative_OTUs_perLocation$percent_taxa <- (relative_OTUs_perLocation$relCounts / relative_OTUs_perLocation$TotalCounts * 100)

##1b. For each unique bat species, grouping all sites (ie. how many times did we see an order of insect get eaten by a given bat, regardless of site?)
all_sites_orderCounts <- counts_table %>% 
  group_by(BatSpecies, order_name) %>% 
  summarise(Freq = sum(Freq))
names(all_sites_orderCounts) <- c("BatSpecies", "order_name", "relCounts")

allOTUs_per_Bat <- counts_table %>% 
  group_by(BatSpecies) %>% 
  summarise(Freq = sum(Freq))
names(allOTUs_per_Bat) <- c("BatSpecies", "TotalCounts")

#then merge the two tables by "LocationName"
relative_OTUs_perBatSpecies <- merge(all_sites_orderCounts, allOTUs_per_Bat, by = "BatSpecies")
#calculate the relative frequency per insect order per Location
relative_OTUs_perBatSpecies$percent_taxa <- (relative_OTUs_perBatSpecies$relCounts / relative_OTUs_perBatSpecies$TotalCounts * 100)


##2. For each unique bat species at each unique site
  ## this already exists - it's the "counts_table" object created initially!!...

allOTUs_by_BatAndSite <- counts_table %>% 
  group_by(BatSpecies, LocationName) %>% 
  summarise(Freq = sum(Freq))
names(allOTUs_by_BatAndSite) <- c("BatSpecies", "LocationName", "TotalCounts")

#merge the original 'counts_table' object with the "TotalCounts" field in 'allOTUs_by_BatAndSite' object with R package 'dplyr'
library(dplyr)
relative_OTUs_perBatAndSite <- counts_table %>% inner_join(allOTUs_by_BatAndSite)
names(relative_OTUs_perBatAndSite) <- c("BatSpecies", "LocationName", "order_name", "relCounts", "TotalCounts")
#calculate relative frequency per insect order per Location AND BatSpecies
relative_OTUs_perBatAndSite$percent_taxa <- (relative_OTUs_perBatAndSite$relCounts / relative_OTUs_perBatAndSite$TotalCounts * 100)


##3a. this calculates the number of OTUs called per for each BatSpecies at each Location (taxonomic level not specified)
library(plyr)
OTUcounts_perBatAndSite <- ddply(metaPlusOTUdf, .(metaPlusOTUdf$BatSpecies, metaPlusOTUdf$LocationName), nrow)
names(OTUcounts_perBatAndSite) <- c("BatSpecies", "LocationName", "Freq")

##3b. this calculates the number of OTUs called per sample of guano while also retaining the BatSpecies and LocationName information
library(plyr)
OTUcounts_perSample <- ddply(metaPlusOTUdf, .(metaPlusOTUdf$tagID, metaPlusOTUdf$BatSpecies, metaPlusOTUdf$LocationName), nrow)
names(OTUcounts_perSample) <- c("tagID", "BatSpecies", "LocationName", "Freq")

#CLEANUP
rm(counts_table, all_bats_orderCounts, allOTUs_per_LocationName, all_sites_orderCounts, allOTUs_per_Bat, allOTUs_by_BatAndSite)


################## a few additional subsetting actions you could take   ##################

#alternatively, you could subset by 'LocationName'. We're going to pull out just the Acadia ME and the Ely Mine VT samples for this instance:
ACA_samples <- subset(metaPlusOTUdf, metaPlusOTUdf$LocationName == "ACA")
ELY_samples <- subset(metaPlusOTUdf, metaPlusOTUdf$LocationName == "ELY")
MSH_samples <- subset(metaPlusOTUdf, metaPlusOTUdf$LocationName == "MSH")
SAG_samples <- subset(metaPlusOTUdf, metaPlusOTUdf$LocationName == "SAG")
SRY_samples <- subset(metaPlusOTUdf, metaPlusOTUdf$LocationName == "SRY")
VRA_samples <- subset(metaPlusOTUdf, metaPlusOTUdf$LocationName == "VRA")

#calculate how many of each sample are represented within those subsets:
#for Acadia, ME
length(unique(ACA_samples$tagID))
  #10 individuals

#for Ely Mine, VT
length(unique(ELY_samples$tagID))
  #133 individuals

#for Mashpee, MA
length(unique(MSH_samples$tagID))
  #21 individuals

#for Sag Harbor, NY
length(unique(SAG_samples$tagID))
  #3 individuals

#for Shriley, NY
length(unique(SRY_samples$tagID))
  #12 individuals  

#for Yorktown Naval Base, VA
length(unique(VRA_samples$tagID))
  #1 individual

#you can also subset data by bat species if you want to assume that bats have an opportunity to eat the same things at all these locations... likely?? no.
MYLE_df <- subset(metaPlusOTUdf, metaPlusOTUdf$BatSpecies == "MYLE")
LABO_df <- subset(metaPlusOTUdf, metaPlusOTUdf$BatSpecies == "LABO")
EPFU_df <- subset(metaPlusOTUdf, metaPlusOTUdf$BatSpecies == "EPFU")
MYLU_df <- subset(metaPlusOTUdf, metaPlusOTUdf$BatSpecies == "MYLU")
MYSE_df <- subset(metaPlusOTUdf, metaPlusOTUdf$BatSpecies == "MYSE")

#build species-specific frequency tables
MYLE_freq <- as.data.frame(table(MYLE_df$species_name))
LABO_freq <- as.data.frame(table(LABO_df$species_name))
EPFU_freq <- as.data.frame(table(EPFU_df$species_name))
MYLU_freq <- as.data.frame(table(MYLU_df$species_name))
MYSE_freq <- as.data.frame(table(MYSE_df$species_name))



########################################################################################
############################    SECTION 4: data analysis  ##############################
########################################################################################


########################################################################################
################# PART A - if you want to get very species specific ####################
### will generate:
  ### 1. a data frame containing only those rows with taxa classified completely to species name ('ufits_species_df')
  ### 2a. a list of those species
  ### 2b. write that list to disk
  ### 3a. a table listing the frequency of each species (data set wide, not site or species or sex, etc. specific)
  ### 3b. write that table to disk

#subset out only cases in which species is known
ufits_species_df <- metaPlusOTUdf[!grepl("unknown", metaPlusOTUdf$species_name),]

#generate a df of just the 'species_name' and 'sampleID' fields
ufits_speciesLIST <- subset(ufits_species_df, select = c("tagID", "species_name"))

#remove all duplicates to generate a list of the species
ufits_species_list <- ufits_speciesLIST[!duplicated(ufits_speciesLIST$species_name),]
ufits_species_list$tagID <- NULL
#save as text file
write.table(ufits_species_list, "updatedbriRound1_specieslist.txt", row.names = F, quote = F, col.names = F)

#create a frequency table to determine the occurence of each species
freq_species <- as.data.frame(table(ufits_species_df$species_name))
  #add in the underscore if you want to: freq_species <- as.data.frame(sapply(freq_species,gsub,pattern=" ",replacement="_"))

#save as text file
write.table(freq_species, "updatedbriRound1_sp_freqtable.txt", row.names = F, quote = F, col.names = F, sep = "\t")

#CLEANUP
rm(ufits_species_df, ufits_species_list, ufits_speciesLIST)
rm(freq_species)
###############################   END PART A    #########################################
#########################################################################################


#########################################################################################
########### PART B - Broad analyses - OTUs per Site, BatSpecies, or Sample ##############

#very basic: how many OTUs are called per bat species per site?
#see object: "relative_OTUs_perBatAndSite" from SECTION 3, #1
head(relative_OTUs_perBatAndSite)
  #if bat species is unknown:
  #relative_OTUs_perBatAndSite$BatSpecies <- as.character(relative_OTUs_perBatAndSite$BatSpecies)    
    #unnecessary if object imported with 'stringsAsFactors = F'
  #relative_OTUs_perBatAndSite[is.na(relative_OTUs_perBatAndSite)] <- "undetermined"

#how many times is an OTU/species called to a given bat species for a given order?
#see object: "relative_OTUs_perBatAndSite" from SECTION 3, #2
head(relative_OTUs_perBatSpecies)

#how many total OTUs were identified for each BatSpecies at each LocationName?
#see object: "OTUcounts_perSample" from SECTION 3, #3a.
head(OTUcounts_perBatAndSite)

#what kind of range existed for OTU detection on a per-sample basis across every BatSpecies at each Location? 
#see object: "OTUcounts_perSample" from SECTION 3, #3b.
head(OTUcounts_perSample)

###############################   END PART B    #########################################
#########################################################################################


#########################################################################################
######### PART C - Site-specific analyses: Sex, Age, and Reproductive Status ############

#working with ELY_sample object only... (has by far most individual samples)
#first thing to do is generate a new df of just unique tagIDs to determine counts for things like Sex, Age, and Reproductive Status
ELY_uniqueSamples <- subset(ELY_samples, !duplicated(tagID))

#how many of each Sex are there?
sex_table = as.data.frame(table(ELY_uniqueSamples$Sex))
  #drop out single column with unknown sex for any sex-based analysis

#how many of each Age are there?
age_table = as.data.frame(table(ELY_uniqueSamples$Age))
  #way more adults than juvenilles... not likely a useful category for analyses

#how many of each reproductive status are there?
## first likely need to separate out males and females
## females only
reproStat_Fdf <- subset(ELY_uniqueSamples, Sex == "F")
  # and the table of counts
  reproStat_Ftable = as.data.frame(table(reproStat_Fdf$ReproStat))

  ## males only
reproStat_Mdf <- subset(ELY_uniqueSamples, Sex == "M")
  # and the table of counts
  reproStat_Mtable = as.data.frame(table(reproStat_Mdf$ReproStat))
    # likely have the power to analyze at this level


###############################   END PART C    #########################################
#########################################################################################



##################################################################################################
########################      SECTION 5: Making the plots     ####################################
##################################################################################################

##################################################################################################
##############################      Part A: Treemaps     #########################################

library(treemap)

treemap_df <- metaPlusOTUdf
treemap_df <- as.data.frame(treemap_df)
treemap_df$counts <- 1

#option1 - look at genus eaten among all bats in study, grouped by taxonomic order
treemap(treemap_df, 
        index=c("order_name", "genus_name"),
        vSize="counts",
        vColor="family_name",
        type="categorical",
        #title = "Numbers of counts OTU at genera level detected",
        #title.legend = "genus name",
        position.legend = "none",
        align.labels=list(c("center", "top"), c("left", "bottom")), 
        fontsize.labels=c(14,6),
        lowerbound.cex.labels=0.2, 
        overlap.labels=0.4,
        palette="RdBu"
        
)


#option2a - look at ORDER eaten among bats grouped by location
treemap(treemap_df, 
        index=c("LocationName", "order_name"),
        vSize="counts",
        vColor="order_name",
        type="categorical",
        #title = "Numbers of counts OTU at order level detected",
        #title.legend = "genus name",
        position.legend = "none",
        align.labels=list(c("center", "top"), c("left", "bottom")), 
        fontsize.labels=c(14,6),
        lowerbound.cex.labels=0.2, 
        overlap.labels=0.4,
        palette="Set1"
        
)

#option2b - look at SPECIES eaten among bats grouped by location
treemap(treemap_df, 
        index=c("LocationName", "species_name"),
        vSize="counts",
        vColor="order_name",
        type="categorical",
        #title = "Numbers of counts OTU at species level detected",
        #title.legend = "genus name",
        position.legend = "bottom",
        align.labels=list(c("center", "top"), c("left", "bottom")), 
        fontsize.labels=c(14,6),
        lowerbound.cex.labels=0.2, 
        overlap.labels=0.4,
        palette="Set1"
        
)


#option3a - look at orders eaten among bats at a single location (say in Acadia NP)
aca_only <- subset(metaPlusOTUdf, metaPlusOTUdf$LocationName == "aca")
aca_only[is.na(aca_only)] <- "unknown"
treemap_aca <- aca_only
treemap_aca <- as.data.frame(treemap_aca)
treemap_aca$counts <- 1

treemap(treemap_aca, 
        index=c("BatSpecies", "order_name"),
        vSize="counts",
        vColor="order_name",
        type="categorical",
        title = "Numbers of counts OTU at species level detected in Acadia NP",
        title.legend = "taxonomic order",
        position.legend = "none",
        align.labels=list(c("center", "top"), c("left", "bottom")), 
        fontsize.labels=c(14,6),
        lowerbound.cex.labels=0.2, 
        #fontsize.legend = 6,
        overlap.labels=0.4,
        #palette=palette14
        palette=palette14,
        
)

###############################   END PART A    #########################################
#########################################################################################



#########################################################################################
##########################           PART B - ggplots        ############################

#load ggplot then play around with formatting
library(ggplot2)
require(graphics)
require(RColorBrewer)

###for most plots you're going to want to tailor the color scheme... there's a bunch of places to read up on choices:
  # see here for examples: http://www.r-bloggers.com/choosing-colour-palettes-part-ii-educated-choices/)
  # if you're using "Paired" as the palette it's important to know how many colors you can pull from each palette...
  ## ... you want the max possible for the qualitative type
####key thing here: specify the number of unique colors you want to pull form that palette.
  # and here's yet another palette option... go here: http://tools.medialab.sciences-po.fr/iwanthue/


tinyPalette = colorRampPalette(brewer.pal(6, "Paired"))(6)




##################################################################################################
#####   Plot1: How many times is an OTU detected per BatSpecies at a given LocationName?   ######
p <- ggplot(data = OTUcounts_perBatAndSite)

#we're going to be filling by LocationName, so how many locations are there? there are 6... see it here:
length(unique(OTUcounts_perBatAndSite$LocationName))

#chose your palette with the right number of colors
p1palette = c("#77e1f0","#58b31e","#aeb6ff","#ffc732","#1f658b","#e49e92")

#plot data 
r <- p +
  aes(x = BatSpecies, y = Freq, fill = LocationName) +
  #geom_bar(stat = "identity", width = 0.4) +
  geom_bar(stat = "identity", color="dark gray", width = 0.7, size = 0.3) +
  scale_fill_manual(values = p1palette) + 
  guides(col = guide_legend(nrow = 10, byrow = TRUE)) + 
  theme(legend.position = "right") + 
  ylab("Number of times an OTU is detected") +
  xlab("Bat species") + 
  theme(axis.text.y = element_text(size = 8)) +
  theme(axis.text.x = element_text(size = 10)) +
  labs(title = "Number of OTUs detected per Bat Species at each Location using UFITS/DADA2 pipeline") +
  theme(legend.key.size = unit(4, "mm"))
  
r

###################################        END Plot 1        #####################################
##################################################################################################



#####################################################################################################################
####  Plot2: a faceted option for keeping the location's lumped together but viewing the bat species separately: ####

#define your palette:
#I selected from the following list:
p2palette = c("#5cdcca", "#dc80a3","#ffad8b","#e09cd7","#538b51","#7b9de3","#dcca81","#2eb0e5",
              "#b67a4a","#72d7fa","#3eb396","#278791","#bb9755","#4a81a4","#add69a","#8a729d","#5e8461","#a76a71","#95b7ca","#e7c4a5")
               
#and plot
s <- ggplot(relative_OTUs_perBatAndSite, aes(x = order_name, y = relCounts, fill = order_name)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = p2palette) + 
  guides(col = guide_legend(nrow = 1)) + 
  #guides(col = guide_legend(override.aes = list(nrow = 1))) + 
  theme(legend.position = "right") + 
  ylab("Number of Observations") +
  xlab("Bat species") + 
  theme(axis.text.y = element_text(size = 8)) +
  theme(axis.text.x=element_blank()) +
  labs(title = "Number of Taxa detected (by Order) per Bat Species at each Location using UFITS/DADA2 pipeline") +
  theme(legend.key.size = unit(7, "mm")) +
  theme(legend.title = element_text(size=16, face="bold")) +
  scale_color_discrete(name="Order") +
  #facet_wrap(LocationName ~ BatSpecies, ncol = 5)
  facet_wrap(LocationName ~ BatSpecies)
  
s

(legend.title = element_text(colour="chocolate", size=16, face="bold"))+
  scale_color_discrete(name="This color is\ncalled chocolate!?")

###################################        END Plot 2        #####################################
##################################################################################################

#a dodged and faceted barchart option:
t <- ggplot(relative_OTUs_perBatAndSite, aes(x = order_name, y = relCounts, fill = order_name)) + 
  geom_bar(stat = "identity") +
  facet_grid(. ~ BatSpecies) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  scale_fill_manual(values = p2palette)
t

###### let's get crazy here... going to make a two-variable faceted figure:

u <- ggplot(relative_OTUs_perBatAndSite, aes(x = order_name, y = relCounts, fill = order_name)) + 
  geom_bar(stat = "identity", width = 0.4) +
  facet_grid(LocationName ~ BatSpecies) +
  scale_fill_manual(values = p2palette) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank()) +
  ylab("Fraction of unique taxonomic orders among all taxa in set") +
  xlab("Bat species") + 
  theme(axis.text.y = element_text(size = 8)) +
  #theme(axis.text.x = element_text(size = 10)) +
  labs(title = "Fraction of orders identified per bat species per site using UFITS/DADA2 pipeline") +
  theme(legend.key.size = unit(4, "mm"))
u


#tiny switch will make this using relative percents:

urel <- ggplot(relative_OTUs_perBatAndSite, aes(x = order_name, y = percent_taxa, fill = order_name)) + 
  geom_bar(stat = "identity", width = 0.4) +
  facet_grid(LocationName ~ BatSpecies) +
  scale_fill_manual(values = p2palette) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank()) +
  ylab("Fraction of unique taxonomic orders among all taxa in set") +
  xlab("Bat species") + 
  theme(axis.text.y = element_text(size = 8)) +
  #theme(axis.text.x = element_text(size = 10)) +
  labs(title = "Fraction of orders identified per bat species per site using UFITS/DADA2 pipeline") +
  theme(legend.key.size = unit(4, "mm"))
urel



elyplot <- ggplot(ely_diets_byorder) + aes(x = order_name, y = fraction, fill = order_name)
v <- elyplot +
  geom_bar(stat = "identity", width = 0.4) +
  facet_grid(. ~ BatSpecies) +
  scale_fill_manual(values = newPalette) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank()) +
  ylab("Percent order detected") +
  xlab("Bat species") + 
  theme(axis.text.y = element_text(size = 8)) +
  #theme(axis.text.x = element_text(size = 10)) +
  labs(title = "Number of orders identified per bat species at Ely Mine, VT, using UFITS/DADA2 pipeline") +
  theme(legend.key.size = unit(4, "mm"))
v


acaplot <- ggplot(aca_diets_byorder) + aes(x = order_name, y = fraction, fill = order_name)
w <- acaplot +
  geom_bar(stat = "identity", width = 0.4) +
  facet_grid(. ~ BatSpecies) +
  scale_fill_manual(values = newPalette) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank()) +
  ylab("Percent order detected") +
  xlab("Bat species") + 
  theme(axis.text.y = element_text(size = 8)) +
  #theme(axis.text.x = element_text(size = 10)) +
  labs(title = "Number of orders identified per bat species at Acadia NP, ME, using UFITS/DADA2 pipeline") +
  theme(legend.key.size = unit(4, "mm"))
w


######### goign the manual way in Excel
#dumped the "bat_sp_freq" table to desktop to edit in Excel by manually summing up the unique LocationName and BatSpecies combinations. 
#I think we can use 'data.table' for this but I didn't have the time'
write.table(bat_sp_freq, "~/Desktop/fixthis.txt", quote = F, sep = "\t", row.names = F, col.names = T)


palette14 = c("#a6cee3", 
              "#1f78b4", 
              "#b2df8a",
              "#33a02c",
              "#fb9a99",
              "#e31a1c", 
              "#fdbf6f",
              "#ff7f00",
              "#cab2d6",
              "#6a3d9a",
              "#ffff99",
              "#b15928",
              "#87d2b4",
              "#e0b792")
  