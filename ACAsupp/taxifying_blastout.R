#formatting NCBI BLAST output (GI and Accession numbers) to include taxonmy information with 'taxize'
#written 1.6.16
#modified 4.3.17
#Devon O'Rourke
#from BRI wingpunch samples submitted by Dave Yates

########################################################################################
############################    part 1 - data wrangling  ###############################
########################################################################################

#import sorted blast output .txt file as new df
setwd("~/Documents/Lab.Foster/bri-wingpunch/")
blastdf1 <- read.table("BRIwp.blastCleaned.txt", quote="\"", comment.char="")

#add colnames
colnames(blastdf1) <- c("SampleID", "Accession", "pid", "length", "bit", "eval", "TAXid")
library(stringr)
newcol <- as.data.frame(str_split_fixed(blastdf1$TAXid, ";", 2))
  # adjust the interger value to however many the longest value of delimited values are within the 'TAXid' field
blastdf2 <- cbind(blastdf1, newcol)
rm(blastdf1, newcol)

blastdf2$TAXid <- NULL
colnames(blastdf2) <- c("SampleID", "Accession", "pid", "length", "bit", "eval", "TAXid1", "TAXid2")
blastdf2$TAXid2[blastdf2$TAXid2 == ""] <- NA

#dump all fields after initial "TAXid" column
mrg.df1 <- blastdf2[,c(1:7)]
colnames(mrg.df1)[7] <- "TAXid"

#create separate data frames with each additional "TAXid" column (other than "TAXid1")
mrg.df2 <- blastdf2[,c(1:6,8)]
#remove entries where values in "TAXid1" are present in "TAXid2" (ie. keep just the unique entries for the "TAXid2" column)
mrg.df2 <- mrg.df2[!is.na(mrg.df2$TAXid2),]
colnames(mrg.df2)[7] <- "TAXid"

#now merge those two data frames vertically:
BOdf.clean <- rbind(mrg.df1, mrg.df2)
rm(mrg.df1, mrg.df2, blastdf2)

#remove instances where the SampleID and TAXid are identical
blastdf <- BOdf.clean[!duplicated(BOdf.clean[c("SampleID", "TAXid")]),]
rm(BOdf.clean)

#reorder by SampleID, then BitScore ("bit"), then PercentIdentity ("pid")
blastdf <- blastdf[with(blastdf, order(SampleID, -bit, -pid)), ] 


########################################################################################
#########################    part 2 - taxonomy classifying  ############################
########################################################################################

#load in necessary package
library(taxize)

#make a list of TAXid's
taxids_toquery <- as.character(unique(blastdf$TAXid))

#make sure list isn't particularly long (keep less than about 200 characters)
#See notes in section "Part X" below for example if you have more than ~200 queries and need to split up your data frame...
taxinfo <- (classification(taxids_toquery, db = 'ncbi'))

#bind info together:
taxa_chunk <- cbind(taxinfo)

#specify columns you want:
taxa_chunk <- taxa_chunk[,c("kingdom", "phylum", "class", "order", "family", "genus", "species", "query")]

#and rename them for an incoming merging with the 'blastdf' object
colnames(taxa_chunk) <- c("ncbi-kingdom", 
                            "ncbi-phylum", 
                            "ncbi-class", 
                            "ncbi-order", 
                            "ncbi-family", 
                            "ncbi-genus", 
                            "ncbi-species",
                            "TAXid")

#and clean up your mess
rm(taxids_toquery, taxinfo)

########################################################################################
##################    part 3 - merge NCBI data with BLAST results  #####################
########################################################################################

blastdf.withTaxa <- merge(blastdf, taxa_chunk, by = "TAXid")
rm(blastdf, taxa_chunk)

#drop unnecessary taxonomic info:
blastdf.withTaxa$`ncbi-kingdom` <- NULL
blastdf.withTaxa$`ncbi-phylum` <- NULL
blastdf.withTaxa$`ncbi-class` <- NULL
blastdf.withTaxa$`ncbi-order` <- NULL
blastdf.withTaxa$`ncbi-family` <- NULL


#reorder by SampleID, then BitScore ("bit"), then PercentIdentity ("pid")
blastdf.withTaxa <- blastdf.withTaxa[with(blastdf.withTaxa, order(SampleID, -bit, -pid)), ] 

#split the SampleID column into "WellNumber" and "BATid" fields
newcols <- as.data.frame(str_split_fixed(blastdf.withTaxa$SampleID, "_", 2))
colnames(newcols) <- c("WellNumber", "BATid")

#paste back into 'blastdf.withTaxa'
blastdf.withTaxa <- cbind(blastdf.withTaxa, newcols)
blastdf.withTaxa$SampleID <- NULL
rm(newcols)

colnames(blastdf.withTaxa) <- c("taxID", "GenbankIDs", "PercID", "Length", "BitScore", "Evalue", "Genus", "Species", "WellNumber", "BatID")

#let's make an assumption that the results can't include a few species. Specifically the two western Myotis species:
##"Myotis californicus" (California myotis), or
##"Myotis ciliolabrum" (Western small footed)

blastdf.lessWest <- blastdf.withTaxa[!grepl("257884|257882", blastdf.withTaxa$taxID),]
  #this produces a single-species call for every sample. nice!

#drop unnecessary fields
blastdf.lessWest$Genus <- NULL
blastdf.lessWest$Length <- NULL
blastdf.lessWest$Evalue <- NULL

#edit the BatID values to remove the primer name attached by the Sanger sequencing center output
blastdf.lessWest$BatID <- sub("_\\S*", "", blastdf.lessWest$BatID)

#reorder by SampleID, then BitScore ("bit"), then PercentIdentity ("pid")
blastdf.lessWest <- blastdf.lessWest[with(blastdf.lessWest, order(BatID, Species)), ] 

#write file to disk:
write.table(blastdf.lessWest, file = "PerSampleTaxaClassifications.txt", row.names = F, quote = F, col.names = T)
write.csv(blastdf.lessWest, file = "PerSampleTaxaClassifications.csv", row.names = F, quote = F)

########################################################################################
###########    part X - merge NCBI data with BLAST results with many samples ###########
########################################################################################

#make a list of TAXid's
taxids_toquery <- as.character(unique(blastdf$TAXid))
#split into smaller lists of ~200
query_chunk1 <- taxids_toquery[1:200]
query_chunk2 <- taxids_toquery[201:400]
query_chunk3 <- taxids_toquery[401:583]

#get classification IDs from ncbi nr database
taxinfo_chunk1 <- (classification(query_chunk1, db = 'ncbi'))
taxinfo_chunk2 <- (classification(query_chunk2, db = 'ncbi'))
taxinfo_chunk3 <- (classification(query_chunk3, db = 'ncbi'))

#bind it together
taxa_chunk1_raw <- cbind(taxinfo_chunk1)
taxa_chunk2_raw <- cbind(taxinfo_chunk2)
taxa_chunk3_raw <- cbind(taxinfo_chunk3)

#get the columns you need:
taxa_chunk1 <- taxa_chunk1_raw[,c("kingdom", "phylum", "class", "order", "family", "genus", "species", "query")]
taxa_chunk2 <- taxa_chunk2_raw[,c("kingdom", "phylum", "class", "order", "family", "genus", "species", "query")]
taxa_chunk3 <- taxa_chunk3_raw[,c("kingdom", "phylum", "class", "order", "family", "genus", "species", "query")]

#then bind them all together
ncbi_taxa_df = rbind(taxa_chunk1, taxa_chunk2, taxa_chunk3)
#and rename them for an incoming merging with the 'blastdf' object
colnames(ncbi_taxa_df) <- c("ncbi-kingdom", 
                            "ncbi-phylum", 
                            "ncbi-class", 
                            "ncbi-order", 
                            "ncbi-family", 
                            "ncbi-genus", 
                            "ncbi-species",
                            "TAXid")

#and clean up your mess
rm(taxa_chunk1, taxa_chunk2, taxa_chunk3, taxa_chunk1_raw, taxa_chunk2_raw, taxa_chunk3_raw, blastdf, ncbi_taxa_df)
rm(query_chunk1, query_chunk2, query_chunk3, taxinfo_chunk1, taxinfo_chunk2, taxinfo_chunk3, taxids_toquery)