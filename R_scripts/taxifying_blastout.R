#updating SINTAX and UTAX taxonmy information with 'taxize'
#written 1.6.16
#Devon O'Rourke
#for ufits analysis from cleaned blast output file to determine taxonomy information

########################################################################################
############################    part 1 - data wrangling  ###############################
########################################################################################

#import sorted blast output .txt file as new df
setwd("~/Documents/Lab.Foster/guano/guano_bri/bri-chunk1/")
blastdf1 <- read.table("cln3blastout.txt", quote="\"", comment.char="")

#add colnames
colnames(blastdf1) <- c("OTUid", "Accession", "pid", "length", "bit", "eval", "TAXid")
library(stringr)
newcol <- as.data.frame(str_split_fixed(blastdf1$TAXid, ";", 4))
blastdf2 <- cbind(blastdf1, newcol)
blastdf2$TAXid <- NULL
colnames(blastdf2) <- c("OTUid", "Accession", "pid", "length", "bit", "eval", "TAXid1", "TAXid2", "TAXid3", "drop")
blastdf2$drop <- NULL
blastdf2$TAXid2[blastdf2$TAXid2 == ""] <- NA
blastdf2$TAXid3[blastdf2$TAXid3 == ""] <- NA

#now fix those weird cases where multiple unique TAXids have the same variable values
library(reshape)
blastdf2 <- melt(blastdf2, id=c("OTUid", "Accession", "pid", "length", "bit", "eval"))
blastdf2$variable <- NULL
colnames(blastdf2)[7] <- "TAXid"
blastdf <- na.omit(blastdf2)

#cleanup
rm(blastdf1, blastdf2, newcol)


########################################################################################
#########################    part 2 - taxonomy classifying  ############################
########################################################################################

#load in necessary package
library(taxize)

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

########################################################################################
##################    part 3 - merge NCBI data with BLAST results  #####################
########################################################################################

blastTaxa_df <- merge(blastdf, ncbi_taxa_df, by = "TAXid")

#I'm not convinced it's worthwhile to compare BLAST outputs with anything less than 95% alignment... 
#what we're really trying to do is replace the obvious things - the OTUs which were assigned as "Animalia" only...
##...but really could be assigned confidently to a genus/species because of really high BLAST alignemnt...
### we could have avoided this entirely if we switched our 'blastn' parameters to query for a '-perc_identity' that was...
### much higher than what we did for this data set (it was set to 79.9); we did this because...
### I wanted to see how many OTUs would have reads higher or lower than a certain bit score / pid / eval...

#subset for high pid
min95blastTaxa_df <- subset(blastTaxa_df, pid >= 95)

#get rid of the unnecessary "ncbi-kingdom" info
min95blastTaxa_df$`ncbi-kingdom` <- NULL

#sort this df by the 'OTUid' first and the 'bit' score second
min95blastTaxa_df <- with(min95blastTaxa_df, min95blastTaxa_df[order(OTUid, bit),])


########################################################################################
####    part 4 - generate data frame for comparison to original UTAX/SINTAX calls  #####
########################################################################################

#import the table of OTUs originally assigned taxonomic information from SINTAX and UTAX
clntaxheads <- read.delim("~/Documents/Lab.Foster/guano/guano_bri/bri-chunk1/clntaxheads.txt", header=FALSE, stringsAsFactors=FALSE)

#rename headers
colnames(clntaxheads) <- c("OTUid", "classifier", "kingdom", "phylum", "class", "order", "family", "genus", "species")

#remove unwanted characters
clntaxheads$OTUid <- substring(clntaxheads$OTUid, 2)
clntaxheads$kingdom <- substring(clntaxheads$kingdom, 3)
clntaxheads$phylum <- substring(clntaxheads$phylum, 3)
clntaxheads$class <- substring(clntaxheads$class, 3)
clntaxheads$order <- substring(clntaxheads$order, 3)
clntaxheads$family <- substring(clntaxheads$family, 3)
clntaxheads$genus <- substring(clntaxheads$genus, 3)
clntaxheads$species <- substring(clntaxheads$species, 3)

#insert "NA" for empty spaces
clntaxheads[clntaxheads == ""] <- NA

#make list of OTUs where SINTAX or UTAX has fully assigned taxonomy
##we're assuming SINTAX gets it right more than NCBI's nr database in case of conflict
suTAX_otus_df <- na.omit(clntaxheads)
suTAXotu_list <- suTAX_otus_df$OTUid

#make a df of just the incompletely assigned values
suTAX_incompleteOTUs_df <- subset(clntaxheads, is.na(clntaxheads$species))

########################################################################################
#######    part 5 - compare filtered BLAST results with SINTAX/UTAX results    #########
########################################################################################


#remove the OTUs which aren't going to be compared from the blastTaxa_df (the fully classified ones form SINTAX/UTAX)
blast_df_Compare2suTAX <- min95blastTaxa_df[!min95blastTaxa_df$OTUid %in% suTAXotu_list, ]

write.table(blast_df_Compare2suTAX, file = "blastout_toparse.txt", sep = "\t", quote = F, col.names = T, row.names = F)
