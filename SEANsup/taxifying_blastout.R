#formatting NCBI BLAST output (GI and Accession numbers) to include taxonmy information with 'taxize'
#written 1.6.16
#modified 4.17.17
#Devon O'Rourke
#from Sean's samples 

########################################################################################
############################    part 1 - data wrangling  ###############################
########################################################################################

#import sorted blast output .txt file as new df
setwd("~/Documents/Lab.Foster/guano/guano_sean/amptk_output/")
blastdf1 <- read.table("cleanBlastout.txt", quote="\"", comment.char="")

#add colnames
colnames(blastdf1) <- c("OTUid", "Accession", "pid", "length", "bit", "eval", "TAXid")
library(stringr)
newcol <- as.data.frame(str_split_fixed(blastdf1$TAXid, ";", 2))
  # adjust the interger value to however many the longest value of delimited values are within the 'TAXid' field
blastdf2 <- cbind(blastdf1, newcol)
rm(blastdf1, newcol)

blastdf2$TAXid <- NULL
colnames(blastdf2) <- c("OTUid", "Accession", "pid", "length", "bit", "eval", "TAXid1", "TAXid2")
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

#remove instances where the OTUid and TAXid are identical
blastdf <- BOdf.clean[!duplicated(BOdf.clean[c("OTUid", "TAXid")]),]
rm(BOdf.clean)

#reorder by OTUid, then BitScore ("bit"), then PercentIdentity ("pid")
blastdf <- blastdf[with(blastdf, order(OTUid, -bit, -pid)), ] 


########################################################################################
#########################    part 2 - taxonomy classifying  ############################
########################################################################################

#load in necessary package
library(taxize)

#make a list of TAXid's
taxids_toquery <- as.character(unique(blastdf$TAXid))

#assign classifications by searching NCBI database for taxonomy descriptions given taxid values
  #(this takes a few minutes)
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

#reorder by OTUid, then BitScore ("bit"), then PercentIdentity ("pid")
blastdf.withTaxa <- blastdf.withTaxa[with(blastdf.withTaxa, order(OTUid, -bit, -pid)), ] 

#write file to disk:
setwd("~/Documents/Lab.Foster/guano/guano_sean/amptk_output/")
write.csv(blastdf.withTaxa, file = "OTUidTaxaClassifications.csv", row.names = F, quote = F)