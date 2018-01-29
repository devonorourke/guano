## purpose: amptk OTU analysis
## written: 28-jan-2018
## author: devon orourke

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # 
               ######     Part 1 - data wrangling     ######     
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #

## load libraries:
library(data.table)
library(reshape2)
library(tidyr)

## read in data:
w1.df <- fread('https://raw.githubusercontent.com/devonorourke/guano/master/Rutgers/rut16_h.otu_table.taxonomy.txt')

## reformat matrix into data.frame, split strings into appropriate columns
w1.df[w1.df == 0] <- NA       # we loaded a binary matrix; all zero's indicate absence of data, all 1's indicate presence of an amplicon
tmp1.df <- melt(w1.df)           # converting from wide format to long format (useful for plots)
rm(w1.df)
colnames(tmp1.df) <- c("OTUid", "Taxonomy", "SampleID", "Presence")
tmp2.df <- tmp1.df[complete.cases(tmp1.df),]
rm(tmp1.df)
tmp2.df$Presence <- NULL
tmp2.df <- separate(data = tmp2.df, 
                    col = Taxonomy, 
                    into = c("TaxMethod", 
                             "AlignScore", 
                             "BOLDid", 
                             "kingdom_name", 
                             "phylum_name", 
                             "class_name", 
                             "order_name", 
                             "family_name", 
                             "genus_name", 
                             "species_name"), 
                    sep = "\\;|,|\\|")

## reformat taxonomy columns to discard given prefix
tmp2.df$kingdom_name <- sub("k:", "", tmp2.df$kingdom_name)
tmp2.df$phylum_name <- sub("p:", "", tmp2.df$phylum_name)
tmp2.df$class_name <- sub("c:", "", tmp2.df$class_name)
tmp2.df$order_name <- sub("o:", "", tmp2.df$order_name)
tmp2.df$family_name <- sub("f:", "", tmp2.df$family_name)
tmp2.df$genus_name <- sub("g:", "", tmp2.df$genus_name)
tmp2.df$species_name <- sub("s:", "", tmp2.df$species_name)

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ # 
               ######     Part 2 - frequency tables     ######     
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #