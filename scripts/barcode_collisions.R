## calculating barcode collisions
## written: 05.11.2017
## author: devon o'rourke

## install package - do once (commented out)
    # install.packages("stringdist")

## load library
library(stringdist)

## load data, create workable data.frame
setwd("~/Desktop/guano/")
barcodes <- read.delim(file = "8bp_barcode_seqs_names.txt", header = T, stringsAsFactors = F)
colnames(barcodes) <- c("i5index_sequence", "i5index_name", "i7index_sequence", "i7index_name")
seqnamesi5 <- barcodes[,1:2]
colnames(seqnamesi5) <- c("index_sequence", "index_name")
seqnamesi7 <- barcodes[,3:4]
colnames(seqnamesi7) <- c("index_sequence", "index_name")
all_barcodes <- rbind(seqnamesi5, seqnamesi7)
rm(seqnamesi5, seqnamesi7)

## calculate i5 edit distances and convert to data.frame
i5dm <- stringdistmatrix(barcodes$i5index_sequence, method = "hamming")
i5dm_df <- as.data.frame(as.matrix(i5dm))
colnames(i5dm_df) <- barcodes$i5index_name
row.names(i5dm_df) <- barcodes$i5index_name

## exclude self:self calculations within distance.matrix
i5dm_df[i5dm_df == 0] <- NA

## exclude self:non-self hamming distances greater than 2 (so 3, 4, 5, 6, 7, or 8 are not analyzed)
i5dm_df[i5dm_df > 2] <- NA

## get a list of pairwise comparisons that passed filters
## note that the only value to identify is `==2` as all values greater were filtered, and no barcodes have a Hamming distance of 1
i5matches <- data.frame(as.matrix(which(i5dm_df==2,arr.ind=T)))   
## because you're looking through a pairwise matrix, output will have a reciprocal entry which is a duplicate collision...
## this is filtered in the last step
## I can't for the life of me figure out how to print out these values short of this stupid merge script below... 
## ...but there has to be a way to get a matrix index and pull down the row/column names.

## let's convert the matrix index values to actual row/column names:
i5hitlist <- data.frame(barcodes[,2])
i5row <- data.frame(rownames(i5hitlist))
i5rowlist <- cbind(i5hitlist, i5row)
i5collist <- i5rowlist
colnames(i5rowlist) <- c("barcode_name", "row")
colnames(i5collist) <- c("barcode_name", "col")

i5match1 <- merge(i5matches, i5rowlist)
i5collisions <- merge(i5match1, i5collist, by="col")
i5collisions$col <- NULL
i5collisions$row <- NULL
colnames(i5collisions) <- c("barcode1", "barcode2")

## and let's sort out the duplicate entries
i5.u.collisions <- data.frame(unique(t(apply(i5collisions, 1, sort))))
colnames(i5.u.collisions) <- c("barcode1", "barcode2")

## finally, if you want to see what the actual sequences are... another dumb merge script:
i5tmp1 <- merge(i5.u.collisions, all_barcodes, by.x = "barcode1", by.y="index_name")
colnames(i5tmp1) <- c("barcode1", "barcode2", "barcode1_sequence")
i5.u.collisions.wSeqs <- merge(i5tmp1, all_barcodes, by.x = "barcode2", by.y="index_name")
colnames(i5.u.collisions.wSeqs) <- c("barcode2_name", "barcode1_name", "barcode1_seq", "barcode2_seq")
u.i5.collisions.wSeqs <- i5.u.collisions.wSeqs[,c(2,1,3,4)]
write.table(u.i5.collisions.wSeqs, file = "i5collisions.txt", row.names = F, quote = F)


##############################

## calculate i7 edit distances and convert to data.frame
i7dm <- stringdistmatrix(barcodes$i7index_sequence, method = "hamming")
i7dm_df <- as.data.frame(as.matrix(i7dm))
colnames(i7dm_df) <- barcodes$i7index_name
row.names(i7dm_df) <- barcodes$i7index_name

## exclude self:self calculations within distance.matrix
i7dm_df[i7dm_df == 0] <- NA

## exclude self:non-self hamming distances greater than 2 (so 3, 4, 5, 6, 7, or 8 are not analyzed)
i7dm_df[i7dm_df > 2] <- NA

## get a list of pairwise comparisons that passed filters
## note that the only value to identify is `==2` as all values greater were filtered, and no barcodes have a Hamming distance of 1
i7matches <- data.frame(as.matrix(which(i7dm_df==2,arr.ind=T)))   
## because you're looking through a pairwise matrix, output will have a reciprocal entry which is a duplicate collision...
## this is filtered in the last step
## I can't for the life of me figure out how to print out these values short of this stupid merge script below... 
## ...but there has to be a way to get a matrix index and pull down the row/column names.

## let's convert the matrix index values to actual row/column names:
i7hitlist <- data.frame(barcodes[,4])
i7row <- data.frame(rownames(i7hitlist))
i7rowlist <- cbind(i7hitlist, i7row)
i7collist <- i7rowlist
colnames(i7rowlist) <- c("barcode_name", "row")
colnames(i7collist) <- c("barcode_name", "col")

i7match1 <- merge(i7matches, i7rowlist)
i7collisions <- merge(i7match1, i7collist, by="col")
i7collisions$col <- NULL
i7collisions$row <- NULL
colnames(i7collisions) <- c("barcode1", "barcode2")

## and let's sort out the duplicate entries
i7.u.collisions <- data.frame(unique(t(apply(i7collisions, 1, sort))))
colnames(i7.u.collisions) <- c("barcode1", "barcode2")

## finally, if you want to see what the actual sequences are... another dumb merge script:
i7tmp1 <- merge(i7.u.collisions, all_barcodes, by.x = "barcode1", by.y="index_name")
colnames(i7tmp1) <- c("barcode1", "barcode2", "barcode1_sequence")
i7.u.collisions.wSeqs <- merge(i7tmp1, all_barcodes, by.x = "barcode2", by.y="index_name")
colnames(i7.u.collisions.wSeqs) <- c("barcode2_name", "barcode1_name", "barcode1_seq", "barcode2_seq")
u.i7.collisions.wSeqs <- i7.u.collisions.wSeqs[,c(2,1,3,4)]
write.table(u.i7.collisions.wSeqs, file = "i7collisions.txt", row.names = F, quote = F)

########## ########## ########## ########## ########## ########## 

## or just stack both data.frames together:
write.table(rbind(u.i5.collisions.wSeqs, u.i7.collisions.wSeqs), file = "all.collisions.txt", row.names = F, quote = F)
