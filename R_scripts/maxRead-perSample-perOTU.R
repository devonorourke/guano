## script to identify the OTU from each sample has the highest proportion of reads in amptk-filter output
## written: 28-jan-2018
## author: devon orourke

## read in data:
library(readr)
reads.df <- read_csv("~/Desktop/guano/Rutgers/mockIn_noNorm.sorted.csv")

## create matrix from data.frame, creating row.names from first column of data.frame:
reads.mat <- as.matrix(reads.df[,-1])
namelist <- reads.df$otuID
rownames(reads.mat) <- namelist

## create output
maxRead.df <- data.frame(row.names = colnames(reads.mat),
                MaxVal = apply(reads.mat, 2, max),
                WhichMax = apply(reads.mat, 2, which.max))

## this creates a "WhichMax" value which is the row number, not row name. Working a bit more to fix that:
maxRead.df$SampleID <- rownames(maxRead.df)
counter = (1:953)
swap.df <- data.frame(counter, namelist)
colnames(swap.df) <- c("WhichMax", "OTUid")

Final.df <- merge(maxRead.df, swap.df)
Final.df$WhichMax <- NULL

write.table(Final.df, file = "~/Desktop/guano/Rutgers/tableS2.txt", row.names = F, col.names = T, sep = "\t", quote = F)
