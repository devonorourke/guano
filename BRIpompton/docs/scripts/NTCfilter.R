setwd("~/Repos/guano/BRIpompton/data/amptk/")
## load in noFilt.final.csv from `amptk filter` output:
df <- read.csv(file = "noFilt.final.csv")

## pull out only negative samples (my samples always have the prefix "NTC"):
ntc.df <- df[grepl("NTC", names(df))]
row.names(ntc.df) <- df$OTUid   ## remember what the OTUids are; use rowname so you can do math on data.frame

## find the highest value among NTC elements per row, then subtract that value from original df:
## first, find the maximum value in any NTC sample, per OTU:
ntcRowmax <- data.frame(apply(ntc.df, 1, max))
colnames(ntcRowmax) <- "counts"

## drop the OTUid column and keep as rowname
df1 <- df[,-1]  
rownames(df1) <- df$OTUid

## subtract the max value of reads in any NTC sample observed among for each OTU across all samples in the original matrix
filt.df <- sweep(df1,1,ntcRowmax$counts,"-")  ## for 'sweep' details see: https://bioinfomagician.wordpress.com/2014/08/12/my-favorite-commands-part3-sweep-function-in-r/
rm(df1)

## and reduce values less than 0 to 0:
filt.df[filt.df<0] <- 0

##notrun: now let's repeat the same process but focus on the remaining values in the `mock` vector:
#notrun: filt.df <- sweep(filt.df,1,filt.df$mock.IM4p11.2,"-")
#notrun: filt.df[filt.df<0] <- 0   # set neg values back to zero again...

## now get format back into what Jon expects:
filt.df$OTUid <- df$X.OTU.ID  # adding back in the OTU column
filt.df <- filt.df[,c(104,1:103)] # positioning OTUid column back to first position
  # note we'll need to rename the OTUid column name in bash once imported back
  ## this was done after writing the file: 
  ## ```
  ## sed -i 's/OTUid/#OTU ID/' NTCreduced.otu_table.csv
  ## ````

rm(df, ntc.df, ntcRowmax)
## wand write to disk:
setwd("~/Repos/guano/BRIpompton/data/amptk/")
write.table(filt.df, file = "NTCreduced.otu_table.txt", quote = F, row.names = FALSE, sep = '\t')
  ## see above note about `sed` command to rename the OTUid header in the first field
          