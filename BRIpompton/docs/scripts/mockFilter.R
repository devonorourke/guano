## subset new dataframe containing all OTUid's not assigned to a mock community:
nonmock.df <- filt.df[! grepl("pident", filt.df$OTUid), ]

## subset that data.frame to include only values where the mock community have > 0 reads:
inmock <- subset(nonmock.df, mock.IM4p11.2 > 0)
rownames(inmock) <- inmock$OTUid
inmock$OTUid <- NULL


## find the max value among any sample for that OTU
inmocksummary <- data.frame(apply(inmock, 1, max))
colnames(inmocksummary) <- "rowMax"

## find the sum of each of these rows:
inmocksummary$rowSums <- (apply(inmock, MARGIN = 1, FUN = function(x) sum(x)))

## find number of samples with non-zero number of reads per row (per OTU)
inmocksummary$numHits <- rowSums(inmock != 0)

## get an average of the number of reads per OTU
inmocksummary$OTUavg <- (inmocksummary$rowSums) / (inmocksummary$numHits)

## and add in the mock column only into the summary data.frame
inmocksummary$mock.IM4p11.2 <- inmock$mock.IM4p11.2

## remove instance where there are <= 2 samples (these OTUs only present in one sample that will be dropped regardless)
inmocksummary <- subset(inmocksummary, numHits > 1)

