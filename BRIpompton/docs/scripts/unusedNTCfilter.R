setwd("~/Repos/guano/BRIpompton/data/amptk/")
## load in *.sorted.csv from `amptk filter` output:
df <- read.csv(file = "mockIn.sorted.csv")

## pull out only negative samples (my samples always have the prefix "NTC"):
ntc.df <- df[grepl("NTC", names(df))]
row.names(ntc.df) <- df$OTUid   ## remember what the OTUids are; use rowname so you can do math on data.frame

## removes rows where no elements > 9, as we'll set `amptk --subtract` filter to == 10
ntc1.df <- ntc.df[apply(ntc.df, MARGIN = 1, function(x) any(x > 9)),]
is.na(ntc1.df) <- ntc1.df == 0
## get a data.frame of the average of the rows, sorted ascending order:
OTUavg <- data.frame(sort(rowMeans(ntc1.df, na.rm = TRUE)))
colnames(OTUavg) <- "counts"
