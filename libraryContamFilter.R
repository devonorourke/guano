## title: Filtering OTUtable for duplicates in different libraries
## author: devon o'rourke via stack exchange help from user 'caw5cv'
## written: 18-March-2018
## see: https://stackoverflow.com/questions/49351179/filtering-single-matrix-conditional-upon-column-name-regex-and-element-value/49351520#49351520

library(data.table)
ra <- c(1,2,0,0,0,0)
rb <- c(3,0,4,0,0,0)
rc <- c(0,0,5,6,7,0)
rd <- c(0,0,8,0,0,0)
re <- c(0,0,0,0,0,99)
mat <- rbind(ra,rb,rc,rd,re)
colnames(mat) <- c("s1-G1", "s2-G1", "s3-G2", "s4-G2", "s5-G3", "s6-G3")
rownames(mat) <- c("OTU1", "OTU2", "OTU3", "OTU4", "OTU5")

# make a list of row names to add in as new column in new data frame about to be created:
myrownames <- rownames(mat)

# Convert to data.table
dat.wide <- data.table(mat)
dat.wide[, OTU := myrownames]

# Convert to long format
dat.long <- melt(dat.wide)

# split sample and group name into separate columns, creating new column "value"
dat.long[, c("sample","group") := tstrsplit(variable, split="-")]

# remove original joint sample-group column
dat.long[, variable := NULL]

# get unique OTUs
unique.OTUs <- dat.long[, list(N=sum(value)), by=list(group, OTU)][, list(Ngroups=sum(N>0)), by=OTU][Ngroups==1]$OTU
df_uniq <- dat.wide[OTU %in% unique.OTUs]

# get duplicated OTUs
df_dupd <- dat.wide[! (OTU %in% unique.OTUs)]