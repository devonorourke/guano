## plotting various OTU tables with NMDS:
## generated from `collapsingSWYdataBySharedTaxa.R`

library(vegan)
library(MASS)
library(dplyr)
library(tidyr)
library(ggplot2)


## read in data:
setwd("~/Desktop/guano/oroAnalyses/qiime/data/diversity/dmats/")
mat.O <- read.csv("mat.c.sharedO.csv")
row.names(mat.O) <- mat.O$X
mat.O$X <- NULL

mat.F <- read.csv("mat.c.sharedF.csv")
row.names(mat.F) <- mat.F$X
mat.F$X <- NULL

mat.G <- read.csv("mat.c.sharedG.csv")
row.names(mat.G) <- mat.G$X
mat.G$X <- NULL

## calculate distances in community composition between sampling units using 4 distinct metrics...
## "brayB" with `binary=TRUE` uses presence/absence data, and is equivalent to Sorensen/Whittaker index
## "brayR" with `binary=FALSE` uses relative abundance; using counts of detections for common taxa sampled
## "raup1" will use presence/absence data, but assumes sampling probabilities proportional to species frequencies
## "raup0" will use presence/absence data, but uses equal sampling proabilities among all species...
## ..this is equivalent to `vegdist(method="raup")` in vegan

NMDS.F.raup0 <- raupcrick(mat.F, null = "r0")
NMDS.F.raup1 <- raupcrick(mat.F, null = "r1")
NMDS.F.brayB <- vegdist(mat.F, method = "bray", binary = TRUE)
NMDS.F.brayR <- vegdist(mat.F, method = "bray", binary = FALSE)


## stress plots
# define a function that performs an NMDS for `k` dimensions and plots the relationship of dimensions vs the stress
# adapted from `goeveg` package (https://rdrr.io/cran/goeveg/), modifying their function to allow for other non-vegdist functions to be added..
# ..in metaMDS call; rather than using their default plot, we just create a vector of the N values of k tested, then use ggplot on the output
# also turned the default metaMDS `autotransformation` to `FALSE`, and increased the `trymax` to 40, and the `k` to 8, and added in `try` variable

# there are 2 separate functions: 
#     - NMDSscreeR: works for distance matricies created with `raupcrick` function specifically
#     - NMDSscreeV: works for any distance matrix created using `vegdist` function (and associated methods)

## define both functions:
## NMDSscreeR:
NMDS.screeR <- function (dmat, distfun = raupcrick, k = 8, trymax = 30, try = 10, autotransform = FALSE) 
{
  stress <- 0
  for (i in 1:k) {
    nmds_i <- metaMDS(dmat, distfun = distfun, k = i, try = try, trymax = trymax, engine = "monoMDS", autotransform = autotransform)
    stress[i] <- nmds_i$stress
  }
  print(stress)
}
## NMDSscreeV:
NMDS.screeV <- function (dmat, distance = "bray", k = 8, trymax = 30, try = 10, autotransform = FALSE) 
{
  stress <- 0
  for (i in 1:k) {
    nmds_i <- metaMDS(dmat, distance = distance, k = i, try = try, trymax = trymax, engine = "monoMDS", autotransform = autotransform)
    stress[i] <- nmds_i$stress
  }
  print(stress)
}


## We can now collect the results for each of the four distance methods
## we first create a data.frame from the output of the 8 dimensions of stress data..
## ..then modify the column name and create a 2-column object with $Stress and $Dimensions for our plots

F.brayR.scree <- as.data.frame(NMDS.screeV(NMDS.F.brayR))
colnames(F.brayR.scree)[1] <- "Stress"
F.brayR.scree$Dimensions <- seq(1:(nrow(F.brayR.scree)))
F.brayR.scree %>% as_tibble()

F.brayB.scree <- as.data.frame(NMDS.screeV(NMDS.F.brayB))
colnames(F.brayB.scree)[1] <- "Stress"
F.brayB.scree$Dimensions <- seq(1:(nrow(F.brayB.scree)))

F.raup0.scree <- as.data.frame(NMDS.screeV(NMDS.F.raup0, distance = "raup"))
colnames(F.raup0.scree)[1] <- "Stress"
F.raup0.scree$Dimensions <- seq(1:(nrow(F.raup0.scree)))

F.raup1.scree <- as.data.frame(NMDS.screeR(NMDS.F.raup1))
colnames(F.raup1.scree)[1] <- "Stress"
F.raup1.scree$Dimensions <- seq(1:(nrow(F.raup1.scree)))


## run the plot:

 
## display results of any of the distance metrics

