## plotting various OTU tables with NMDS:
## generated from `collapsingSWYdataBySharedTaxa.R`

library(vegan)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

## read in data:
setwd("~/Desktop/guano/oroAnalyses/qiime/data/diversity/dmats/")

mat.F <- read.csv("mat.c.sharedF.csv")
row.names(mat.F) <- mat.F$X
mat.F$X <- NULL

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
NMDS.screeR <- function (dmat, distfun = raupcrick, k = 8, trymax = 30, try = 20, autotransform = FALSE) 
{
  stress <- 0
  for (i in 1:k) {
    nmds_i <- metaMDS(dmat, distfun = distfun, k = i, try = try, trymax = trymax, engine = "monoMDS", autotransform = autotransform)
    stress[i] <- nmds_i$stress
  }
  print(stress)
}
## NMDSscreeV:
NMDS.screeV <- function (dmat, distance = "bray", k = 8, trymax = 30, try = 20, autotransform = FALSE) 
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
## ..then modify the column name and create a 3-column object with $Dimensions (x-axis), $Stress (y-axis) and $Matrix (distance matrix input)

F.brayR.scree <- as.data.frame(NMDS.screeV(NMDS.F.brayR))
colnames(F.brayR.scree)[1] <- "Stress"
F.brayR.scree$Dimensions <- seq(1:(nrow(F.brayR.scree)))
F.brayR.scree$Matrix <- "F.brayR.scree"
F.brayR.scree %>% as_tibble()

F.brayB.scree <- as.data.frame(NMDS.screeV(NMDS.F.brayB))
colnames(F.brayB.scree)[1] <- "Stress"
F.brayB.scree$Dimensions <- seq(1:(nrow(F.brayB.scree)))
F.brayB.scree$Matrix <- "F.brayB.scree"

F.raup0.scree <- as.data.frame(NMDS.screeV(NMDS.F.raup0, distance = "raup"))
colnames(F.raup0.scree)[1] <- "Stress"
F.raup0.scree$Dimensions <- seq(1:(nrow(F.raup0.scree)))
F.raup0.scree$Matrix <- "F.raup0.scree"

F.raup1.scree <- as.data.frame(NMDS.screeR(NMDS.F.raup1))
colnames(F.raup1.scree)[1] <- "Stress"
F.raup1.scree$Dimensions <- seq(1:(nrow(F.raup1.scree)))
F.raup1.scree$Matrix <- "F.raup1.scree"


## combine the four individual data frames:
F.scree <- rbind(F.brayB.scree, F.brayR.scree, F.raup0.scree, F.raup1.scree)
## and plot...
ggplot(F.scree, aes(Dimensions, Stress, colour=Matrix)) + 
  geom_point() + 
  scale_colour_manual(values = c("red", "orange", "blue", "purple"),
                      labels = c("Sorensen", "Bray-Curtis", "Raup-Crick {r=0}", "Raup-Crick {r=1}")) +
  #scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +    ## useful for specifying number of decimals presented on y-axis
  geom_line() +
  labs(title = "Ordination stress reduction with K-dimensional space",
       subtitle = "NMDS ordination with shared Family-level taxa \n Distances calculated with Bray or Raup-Crick methods",
       colour = "Distance method used") +
  theme_light()
 

# ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ #

# We can go through the whole process again, but substitute in a different table (say Genus-collapsed data):

# read data
mat.G <- read.csv("mat.c.sharedG.csv")
row.names(mat.G) <- mat.G$X
mat.G$X <- NULL

# create distance matricies
NMDS.G.raup0 <- raupcrick(mat.G, null = "r0")
NMDS.G.raup1 <- raupcrick(mat.G, null = "r1")
NMDS.G.brayB <- vegdist(mat.G, method = "bray", binary = TRUE)
NMDS.G.brayR <- vegdist(mat.G, method = "bray", binary = FALSE)

# create the stress data:
G.brayR.scree <- as.data.frame(NMDS.screeV(NMDS.G.brayR))
colnames(G.brayR.scree)[1] <- "Stress"
G.brayR.scree$Dimensions <- seq(1:(nrow(G.brayR.scree)))
G.brayR.scree$Matrix <- "G.brayR.scree"
G.brayR.scree %>% as_tibble()

G.brayB.scree <- as.data.frame(NMDS.screeV(NMDS.G.brayB))
colnames(G.brayB.scree)[1] <- "Stress"
G.brayB.scree$Dimensions <- seq(1:(nrow(G.brayB.scree)))
G.brayB.scree$Matrix <- "G.brayB.scree"

G.raup0.scree <- as.data.frame(NMDS.screeV(NMDS.G.raup0, distance = "raup"))
colnames(G.raup0.scree)[1] <- "Stress"
G.raup0.scree$Dimensions <- seq(1:(nrow(G.raup0.scree)))
G.raup0.scree$Matrix <- "G.raup0.scree"

G.raup1.scree <- as.data.frame(NMDS.screeR(NMDS.G.raup1))
colnames(G.raup1.scree)[1] <- "Stress"
G.raup1.scree$Dimensions <- seq(1:(nrow(G.raup1.scree)))
G.raup1.scree$Matrix <- "G.raup1.scree"

# combine four datasets and plot:
G.scree <- rbind(G.brayB.scree, G.brayR.scree, G.raup0.scree, G.raup1.scree)
## and plot...
ggplot(G.scree, aes(Dimensions, Stress, colour=Matrix)) + 
  geom_point() + 
  scale_colour_manual(values = c("red", "orange", "blue", "purple"),
                      labels = c("Sorensen", "Bray-Curtis", "Raup-Crick {r=0}", "Raup-Crick {r=1}")) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.0001)) +
  geom_line() +
  labs(title = "Ordination stress reduction with K-dimensional space",
       subtitle = "NMDS ordination with shared Genus-level taxa \n Distances calculated with Bray or Raup-Crick methods",
       colour = "Distance method used") +
  theme_light()

# ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ # ~ #

# and again for Order-collapsed data:
# read data
mat.O <- read.csv("mat.c.sharedO.csv")
row.names(mat.O) <- mat.O$X
mat.O$X <- NULL

# create distance matricies
NMDS.O.raup0 <- raupcrick(mat.O, null = "r0")
NMDS.O.raup1 <- raupcrick(mat.O, null = "r1")
NMDS.O.brayB <- vegdist(mat.O, method = "bray", binary = TRUE)
NMDS.O.brayR <- vegdist(mat.O, method = "bray", binary = FALSE)

# create the stress data:
O.brayR.scree <- as.data.frame(NMDS.screeV(NMDS.O.brayR))
colnames(O.brayR.scree)[1] <- "Stress"
O.brayR.scree$Dimensions <- seq(1:(nrow(O.brayR.scree)))
O.brayR.scree$Matrix <- "O.brayR.scree"
O.brayR.scree %>% as_tibble()

O.brayB.scree <- as.data.frame(NMDS.screeV(NMDS.O.brayB))
colnames(O.brayB.scree)[1] <- "Stress"
O.brayB.scree$Dimensions <- seq(1:(nrow(O.brayB.scree)))
O.brayB.scree$Matrix <- "O.brayB.scree"

O.raup0.scree <- as.data.frame(NMDS.screeV(NMDS.O.raup0, distance = "raup"))
colnames(O.raup0.scree)[1] <- "Stress"
O.raup0.scree$Dimensions <- seq(1:(nrow(O.raup0.scree)))
O.raup0.scree$Matrix <- "O.raup0.scree"

O.raup1.scree <- as.data.frame(NMDS.screeR(NMDS.O.raup1))
colnames(O.raup1.scree)[1] <- "Stress"
O.raup1.scree$Dimensions <- seq(1:(nrow(O.raup1.scree)))
O.raup1.scree$Matrix <- "O.raup1.scree"

# combine four datasets and plot:
O.scree <- rbind(O.brayB.scree, O.brayR.scree, O.raup0.scree, O.raup1.scree)
## and plot...
ggplot(O.scree, aes(Dimensions, Stress, colour=Matrix)) + 
  geom_point() + 
  scale_colour_manual(values = c("red", "orange", "blue", "purple"),
                      labels = c("Sorensen", "Bray-Curtis", "Raup-Crick {r=0}", "Raup-Crick {r=1}")) +
  #scale_y_continuous(labels = scales::number_format(accuracy = 0.0001)) +
  geom_line() +
  labs(title = "Ordination stress reduction with K-dimensional space",
       subtitle = "NMDS ordination with shared Order-level taxa \n Distances calculated with Bray or Raup-Crick methods",
       colour = "Distance method used") +
  theme_light()
