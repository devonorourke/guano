## purpose: Sankey plot for visualizing taxonomic proportions across all levels
## written: 28-jan-2018
## author: devon orourke
## resource: see http://corybrunson.github.io/ggalluvial/index.html and https://cran.r-project.org/web/packages/ggalluvial/vignettes/ggalluvial.html

## load libraries
library(ggalluvial)
library(dplyr)
library(data.table)

## load data
  ## not run: master.df <- fread('https://raw.githubusercontent.com/devonorourke/guano/master/Rutgers/master.csv')
tmp.df <- master.df
tmp.df[is.na(tmp.df)] <- "undefined"
tmp.df$Location <- sub("^$", "-control", tmp.df$Location)
tmp.df$CollectDate <- sub("^$", "-control", tmp.df$CollectDate)

tmp.table <- tmp.df %>%
  count(Location, phylum_name, class_name, order_name, family_name, genus_name, species_name)

## coloring schemes can be altered depending on what part of plot desired to be showed.
## see http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf for one color list to choose from

## full shenanigans, coloring through taxonomic order , grouping color hues by taxonomic class
## (I think this is just too busy)...
all25colored <- c("#143b44","#76553f","#f8c7bb","#e03b2d","#442c61",
                  "#843e27","#b5d632","#fbaa72","#cb6e38","#b46569",
                  "#c1723d","#5c2ba4","#da9395","#a25325","#ffaa9a",
                  "#68a6c9","#ab4638","#bf9584","#e2835b","#97d54c",
                  "#e07e83","#3f646e","#68823f","#cf8957","#1c5872")

ggplot(tmp.table,
       aes(weight = n, axis1 = phylum_name, axis2 = class_name, axis3 = order_name, 
           axis4 = family_name, axis5 = genus_name, axis6 = species_name)) +
  geom_alluvium(aes(fill = order_name), width = 1/12) +
  scale_fill_manual(values = all25stark) +
  geom_stratum(width = 1/8, color = "darkgray") +   ## width = node/box chunk size, fill = box fill, color = box border
  geom_label(stat = "stratum", label.strata = TRUE, size = 1, fill = "honeydew") +             ## inserts labels for all values in a node; too many to use
  scale_x_continuous(breaks = 1:6, labels = c("phylum", "class", "order", "family", "genus", "species")) +
  ggtitle("What's for dinner or breakfast...") +
  labs(y = "number of detections")


## same plot structure, but colors only show unique insect orders, with all other arthropods in gray shade
## less busy, but still bad...

onlyInsectPal <- c("gray60","#be3c54","#83c955","gray60","gray60",
                   "#7554bf","gray60","#cab542","gray60","#b84c9e",
                   "#6ad29b","gray60","#d06e3a","#68b8c1","#7e5137",
                   "gray60","#9f9dda","#4c6f3a","#575483","gray60",
                   "#c7c18a","gray60","gray60","#d18a9b","gray60")

ggplot(tmp.table,
       aes(weight = n, axis1 = phylum_name, axis2 = class_name, axis3 = order_name, 
           axis4 = family_name, axis5 = genus_name, axis6 = species_name)) +
  geom_alluvium(aes(fill = order_name), width = 1/12) +
  scale_fill_manual(values = onlyInsectPal) +
  geom_stratum(width = 1/8, color = "darkgray") +   ## width = node/box chunk size, fill = box fill, color = box border
  geom_label(stat = "stratum", label.strata = TRUE, size = 1.8, fill = "honeydew") +             ## inserts labels for all values in a node; too many to use
  scale_x_continuous(breaks = 1:6, labels = c("phylum", "class", "order", "family", "genus", "species")) +
  ggtitle("What's for dinner or breakfast...") +
  labs(y = "number of detections")

## same plot structure, but now coloring just by taxonomic class (not order)
## simpler legend and color comprehension

septaPal <- c("#e56531", "#017ded", "#b99b00", "#d86edf", "#005e11",
              "#79003e", "#3ceebc")

ggplot(tmp.table,
       aes(weight = n, axis1 = phylum_name, axis2 = class_name, axis3 = order_name, 
           axis4 = family_name, axis5 = genus_name, axis6 = species_name)) +
  geom_alluvium(aes(fill = class_name), width = 1/12) +
  scale_fill_manual(values = septaPal) +
  geom_stratum(width = 1/8, color = "darkgray") +   ## width = node/box chunk size, fill = box fill, color = box border
  geom_label(stat = "stratum", label.strata = TRUE, size = 2, fill = "honeydew") +             ## inserts labels for all values in a node; too many to use
  scale_x_continuous(breaks = 1:6, labels = c("phylum", "class", "order", "family", "genus", "species")) +
  ggtitle("What's for dinner or breakfast...") +
  labs(y = "number of detections")

# --------------------------------------------------------------------------- #
