## attempt at using ggiraph package for Github Pages markdown posts

## written by: devon o'rourke
## created on: 2018-01-14
## modified on: 13-mar-2018


# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #

# run once - commented out by default
# install.packages('data.table')
# install.packages('ggplot2')
# install.packages('ggiraph')

library(data.table)
library(ggplot2)
library(ggiraph)  # see: https://davidgohel.github.io/ggiraph/reference/index.html

# --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- # --- #

# load initial data:
master.df <- fread('https://raw.githubusercontent.com/devonorourke/guano/master/Perlut/data/Routput/master.csv')

# how many unique taxa are there for each taxonomic level?
## for example: taxonomic order?
length(unique(master.df$order_name))
  # there are 35 unique taxonomic levels; that's a lot to try to visualize at once... 
## try more inclusive taxonomic class
length(unique(master.df$class_name))
  # there are just 9 ... much easier to visualize

# make a color palette to reflect 9 unique groups for plotting
## Qualitative color schemes by Paul Tol - see this post: https://tradeblotter.wordpress.com/2013/02/28/the-paul-tol-21-color-salute/
tol9 = c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#CC6677", "#882255", "#AA4499")

# order the items by the taxonomic level to be plotted (rather than by species name)


#figure with my data for bar chart:
b <- ggplot(master.df, aes(x = NestBox, fill = class_name)) +
  geom_bar_interactive(aes(data_id = species_name,
                           tooltip = species_name,
                           onclick = onclick), size = 2) +
  labs(title = "Frequency of arthropods detected in guano \n from New Hampshire birds \n grouped by taxonomic class",
       subtitle = "guano and metadata collected by Noah Perlut;\nmolecular and bioinformatic work performed by Devon O'Rourke",
       x = "Nest box location",
       y = "Number of times observed") +
  theme(legend.position = "right") +
  guides(fill=guide_legend(title="Taxonomic Class")) +
  theme(plot.subtitle=element_text(size=8, hjust=0.5, face="italic", color="black")) +
  theme(plot.caption=element_text(size=8, hjust=0.5, face="italic", color="black")) +
  scale_fill_manual(values=tol9)

ggiraph(code = print(b), zoom_max = 5)
