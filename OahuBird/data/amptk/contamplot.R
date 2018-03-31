## info: this workflow was derived from the 'OTUanalysis.R' script for that OahuBird project

require(scales)
require(ggplot2)
require(ggrepel)
require(ggiraph)

Match_summary$onclick <- sprintf("window.open(\"%s%s\")",
                                 "http://v4.boldsystems.org/index.php/Public_BarcodeCluster?clusteruri=",
                                 as.character(Match_summary$BOLDalt))

## plot without interactive information
p <- ggplot(Match_summary, aes(x = TotalCounts, y = n))
q <- p + geom_point() +
  geom_text_repel(data=subset(Match_summary, TotalCounts > 200000 | n > 100),
                  aes(x = TotalCounts, y = n, label = BOLDalt)) +
  geom_point(data=subset(Match_summary, TotalCounts > 200000 | n > 100),
             aes(x = TotalCounts, y = n), color = 'red') +
  scale_x_continuous(labels = comma)
q

## plot with interactive info
pi <- ggplot(Match_summary, aes(x = TotalCounts, y = n)) + scale_x_continuous(labels = comma) +
  #geom_text_repel(data=subset(Match_summary, TotalCounts > 200000 | n > 100), aes(x = TotalCounts, y = n, label = BOLDalt)) +
    geom_point_interactive(
    aes( data_id = BOLDalt, tooltip = BOLDalt, onclick = onclick), size = 2 )
ggiraph(code = print(pi), hover_css = "fill:blue;r:7pt")
