#' Plots connected scaffolds in a network graph
#'
#' A nice long description
#'
#' @usage mmplot_network(data, network, nconnections, ewidth, nsize, labels, color)
#'
#' @param data The dataframe with the network to be plotted.
#' @param network A network file with connections between scaffolds.
#' @param nconnections The minimum number of connections.
#' @param ewidth Scale the width of edges by connections.
#' @param nsize Scale the size of nodes by length.
#' @param labels If scaffold names are to be plotted.
#' @param color Color by gc, phylum or coverage.
#' 
#' @return An igraph network plot.
#' 
#' @export
#' @author Soren M. Karst \email{smk@@bio.aau.dk}
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}

mmplot_network <- function(data, network, nconnections = 2, ewidth = 1, nsize = 1, labels = F, color = "phylum"){

  snetwork <- subset(network, (network$scaffold1 %in% data$scaffolds$scaffold & network$scaffold2 %in% data$scaffolds$scaffold) & network$connections >= nconnections)
  
  g <- graph.data.frame(snetwork, directed = F)
  
  g.stats <- subset(data$scaffolds, data$scaffolds$scaffold %in% V(g)$name)
  g.stats <- g.stats[match(V(g)$name,g.stats$scaffold),]
  
  V(g)$size <- sqrt(g.stats$length)*0.1*nsize
  if (labels == F){V(g)$name <- ""}
  V(g)$frame.color <- "black"
  E(g)$width <- log(E(g)$connections)*ewidth 
  
  if(color == "gc"){
    rgb.c <- colorRampPalette(c("red", "green", "blue"))
    rgb.a <- adjustcolor(rgb.c(max(g.stats$gc)-min(g.stats$gc)+1), alpha.f = 1)
    palette(rgb.a)
    V(g)$color <- g.stats$gc-min(g.stats$gc)+1
  }
  if (color == "phylum"){
    gg_color_hue <- function(n) {
      hues = seq(15, 375, length=n+1)
      hcl(h=hues, l=65, c=100)[1:n]
    }
    palette(gg_color_hue(length(levels(data$scaffolds$phylum))))
    V(g)$color <- g.stats$phylum
  }
  
  if (color != "gc" & color != "phylum"){
    rgb.c <- colorRampPalette(c("red", "green", "blue"))
    rgb.a <- adjustcolor(rgb.c(100), alpha.f = 1)
    palette(rgb.a)
  
    diff <- log10(max(g.stats[,color])+1) - log10(min(g.stats[,color])+1)
    V(g)$color <- (log10(g.stats[,color]+1) -log10(min(g.stats[,color])+1)) * (99 / diff) +1
  }
  
  p <- plot(g)
  return(p)
}
