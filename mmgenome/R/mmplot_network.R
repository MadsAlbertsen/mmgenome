#' Plots connected scaffolds in a network graph
#'
#' A nice long description
#'
#' @usage mmplot_network(data, network)
#'
#' @param data (required) The dataframe with the network to be plotted.
#' @param network (required) A network file with connections between scaffolds.
#' @param nconnections The minimum number of connections (default: 2).
#' @param color Color by gc, phylum or a specific coverage dataset (default: "gc").
#' @param labels If scaffold names are to be plotted (default: F).
#' 
#' @return An igraph network plot.
#' 
#' @export
#' @author Soren M. Karst \email{smk@@bio.aau.dk}
#' @author Mads Albertsen \email{MadsAlbertsen85@@gmail.com}

mmplot_network <- function(data, network, nconnections = 2, labels = F, color = "gc"){

  ## Subset the network
  
  snetwork <- subset(network, network$scaffold1 %in% data$scaffolds$scaffold & network$scaffold2 %in% data$scaffolds$scaffold & network$connections >= nconnections)
  
  ## Convert to graph 
  
  g <- graph.data.frame(snetwork, directed = F)
  
  ## Calculate a layout   
  
  t <- layout.fruchterman.reingold(g)
  
  ## Extract layout coordinates
  
  gpoints <- data.frame( "scaffold" = V(g)$name, "x" = t[,1], "y" = t[,2])
  gpoints1 <- merge(gpoints, data$scaffolds)
  
  ## Extract link coordinates 
  
  links <- merge(snetwork, gpoints1[,1:3], by.x = "scaffold1", by.y = "scaffold")
  links1 <- merge(links, gpoints1[,1:3], by.x = "scaffold2", by.y = "scaffold")
  colnames(links1)[4:7] <- c("x", "y", "xend", "yend")
  
  ## Plot the data
  
  if (color == "phylum"){
    p <- ggplot(data = gpoints1, aes(x = x, y = y, size = length, color = phylum)) +
      geom_segment(data=links1, aes(x=x, y=y, xend=xend, yend=yend), color = "darkgrey", size = log10(links1$connections)) +
      geom_point(alpha=0.1, color = 'black') +
      geom_point(data=subset(gpoints1, phylum != "NA"), shape = 1) +
      scale_size_area(name= "Scaffold length", max_size=20) +
      guides(colour = guide_legend(override.aes = list(alpha = 1, size = 5, shape = 19)))
  }
  
  if (color == "gc"){
    p <- ggplot(data = gpoints1, aes(x = x, y = y, size = length, color = gc)) +
      geom_segment(data=links1, aes(x=x, y=y, xend=xend, yend=yend), color = "darkgrey", size = log10(links1$connections)) +
      geom_point(alpha=0.7) +
      scale_size_area(name= "Scaffold length", max_size=20) +
      scale_colour_gradientn(colours = c("red", "green", "blue"))
  }
  
  if(color != "gc" & color != "phylum"){
    options(digits=2)
    p <- ggplot(data=gpoints1, aes_string(x = "x", y = "y", size = "length", color = color)) +       
      geom_segment(data=links1, aes(x=x, y=y, xend=xend, yend=yend), color = "darkgrey", size = log10(links1$connections)) +
      geom_point(alpha = 0.7) +
      scale_size_area(name = "Scaffold length", max_size = 20) +
      scale_colour_gradientn(name = paste("Coverage ", color, sep=""), colours = c("red", "green", "blue"), trans = "log")
  }
  
  if (labels == T){
    p <- p + geom_text(label = gpoints1$scaffold, size = 4, color = "black")
  }
  
  p <- p + 
    #theme_bw() + 
    theme(panel.grid.major=element_blank(), 
          panel.grid.minor=element_blank(), 
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.ticks=element_blank(),
          panel.border=element_blank())
  
  return(p)
}
